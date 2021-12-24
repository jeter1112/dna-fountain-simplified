# encoding: utf-8

from utils.robust_solition import PRNG
from utils.droplet import Droplet
from reedsolo import RSCodec
import utils.scr_rept as sr
import utils.file_process as fp
import utils.Colorer
import logging,sys,json,numpy,random,operator,math,tqdm


class LFSR():
    def __init__(self):
        pass
    def lfsr(self, state, mask):
        #Galois lfsr:
        result = state
        nbits = mask.bit_length()-1
        while True:
            result = result << 1
            xor = result >> nbits
            if xor != 0:
                result ^= mask
            yield result

    def lfsr32p(self):
        #this function returns a hard coded polynomial (0b100000000000000000000000011000101).
        #The polynomial corresponds to 1 + x^25 + x^26 + x^30 + x^32, which is known 
        #to repeat only after 32^2-1 tries. Don't change unless you know what you are doing.
        return 0b100000000000000000000000011000101

    def lfsr32s(self):
        #this function returns a hard coded state for the lfsr (0b001010101)
        #this state is the inital position in the register. You can change it without a major implication.
        return 0b001010101

    def lfsr_s_p(self):
        return self.lfsr(self.lfsr32s(),self.lfsr32p())

class DNAFountain():

    def __init__(self, file_in, file_size, chunk_size, alpha, max_repeat,
            stop = None,rs = 0, c_dist = 0.1, delta = 0.5, 
            np = False,max_homopolymer = 3,gc = 0.05):

        #alpha is the redundency level
        #stop is whether we have a limit on the number of oligos
        #chunk_size and file_size are in bytes
        #rs is the number of bytes for reed-solomon error correcting code over gf(2^8).
        #c_dist is a parameter of the degree distribution
        #delta is a parameter of the degree distribution
        #np: should we use numpy random number generator? Faster, but incompatible for previous versions
        #max_homopolymer: the largest homopolymer allowed
        #gc: the allowable range of gc +- 50%

        #things realted to data:
        self.file_in = file_in
        self.chunk_size = chunk_size
        self.num_chunks = int(math.ceil(file_size / float(chunk_size)))
        self.file_size = file_size
        self.alpha = alpha
        self.stop = stop
        self.final = self._calc_stop()

        #things related to random mnumber generator
        self.lfsr = LFSR().lfsr_s_p() #starting an lfsr with a certain state and a polynomial for 32bits.
        self.lfsr_l = len('{:b}'.format(LFSR().lfsr32p())) - 1 #calculate the length of lsfr in bits 
        self.seed = next(self.lfsr)


        self.PRNG = PRNG(K = self.num_chunks, delta = delta, c = c_dist, np = np) #creating the solition distribution object
        self.PRNG.set_seed(self.seed)

        #things related to error correcting code:
        self.rs = rs #the number of symbols (bytes) to add
        self.rs_obj = RSCodec(self.rs)#initalizing an reed solomon object

        #things related to biological screens:
        self.gc = gc
        self.max_homopolymer = max_homopolymer
        self.tries = 0 #number of times we tried to create a droplet
        self.good = 0 #droplets that were screened successfully.
        self.oligo_l = self._calc_oligo_length()

        sr.prepare(max_repeat)


    def _calc_oligo_length(self):
        #return the number of nucleotides in an oligo:
        bits = self.chunk_size * 8 + self.lfsr_l + self.rs * 8
        return bits/4
    def _calc_stop(self):
        if self.stop is not None:
            return self.stop
        stop = int(self.num_chunks*(1+self.alpha))+1
        return stop

    def _rand_chunk_nums(self):
        #This funcation returns a subset of segments based on the solition distribution.
        #It updates the lfsr to generates a new seed.
        #This function creates a fresh seed for the droplet and primes the solition inverse cdf sampler
        self.seed = next(self.lfsr) #deploy one round of lfsr, and read the register.
        self.PRNG.set_seed(self.seed) #update the seed with the register
        d, ix_samples = self.PRNG.get_src_blocks_wrap()
        return d, ix_samples #return a list of segments.
    def droplet(self):
        #creating a droplet.
        data = None
        d, num_chunks = self._rand_chunk_nums() #creating a random list of segments.
        for num in num_chunks: #iterating over each segment
            if data is None: #first round. data payload is empty.
                data = self.file_in[num] #just copy the segment to the payload.
            else: #more rounds. Starting xoring the new segments with the payload.
                data = list(map(operator.xor, data, self.file_in[num]))
        self.tries += 1 #upadte counter.
        
        return Droplet(data = data, seed = self.seed, rs = self.rs,
                    rs_obj = self.rs_obj, num_chunks = num_chunks, degree = d)
    
    def screen(self, droplet):
        if sr.screen_repeat(droplet, self.max_homopolymer, self.gc):
        #if self.screen_obj.screen(droplet.toDNA(), self.oligo_l):
            self.good += 1
            return 1
        return 0





class Encode():
    def __init__(self, file_in, out = None, size = 128, max_homopolymer = 4, gc = 0.2, 
                 rs = 0, delta = 0.05, c_dist = 0.1, stop = None, alpha = 0.07, no_fasta = False, rand_numpy = False):
        '''file_in: file to encode,
        out: File with DNA oligos,
        size: number of information bytes per message,type = int
        max_homopolymer: the largest number of nt in a homopolymer,type = int
        gc: the fraction of gc content above/below 0.5 (example:0.1 means 0.4-0.6),type = restricted_float
        rs: Number of bytes for rs codes,type = int
        delta: Degree distribution tuning parameter,type = float
        c_dist: Degree distribution tuning parameter,type = float
        stop: Maximal number of oligos,type = int
        alpha: How many more fragments to generate on top of frst k (example: 0.1 will generate 10 percent more fragments),type = float
        no_fasta: Print oligo without a fasta header, type = bure
        rand_numpy: Uses numpy random generator. Faster but not compatible with older versions, type = bure'''

        logging.basicConfig(level=logging.DEBUG)
        
        self.file_in = file_in
        if out:
            self.out = out
        else:
            self.out = file_in+'.dna'
        self.size = size
        self.max_homopolymer = max_homopolymer
        self.gc = float(gc)
        if self.gc < 0.0 or self.gc > 1.0:
            logging.error("%s not in range [0.0, 1.0]", self.gc)
            sys.exit(0)
        self.rs = rs
        self.delta = delta
        self.c_dist = c_dist
        self.stop = stop
        #self.alpha = alpha
        #pirnt(11123123)
        self.alpha = 0.01
        self.no_fasta = no_fasta
        self.rand_numpy = rand_numpy

    def main(self):
        logging.info("Reading the file. This may take a few mintues")
        file_in = fp.write_tar(self.file_in)
        with open(self.file_in,'rb') as f:
            file_in = f.read()
        f_in, file_size = fp.read_file(file_in, self.size)

        f = DNAFountain(file_in = f_in, file_size = file_size, chunk_size = self.size, max_repeat = self.max_homopolymer,
                        rs = self.rs, max_homopolymer = self.max_homopolymer,
                        gc = self.gc, delta = self.delta, c_dist = self.c_dist, 
                        np = self.rand_numpy, alpha = self.alpha, stop = self.stop)

        logging.info("Upper bounds on packets for decoding is %d (x%f)  with %f probability\n", int(json.loads(f.PRNG.debug())['K_prime']), 
                                                                   json.loads(f.PRNG.debug())['Z'],
                                                                   json.loads(f.PRNG.debug())['delta'])
        if (self.out == '-'):
            out = sys.stdout
        else: 
            out = open(self.out, 'w')
            pbar = tqdm.tqdm(total= f.final, desc = "Valid oligos")

        used_bc = dict()
        while f.good < f.final:
            d = f.droplet()

            if f.screen(d):
                if not self.no_fasta:
                    out.write(">packet {}_{}\n".format(f.good, d.degree))
                out.write("{}\n".format(d.to_human_readable_DNA()))

                if d.seed in used_bc:
                    logging.error("Seed %d has been seen before\nDone", d.seed)
                    sys.exit(1)

                used_bc[d.seed] = 1

                if (self.out != '-'):
                    pbar.update()

        if (self.out != '-'):
            pbar.close()
        logging.info("Finished. Generated %d packets out of %d tries (%.3f)", f.good, f.tries, (f.good+0.0)/f.tries)
        logging.info("Out file's name is '%s'", self.out)
        out.close()
    
if __name__ == '__main__':
    f = 'mengnalisha.tar.gz'
    o = 'mengnalisha.tar.gz.dna'
    Encode(file_in=f,out=o,size=32,rs=2,max_homopolymer=3,gc=0.05,delta=0.001,c_dist=0.025,stop=2000,no_fasta=True).main()
	
	
	
	