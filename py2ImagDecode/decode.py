"""
Copyright (C) 2016 Yaniv Erlich
License: GPLv3-or-later. See COPYING file for details.
"""


from utils.robust_solition import PRNG
from utils.droplet import Droplet
from reedsolo import RSCodec
from collections import defaultdict
from shutil import copyfile
from tqdm import tqdm
import struct

import utils.Colorer
import utils.file_process as fp
import json,re,sys,os,logging,operator,numpy,random






import utils.scr_rept as sr




class Aggressive:
    def __init__(self, g, file_in, times, min_coverage = 1):
        
        self.entries = []
        self.glass = g
        self.file_in = file_in
        self.times = times
        self.on = False
        self.min_coverage = min_coverage
        self.glass_file = None
        self.md5_dict = defaultdict(list)

    def turn_on(self, coverage, seen_seed):
        self.seen_seed = seen_seed
        #coverage is too high
        if coverage > self.min_coverage:
            return 0

        return 1

    def start(self):

        #turnning on!
        logging.debug("We are in aggressive mode. Hold tight")
        self.glass_file = self.glass.save()
        logging.debug("Frozen glass at %s", self.glass_file)
        self.get_all_remaining_lines()
        logging.debug("Read all remaining lines (%d)", len(self.lines))
        self.saveme('del.agg.1.tmp')
        return self.loop()


    def loop(self):
        random.seed(1)
        for i in range(self.times):

            logging.debug("Try %d out of %d", i+1, self.times)
            g = self.load_glass(self.glass_file)
            logging.debug("Loaded glass successfully from %s", self.glass_file)
            logging.debug("Glass status is alive? %r", g.alive())
            lines = self.shuffle_lines()
            logging.debug("Shuffled lines [we have %d lines]", len(lines))
            outstring, errors, n_line = self.reciever(lines, g)
            logging.debug("Finished. Read %d additional lines. %d were rejected", n_line, errors)

            if outstring is not None:
                self.save(i, outstring)
            else:
                logging.debug("Can't decode file in this try")

        logging.debug("Finished aggressive decoding. Let's evaluate...")
        return self.find_best()

    def find_best(self):

        best = 0
        best_file = None
        for md5_str in self.md5_dict.keys():
            value = len(self.md5_dict[md5_str])
            logging.debug("For MD5 %s, we have %d successes. For example: %s", md5_str, value, self.md5_dict[md5_str][0])
            if best < value:
                best = value
                best_file = self.md5_dict[md5_str][0]

        logging.debug("Best is %s with %d successes", best_file, best)
        return best_file, best

    def save(self, index, outstring):

        outname = ''.join([self.glass_file, str(i)])
        with open(outname, "w") as o:
            o.write(outstring)
        logging.debug("Results of decoding are at %s", outname)
        o.close()
        
        md5_str = md5.new(outstring).hexdigest()
        logging.debug("Results of decoding are at %s with MD5: %s", outname, md5_str)
        self.md5_dict[md5_str].append(outname)
        return 1

    def get_all_remaining_lines(self):
        self.lines = self.file_in.readlines()

    def shuffle_lines(self):
        lines = self.lines
        shuffle(lines)
        return lines

    def load_glass(self, name):
        with open(name, 'rb') as input:
            return pickle.load(input)
        
    def reciever(self, lines, g):

        errors = 0
        n_line = 0

        logging.debug("Starting reciever. Already %d chunks are done", g.chunksDone())
        for dna in lines:

            if 'N' in dna:
                continue
            coverage, dna = dna.rstrip('\n').split()
            seed, data = g.add_dna(dna)
            n_line += 1


            if seed == -1: #reed-solomon error!
                errors += 1

            if n_line % 100 == 0:
                logging.info("After reading %d additional lines, %d chunks are done. So far: %d rejections %d barcodes", 
                             n_line, g.chunksDone(), errors, g.len_seen_seed())
            
            if g.isDone():
                logging.debug("Done![don't get too excited...]")
                break

        if not g.isDone():
             return None, errors, n_line

        return g.getString(), errors, n_line
       
    def saveme(self, name):
        with open(name, 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)
        return name


class Glass:
    def __init__(self, num_chunks, out, header_size = 4, 
                 rs = 0, c_dist = 0.1, delta = 0.05, 
                flag_correct = True, gc = 0.2, max_homopolymer = 4, 
                max_hamming = 100, decode = True, chunk_size = 32, exDNA = False, np = False, truth = None):
        
        self.entries = []
        self.droplets = set()
        self.num_chunks = num_chunks
        self.chunks = [None] * num_chunks
        self.header_size = header_size
        self.decode = decode
        self.chunk_size = chunk_size
        self.exDNA = exDNA
        self.np = np
        self.chunk_to_droplets = defaultdict(set)
        self.done_segments = set()
        self.truth = truth
        self.out = out
        self.PRNG = PRNG(K = self.num_chunks, delta = delta, c = c_dist, np = np)
        self.max_homopolymer = max_homopolymer
        self.gc = gc
        sr.prepare(self.max_homopolymer)
        self.max_hamming = max_hamming
        self.rs = rs
        self.correct = flag_correct
        self.seen_seeds = set()
        
        if self.rs > 0:
            self.RSCodec = RSCodec(rs)
        else:
            self.RSCodec = None 

    def add_dna(self, dna_string):
        #header_size is in bytes
        # data = self.dna_to_byte(dna_string)
        if self.exDNA:
            data = self._dexpandable_alphabet(dna_string, len(dna_string), n_symbols = 65, n_bytes = 21, alphabet_size = 6)
        else:
            data = self._dna_to_int_array(dna_string)

        #error correcting:
        if self.rs > 0:
            #there is an error correcting code
            if self.correct: #we want to evaluate the error correcting code
                try:

                    data_corrected = list(self.RSCodec.decode(data)[0])

                except:
                    return -1, None #could not correct the code
                #we will encode the data again to evaluate the correctness of the decoding
                data_again = list(self.RSCodec.encode(data_corrected)) #list is to convert byte array to int
                #measuring hamming distance between raw input and expected raw input
                if numpy.count_nonzero(data != list(data_again)) > self.max_hamming: 
                    #too many errors to correct in decoding                    
                    return -1, None
            else: #we don't want to evaluate the error correcting code (e.g. speed)
                data_corrected  = data[0:len(data) - self.rs] #just parse out the error correcting part
        else:
            data_corrected = data

        #seed, data = split_header(data, self.header_size)
        seed_array = data_corrected[:self.header_size]
        seed = sum([int(x)*256**i for i, x in enumerate(seed_array[::-1])])
        payload = data_corrected[self.header_size:]

        #more error detection (filter seen seeds)
        if seed in self.seen_seeds:
            return -1, None
        self.add_seed(seed)
        if self.decode:
            #create droplet from DNA
            self.PRNG.set_seed(seed)
            ix_samples = self.PRNG.get_src_blocks_wrap()[1]


            d = Droplet(payload, seed, ix_samples)
            #more error detection (filter DNA that does not make sense)
            if not sr.screen_repeat(d, self.max_homopolymer, self.gc):
                return -1, None
            self.addDroplet(d)
        return seed, data
    def _dna_to_int_array(self, dna_str):
    #convert a string like ACTCA to an array of ints like [10, 2, 4]
        num = dna_str.replace('A','0').replace('C','1').replace('G','2').replace('T','3')
        s = ''.join('{0:02b}'.format(int(num[t])) for t in range(0, len(num),1))
        data = [int(s[t:t+8],2) for t in range(0,len(s), 8)]
        return data

    def _dexpandable_alphabet(self, dna, l, n_symbols, n_bytes, alphabet_size = 6):

        res = list()
        for i in range(0,l, n_symbols):
            slice_array = dna[i:i+n_symbols]

            number_int = int(0)
            for pos,val in enumerate(slice_array):
                number_int += int(val) * 6**(n_symbols-(pos+1)) #little endian

            digits = list()
            while number_int > 0:
                digits.insert(0,int(n % 256)) #little endian.
                number_int = number_int // 256

            digits = [0] * (n_bytes-len(digits)) + digits

            res += digits

        return res
    def addDroplet(self, droplet):
        self.droplets.add(droplet)
        for chunk_num in droplet.num_chunks:
            self.chunk_to_droplets[chunk_num].add(droplet) #we document for each chunk all connected droplets        
        
        self.updateEntry(droplet) #one round of message passing

        
    def updateEntry(self, droplet):

        #removing solved segments from droplets

        for chunk_num in (droplet.num_chunks & self.done_segments):
            #if self.chunks[chunk_num] is not None:
                #we solved already this input segment. 
                
            droplet.data = list(map(operator.xor, droplet.data, self.chunks[chunk_num]))
            #subtract (ie. xor) the value of the solved segment from the droplet.
            droplet.num_chunks.remove(chunk_num)
            #cut the edge between droplet and input segment.
            self.chunk_to_droplets[chunk_num].discard(droplet)
            #cut the edge between the input segment to the droplet               

        #solving segments when the droplet have exactly 1 segment
        if len(droplet.num_chunks) == 1: #the droplet has only one input segment
            lone_chunk = droplet.num_chunks.pop() 

            self.chunks[lone_chunk] = droplet.data #assign the droplet value to the input segment (=entry[0][0])

            self.done_segments.add(lone_chunk) #add the lone_chunk to a data structure of done segments.
            if self.truth:
                self.check_truth(droplet, lone_chunk)
            self.droplets.discard(droplet) #cut the edge between the droplet and input segment
            self.chunk_to_droplets[lone_chunk].discard(droplet) #cut the edge between the input segment and the droplet
            
            #update other droplets
            for other_droplet in self.chunk_to_droplets[lone_chunk].copy():
                self.updateEntry(other_droplet)

    def getString(self):
        #return ''.join(x or ' _ ' for x in self.chunks)
        res = ''
        for x in self.chunks:
            res += ''.join(map(chr, x))

        return res



    def check_truth(self, droplet, chunk_num):
        try:
            truth_data = self.truth[chunk_num]
        except:
            logging.error("chunk: %s does not exist.", chunk_num)
            exit(1)

        if not droplet.data == truth_data:
            #error
            logging.error("Decoding error in %s.\nInput is: %s\nOutput is: %s\nDNA: %s",
                          chunk_num, truth_data, droplet.data, droplet.to_human_readable_DNA(flag_exDNA = False))
            exit(1)
        else:

            return 1
        
    def save(self):
        '''name = self.out + '.glass.tmp'
        with open(name, 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)
        return name'''
        pass
    def add_seed(self, seed):
        self.seen_seeds.add(seed)
        
    def len_seen_seed(self):
        return len(self.seen_seeds)
    def isDone(self):
        if self.num_chunks - len(self.done_segments) > 0:
            return None 
        return True
    def chunksDone(self):
        return len(self.done_segments)
    def alive(self):
        return True




class Decode():
    def __init__(self, file_in, out = None, chunk_num = 128, header_size = 4, rs = 0, delta = 0.05, 
                 c_dist = 0.1, fasta = False, no_correction = False, debug_barcodes = None, 
                 gc = 0.5, max_homopolymer = 4, mock = False, max_hamming = 100, max_line = None, 
                 expand_nt = False, size = 32, rand_numpy = False, truth = None, aggressive = None):
        '''
        file_in: file to decode
        header_size: number of bytes for the header; type = int
        chunk_num: the total number of chunks in the file; type = int
        rs: number of bytes for rs codes; type = int
        delta: Degree distribution tuning parameter; type = float
        c_dist: Degree distribution tuning parameter; type = float
        out: Output file; type = str
        fasta: Input file is FASTA
        no_correction: Skip error correcting
        debug_barcodes: Compare input barcodes to output; type = str
        gc: range of gc content; type = float
        max_homopolymer: the largest number of nt in a homopolymer; type = int
        mock: Don't decode droplets. Just evaluate correctness of data
        max_hamming: How many differences between sequenced DNA and corrected DNA to tolerate; type = int
        max_line: If defined; type = int
        expand_nt: Use a 6-nucleotide version
        size: The number of bytes of the data payload in each DNA string; type = int
        rand_numpy: Uses numpy random generator. Faster but not compatible with older versions
        truth: Reading the `true` input file. Good for debuging; type = str
        aggressive: Aggressive correction of errors using consensus in file building; type = int
        '''
        self.file_in = file_in
        if out:
            self.out = out
        else:
            self.out = temp_name('final_file')
        self.header_size = header_size
        self.chunk_num = chunk_num
        self.rs = rs
        self.delta = delta
        self.c_dist = c_dist
        self.fasta = fasta
        self.no_correction = no_correction
        self.debug_barcodes = debug_barcodes
        self.gc = gc
        self.max_homopolymer = max_homopolymer
        self.mock = mock
        self.max_hamming = max_hamming
        self.max_line = max_line
        self.expand_nt = expand_nt
        self.size = size
        self.rand_numpy = rand_numpy
        logging.basicConfig(level=logging.DEBUG)
        sys.setrecursionlimit(10000000)
        if truth is not None:
            truth = fp.write_tar(truth)
            self.truth = fp.read_file(truth, self.size)[0]
        else:
            self.truth = None

        if debug_barcodes:
            self.valid_barcodes = self._load_barcodes()
        else:
            self.valid_barcodes = None

        if aggressive:
            self.aggressive = Aggressive(g = g, file_in = f, times = aggressive)
        else:
            self.aggressive = None
            
    
    
    def _load_barcodes(self):
        valid_barcodes = dict()
        try:
            f = open(self.debug_barcodes, 'r')
        except:
            logging.error("%s file not found", self.debug_barcodes)
            sys.exit(0)
        for dna in f:
            if (re.search(r"^>", dna)):
                continue
            valid_barcodes[dna.rstrip("\n")] = 1
        return valid_barcodes

    def _aggressive_mode(self,coverage, seen_seeds):
        ####Aggresive mode
        if self.aggressive is not None and self.aggressive.turn_on(int(coverage), seen_seeds):
            best_file, value = self.aggressive.start()
            if best_file is not None:
                copyfile(best_file, self.out)
                logging.info("Done!")
            else:
                logging.error("Could not decode all file...")

            sys.exit(1)
        ### End of aggressive mode

    def _link_glass(self):
        return Glass(self.chunk_num, header_size = self.header_size, rs = self.rs, 
                     c_dist = self.c_dist, delta = self.delta, flag_correct = not (self.no_correction), 
                     gc = self.gc, max_homopolymer = self.max_homopolymer, max_hamming = self.max_hamming, 
                     decode = not (self.mock), exDNA = self.expand_nt, chunk_size = self.size, 
                     np = self.rand_numpy, truth = self.truth, out = self.out)
    def _read_file(self):
        if self.file_in == '-':
            f = sys.stdin
        else:    
            try: 
                f = open(self.file_in, 'r')
            except:
                logging.error("%s file not found", self.file_in)
                sys.exit(0)
        return f
    def main(self):
        glass = self._link_glass()
        f = self._read_file()
        line = 0
        errors = 0
        seen_seeds = defaultdict(int)

        while True:
            try:     
                dna = f.readline().rstrip('\n')
                
                if len(dna) == 0:
                    logging.info("Finished reading input file!")
                    break
                if ('N' in dna) or (self.fasta and re.search(r"^>", dna)):
                    continue
            except:
                logging.info("Finished reading input file!")
                break

            #when the file is in the format of coverage \t DNA
            if (len(dna.split()) == 2):

                coverage, dna = dna.split()
                self._aggressive_mode(coverage, seen_seeds)

            line += 1

            seed, data = glass.add_dna(dna)
            if line <3000 :
                pass
            else:
                exit()
            if seed == -1: #reed-solomon error!
                errors += 1
            else:
                if self.valid_barcodes:
                    if not dna in self.valid_barcodes:
                        logging.error("Seed or data %d in line %d are not valid:%s", seed, line, dna)
                    else:
                        seen_seeds[dna] += 1
                else:
                    seen_seeds[seed] += 1       

            if line % 1000 == 0:
                logging.info("After reading %d lines, %d chunks are done. So far: %d rejections (%f) %d barcodes", 
                             line, glass.chunksDone(), errors, errors/(line+0.0), glass.len_seen_seed())

            if line == self.max_line:
                logging.info("Finished reading maximal number of lines that is %d",self.max_line)
                break

            if glass.isDone():
                logging.info("After reading %d lines, %d chunks are done. So far: %d rejections (%f) %d barcodes", 
                             line, glass.chunksDone(), errors, errors/(line+0.0), glass.len_seen_seed())
                logging.info("Done!")
                break
        f.close()
        if not glass.isDone():
            logging.error("Could not decode all file...")
            sys.exit(1)

        outstring = glass.getString()

        with open(self.out, 'wb') as f:
            f.write(outstring)

        logging.info("MD5 is %s", fp.get_md5(outstring))
        logging.info("Out file's name is '%s',that is type of '.tar.gz'", self.out)
        json.dump(seen_seeds, open("seen_barocdes.json",'w'), sort_keys = True, indent = 4)


if __name__=='__main__':
    f = '50-SF.txt'
    o = '50-SF.jpg'

    Decode(header_size = 4, rs = 5, delta = 0.05, c_dist = 0.1, chunk_num = 1494, max_homopolymer = 3, gc = 0.05, file_in = f, out = o).main()
