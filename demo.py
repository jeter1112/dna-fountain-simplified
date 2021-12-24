# encoding: utf-8

from time import sleep
from os import path
import math

from encode import Encode
from decode import Decode

source = 'test.txt'
encoded = 'test.encoded'
decoded = 'test.decoded'
sourcesize = (path.getsize(source))
size = 4 
chunk = math.ceil(sourcesize / size)
print('encoding................')
sleep(0.5)


Encode(file_in=source,out=encoded,
       size=size,
       rs=2,  ## reedsolo code byte number
       
       ## biological parameter
       max_homopolymer=3,
       gc=0.05,
       delta=0.001,
       c_dist=0.025,
       
       ## ensure stop > chunk size
       stop = chunk*5,
       no_fasta=True).main()
print('decoding................')
sleep(0.5)
Decode(file_in = encoded, out = decoded,
       header_size = 4,  ## seed size 
       chunk_num = chunk, ## source file chunk number
       rs = 2,  ## reedsolo code byte number
       
       delta = 0.001, 
       c_dist = 0.025,  
       max_homopolymer = 3, 
       gc = 0.05, 


       max_hamming = 0).main()
print('complete!')
print('source file:', source ,'\nencoded file:', encoded ,'\nrecover file:', decoded)