# -*- coding: utf-8 -*-  
import hashlib
import logging
import numpy
import sys
import tarfile
import os

def read_file(file_in, chunk_size, np = False):
    '''
    读取文件，并进行矩阵转化
    input{
            file_in: file name; type: str.
            chunk_size: number of information bytes per message; type: int.
            np: numpy request; type: hash
    };
    return{
            data_array: data array; type: list/numpy.ndarray.
            len_data: file size after zero padding; type: int.
    }
    '''

    data = file_in


    # 原始文件大小 # 文件二进制长度
    len_data = len(data)
    logging.debug("Input file has %d bytes", len_data)

    # 文件大小补齐为 chunk_size倍数，有待调试
    pad = -len_data % chunk_size
    if pad > 0:
        logging.debug("Padded the file with %d zero to have a round number of blocks of data", pad)    
        data += b"\x00" * pad #zero padding.
        len_data = len(data)

    logging.info("File MD5 is %s", get_md5(data))
    segments = len_data//chunk_size
    logging.info("There are %d input segments", segments)

    if np:
        data_array = numpy.fromstring(data, dtype = numpy.uint8)
        data_array.shape = (segments, chunk_size)
    else:
        data_array = [None] * (segments)
        for num in range(segments):
            start = chunk_size * num
            end = chunk_size * (num+1)
            chunk_binary = data[start:end]
            chunk_ords = [None] * chunk_size
            for pos in range(chunk_size):
                chunk_ords[pos] = chunk_binary[pos]
            data_array[num] = chunk_ords

    return (data_array, len_data)
def get_md5(data):
    m = hashlib.md5()
    m.update(data)
    return m.hexdigest()


def temp_name(a):
    num = 0
    name = a
    while os.path.exists(name):
        num += 1
        name = a
        name = name + '_' + str(num)
    return name

def write_tar(file_name):
    name = temp_name('temp')
    #if ' ' in file_name:
        #file_name = file_name.replace(' ','_')
    try:
        with tarfile.open(name, "w:gz") as write:
            write.add(file_name)
    except: 
        logging.error("%s file not found", file_name)
        sys.exit(0)
    with open(name,'rb') as f:
        t = f.read()
    os.remove(name)
    #pad = -len(t) % 512
    #if pad == 0:
    return t
    ##else:
    #    return t + b'\x00' * pad
    
def read_tar(file_bin):#,outfiledir = None):
    #name = temp_name('temp')
    #if outfiledir:
    #    if not os.path.exists(outfiledir):
    #        os.mkdir(outfiledir)
    #else:
    #    outfiledir = temp_name('fountain_out_filedir')
    #    os.mkdir(outfiledir)
    with open(name,'wb') as f:
        f.write(file_bin)
        
    #with tarfile.open(name, "r:gz") as read:
    #    filename = str(read.getnames()[0])
        
    #    read.extractall(outfiledir)
    #os.remove(name)
    #return outfiledir + '/' + filename
    return name
    
    
    
    