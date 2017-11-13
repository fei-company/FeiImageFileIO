'''
    The header and extended header definitions of the MRC format for tomography

    Usage :

    from MrcHeader import parse_header, parse_ext_header, parse_mrc_file
    # parse a whole file
    parse_mrc_file(filename)
    # to parse a normal header
    parse_header(header_str)
    # to parse an extended header
    parse_ext_header(extended_header_str)

'''
import numpy as np
from numpy import ndarray

def parse_imagic_header(headerfilename):
    raw_hed = np.fromfile(headerfilename, dtype='uint32', count = 256)
    dtype, sizeOfType = ImagicTypeToNumpyType(raw_hed[14])
    zsize = max(raw_hed[1], max(raw_hed[60],raw_hed[61])) # imagic has seperate size attribute for volumes and nr of images. I don;t care for now. so I take the max.
    print 'imagic z-size attributes tells me ',raw_hed[1], '  ',  raw_hed[60], '   ',raw_hed[61], ' so I choose as z-size : ', zsize
    hed = {'size': (zsize, raw_hed[12], raw_hed[13]), 'dtype': dtype, 'sizeOfType': sizeOfType}
    # todo: more; 'pixelSpacing'
    return hed

def write_imagic_header(headerfilename, header):
    ndarray.tofile(header, headerfilename)


def CreateDefaultImagicHeader(xs, ys, zs, dtype=None, isVolumeData=False): #todo pixelsize, other attribs
    header = np.zeros(zs*256, dtype='uint32')
    imgTypeCode = NumpyTypeToImagicType(dtype)
    for i in range(zs):
        header[1+i*256] = zs-i-1
        header[3+i*256] = 1 # NBLOCKS
        header[12+i*256] = ys
        header[13+i*256] = xs
        header[14+i*256] = imgTypeCode
        if isVolumeData:
            header[60+i*256] = zs # IZLP
            header[61+i*256] = 1  # I4LP
        else:
            header[60+i*256] = 1
            header[61+i*256] = zs
        #header[67+i*256] imavers??
        header[68+i*256] = 33686018 # machine stamp
    return header

    
# def imagicstringconvert(str):
    # r = 0
    # for i in range(4):
        # r = r + ord(str[i])*(256**i)
    # return r
    
def NumpyTypeToImagicType(numpyDtype):
    if (numpyDtype == 'int16'):
        return 1196707401 # imagicstringconvert('INTG')
    if (numpyDtype == 'uint16'): # HACK, this is not the right one!
        return 1196707401 # imagicstringconvert('INTG')
    if (numpyDtype == 'float32'):
        return 1279346002 # imagicstringconvert('REAL')
    if (numpyDtype == 'complex64'):
        return 1347243843 # imagicstringconvert('COMP')
    else:
        print "unknown Imagic type! ", numpyDtype   
        return 0
    
def ImagicTypeToNumpyType(imgType):
    if (imgType == 1196707401):
        size_of = 2
        type = 'int16'
    elif (imgType == 1279346002):
        size_of = 4
        type = 'float32'
    elif (imgType == 1347243843):
        size_of = 8
        type = 'complex64'
    else:
        size_of = 0
        type = None
        print 'unsupported Imagic type!'
    return type, size_of
    