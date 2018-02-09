'''
    Implementation of FeiImageFile class, which reads and writes MRC + imagic, and read ser.

    Basic usage:
    reading images:
        a = FeiImageFile('filename.mrc','r')
        nImg = a.n_images()
        for i in range(nImg):
            image = a.read_image_at(i) # or read_image-> using an internal counter that is incremented
            imageMeta = a.get_metadata_at(i) # returns dictionairy with the headers with metadata, see FeiMrcHeader.
            # do something with the image
        a.close()
    reading volume:
        a.read_volume() # assumes the MRC is not a stack of images but a volume.


    writing:
        a = FeiImageFile('filename.mrc', 'w')
        a.write_image(image1)
        a.write_image(image2)
        #or
        a.write_volume(vol)
        a.close() # important! --> only then header will be filled with correct min/max/mean and size.
    
    TO DOCUMENT new features:  appending; writing/reading FEI extended header; with construction is also allowed now.
'''

import numpy
from numpy import pi, fromfile, ndarray, memmap, array, diff, unique
import os
from FeiMrcHeader import *
from FeiUtilFuncs import *
import struct


class FeiImageFile:
    # note: mirrorY for MRC is automatically determined so overridden...
    def __init__(self, filename, accessatr, writeFeiExtHeader=False, useMemMap=True, mirrorY=False):
        fileext = getFileExtension(filename)
        if (fileext == 'ser'):
            self.__filetype = 'ser'
            self.__datafilename = filename
	    useMemMap = False
        elif (fileext == 'hed' or fileext == 'img'):
            self.__filetype = 'imagic'
            self.__datafilename = filename[:-4]+'.img'
            self.__headerfilename = filename[:-4]+'.hed'
        elif (fileext == 'raw'):
            self.__filetype = 'feiraw'
            self.__datafilename = filename
        else:
            if (fileext in ['mrc','mrcs', 'rec', 'ccp4']):
                self.__datafilename = filename
            else:
                self.__datafilename = filename + '.mrc' # default extension auto add
            self.__filetype = 'mrc'
        self._openatr = accessatr
        if accessatr =='a':
            if not os.path.isfile(self.__datafilename):
                self._openatr = 'w'
                accessatr ='w'
            else:
                self._openatr = 'r+'

        if self.__filetype != 'mrc' and accessatr != 'r':
            raise ValueError('FeiImageFile can only write MRC files.')

        self.__filehandle = open(self.__datafilename, self._openatr+'b')

        self.__imcount = 0
        self.__writeCount = 0
        self.size = (0,0,0)
        self.__metadata = None
        self.__swapEndianness = False
        self.__mirrorY = mirrorY
        self.__useMemMap = useMemMap
        self.__mmap = None
        self.__min = None
        self.__max = None
        self.__sum = None
        self.__commonImgMetadata = {}
        
        # Used to read SER sequences
        self.__framestride = None

        if accessatr != 'w' and accessatr != 'w+': 
            # read or append mode
            if self.__filetype == 'mrc':
                mrc_header = read_mrc_header(self.__filehandle)
                if mrc_header:
                    self.size = (mrc_header['header']['nz'], mrc_header['header']['ny'], mrc_header['header']['nx'])
                    self.voxelType, self.sizeOfType = get_voxel_type_from_mrc_header(mrc_header)
                    self.__swapEndianness = mrc_header.get('swap_endianness', False)
                    #print 'self.__swapEndianness = ', self.__swapEndianness
                    self.__mirrorY = mrc_header.get('mirror_y', False)
                    #print 'self.__mirrorY = ', self.__mirrorY
                    if mrc_header.has_key('ext_header'):
                        # case with 'old' FEI ext header. New FEI ext header not (yet) supported
                        self.__metadata = mrc_header['ext_header']
                    elif mrc_header['header']['nversion']>=20140 and mrc_header['header']['xlen']>0:
                        # no ext header. if it is a 2014 file we can get the pixel size from the normal header, we assume
                        self.__metadata = {'pixelSpacing': 1e-10*mrc_header['header']['xlen']/mrc_header['header']['nx']}
                    self.__memmapoffset = int(mrc_header_size + mrc_header['header']['next'])
                    self.__mrc_header = mrc_header
                elif accessatr=='a':
                    # situation in which a file existed but it was empty or a corrupt noncomplete header. assume it is empty and go to write mode.
                    accessatr = 'w'
                    self.__filehandle.seek(0)
                else:
                    raise IOError('','Corrupt MRC file '+filename)
            else:
                if self.__filetype == 'ser':
                    theHeader = read_ser_header(self.__filehandle)
                elif self.__filetype == 'imagic': 
                    theHeader = parse_imagic_header(self.__headerfilename)
                else: # self.__filetype == 'feiraw':
                    theHeader = read_feiraw_header(self.__filehandle)
                self.__memmapoffset = theHeader['dataStartPosition']
                self.size = theHeader['size']
                self.voxelType = theHeader['dtype']
                self.sizeOfType = theHeader['sizeOfType']
                self.__memmapoffset = theHeader['dataStartPosition']
                self.__framestride = theHeader.get('framestride')

        if accessatr != 'w' and accessatr != 'w+': #recheck for non-write mode, accessatr can be changed.
            actualFileSize = os.path.getsize(self.__datafilename)
            if self.__memmapoffset + self.size[0]*self.size[1]*self.size[2]*self.sizeOfType > actualFileSize:
                print self.__memmapoffset + self.size[0]*self.size[1]*self.size[2]*self.sizeOfType, actualFileSize
                newsize = (actualFileSize - self.__memmapoffset) / (self.size[1]*self.size[2]*self.sizeOfType)
                print 'Warning: header of "',filename,'" indicated file containing ',self.size[0], ' images, while file size suggests it only contains ', newsize,' images.'
                self.size = (newsize, self.size[1], self.size[2])
            self.__header_to_be_written = False
        else:
            # write mode; the file is new
            self.__header_to_be_written = True
            self.size = None
            self.voxelType = None
        if accessatr=='a':
            self.__filehandle.seek(0, 2) #go to files end
            self.__imcount = self.size[0]
            self.__min = self.__mrc_header['header']['amin'] 
            self.__max = self.__mrc_header['header']['amax'] 
            self.__sum = self.__mrc_header['header']['amean'] * (self.size[0] * self.size[1] * self.size[2])
            self.__writeCount = self.size[0]
        self.__writemode = accessatr!='r'
        if writeFeiExtHeader == True:
            self.__writeFeiExtHeader = 128 # the default
        else:
            self.__writeFeiExtHeader = writeFeiExtHeader 
        if self.__framestride is None:
            #easier to precompute
            if not self.voxelType is None:
                self.__framestride = self.size[1]*self.size[2]*self.sizeOfType

    def read_image(self):
        self.__imcount  += 1
        if self.__imcount <= self.size[0]:
            return self.read_image_at(self.__imcount)
        return None

    def read_image_at(self, p):
        if p < 1 or p>self.size[0]:
            raise ValueError('FeiImageFile::read_image_at: index of %i is outside the total stack containing %i images' %(p, self.size[0]))
        if self.__useMemMap:
            if self.__mmap is None:
                self.__openmemmap()
            self.__imcount = p;
            r = self.__mmap[p-1]
            if not self.__swapEndianness and not self.__mirrorY:
                r = numpy.array(r) # make a copy of the memmap to a normal numpy array to prevent all kinds of complications later.
        else:
            self.__filehandle.seek((p-1)*self.__framestride + self.__memmapoffset)
            r = numpy.fromfile(self.__filehandle, dtype=self.voxelType, count=self.size[1]*self.size[2]).reshape((self.size[1], self.size[2]))
        if self.__swapEndianness:
            r = numpy.ndarray.byteswap(r)
        if self.__mirrorY:
            r = numpy.flipud(r)
        return r

    def get_metadata_at(self, p):
        if p < 1 or p>self.size[0]:
            raise ValueError('FeiImageFile::get_metadata_at: index of %i is outside the total stack containing %i images' %(p, self.size[0]))
        if self.__metadata is None:
            return {}
        elif type(self.__metadata)==dict:
            return self.__metadata # common for all the images
        else: # it will be a list
            return self.__metadata[p-1]

    def read_volume(self):
        if self.__useMemMap:
            if self.__mmap is None:
                self.__openmemmap()
            r = self.__mmap
            if self.__mirrorY:
                r = numpy.array(r) # make a copy of the memmap to a normal numpy array to prevent all kinds of complications later.
        else:
            self.__filehandle.seek(self.__memmapoffset)
            if (self.__framestride != self.size[1]*self.size[2]*self.sizeOfType):
                r = numpy.empty([self.size[0], self.size[1], self.size[2]], dtype=self.voxelType)
                for i in range(self.size[0]):
                    self.__filehandle.seek(i*self.__framestride + self.__memmapoffset)
                    r[i] = numpy.fromfile(self.__filehandle, dtype=self.voxelType, count=self.size[1]*self.size[2]).reshape((self.size[1], self.size[2]))
            else:
                r = numpy.fromfile(self.__filehandle, dtype=self.voxelType, count=self.size[0]*self.size[1]*self.size[2]).reshape((self.size[0], self.size[1], self.size[2]))
        if self.__mirrorY:
            for i in range(len(r)):
                r[i] = numpy.flipud(r[i])
        if self.__swapEndianness:
            r = numpy.ndarray.byteswap(r)
        return r


    def get_mmap(self):
        return self.__mmap

    def get_size(self):
        return self.size

    def get_type(self):
        return self.voxelType

    def n_images(self):
        return self.size[0]

    # optionally call prior to the first write (not needed anymore, was neccesary to get correct header)
    def set_size(self, size):
        self.size = size

    def get_mrcHeader(self):
        return self.__mrc_header

    # EXPERIMENTAL STUFF - not fool proof. only call it before writing any image
    def set_mrcHeader(self, header):
        write_mrc_file_header(header, self.__filehandle)
        self.__mrc_header = header
        self.__header_to_be_written = False

    def __fillExtHeader(self, mrcExtHeader, metadata):
        for key in  ['pixelSpacing', 'acceleratingVoltage','exposureTime','roughDefocus','alphaTilt','betaTilt', 'appliedDefocus']:
            if key in metadata:
                mrcExtHeader[key] = metadata[key]

    def __write_header(self, imidx, metadata):
        if self.__header_to_be_written:
            if self.__filetype in ['mrc']:
                pixelSpacing = metadata.get('pixelSpacing') if metadata else None
                self.__mrc_header = CreateDefaultMrcHeader(
                    self.size[2], self.size[1], self.size[0],
                    dtype=self.voxelType, amin=self.__min.real, amax=self.__max.real, 
                    amean=self.__sum.real/float(self.size[2]*self.size[1]*self.size[0]),
                    pixelSpacing=pixelSpacing)
                if metadata is not None and self.__writeFeiExtHeader:
                    self.__mrc_header['ext_header'] = [CreateExtHeaderEntry() for i in range (self.__writeFeiExtHeader)] #128 is FEI DEFAULT
                    self.__fillExtHeader(self.__mrc_header['ext_header'][imidx], metadata)   
                write_mrc_file_header(self.__mrc_header, self.__filehandle)
            else: #imagic
                imagic_header = CreateDefaultImagicHeader(self.size[2], self.size[1], self.size[0], dtype=self.voxelType)
                write_imagic_header(self.__headerfilename, imagic_header)
            self.__header_to_be_written = False
        elif metadata is not None and self.__writeFeiExtHeader:
            if imidx<len(self.__mrc_header['ext_header']):
                self.__fillExtHeader(self.__mrc_header['ext_header'][imidx], metadata)
            else:
                print 'FeiImageFileIO warning: can\'t write ext header for image at index',imidx,'; #extended headers is restricted to', len(self.__mrc_header['ext_header'])

    def write_image(self, im, metadata=None):
        if self.size is None:
            self.size = (1, im.shape[0], im.shape[1])
        else:
            self.size = (self.size[0] + 1, self.size[1], self.size[2])
        if self.voxelType is None:
            self.voxelType = im.dtype.name
            self.sizeOfType = im.dtype.itemsize
            self.__framestride = self.size[1] * self.size[2] * self.sizeOfType
        if self.__min is None:
            self.__min = numpy.min(im)
            self.__max = numpy.max(im)
            self.__sum = float(numpy.sum(im))
        else:
            self.__min = min(self.__min, numpy.min(im))
            self.__max = max(self.__max, numpy.max(im))
            self.__sum = self.__sum + float(numpy.sum(im))
        self.__write_header(self.__writeCount, metadata)
        if self.__mirrorY:
            ndarray.tofile(numpy.flipud(im), self.__filehandle)
        else:
            ndarray.tofile(im, self.__filehandle)
        self.__imcount = self.__imcount + 1
        self.__writeCount += 1


    # TODO issue: mean,min,max get erroneous when overwriting files!
    def write_image_at(self, index, im, metadata=None):
        if not self.size or index==self.size[0]+1:
            self.write_image(im, metadata=metadata)
        else:
            seekPosFromBack = (index-1-self.size[0]) * self.size[1]*self.size[2] * self.sizeOfType
            print seekPosFromBack
            if seekPosFromBack>0:
                 raise ValueError('FeiImageFile::write_image_at specified index %i is beyond current file containg %i images'%(index, self.size[0]))
            self.__filehandle.seek(seekPosFromBack, 2)
            if self.__min is None:
                self.__min = numpy.min(im)
                self.__max = numpy.max(im)
                self.__sum = float(numpy.sum(im))
            else:
                self.__min = min(self.__min, numpy.min(im))
                self.__max = max(self.__max, numpy.max(im))
                self.__sum = self.__sum + float(numpy.sum(im))
            self.__write_header(index-1, metadata)
            if self.__mirrorY:
                ndarray.tofile(numpy.flipud(im), self.__filehandle)
            else:
                ndarray.tofile(im, self.__filehandle)


    def write_volume(self, vol, metadata=None):
        if self.size is None:
            self.size = vol.shape
        if self.voxelType is None:
            self.voxelType = vol.dtype.name
            self.sizeOfType = vol.dtype.itemsize
            self.__framestride = self.size[1] * self.size[2] * self.sizeOfType
        if self.__min is None:
            self.__min = numpy.min(vol)
            self.__max = numpy.max(vol)
            self.__sum = numpy.sum(vol)
        else:
            self.__min = min(self.__min, numpy.min(vol))
            self.__max = max(self.__max, numpy.max(vol))
            self.__sum = self.__sum + numpy.sum(vol)
        self.__write_header(0, metadata)
        if self.__mirrorY:
            for im in vol:
                ndarray.tofile(numpy.flipud(im), self.__filehandle)
        else:
            ndarray.tofile(vol, self.__filehandle)
        self.__imcount = self.__imcount + len(vol)
        self.__writeCount +=  len(vol)

    def get_last_image_name(self):
        return self.__datafilename+','+str(self.__imcount)

    def get_filename(self):
        return self.__datafilename

    def close(self, autoCorrectHeader=True):
        if autoCorrectHeader and self.__writemode and self.__filetype == 'mrc':
            self.__mrc_header['header']['nz'] = self.__imcount
            self.__mrc_header['header']['mz'] = self.__imcount
            self.__mrc_header['header']['zlen'] = self.__imcount
            self.__mrc_header['header']['amin'] = self.__min.real
            self.__mrc_header['header']['amax'] = self.__max.real
            self.__mrc_header['header']['amean'] = self.__sum.real / float(self.__writeCount * self.size[1] * self.size[2])
            self.__filehandle.seek(0)
            write_mrc_file_header(self.__mrc_header, self.__filehandle, skipExtendedHeader = not self.__writeFeiExtHeader)
        if self.__mmap is not None:
            del self.__mmap  #close did not work anymore. but this ensure garbage collection I suppose
        self.__filehandle.close()

    def __openmemmap(self):
        if self._openatr != 'r':
            raise ValueError('Reading images/volumes (and opening memory maps) currently only allowed in read mode!')
        self.__mmap = memmap(self.__filehandle, offset=self.__memmapoffset, dtype=self.voxelType, mode=self._openatr, shape=self.size)
        if (not self.__framestride is None) and (self.__framestride != self.size[1]*self.size[2]*self.sizeOfType):
            strides = array(self.__mmap.strides)
            strides[0] = self.__framestride
            self.__mmap = numpy.lib.stride_tricks.as_strided(self.__mmap, strides=strides)

    def __enter__ (self):
        return self

    def __exit__(self, exc, value, tb):
        self.close()


# SER helper functions

ser_types =  ["uint8", "uint16", "uint32", 
  "int8", "int16", "int32", "float32", "float64", 
  "complex64", "complex128"]

_sizeof_dtypes =  {"uint8":1, "uint16":2, "uint32":4, "int8":1, "int16":2, "int32":4, "float32":4, "float64":8, "complex64":8, "complex128":16}

def get_voxel_type_from_ser_header(dataType):
    result = ser_types[dataType-1]
    print '##get_voxel_type_from_ser_header ', result
    return result, _sizeof_dtypes[result]


def read_ser_header(filehandle):
    #filehandle.seek(22);
    byteOrder = struct.unpack('H', filehandle.read(2))[0]
    seriesID = struct.unpack('H', filehandle.read(2))[0]
    seriesVersion = struct.unpack('H', filehandle.read(2))[0]
    dataTypeID = struct.unpack('I', filehandle.read(4))[0] # 	0x4120 1D - 0x4122 2D 
    tagTypeID = struct.unpack('I', filehandle.read(4))[0]
    totalNumberElements = struct.unpack('I', filehandle.read(4))[0]
    validNumberElements = struct.unpack('I', filehandle.read(4))[0]

    if seriesVersion == 0x210:
        offsetArrayOffset = struct.unpack('I', filehandle.read(4))[0]
    elif seriesVersion == 0x220:
        offsetArrayOffset = struct.unpack('Q', filehandle.read(8))[0]
    else:
        raise ValueError('Unknown SER series version 0x%x'%seriesVersion)

    filehandle.seek(offsetArrayOffset)
    #dataOffset = struct.unpack('I', filehandle.read(4))[0]
    if seriesVersion == 0x210:
        dataOffsetArray = struct.unpack(totalNumberElements * 'I', filehandle.read(4 * totalNumberElements))
    else:
        dataOffsetArray = struct.unpack(totalNumberElements * 'Q', filehandle.read(8 * totalNumberElements))  

    strides = diff(dataOffsetArray)
    validNumberElements = len(numpy.where(diff(dataOffsetArray)==strides[0])[0])
    #strides = unique(diff(dataOffsetArray))
    #if len(strides) > 1:
    #    raise ValueError('Only contiguously stored data supported for SER files')
    #elif len(strides) < 1:
    #	strides = unique(dataOffsetArray)

    dataOffset = dataOffsetArray[0]
    filehandle.seek(dataOffset + 40)
    dataType = struct.unpack('H', filehandle.read(2))[0]
    dataSizeX = struct.unpack('I', filehandle.read(4))[0]
    dataSizeY = struct.unpack('I', filehandle.read(4))[0]
    dType, sizeOfType = get_voxel_type_from_ser_header(dataType)
    hed = {'size': (validNumberElements, dataSizeY, dataSizeX), 'dtype': dType, 'sizeOfType': sizeOfType, 'dataStartPosition': dataOffset + 50, 'framestride': strides[0]}
    return hed

#FEI raw helping function
# TODO check the first byes for "FEI RawImage", support other encodings.
def read_feiraw_header(filehandle):
    filehandle.seek(17)
    width = struct.unpack('I', filehandle.read(4))[0]
    height = struct.unpack('I', filehandle.read(4))[0]
    dType = "int32"
    filehandle.seek(49)
    hed = {'size': (1, height, width), 'dtype': dType, 'sizeOfType': _sizeof_dtypes[dType], 'dataStartPosition': 49}
    return hed


# imagic helper functions

def parse_imagic_header(headerfilename):
    raw_hed = numpy.fromfile(headerfilename, dtype='uint32', count = 256)
    dtype, sizeOfType = ImagicTypeToNumpyType(raw_hed[14])
    zsize = max(raw_hed[1], max(raw_hed[60],raw_hed[61])) # imagic has seperate size attribute for volumes and nr of images. I don;t care for now. so I take the max.
    print 'imagic z-size attributes tells me ',raw_hed[1], '  ',  raw_hed[60], '   ',raw_hed[61], ' so I choose as z-size : ', zsize
    hed = {'size': (zsize, raw_hed[12], raw_hed[13]), 'dtype': dtype, 'sizeOfType': sizeOfType}
    # todo: more; 'pixelSpacing'
    return hed

def write_imagic_header(headerfilename, header):
    ndarray.tofile(header, headerfilename)


def CreateDefaultImagicHeader(xs, ys, zs, dtype=None, isVolumeData=False): #todo pixelsize, other attribs
    header = numpy.zeros(zs*256, dtype='uint32')
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



# CONVENIENCE functions
def readImageFile(fn):
    fh = FeiImageFile(fn,'r')
    im = fh.read_image()
    fh.close()
    return im

def writeImageFile(fn, im):
    fh = FeiImageFile(fn,'w')
    im = fh.write_image(im)
    fh.close()
