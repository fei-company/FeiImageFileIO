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

import struct

''' header sizes '''
mrc_header_size = 1024
mrc_ext_header_count = 1024
mrc_ext_header_size = 128

_sizeof_dtypes =  {"uint8":1, "uint16":2, "uint32":4, "int8":1, "int16":2, "int32":4, "float32":4, "float64":8, "complex64":8, "complex128":16}


''' python struct formats for different used types '''
_2_int = 'h' # 2 bytes int
_4_int = 'i' # 4 bytes int
_4_flt = 'f' # 4 bytes float
_byte = 'B'  # 1 byte
_str = 'c'  # string

''' The MRC Header definition as list
    (attribute name, length, type)
'''
header_def = [
    ('nx',          1, _4_int),  # 4-byte integer        0   The number of pixels in the x direction of the image.
    ('ny',          1, _4_int),  # 4-byte integer        4   The number of pixels in the y direction of the image.
    ('nz',          1, _4_int),  # 4-byte integer        8   The number of pixels in the z direction of the image. Effectively in the tomography tilt series this means the number of images in the tilt series.
    ('mode',        1, _4_int),  # 4-byte integer        12  Defines the data type. Should always be 1 (2-byte integer) in the case of tomography.
    ('nxstart',     1, _4_int),  # 4-byte integer        16  set to 0 : not used; lower bound of colums
    ('nystart',     1, _4_int),  # 4-byte integer        20  set to 0 : not used; lower bound of rows
    ('nzstart',     1, _4_int),  # 4-byte integer        24  set to 0 : not used; lower bound of sections
    ('mx',          1, _4_int),  # 4-byte integer        28  set to nx : not used; grid size x
    ('my',          1, _4_int),  # 4-byte integer        32  set to ny : not used; grid size y
    ('mz',          1, _4_int),  # 4-byte integer        36  set to nz : not used; grid size z
    ('xlen',        1, _4_flt),  # 4-byte floating point 40  set to mx : not used; cell size in x Angstroms (pixel spacing=xlen/mx)
    ('ylen',        1, _4_flt),  # 4-byte floating point 44  set to my : not used; cell size in y Angstroms (pixel spacing=ylen/my)
    ('zlen',        1, _4_flt),  # 4-byte floating point 48  set to mz : not used; cell size in z Angstroms (pixel spacing=zlen/mz)
    ('alpha',       1, _4_flt),  # 4-byte floating point 52  set to 90 : not used; cell angle in degrees
    ('beta',        1, _4_flt),  # 4-byte floating point 56  set to 90 : not used; cell angle in degrees
    ('gamma',       1, _4_flt),  # 4-byte floating point 60  set to 90 : not used; cell angle in degrees
    ('mapc',        1, _4_int),  # 4-byte integer        64  set to 1 : not used; mapping colums, rows, sections on axis (x=1,y=2,z=3)
    ('mapr',        1, _4_int),  # 4-byte integer        68  set to 2 : not used; mapping colums, rows, sections on axis (x=1,y=2,z=3)
    ('maps',        1, _4_int),  # 4-byte integer        72  set to 3 : not used; mapping colums, rows, sections on axis (x=1,y=2,z=3)
    ('amin',        1, _4_flt),  # 4-byte floating point 76  minimum pixel value of all images in file
    ('amax',        1, _4_flt),  # 4-byte floating point 80  maximum pixel value of all images in file
    ('amean',       1, _4_flt),  # 4-byte floating point 84  mean pixel value of all images in file
    ('ispg',        1, _4_int),  # 2-byte integer        88  set to 0 : not used; space group number (0 for images)
    ('next',        1, _4_int),  # 4-byte integer        92  This value gives the offset (in bytes) from the end of the file header to the first dataset (image). Thus you will find the first image at 1024 + next bytes.
    ('extra',       8, _byte),  # array of 8 bytes       98  set to 0 : not used; extra 8 bytes data
    ('exttype',     4, _str),  # extended header type   104  'FEI1' in case of FEI image.
    ('nversion',    1, _4_int),  # version of MRC file.  108  year*10 + version within the year (base 0)
    ('extra2',      84, _byte),  # array of 84 bytes      112  set to 0 : not used; extra 84 bytes data
    ('origin',      3, _4_int),  # origin of fourier transform spec   196  set to 0 : not used; extra 8 bytes data
    ('map',         4, _str),  # string of 4 bytes;     208 'MAP ' for indentifiying
    ('machst',      1, _4_int),  # machine type           212 
    ('rms',         1, _4_flt), # rms deviation          216
    ('nlabl',       1, _4_int),  # 4-byte integer        220 set to 1; number of labels
    ('labl',        800, _str),  # 10 arrays of 80 characters
                                 #                       224 Arrays of characters that can be used for description.
                                 #                           Label 0 is used for copyright information (FEI)
]

''' The OLD MRC Extended Header definition as list
    (attribute name, length, type)
    TODO: the new FEI ext header if we want it
'''
ext_header_def = [
    ('alphaTilt',         1, _4_flt),  # 4-byte floating point   Alpha tilt, in degrees.
    ('betaTilt',         1, _4_flt),  # 4-byte floating point   Beta tilt, in degrees
    ('x_stage',        1, _4_flt),  # 4-byte floating point   Stage x position. Normally in SI units (meters), but some older files may be in micrometers. Check by looking at values for x,y,z. If one of these exceeds 1, it will be micrometers.
    ('y_stage',        1, _4_flt),  # 4-byte floating point   Stage y position. For testing of units see x_stage.
    ('z_stage',        1, _4_flt),  # 4-byte floating point   Stage z position. For testing of units see x_stage.
    ('x_shift',        1, _4_flt),  # 4-byte floating point   Image shift x. For testing of units see x_stage.
    ('y_shift',        1, _4_flt),  # 4-byte floating point   Image shift y. For testing of units see x_stage.
    ('roughDefocus',        1, _4_flt),  # 4-byte floating point   Defocus as read from microscope. For testing of units see x_stage.
    ('exposureTime',       1, _4_flt),  # 4-byte floating point   Exposure time in seconds.
    ('mean_int',       1, _4_flt),  # 4-byte floating point   Mean value of image.
    ('tilt_axis',      1, _4_flt),  # 4-byte floating point   The orientation of the tilt axis in the image in degrees. Vertical to the top is 0 [deg], the direction of positive rotation is anti-clockwise.
    ('pixelSpacing',     1, _4_flt),  # 4-byte floating point   The pixel size of the images in SI units (meters).
    ('magnification',  1, _4_flt),  # 4-byte floating point   The magnification used for recording the images.
    ('acceleratingVoltage', 1, _4_flt),  # 4-byte floating point   Value of the high tension in SI units (volts).
    ('binning',        1, _4_flt),  # 4-byte floating point   The binning of the CCD or STEM acquisition
    ('appliedDefocus', 1, _4_flt),  # 4-byte floating point   The intended application defocus in SI units (meters), as defined for example in the tomography parameters view.
    ('remainder',      16,_4_flt),  # 4-byte floating point values, filling up to 128 bytes
                                    #                         Not used
]

        
def get_format_string(header, endianness = '<'):
    ''' return the struct format string given the header list '''
    #return ''.join(elem[1]*elem[2] for elem in header)
    result = ''.join('%d%s' % (elem[1], elem[2]) for elem in header)
    return endianness+result

def _parse_header(header_def, header_str, endianness = '<'):
    ''' parse the header string and return as dict
        using the given header definition (list)
    '''
    ret = {}
    elems = struct.unpack(get_format_string(header_def, endianness = endianness), header_str)
    index = 0
    for attr, size, _type in header_def :
        ret[attr] = elems[index:index+size]
        # for size 1 just take take only element
        if _type == _str:
            ret[attr] = ''.join(ret[attr])
        elif( size == 1) :
            ret[attr] = ret[attr][0]
        index += size
    ret['skip' ] = False
    return ret

def parse_mrc_header(header_str):
    ''' parse the header string and return as dict '''
    std_header = _parse_header(header_def, header_str)
    extra_prop = {}
    if std_header['mode']>20:
        print 'Warning: found out by weird values that MRC is probably big endian!'
        std_header = _parse_header(header_def, header_str, endianness = '>')
        extra_prop['swap_endianness'] = True
    extra_prop['mirror_y'] = False # when exactly? TODOOOO
    return std_header, extra_prop

def parse_mrc_ext_header(header_str):
    ''' parse the header extended string and return as dict '''
    return _parse_header(ext_header_def, header_str)

def read_mrc_header(filehandle, skipExtendedHeader = False):
    ''' parse a mrc file, return headers '''
    #headerlength = mrc_header_size + mrc_ext_header_count * mrc_ext_header_size
    content = filehandle.read(mrc_header_size)
    if len(content) < mrc_header_size:
        return None
    std_header, result = parse_mrc_header(content)
    next = std_header['next']
    # print '*****writing mrc info'
    # outpf = open("d:\\mrc.txt",'w')
    # outpf.write('---------normal header---------')
    # for a in std_header:
        # outpf.write('%s: %s\n'%(a, str(std_header[a])))
    # outpf.close()
    
    result['status'] = 'ok'
    result['header'] = std_header
    
    if next>0 and not skipExtendedHeader and std_header['nversion']<20140:  
       # note:last check is a fishy way to check we have an old MRC file
       # in that case we ASSUME an old-fashioned FEI ext header
       cnt = int(next / mrc_ext_header_size)
       content = filehandle.read(next)    
       ext_header = []
       for count in range(cnt):
           ext_header.append(parse_mrc_ext_header(content[count*mrc_ext_header_size:(1+count)*mrc_ext_header_size]))
       result['ext_header'] = ext_header
    # print '*****writing mrc info'
    # outpf = open("c:\\mrc.txt",'w')
    # outpf.write('---------normal header---------')
    # for a in result['header']:
        # outpf.write('%s: %s\n'%(a, str(result['header'][a])))
    # if result.has_key('ext_header'):
        # outpf.write('---------extended header---------')
        # for a in result['ext_header'][0]:
            # print 'a = ', a
            # outpf.write('%s: %s\n'%(a, str(result['ext_header'][0][a])))
    # outpf.close()
    return result
                

def write_mrc_header(header_values, header_def, f):
    for item in header_def:
       itemkey = item[0]
       itemlength = item[1]
       itemtype = item[2]
       if (header_values.has_key(itemkey)):
          data = header_values[itemkey]
       elif itemtype == _str:
          data = ''
       else:
          data = 0
       if itemtype == _str:
          format = '%ds' % (itemlength)
          binary = struct.pack(format, data)
          f.write(binary)
       elif itemlength == 1:
          format = '%s' % (itemtype)
          binary = struct.pack(format, data)
          f.write(binary)
       else:
          format = '%s' % (itemtype)
          for i in range(0, itemlength):
              binary = struct.pack(format,  data if type(data) not in [list,tuple] else data[i])
              f.write(binary)

# def write_mrc_padding_bytes(nr_padding, f):
    # for i in range(nr_padding):
        # binary = struct.pack("B",0)
        # f.write(binary)
       
def write_mrc_file_header(mrc_header_values, f, skipExtendedHeader=False):
    header_values = mrc_header_values['header']
    extH = not skipExtendedHeader and 'ext_header' in mrc_header_values
    if extH:
        ext_header_values = mrc_header_values['ext_header']
        nextParam = len(ext_header_values) * mrc_ext_header_size
    else:
        nextParam = 0
    mrc_header_values['header']['next'] = nextParam
    write_mrc_header(header_values, header_def, f)
    if extH:
        written_ext_headers = 0
        #   print 'write_mrc_file_header: number of ext headers : ', len(ext_header_values)
        for count in range(len(ext_header_values)):
            write_mrc_header(ext_header_values[count], ext_header_def, f)
            written_ext_headers = written_ext_headers + 1
        # padding_headers = mrc_ext_header_count - written_ext_headers
        # if padding_headers > 0:
            # padding_bytes = padding_headers * mrc_ext_header_size
            # write_mrc_padding_bytes(padding_bytes, f)
    
def get_voxel_type_from_mrc_header(mrc_header):
    return MrcModeToNumpyType(mrc_header['header']['mode'])
    
def CreateDefaultMrcHeader(xs, ys, zs, dtype=None, amin=-32768, amax=32767, amean=0, pixelSpacing=None):
    mrc_header = {} #header.mrc_header
    mrc_header['header'] = {}
    mrc_header['header']['map'] = "MAP "
    mrc_header['header']['exttype'] = "\0\0\0\0"
    mrc_header['header']['nversion'] = 20140
    mrc_header['header']['machst'] = 17476
    mrc_header['header']['nx'] = xs
    mrc_header['header']['ny'] = ys
    mrc_header['header']['nz'] = zs
    mrc_header['header']['mx'] = xs
    mrc_header['header']['my'] = ys
    mrc_header['header']['mz'] = zs
    if not pixelSpacing:
        pixelSpacing = 1
    else:
        pixelSpacing *= 1e10 # meters to angstroms
    mrc_header['header']['xlen'] = xs * pixelSpacing
    mrc_header['header']['ylen'] = ys * pixelSpacing
    mrc_header['header']['zlen'] = zs * pixelSpacing
    if dtype is not None:
        mrc_header['header']['mode'] = NumpyTypeToMrcMode(dtype)
    mrc_header['header']['amin'] = amin
    mrc_header['header']['amax'] = amax
    mrc_header['header']['amean'] = amean
    mrc_header['ext_header'] = []
    return mrc_header

def CreateExtHeaderEntry():
    extHeader = {}
    extHeader['alphaTilt'] = 0
    extHeader['betaTilt'] = 0
    extHeader['x_stage'] = 0
    extHeader['y_stage'] = 0
    extHeader['z_stage'] = 0
    extHeader['x_shift'] = 0
    extHeader['y_shift'] = 0
    extHeader['roughDefocus'] = 0
    extHeader['exposureTime'] = 0
    extHeader['mean_int'] = 0
    extHeader['tilt_axis'] = 0
    extHeader['pixelSpacing'] = 0
    extHeader['magnification'] = 0
    extHeader['acceleratingVoltage'] = 0
    extHeader['binning'] = 0
    extHeader['appliedDefocus'] = 0
    extHeader['remainder'] = [0.0]*16
    return extHeader


def NumpyTypeToMrcMode(numpyDtype):
    if (numpyDtype == 'uint8'):
        return 0
    if (numpyDtype == 'int16'):
        return 1
    if (numpyDtype == 'float32'):
        return 2
    if (numpyDtype == 'complex64'):
        return 4
    elif (numpyDtype == 'uint16'):
        return 6
    else:
        print "unknown MRC mode"
        return 0
    
def MrcModeToNumpyType(mrcMode):
    type = 0
    if (mrcMode == 0):
        type = 'uint8'
    elif (mrcMode == 1):
        type = 'int16'
    elif (mrcMode == 2):
        type = 'float32'
    elif (mrcMode == 4):
        type = 'complex64'
    elif (mrcMode == 6):
        type = 'uint16'
    else:
        print 'unsupported MRC type!'
    return type, _sizeof_dtypes[type]
