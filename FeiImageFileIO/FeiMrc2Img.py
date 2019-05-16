#!/usr/bin/python

from __future__ import division
from FeiImageFileIO import FeiImageFile
import sys
import os
import getopt
import datetime
import glob
import struct
from numpy import where, sqrt, average, clip, empty, diff, array


ext_header_def = [
    ('alphaTilt',         1, 'f'),  # 4-byte floating point   Alpha tilt, in degrees.
    ('integrationTime',       1, 'f'),  # 4-byte floating point   Exposure time in seconds.
    ('tilt_axis',      1, 'f'),  # 4-byte floating point   The orientation of the tilt axis in the image in degrees. Vertical to the top is 0 [deg], the direction of positive rotation is anti-clockwise.
    ('pixelSpacing',     1, 'f'),  # 4-byte floating point   The pixel size of the images in SI units (meters).
    ('acceleratingVoltage', 1, 'f'),  # 4-byte floating point   Value of the high tension in SI units (volts).
    ('cameraLength', 1, 'f'), #4-byte floating point The calibrated camera length
    ('camera', 16, 'c'),
    ('physicalPixel', 1, 'f' ),
    ('dim', 1, 'i'),
    ('binning', 1, 'i'),
    ('wavelength', 1, 'f'),
    ('noiseReduction',1,'?')
]
_sizeof_dtypes =  {
    "i":4, "f":4, "d":8, "?":1,"s":16}
ext_header_offset = {
    'alphaTilt':(100,'d'),
    'integrationTime': (419,'d'),
    'tilt_axis':(140, 'd'),
    'pixelSpacing':(156, 'd'),
    'acceleratingVoltage':(84, 'd'),
    'camera':(435, 's'),
    'binning':(427, 'i'),
    'noiseReduction':(467,'?')}

def cal_wavelength(V0):
    h=6.626e-34 #Js, Plack's constant
    m=9.109e-31 #kg, electron mass
    e=1.6021766208e-19 #C, electron charge
    c=3e8 #m/s^2, speed

    return h/sqrt(2*m*e*V0*(1+e*V0/(2*m*c*c)))*1e10 #return wavelength in Angstrom

def readMrc(infolder, use_metadata):
    # find mrc files
    mrcList = sorted(glob.glob("{}/*.mrc".format(infolder)))

    nz = len(mrcList)
    # read input
    for k in range(nz):
        a = FeiImageFile(mrcList[k],'r')
        
        if k==0:
            (n,ny,nx) = a.get_size()
            img = empty([nz, ny, nx])
            meta = []
		
        img[k] = a.read_image()
        if use_metadata:
			meta.append(read_ext_header(mrcList[k]))
        a.close()
    
    print "\tNumber of images:", nz
    return img,meta

def read_ext_header(fileName):
    ext_header = {}
    with open(fileName,'rb') as a:
        ext_header['dim'] = struct.unpack('i',a.read(4))[0]
        for key, offset in ext_header_offset.iteritems():
            a.seek(1024+offset[0])
            if 's' not in offset[1] :
                ext_header[key] = struct.unpack(offset[1],a.read(_sizeof_dtypes[offset[1]]))[0]
            else:
                ext_header[key] = ''.join(struct.unpack(offset[1]*_sizeof_dtypes[offset[1]], a.read(_sizeof_dtypes[offset[1]]))).rstrip('\x00')
        if 'Ceta' in ext_header['camera']:
            ext_header['binning'] = 4096/ext_header['dim']
            ext_header['physicalPixel'] = 14e-6
        ext_header['wavelength'] = cal_wavelength(ext_header['acceleratingVoltage'])
        ext_header['cameraLength'] = (ext_header['physicalPixel']*ext_header['binning'])/(ext_header['pixelSpacing']*ext_header['wavelength'])
        #print ext_header
    
    return ext_header

def printHelp():
    sys.stdout.write('FeiMrc2Img.py: convert .mrc to .img\n')
    sys.stdout.write('options:\n')
    sys.stdout.write('-i, --input \t<input folder> (required)\n')
    sys.stdout.write('-o, --output \t<output file pattern> (required)\n')
    sys.stdout.write('-x, --cx \t<center x(px)> (required)\n')
    sys.stdout.write('-y, --cy \t<center y(px)> (required)\n')
    sys.stdout.write('-k, --HT \t<voltage (kV)> (required if no FEI extended header)\n')
    sys.stdout.write('-c, --camera \t<camera length (m)> (required if no FEI extended header)\n')
    sys.stdout.write('-a, --osc \t<rotation angle per frame (deg)> (required if no FEI extended header)\n')
    sys.stdout.write('-F, --FEI \t<Use FEI meta data>\n')
    sys.stdout.flush()

#printHelp()

def getUserInput(argv):
    inputfile = ''
    outputfile = ''
    HT = 200
    cx = -1
    cy = -1
    cameraLength = -1
    osc_range = -1
    use_metadata = False

    try:
        opts, args = getopt.getopt(argv,"hi:o:x:y:k:c:a:F",["help","input=","output=","cx=","cy=","HT=","camera=","osc=","FEI"])
    except getopt.GetoptError:
        printHelp()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h','--help'):
            printHelp()
            sys.exit()
        elif opt in ("-i", "--input"):
            inputfile = arg
            print inputfile
        elif opt in ("-o", "--output"):
            outputfile = arg
        elif opt in ('-x', '--cx'):
            cx = float(arg)
        elif opt in ('-y','--cy'):
            cy = float(arg)
        elif opt in ('-k', '--HT'):
            HT = int(arg) * 1e3
        elif opt in ('-c', '--camera'):
            cameraLength = float(arg)
        elif opt in ('-a', '--osc='):
            osc_range = float(arg)
            gain = float(arg)
        elif opt in ('-F', '--FEI'):
            use_metadata = True
        else:
            print opt,'unrecognized option'

    # check if mandatory values are valid
    if not use_metadata:
        if cameraLength < 0:
            print 'Camera length missing'
            printHelp()
            sys.exit(1)
        if osc_range < 0:
            print 'rotation angle missing'
            printHelp()
            sys.exit(1)
    if inputfile == '':
        print 'Input file missing'
        printHelp()
        sys.exit(1)
    if outputfile == '':
        print 'Output file missing'
        printHelp()
        sys.exit(1)
    if cx < 0 or cy < 0:
        print 'beam center missing, using the center of the image instead'
    return inputfile, outputfile, cx, cy, HT, cameraLength, osc_range, use_metadata

def saveImg(img, file_handler, cxg, cyg, osc_start,osc_range, expTime, pixelSize, binning, wavelength, cameraLength):

    # saving
	header = "BEAM_CENTER_X=%-.9g;\nBEAM_CENTER_Y=%-.9g;\n" % (cyg*pixelSize,cxg*pixelSize)
	# header = "BEAM_CENTER_X=%-.9g;\nBEAM_CENTER_Y=%-.9g;\n"%((nx-cxg)*pixelSize,cyg*pixelSize)
	header += "BIN=%dx%d;\n" % (binning,binning)
	header += "BYTE_ORDER=little_endian;\n"
	header += "DATE=%s;\n" % (datetime.datetime.now().strftime("%a %b %d %H:%M:%S %Y"))
	header += "DETECTOR_SN=unknown;\n"
	header += "DIM=2;\nDISTANCE=%-.9g;\n" % (cameraLength*1000)
	header += "OSC_RANGE=%-.9g;\nOSC_START=%-.9g;\nPHI=%-.9g;\n" % (osc_range, osc_start, osc_start)
	header += "SIZE1=%d;\nSIZE2=%d;\n" % (nx,ny)
	header += "PIXEL_SIZE=%-.9g;\nTIME=%-.9g;\nWAVELENGTH=%-.9g;\n" % (pixelSize,expTime,wavelength)
	header += "TWOTHETA=0;\nTYPE=unsigned_short;\n}\n"
	if len(header)<492:
		header = "{\nHEADER_BYTES=512;\n"+header
		header = "{:<512}".format(header)
		# print k,header
	#with open("%s//%s_%03d.img" % (outdir, filename, k+1), 'wb') as f0:
	file_handler.write(header)
	file_handler.write((img+0.5).astype('uint16'))
	return
	
def most_common(lst):
    return max(set(lst), key=lst.count)

if __name__ == '__main__':
	infolder, outputfile, cx, cy, HT, cameraLength, osc_range, use_metadata = getUserInput(sys.argv[1:])
	
	
	try:
		print "Input folder:", infolder
		img,meta = readMrc(infolder, use_metadata)
	except Exception:
		print 'Reading', infolder, 'failed'
		sys.exit(1)
		

	# Take care of the negative values
	dark_noise = img[where(img<0)]
	mm = sqrt(average(dark_noise*dark_noise))
	print '\tDark nosie sigma:',mm,
	mm = min(mm * 5, abs(dark_noise.min()))
	print ', add:', mm
	clip(img+mm, 0, 8000, out=img)
	
	n, ny, nx = img.shape

	if cx < 0 or cy < 0:
		cx = nx/2
		cy = ny/2
	nz, ny, nx = img.shape

    # make a folder
	outdir,filename=os.path.split(outputfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	
	print "Saving to %s###.img" % outputfile
	if not use_metadata:
		binning = 4096/nx
		pixelSize = 0.014*binning
		wavelength = cal_wavelength(HT)
		osc_start = 0
		expTime = 1.0
	else:
		binning = meta[0]['binning']
		pixelSize = meta[0]['physicalPixel']*1e3*binning
		wavelength = meta[0]['wavelength']
		osc_start = round(meta[0]['alphaTilt'],1)
		osc_range = []
		for m in meta:
			osc_range.append(round(m['alphaTilt'],1))
		osc_range = most_common(list(diff(array(osc_range))))
		expTime = meta[0]['integrationTime']
		cameraLength = meta[0]['cameraLength']
	for k in range(n):
		with open("%s//%s%03d.img" % (outdir, filename, k+1), 'wb') as f0:
			saveImg(img[k], f0, cx, cy, osc_start+k*osc_range, osc_range, expTime, pixelSize, binning, wavelength, cameraLength)
			
	#except Exception:
	#    print 'Saving', outfolder, 'failed'
	#    sys.exit(1)


else:
    print "Please do not import this program" 
