#!/usr/bin/python

from __future__ import division
import numpy as np
from FeiImageFileIO import FeiImageFile
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.measurements import center_of_mass
from scipy.signal import correlate2d
import sys
import os
import getopt
import datetime
import glob


def readMrc(infolder):
    # find mrc files
    mrcList = sorted(glob.glob("{}/*.mrc".format(infolder)))
    nz = len(mrcList)
    # read input
    for k in range(nz):
        a = FeiImageFile(mrcList[k],'r')
        if k==0:
            (n,ny,nx) = a.get_size()
            img = np.empty([nz, ny, nx])

        img[k] = a.read_image()
        a.close()
    
    print "\tNumber of images:", nz
    #for i in range(nz):
    #     print i, img[i].mean(), img[i].std(),img[i].min(),img[i].max()
    return img

def readTxt(infolder, use_stage_logging = False): 
    # find mrc files
    txtList = sorted(glob.glob("{}/*.txt".format(infolder)))
    nz = len(txtList)
    
    # read input
    paraList = np.zeros([nz,6])
    tilt_before = 0
    if use_stage_logging:
        for k in range(nz):
            with open(txtList[k],'r') as f:
                for line in f.readlines():
                    if 'CETA Binning' in line:
                        paraList[k,0] = int(line.split(':')[-1])
                        #binning = int(line.split(':')[-1])
                    if 'Current Stage Tilt' in line:
                        #osc_current = float(line.split(':')[-1])
                        paraList[k,3] = float(line.split(':')[-1])
                        paraList[k-1,2] = paraList[k,3] - tilt_before
                        tilt_before = paraList[k,3]
                    if 'AccelerationVoltage:200000,' in line:
                        #wavelength = 0.02508
                        paraList[k,5] = 0.02508
                    elif 'AccelerationVoltage:300000,' in line:
                        #wavelength = 0.0197
                        paraList[k,5] = 0.0197
                    if 'Cameralength' in line:
                        #camera_length = float(line[:-3].split(':'))
                        paraList[k,1] = float(line[:-3].split(':')[-1])*1000
                    if 'Frame Time' in line:
                        #expTime = float(line.split(':')[-1])
                        paraList[k,4] = float(line.split(':')[-1])
        paraList[nz-1,2] = paraList[nz-2,2] 
    else:
        for k in range(nz):
            with open(txtList[k],'r') as f:
                for line in f.readlines():
                    if 'CETA Binning' in line:
                        paraList[k,0] = int(line.split(':')[-1])
                        #binning 
                    if 'Tilt start angle' in line:
                        #osc_current 
                        osc_start = float(line.split(':')[-1])
                    if 'Tilt Per Image' in line:
                        #osc_range 
                        paraList[k,2] = float(line.split(':')[-1])
                    if 'AccelerationVoltage:200000,' in line:
                        #wavelength 
                        paraList[k,5] = 0.02508
                    elif 'AccelerationVoltage:300000,' in line:
                        #wavelength 
                        paraList[k,5] = 0.0197
                    if 'Cameralength' in line:
                        #camera_length 
                        paraList[k,1] = float(line[:-3].split(':')[-1])*1000
                    if 'Frame Time' in line:
                        #expTime 
                        paraList[k,4] = float(line.split(':')[-1])
                paraList[k,3] = osc_start+k*paraList[k,2] 
    #for k in range(nz):
    #    print txtList[k], paraList[k,3],paraList[k,2]
    #print 'binning, cameralength, osc_range, osc_current, expTime, wavelength'
    return paraList



def darkCorrecion(img, direction='Y'):
    nz, ny, nx = img.shape
    img_dark = np.empty([nz,ny,nx])
    # first check if bias is drifting
    mm0 = []
    for i in [0,nz-1]:
        ind = np.where(img[i] < 0)
        mm0.append(img[i][ind].mean())
    #print mm0
    mm_diff = abs(np.diff(mm0))
    if mm_diff > 5:
        print '\tDark Correction Enabled: Bias drifted %.2f, adding %.2f ~ %.2f to each ADC block'%(mm_diff,-mm0[0],-mm0[1])
        dk = int(ny/32)
        for i in range(nz):
            for k in range(32):
                if direction == 'X':
                    img_32 = img[i][k*dk:k*dk+dk, :]
                    ind = np.where(img_32 < 0)
                    mm = img_32[ind].mean()
                    img_dark[i][k * dk:k * dk + dk, :] = img[i][k*dk:k*dk+dk, :] - mm
                elif direction == 'Y':
                    img_32 = img[i][:, k*dk:k*dk+dk]
                    ind = np.where(img_32 < 0)
                    mm = img_32[ind].mean()
                img_dark[i][:, k * dk:k * dk + dk] = img[i][:, k*dk:k*dk+dk] - mm
    else:
        print '\tBias drift is not detected, will not correct.'
        img_dark = img
    return img_dark


def findBeamCenter_frame(img_mm, findCenterofMass=False):
    ny,nx = img_mm.shape
    img_mm = gaussian_filter(img_mm, 5.0)
    dc = 128
    patch = img_mm[int(-dc+ny/2):int(dc+ny/2),int(-dc+nx/2):int(dc+nx/2)]
    ind = np.argmax(patch)
    cy = int(0.5 * ind / dc)
    cx = ind % (2 * dc)

    if findCenterofMass:
        dcc = 32
        patch_center=patch[cy-dcc:cy+dcc,cx-dcc:cx+dcc]
        (cyg,cxg)=center_of_mass(patch_center)
        cxg = cxg+cx-dcc+int(-dc+nx/2)
        cyg = cyg+cy-dcc+int(-dc+ny/2)
    else:
        cxg = cx+int(-dc+nx/2)
        cyg = cy+int(-dc+ny/2)
    return cxg, cyg


def findBeamCenter(img,findCenterOfMass=True):
    ## average all frames for beam center
    img_mm = np.mean(img,0)
    cxg, cyg = findBeamCenter_frame(img_mm, findCenterOfMass)


    return cxg, cyg

def register_translation(src_image, target_image,space="real"):
    """
    Copied and modified from skimage
    """
    # images must be the same shape
    if src_image.shape != target_image.shape:
        raise ValueError("Error: images must be same size for "
                         "register_translation")

    # only 2D data makes sense right now
    if src_image.ndim != 2 and upsample_factor > 1:
        raise NotImplementedError("Error: register_translation only supports "
                                  "subpixel registration for 2D images")

    # assume complex data is already in Fourier space
    if space.lower() == 'fourier':
        src_freq = src_image
        target_freq = target_image
    # real data needs to be fft'd.
    elif space.lower() == 'real':
        src_image = np.array(src_image, dtype=np.complex128, copy=False)
        target_image = np.array(target_image, dtype=np.complex128, copy=False)
        src_freq = np.fft.fftn(src_image)
        target_freq = np.fft.fftn(target_image)
    else:
        raise ValueError("Error: register_translation only knows the \"real\" "
                         "and \"fourier\" values for the ``space`` argument.")

    # Whole-pixel shift - Compute cross-correlation by an IFFT
    shape = src_freq.shape
    image_product = src_freq * target_freq.conj()
    cross_correlation = np.fft.ifftn(image_product)

    # Locate maximum
    maxima = np.unravel_index(np.argmax(np.abs(cross_correlation)),
                              cross_correlation.shape)
    midpoints = np.array([np.fix(axis_size / 2) for axis_size in shape])

    shifts = np.array(maxima, dtype=np.float64)
    shifts[shifts > midpoints] -= np.array(shape)[shifts > midpoints]

    # If its only one row or column the shift along that dimension has no
    # effect. We set to zero.
    for dim in range(src_freq.ndim):
        if shape[dim] == 1:
            shifts[dim] = 0

    return shifts


def alignFrames(img):
	nz, ny, nx = img.shape
	img_align = np.empty([nz,ny,nx])
	img_align[0] = img[0].copy()
	dc = 64
	
	
	#for k in range(nz):
	#	print k,findBeamCenter_frame(img[k],True)

	cx0, cy0 = findBeamCenter_frame(img[0],False)
	sigma = 2.0
	if ny == 4096:
		sigma = 6.0
	elif ny == 2048:
		sigma = 3.0

	img_ref = gaussian_filter(img[0,cy0-dc:cy0+dc+1,cx0-dc:cx0+dc+1], sigma)
	shift_max=0
	for k in range(1,nz):
		img_patch = img[k,cy0-dc:cy0+dc+1,cx0-dc:cx0+dc+1]
		shifts = register_translation(img_ref, img_patch, space='real')
		if abs(shifts).max() > shift_max:
			shift_max = abs(shifts).max()
		img_align[k] = np.roll(img[k],shifts.astype('int'))
		#cx, cy = findBeamCenter_frame(img_align[i], False)
		#print i, cx, cy
	print '\tMaximum center beam shift', shift_max,'pixel(s)'
	return img_align


def saveImg(img, para, outfile, cxg, cyg):
	nz, ny, nx = img.shape

	# make a folder
	#outdir,filename=os.path.split(outfile)
	if not os.path.exists(outfile):
		os.makedirs(outfile)

	#parameters: [binning, cameralength, osc_range, osc_current, expTime, wavelength]
	pixelSize = 0.014
	# saving
	for k in range(nz):
		header = "BEAM_CENTER_X=%-.9g;\nBEAM_CENTER_Y=%-.9g;\n" % (cxg*pixelSize*para[k,0],(ny-cyg)*pixelSize*para[k,0])
		header += "BIN=%dx%d;\n" % (para[k,0], para[k,0])
		header += "BYTE_ORDER=little_endian;\n"
		header += "DATE=%s;\n" % (datetime.datetime.now().strftime("%a %b %d %H:%M:%S %Y"))
		header += "DETECTOR_SN=unknown;\n"
		header += "DIM=2;\nDISTANCE=%-.9g;\n" % (para[k,1])
		header += "OSC_RANGE=%-.9g;\nOSC_START=%-.9g;\nPHI=%-.9g;\n" % (para[k,2], para[k,3], para[k,3])
		header += "SIZE1=%d;\nSIZE2=%d;\n" % (nx,ny)
		header += "PIXEL_SIZE=%-.9g;\nTIME=%-.9g;\nWAVELENGTH=%-.9g;\n" % (pixelSize*para[k,0],para[k,4],para[k,5])
		header += "TWOTHETA=0;\nTYPE=unsigned_short;\n}\n"
		if len(header)<492:
			header = "{\nHEADER_BYTES=512;\n"+header
			header = "{:<512}".format(header)
			# print k,header
			with open("%s//Image_%03d.img" % (outfile, k+1), 'wb') as f0:
				f0.write(header)
				f0.write((img[k]+0.5).astype('uint16'))


def printHelp():
    sys.stdout.write('ser2img.py: convert .ser (TIA format) to .img\n')
    sys.stdout.write('options:\n')
    sys.stdout.write('-i, --input \t<inputFolder>\n')
    sys.stdout.write('-o, --output \t<outputfile>\n')
    sys.stdout.write('-x, --cx \t<center x(px)> (required for images with a beam stopper, optional for images without a beam stopper)\n')
    sys.stdout.write('-y, --cy \t<center y(px)> (required for images with a beam stopper, optional for images without a beam stopper)\n')
    sys.stdout.write('-s, --sigma \t<sigma of Gaussian filter> (default 2.0)\n')
    sys.stdout.write('-g, --gain \t<multiply the data by the provided gain> (default 1.0)\n')
    # print '-m, --mode \t '
    sys.stdout.write('-b, --bias_correction\t (Optional processing) Bias correction, sometimes the bias drifted during the ')
    sys.stdout.write('exposure and you observe strips in the image, this may help\n')
    sys.stdout.write('-d, --beam_centering\t (Optional processing) align the beam center, sometimes the beam center drifted during')
    sys.stdout.write(' the exposure, if it is not blocked by the beam stopper and not saturated, this simple routine can fix it\n')
    sys.stdout.write('-l, --lowpass_filtering\t(Optional processing) low-pass filtering, ')
    sys.stdout.write('may help in spots finding and indexing)\n')
    sys.stdout.flush()


def getFilenames(argv):
    inputfile = ''
    outputfile = ''
    cx = -1
    cy = -1
    sigma = 2.0
    gain = 1.0
    opt_flags=[False,False,False,False]

    try:
        opts, args = getopt.getopt(argv,"hi:o:x:y:s:g:ubld",["help","input=","output=","cx=","cy=","sigma=","gain=","use_stage_logging","bias_correction","lowpass_filtering","beam_centering"])
    except getopt.GetoptError:
        printHelp()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h','--help'):
            printHelp()
            sys.exit()
        elif opt in ("-i", "--input pattern"):
            inputfile = arg
        elif opt in ("-o", "--output pattern"):
            outputfile = arg
        elif opt in ('-x', '--cx'):
            cx = float(arg)
        elif opt in ('-y','--cy'):
            cy = float(arg)
        elif opt in ('-k', '--HT'):
            HT = int(arg)
        elif opt in ('-s', '--sigma='):
            sigma = float(arg)
        elif opt in ('-g', '--gain='):
            gain = float(arg)
        
        elif opt in ('-b', '--bias_correction'):
            opt_flags[0] = True
        elif opt in ('-l', '--lowpass_filtering'):
            opt_flags[2] = True
        elif opt in ('-d', '--beam_centering'):
            opt_flags[1] = True
        elif opt in ('-u', '--use_stage_logging'):
            opt_flags[3] = True
        else:
            print opt,'unrecognized option'

    if inputfile == '':
        print 'Input file missing'
        printHelp()
        sys.exit(1)
    if outputfile == '':
        print 'Output file missing'
        printHelp()
        sys.exit(1)
    return inputfile, outputfile, cx, cy, sigma, gain, opt_flags

if __name__ == '__main__':
	infolder, outfolder, cx, cy, sigma, gain, opt_flags = getFilenames(sys.argv[1:])

	try:
		print "Input folder from prototype:", infolder
		img = readMrc(infolder)
		para = readTxt(infolder, opt_flags[3])
		if len(img) <> len(para):
			print "MetaData doesn't match the Data. "
			sys.exit(1)
	except Exception:
		print 'Reading', infolder, 'failed'
		sys.exit(1)

	# Dark correction
	if opt_flags[0]:
		print "Bias Correction"
		img_dark = darkCorrecion(img,'X')
	else:
		img_dark = img

	# Align the center
	if opt_flags[1]:
		print "Align beam center"
		img_align = alignFrames(img_dark)
	else:
		img_align = img_dark
		
	# low_pass filter:
	if opt_flags[2]:
		print "Low-pass filter (Gaussian kernel, sigma=%.1f)" % sigma
		img_lp = gaussian_filter(img_align, sigma)
	else:
		img_lp = img_align
		

	# Take care of the negative values
	if img_lp.min() < 0:
		print 'Negative Values found:'
		img_int16 = img_lp - img.min()
		print '\tAdding all images by  %.2f' % -img_lp.min()
		# gain
		if gain <> 1.0:
			img_int16 = img_int16 * gain
			print '\tScaling all images by  %.2f' % gain
	else:
		img_int16 = img_lp

	if cx < 0 or cy < 0:
		cx,cy = findBeamCenter(img_align)
		print 'Finding Beam center: %.2f, %.2f' % (cx, cy)
	#try:
	print "Saving to %s###.img" % outfolder
	saveImg(img_int16, para, outfolder, cx, cy)
	#except Exception:
	#    print 'Saving', outfolder, 'failed'
	#    sys.exit(1)


else:
    print "Please do not import this program" 
