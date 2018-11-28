#!/usr/bin/env python

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


def readSer(infile):
    # read input
    a = FeiImageFile(infile,'r')
    nz = a.n_images()
    img = a.read_volume()
    a.close()
    print "\tNumber of images:", nz
    # for i in range(nz):
    #     print i, img[i].mean(), img[i].std(),img[i].min(),img[i].max()
    return img


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
    cx0, cy0 = findBeamCenter_frame(img[0],False)
    # print 0, cx0, cy0
    dc = 64
    img_ref = gaussian_filter(img[0,cy0-dc:cy0+dc+1,cx0-dc:cx0+dc+1], 5.0)
    for i in range(1,nz):
        img_patch = img[i,cy0-dc:cy0+dc+1,cx0-dc:cx0+dc+1]
        shifts = register_translation(img_ref, img_patch, space='real')
        # print i, shifts[::-1]
        img_align[i] = np.roll(img[i],shifts.astype('int'))
        # cx, cy = findBeamCenter_frame(img_align[i], False)
        # print i, cx, cy
    return img_align

def processSer(img, gain=1.0):
    nz, ny, nx = img.shape
    # low-pass filter
    # img_lp = np.empty([nz, ny, nx])
    # for i in range(nz):
    #    img_lp[i] = gaussian_filter(img[i].astype('float'), sigma)

    # add bias to get rid of negative values
    img_lp = img - img.min()

    # gain
    if gain <> 1.0:
        img_lp = img_lp * gain
        print '\tScaling all images by  %.2f' % gain

    return img_lp


def saveImg(img, outfile, cxg, cyg, osc_range, wavelength, cameraLength):
    nz, ny, nx = img.shape

    # make a folder
    outdir,filename=os.path.split(outfile)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # set some other parameters (fixed)
    binning = 4096/nx
    pixelSize = 0.014*binning
    expTime = 1
    osc_start = 0
    # saving
    for k in range(nz):
        osc_current = k*osc_range+osc_start
        header = "BEAM_CENTER_X=%-.9g;\nBEAM_CENTER_Y=%-.9g;\n" % (cyg*pixelSize,cxg*pixelSize)
        # header = "BEAM_CENTER_X=%-.9g;\nBEAM_CENTER_Y=%-.9g;\n"%((nx-cxg)*pixelSize,cyg*pixelSize)
        header += "BIN=%dx%d;\n" % (binning,binning)
        header += "BYTE_ORDER=little_endian;\n"
        header += "DATE=%s;\n" % (datetime.datetime.now().strftime("%a %b %d %H:%M:%S %Y"))
        header += "DETECTOR_SN=unknown;\n"
        header += "DIM=2;\nDISTANCE=%-.9g;\n" % (cameraLength)
        header += "OSC_RANGE=%-.9g;\nOSC_START=%-.9g;\nPHI=%-.9g;\n" % (osc_range, osc_current, osc_current)
        header += "SIZE1=%d;\nSIZE2=%d;\n" % (nx,ny)
        header += "PIXEL_SIZE=%-.9g;\nTIME=%-.9g;\nWAVELENGTH=%-.9g;\n" % (pixelSize,expTime,wavelength)
        header += "TWOTHETA=0;\nTYPE=unsigned_short;\n}\n"
        if len(header)<492:
            header = "{\nHEADER_BYTES=512;\n"+header
            header = "{:<512}".format(header)
            # print k,header
            with open("%s//%s_%03d.img" % (outdir, filename, k+1), 'wb') as f0:
                f0.write(header)
                f0.write((img[k]+0.5).astype('uint16'))


def printHelp():
    sys.stdout.write('ser2img.py: convert .ser (TIA format) to .img\n')
    sys.stdout.write('options:\n')
    sys.stdout.write('-i, --input \t<inputfile>\n')
    sys.stdout.write('-o, --output \t<outputfile>\n')
    sys.stdout.write('-x, --cx \t<center x(px)> (required for images with a beam stopper, optional for images without a beam stopper)\n')
    sys.stdout.write('-y, --cy \t<center y(px)> (required for images with a beam stopper, optional for images without a beam stopper)\n')
    sys.stdout.write('-k,- -HT \t<voltage (kV)>\n')
    sys.stdout.write('-c, --camera \t<camera length (m)>\n')
    sys.stdout.write('-a,--osc \t<rotation angle per frame (deg)>\n')
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
    HT = 200
    cx = -1
    cy = -1
    cameraLength = -1
    osc_range = -1
    sigma = 2.0
    gain = 1.0
    opt_flags=[False,False,False]

    try:
        opts, args = getopt.getopt(argv,"hi:o:x:y:k:c:a:s:g:bld",["help","input=","output=","cx=","cy=","HT=","camera=","osc=","sigma=","gain=","bias_correction","lowpass_filtering","beam_centering"])
    except getopt.GetoptError:
        printHelp()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h','--help'):
            printHelp()
            sys.exit()
        elif opt in ("-i", "--input"):
            inputfile = arg
        elif opt in ("-o", "--output"):
            outputfile = arg
        elif opt in ('-x', '--cx'):
            cx = float(arg)
        elif opt in ('-y','--cy'):
            cy = float(arg)
        elif opt in ('-k', '--HT'):
            HT = int(arg)
        elif opt in ('-c', '--camera'):
            cameraLength = float(arg)*1000
        elif opt in ('-a', '--osc='):
            osc_range = float(arg)
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
        else:
            print opt,'unrecognized option'

    # check if mandatory values are valid
    if cameraLength < 0:
        print 'Camera length missing'
        printHelp()
        sys.exit(1)
    elif osc_range < 0:
        print 'rotation angle missing'
        printHelp()
        sys.exit(1)
    elif inputfile == '':
        print 'Input file missing'
        printHelp()
        sys.exit(1)
    elif outputfile == '':
        print 'Output file missing'
        printHelp()
        sys.exit(1)
    return inputfile, outputfile, cx, cy, HT, cameraLength, osc_range, sigma, gain, opt_flags

if __name__ == '__main__':
    
    inputfile, outputfile, cx, cy, HT, cameraLength, osc_range, sigma, gain, opt_flags = getFilenames(sys.argv[1:])

    wavelength=0.02508   #0.0197
    if HT == 300:
        wavelength = 0.0197

    try:
        print "Input .ser file:", inputfile
        img = readSer(inputfile)
    except Exception:
        print 'Reading', inputfile, 'failed'
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
    print "Rounding to unsigned int (gain=%.1f)" % gain
    img_int16 = processSer(img_lp, gain)

    if cx < 0 or cy < 0:
        cx,cy = findBeamCenter(img_align)
        print 'Finding Beam center: %.2f, %.2f' % (cx, cy)
    try:
        print "Saving to %s###.img" % outputfile
        saveImg(img_int16, outputfile, cx, cy, osc_range, wavelength, cameraLength)
    except Exception:
        print 'Saving', outputfile, 'failed'
        sys.exit(1)


else:
    print "Please do not import this program"
