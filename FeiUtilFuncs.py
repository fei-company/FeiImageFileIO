import sys
import time
import glob
import os
import numpy
import math
import fnmatch

# very basic functiosn only in this library. no ref to other Fei py files here please



# ------------------------- file name util funcs ------------------------- 


def getFileExtension(fn):
    return os.path.splitext(fn)[-1][1:].lower()

def is_EMXfile(spec):
    return getFileExtension(spec) == 'emx'

# todo rename. is misleading.
def dataFileNameFromEMX(fn, postfix = '', defaultExtension='mrc'):
    if is_EMXfile(fn):
        return os.path.splitext(fn)[0]+postfix+'.'+defaultExtension
    else:
        if not postfix:
            return fn
        a = os.path.splitext(fn)
        return a[0]+postfix+a[1]

def CreateEMXFileName(fn, postfix = ''):
    return os.path.splitext(fn)[0]+postfix+'.emx'

def CreateFileName(fn, extension, postfix = '', prefix=''):
    return prefix+os.path.splitext(fn)[0]+postfix+'.'+extension

def splitFileNameAndIndex(fn):
    commapos = fn.rfind(',')
    filename = fn
    fileidx = 1
    if commapos>=0:
        filename = fn[:commapos]
        fileidx = int(fn[commapos+1:])
    return filename, fileidx

def copyDictionairyFields(dicOut, dicIn, fields):
    for f in fields:
        if dicIn.has_key(f):
            dicOut[f] = dicIn[f]

def getDictionairyFields(dicIn, fields):
    dicOut = {}
    for f in fields:
        if dicIn.has_key(f):
            dicOut[f] = dicIn[f]
    return dicOut

def isCached(inputFile, cacheFile):
    ### Checks if the cacheFile is still of the inputFile
    ### Return False if the inputFile has been updated since the cacheFile has been updated or no cacheFile exists
    if not os.path.isfile(cacheFile):
        return False
    inputEMX = inputFile[:-4] + '.emx'
    inputMRC = inputFile[:-4] + '.mrc'
    if os.path.isfile(inputMRC):
        return os.path.getmtime(inputMRC) < os.path.getmtime(cacheFile)
    else:
        return os.path.getmtime(inputEMX) < os.path.getmtime(cacheFile)



def find_files(directory, pattern):
    dirLen = len(directory)
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename

# ------------------------- unit conversion funcs ------------------------- 


def ConvertPhysicalSizeToPixels(s, ps):
    if type(s) not in [int,float]:
        if s[-2:] == 'nm':
            return float(s[:-2]) / (ps * 1e9)
        if s[-2:] == 'um':
            return float(s[:-2]) / (ps * 1e6)
        if s[-2:] == 'mm':
            return float(s[:-2]) / (ps * 1e3)
        if s[-1] == 'A':
            return float(s[:-1]) / (ps * 1e10)
        if s[-2:] == 'px':
            return float(s[:-2])
        return float(s)
    return s

# ------------------------- basic math stuff  ------------------------- 

def RotationMatrixToEuler(r):
    theta = numpy.arccos(max(min(r[2,2], 1.0), -1.0))
    if (numpy.abs(numpy.mod(theta + numpy.pi/2, numpy.pi) - numpy.pi/2) < 1e-7):
        phi = 0.0
        psi = numpy.arctan2(-r[1,0], r[0,0])
    else:
        phi = numpy.arctan2(r[2,1], r[2,0])
        psi = numpy.arctan2(r[1,2], -r[0,2])
    return [psi, theta, phi]

# psi theta phi should be the ordering of the angles
def EulerToRotationMatrix(angles): 
    cospsi = numpy.cos(angles[0])
    sinpsi =  numpy.sin(angles[0])
    costheta = numpy.cos(angles[1])
    sintheta = numpy.sin(angles[1])
    cosphi = numpy.cos(angles[2])
    sinphi = numpy.sin(angles[2])
    return numpy.array([
        [cospsi*costheta*cosphi - sinpsi * sinphi, cospsi*costheta*sinphi + sinpsi*cosphi, -cospsi*sintheta],
        [-sinpsi*costheta*cosphi - cospsi*sinphi, -sinpsi*costheta*sinphi + cospsi*cosphi, sinpsi*sintheta],
        [sintheta*cosphi, sintheta*sinphi, costheta]])

