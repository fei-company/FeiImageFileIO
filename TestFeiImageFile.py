import numpy
from _FeiImageFileIO import *

# writes a random image mrc and reads it again
def test_image_write_read():
    print '======= test_image_write_read ======='
    im1 = numpy.random.rand(128,128).astype('float32')
    im2 = numpy.random.rand(128,128).astype('float32')
    
    print 'writing 2 images...'
    with FeiImageFile('testfile.mrc', 'w') as a:
        a.write_image(im1)
        a.write_image(im2)

    print 'reading...'
    with FeiImageFile('testfile.mrc', 'r') as b:
        print '#images',b.n_images()
        im1back = b.read_image()
        im2back = b.read_image()
        
    print 'checking...'
    diffMax1 = numpy.max(numpy.abs((im1back-im1).flatten()))
    diffMax2 = numpy.max(numpy.abs((im2back-im2).flatten()))
    print 'diffMaxs:',diffMax1,diffMax2
    testOK = diffMax1 == 0 and diffMax2 == 0
    print 'TEST OK: ', testOK
    return testOK



def test_image_append():
    print '======= test_image_append ======='
    # at the same time testwithout the useful 'with' construction
    a = FeiImageFile('testfile.mrc', 'a')
    nimg, ysize, xsize = a.get_size()
    print 'appending to file with original dimensions', nimg, ysize, xsize
    im1 = numpy.random.rand(128,128).astype(a.voxelType)
    print im1.dtype
    im2 = numpy.random.rand(128,128).astype(a.voxelType)
    a.write_image_at(1, im2) # writes the image at specified index (overwriting existing one.)
    a.write_image(im1) # always appends the image
    a.close()

    print 'reading...'
    with FeiImageFile('testfile.mrc', 'r') as b:
        print '#images',b.n_images()
        # reversed order for sake of testing image_at
        im1back = b.read_image_at(b.n_images())
        im2back = b.read_image_at(1)
        
    print 'checking...'
    diffMax1 = numpy.max(numpy.abs((im1back-im1).flatten()))
    diffMax2 = numpy.max(numpy.abs((im2back-im2).flatten()))
    print 'diffMaxs:',diffMax1,diffMax2
    testOK = diffMax1 == 0 and diffMax2 == 0
    print 'TEST OK: ', testOK
    return testOK

def test_volume_read_write():
    print '======= test_volume_read_write ======='
    vol = numpy.random.rand(32,32,32).astype('float32')
    
    print 'writing 2 images...'
    with FeiImageFile('testfileVol.mrc', 'w') as a:
        a.write_volume(vol)

    print 'reading...'
    with FeiImageFile('testfileVol.mrc', 'r') as b:
        volback = b.read_volume()
        
    print 'checking...'
    diffMax1 = numpy.max(numpy.abs((volback-vol).flatten()))
    print 'diffMax:',diffMax1
    testOK = diffMax1 == 0
    print 'TEST OK: ', testOK
    return testOK

def test_FEIextHeader():
    pass

def test_unhappy_flows():
    pass


if __name__ == '__main__':
    test_image_write_read()
    test_image_append()
    test_volume_read_write()
    test_FEIextHeader()
    test_unhappy_flows()
    