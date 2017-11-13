from distutils.core import setup

long_descr = '''This package implements a reader and writer of several image formats
used in electron microscopy. Files are memory mapped.'''

setup(name='FeiImageFileIO',
      packages=['FeiImageFileIO'],
      version="0.2",
      description="Package to read and write image files.",
      long_description = long_descr,
      author="FEI Company",
      author_email="spam@fei.com",
      url="http://www.fei.com",
      download_url="http://www.fei.com/",
      license="LGPL",
      platforms=["Any"]
      )
