#!/opt/local/bin/python
from astropy.io import fits as pyfits
import sys


if __name__ == "__main__":
  for fname in sys.argv[1:]:
    header = pyfits.getheader(fname)
    print fname, header["OBJECT"][:15], header["DATE-OBS"], header["UT"], header['EXPTIME']/60.0, header['IMAGETYP']

