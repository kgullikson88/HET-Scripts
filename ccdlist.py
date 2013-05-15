import pyfits
import sys


if __name__ == "__main__":
  for fname in sys.argv[1:]:
    header = pyfits.getheader(fname)
    print fname, header["OBJECT"]

