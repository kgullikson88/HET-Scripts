from astropy.io import fits
import sys


if __name__ == "__main__":
    file_list = []
    for arg in sys.argv[1:]:
        if 1:
            file_list.append(arg)

    for fname in file_list:
        hdulist = fits.open(fname, mode='update')
        object = hdulist[0].header['object']
        new_object = object.split()[0].replace("_", " ")
        print object, new_object
        hdulist[0].header['object'] = new_object
        hdulist.flush()
        hdulist.close()