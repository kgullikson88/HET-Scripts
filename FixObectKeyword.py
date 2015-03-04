import sys

from astropy.io import fits


if __name__ == "__main__":
    file_list = []
    for arg in sys.argv[1:]:
        if 1:
            file_list.append(arg)

    for fname in file_list:
        hdulist = fits.open(fname, mode='update')
        object_name = hdulist[0].header['object']
        new_object = object_name.split()[0].replace("_", " ")
        print object_name, new_object
        hdulist[0].header['object'] = new_object
        hdulist.flush()
        hdulist.close()
        """
        fname2 = fname.split('_telluric')[0] + '.fits'
        real_object = fits.getheader(fname2)['OBJECT']
        hdulist = fits.open(fname, mode='update')
        hdulist[0].header['object'] = real_object
        print(fname, real_object)
        hdulist.flush()
        hdulist.close()
        """
