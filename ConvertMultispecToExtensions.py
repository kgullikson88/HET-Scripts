import FitsUtils
import FindContinuum
from astropy.io import fits as pyfits
import sys
import os
import numpy as np
import pylab
import HelperFunctions

if __name__ == "__main__":
  fileList = []
  blazecorrect = False
  for arg in sys.argv[1:]:
    if "blaze" in arg:
      blazecorrect = True
    else:
      fileList.append(arg)
  for fname in fileList:
    outfilename = "%s-0.fits" %(fname.split(".fits")[0])
    header = pyfits.getheader(fname)
    print header
    orders = FitsUtils.MakeXYpoints(fname)
    orders = orders[::-1]    #Reverse order so the bluest order is first
    if blazecorrect:
      header = pyfits.getheader(fname)
      try:
        blazefile = "%s.fits" %header['BLAZE']
      except KeyError:
        allfiles = os.listdir("./")
        blazefile = [f for f in allfiles if "BLAZE" in f][0]
      try:
        blaze = FitsUtils.MakeXYpoints(blazefile)
        blaze = blaze[::-1]
      except IOError:
        print "Error! blaze file %s does not exist!" %blazefile
        print "Not converting file %s" %fname
        continue
    column_list = []
    for i, order in enumerate(orders):
      #This data is weird. Some parts of the extracted spectra have 0 flux on the edges
      goodindices = np.where(order.y > 1e-4)[0]
      if goodindices.size < 2:
        continue
      left, right = goodindices[0], goodindices[-1]
      

      #Blaze correction
      if blazecorrect:
        order.y /= blaze[i].y
        order.err /= blaze[i].y

      #Trim data
      order.x = order.x[left:right]
      order.y = order.y[left:right]
      order.err = order.err[left:right]
        
      order.cont = FindContinuum.Continuum(order.x, order.y, fitorder=3, lowreject=2, highreject=4)
      columns = columns = {"wavelength": order.x,
                           "flux": order.y,
                           "continuum": order.cont,
                           "error": order.err}
      column_list.append(columns)
      
    HelperFunctions.OutputFitsFileExtensions(column_list, fname, outfilename, mode="new")

      
