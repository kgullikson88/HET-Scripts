import FitsUtils
import FindContinuum
import pyfits
import sys
import numpy
import pylab


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
    orders = FitsUtils.MakeXYpoints(fname)
    orders = orders[::-1]    #Reverse order so the bluest order is first
    if blazecorrect:
      header = pyfits.getheader(fname)
      blazefile = "%s.fits" %header['BLAZE']
      try:
        blaze = FitsUtils.MakeXYpoints(blazefile)
        blaze = blaze[::-1]
      except IOError:
        print "Error! blaze file %s does not exist!" %blazefile
        print "Not converting file %s" %fname
        continue
    for i, order in enumerate(orders):
      #This data is weird. Some parts of the extracted spectra have 0 flux on the edges
      goodindices = numpy.where(order.y > 1e-4)[0]
      left, right = goodindices[0], goodindices[-1]
      #print left, right
      #print order.y
      #print order.y[left], order.y[right]
      

      #Blaze correction
      if blazecorrect:
        #pylab.plot(order.x, order.y/order.y.mean())
        #pylab.plot(order.x, blaze[i].y)
        #pylab.show()
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

      
      if i == 0:
        FitsUtils.OutputFitsFileExtensions(columns, fname, outfilename, mode="new")
      else:
        FitsUtils.OutputFitsFileExtensions(columns, outfilename, outfilename, mode="append")
      
      #pylab.plot(order.x, order.y)
      #pylab.show()
