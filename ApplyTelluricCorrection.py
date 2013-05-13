import pyfits
import FitsUtils
import sys
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import matplotlib.pyplot as plt
import DataStructures
import os
import FindContinuum
import numpy


def ReadCorrectedFile(fname):
  orders = []
  headers = []
  hdulist = pyfits.open(fname)
  numorders = len(hdulist)
  for i in range(1, numorders):
    order = hdulist[i].data
    xypt = DataStructures.xypoint(x=order.field("wavelength"),
                                  y=order.field("model"),
                                  cont=order.field("continuum"),
                                  err=order.field("error"))

    orders.append(xypt)
    headers.append(hdulist[i].header)
  return orders, headers


def Correct(original, corrected, offset=None):
  #Read in the data and model
  original_orders = FitsUtils.MakeXYpoints(original, extensions=True, x="wavelength", y="flux", errors="error", cont="continuum")
  corrected_orders, corrected_headers = ReadCorrectedFile(corrected)
  print len(original_orders), len(corrected_orders)
  if offset == None:
    offset = len(original_orders) - len(corrected_orders)
  print offset
  for i in range(offset, len(original_orders)):
    data = original_orders[i]
    data.cont = FindContinuum.Continuum(data.x, data.y)
    try:
      model = corrected_orders[i-offset]
      header = corrected_headers[i-offset]
      print "Order = %i\nHumidity: %g\nO2 concentration: %g\n" %(i, header['h2oval'], header['o2val'])
    except IndexError:
      model = DataStructures.xypoint(x=data.x, y=numpy.ones(data.x.size))
      print "Warning!!! Telluric Model not found for order %i" %i
    print data.y.shape
    print model.y.shape
    data.y /= model.y
    original_orders[i] = data.copy()
  return original_orders

if __name__ == "__main__":
  if len(sys.argv) > 2:
    original = sys.argv[1]
    corrected = sys.argv[2]
  
    outfilename = "%s_telluric_corrected.fits" %(original.split(".fits")[0])
    print "Outputting to %s" %outfilename

    corrected_orders = Correct(original, corrected, offset=None)
    
    for i, data in enumerate(corrected_orders):
      plt.plot(data.x, data.y/data.cont)
      #Set up data structures for OutputFitsFile
      columns = {"wavelength": data.x,
                 "flux": data.y,
                 "continuum": data.cont,
                 "error": data.err}
      if i == 0:
        FitsUtils.OutputFitsFileExtensions(columns, original, outfilename, mode="new")
      else:
        FitsUtils.OutputFitsFileExtensions(columns, outfilename, outfilename)
    plt.show()
