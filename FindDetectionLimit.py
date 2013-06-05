"""
  THIS FUNCTION IS NOT FOR THE SENSITIVITY ANALYSIS!	
  This function takes a spectrum with a real secondary in it,
    and adds various levels of Gaussian random noise, then attempts
    to still detect the signal.
"""
import sys
import os
import astropy
import DataStructures
import FitsUtils
import matplotlib.pyplot as plt
import Correlate
from astropy import units
import numpy
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import FindContinuum

badregions = [[0, 466],
              [567.5, 575.5],
              [587.5, 593],
              [627, 634.5],
              [686, 706],
              [716, 742],
              [759, 9e9]]

if __name__ == "__main__":
  fileList = []
  model_dir = "%s/School/Research/Models/Sorted/Stellar/Vband/" %(os.environ["HOME"])
  modelfile = "lte57-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted"
  trimsize = 100
  correlation_dir = "DetectionLimitingCCFS/"
  for arg in sys.argv[1:]:
    if 'modelfile' in arg:
      modelfile = arg.split("=")[-1]
    elif "modeldir" in arg:
      model_dir = arg.split("=")[-1]
    elif "trim" in arg:
      trimsize = int(arg.split("=")[-1])
    else:
      fileList.append(arg)

  Correlate.ensure_dir(correlation_dir)
  snrs = [1,2,3,4,5,10,20, 50, 100, 200, 300, 500, 1000]
  
  #Read in the model
  x,y = numpy.loadtxt(model_dir + modelfile, usecols=(0,1), unpack=True)
  model = DataStructures.xypoint(x=x*units.angstrom.to(units.nm)/1.00026, y=10**y)
  
  #Make lists for Correlate.PyCorr2
  models = [model,]
  stars = ["5700",]
  temps = [5700, ]
  gravities = [4.0,]
  metals = [0.0,]

  for fname in fileList:
    orders = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")

    numorders = len(orders)
    for i, order in enumerate(orders[::-1]):
      DATA = interp(order.x, order.y)
      CONT = interp(order.x, order.cont)
      ERROR = interp(order.x, order.err)
      order.x = numpy.linspace(order.x[trimsize], order.x[-trimsize], order.size() - 2*trimsize)
      order.y = DATA(order.x)
      order.cont = CONT(order.x)
      order.err = ERROR(order.x)
      
      #Remove bad regions from the data
      for region in badregions:
        left = numpy.searchsorted(order.x, region[0])
        right = numpy.searchsorted(order.x, region[1])
        order.x = numpy.delete(order.x, numpy.arange(left, right))
        order.y = numpy.delete(order.y, numpy.arange(left, right))
        order.cont = numpy.delete(order.cont, numpy.arange(left, right))
        order.err = numpy.delete(order.err, numpy.arange(left, right))

      #Remove whole order if it is too small
      if order.x.size > 10:
        order.cont = FindContinuum.Continuum(order.x, order.y, lowreject=3, highreject=3)
        orders[numorders -1 -i] = order.copy()
      else:
        print "Removing order %i" %(numorders - 1 - i)
        orders.pop(numorders - 1 - i)
        
    #Add noise to the data
    for snr in snrs:
      orders2 = [order.copy() for order in orders]
      for i, order in enumerate(orders2):
        noise = numpy.random.normal(loc=0.0, scale=1.0/float(snr), size=order.x.size)
        orders2[i].y += order.cont*noise
      
      """
      for order in orders2:
        plt.plot(order.x, order.y/order.cont)
      plt.title("SNR = %i" %snr)
      plt.show()
      """
          
      Correlate.PyCorr2(orders2, resolution=60000, outdir=correlation_dir, models=models, stars=stars, temps=temps, gravities=gravities, metallicities=metals, vsini=20*units.km.to(units.cm), debug=False, outfilebase="%s_SNR%i" %(fname.split(".fits")[0], snr) )
        
        
      
      
      
      
