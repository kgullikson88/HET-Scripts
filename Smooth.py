import numpy
import FitsUtils
import FittingUtilities
import matplotlib.pyplot as plt
import sys
import os
from astropy import units
import DataStructures
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import MakeModel
import HelperFunctions

BadRegions = [[420.7, 421.7],
              [425,429],
              [462,466],
              [481.2, 481.7],
              [482.4, 482.7],
              [573.2, 573.6],
              [687.9, 690.6],
              [703.5, 708],
              [753.1, 755.7]]

def SmoothData(order, windowsize=91, smoothorder=5, lowreject=3, highreject=3, numiters=10, normalize=True):
  denoised = FittingUtilities.Denoise3(order.copy())
  denoised.y = FittingUtilities.Iterative_SV(denoised.y, windowsize, smoothorder, lowreject=lowreject, highreject=highreject, numiters=numiters)
  if normalize:
    denoised.y /= denoised.y.max()
  return denoised
  

def Smooth(fname, window_size=91, numiters=100, lowreject=3, highreject=3, smoothorder=3):
  orders = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
  column_list = []
  for i, order in enumerate(orders):
    #Check for bad regions in the chip (internal reflections I believe)
    for region in BadRegions:
      start, finish = region[0], region[1]
      left = numpy.searchsorted(order.x, start)
      right = numpy.searchsorted(order.x, finish)
      if right - left > 0:
        order.y[left:right] = order.cont[left:right]
    
    print "\n\n"
    done = False
    x = order.x.copy()
    y = order.y.copy()
    indices = numpy.array([True]*x.size)
    iteration = 0
    while not done and iteration < numiters:
      done = True
      iteration += 1
      smoothed = FittingUtilities.savitzky_golay(y, window_size, smoothorder)
      residuals = y/smoothed
      if iteration == 1:
        #Save the first and last points, for interpolation
        first, last = smoothed[0], smoothed[-1]
      mean = residuals.mean()
      std = residuals.std()
      print residuals.size, x.size, y.size
      #plt.plot((residuals - mean)/std)
      #plt.show()
      badindices = numpy.where( numpy.logical_or( (residuals-mean)/std > highreject, (residuals-mean)/std < -lowreject ) )[0]
      if badindices.size > 1 and y.size - badindices.size > 2*window_size and iteration < numiters:
        done = False
        y[badindices] = smoothed[badindices]
        #x = numpy.delete(x, badindices)
        #y = numpy.delete(y, badindices)

if __name__ == "__main__":
  fileList = []
  for arg in sys.argv[1:]:
    fileList.append(arg)
  if len(fileList) == 0:
    fileList = [f for f in os.listdir("./") if f.endswith("telluric_corrected.fits")]
  for fname in fileList:
    orders = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
    column_list = []
    for order in orders:
      #Linearize
      xgrid = numpy.linspace(order.x[0], order.x[-1], order.x.size)
      order = MakeModel.RebinData(order, xgrid)
      
      denoised = SmoothData(order, 61, 5, 2, 2, 10)
      #order2 = order.copy()
      #denoised = FittingUtilities.Denoise3(order2) #, snr=400.0, reduction_factor=0.15)
      #denoised.y = FittingUtilities.Iterative_SV(denoised.y, 91, 5, lowreject=2, highreject=2, numiters=10)

      column = {"wavelength": denoised.x,
                "flux": order.y/denoised.y,
                "continuum": denoised.cont,
                "error": denoised.err}
      column_list.append(column)
      plt.figure(1)
      plt.plot(order.x, order.y/order.y.mean())
      plt.plot(denoised.x, denoised.y/denoised.y.mean())
      plt.figure(2)
      plt.plot(order.x, order.y/denoised.y)
      plt.plot(order.x, (order.y-denoised.y)/numpy.median(order.y))
    plt.show()
    outfilename = "%s_smoothed.fits" %(fname.split(".fits")[0])
    print "Outputting to %s" %outfilename
    HelperFunctions.OutputFitsFileExtensions(column_list, fname, outfilename, mode='new')
