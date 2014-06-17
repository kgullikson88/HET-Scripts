import numpy
#import FitsUtils
import FittingUtilities
import HelperFunctions
import matplotlib.pyplot as plt
import sys
import os
from astropy import units
from astropy.io import fits, ascii
import DataStructures
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import MakeModel
import HelperFunctions
from collections import Counter
from sklearn.gaussian_process import GaussianProcess
import warnings


def SmoothData(order, windowsize=91, smoothorder=5, lowreject=3, highreject=3, numiters=10, expand=0, normalize=True):
  denoised = HelperFunctions.Denoise(order.copy())
  denoised.y = FittingUtilities.Iterative_SV(denoised.y, windowsize, smoothorder, lowreject=lowreject, highreject=highreject, numiters=numiters, expand=expand)
  if normalize:
    denoised.y /= denoised.y.max()
  return denoised




def roundodd(num):
  rounded = round(num)
  if rounded%2 != 0:
    return rounded
  else:
    if rounded > num:
      return rounded - 1
    else:
      return rounded + 1



def GPSmooth(data, low=0.1, high=10, debug=False):
  """
  This will smooth the data using Gaussian processes. It will find the best
  smoothing parameter via cross-validation to be between the low and high.

  The low and high keywords are reasonable bounds for  A and B stars with 
  vsini > 100 km/s.
  """

  smoothed = data.copy()

  # First, find outliers by doing a guess smooth
  smoothed = SmoothData(data, normalize=False)
  temp = smoothed.copy()
  temp.y = data.y/smoothed.y
  temp.cont = FittingUtilities.Continuum(temp.x, temp.y, lowreject=2, highreject=2, fitorder=3)
  outliers = HelperFunctions.FindOutliers(temp, numsiglow=3, expand=5)
  if len(outliers) > 0:
    data.y[outliers] = smoothed.y[outliers]
    
  gp = GaussianProcess(corr='squared_exponential',
                       theta0 = numpy.sqrt(low*high),
                       thetaL = low,
                       thetaU = high,
                       normalize = False,
                       nugget = (data.err / data.y)**2,
                       random_start=1)
  try:
    gp.fit(data.x[:,None], data.y)
  except ValueError:
    #On some orders with large telluric residuals, this will fail.
    # Just fall back to the old smoothing method in that case.
    return SmoothData(data), 91
  if debug:
    print "\tSmoothing parameter theta = ", gp.theta_
  smoothed.y, smoothed.err = gp.predict(data.x[:,None], eval_MSE=True)
  return smoothed, gp.theta_[0][0]


if __name__ == "__main__":
  fileList = []
  plot = False
  vsini_file = "%s/School/Research/Useful_Datafiles/Vsini.csv" %(os.environ["HOME"])
  for arg in sys.argv[1:]:
    if "-p" in arg:
      plot = True
    elif "-vsini" in arg:
      vsini_file = arg.split("=")[-1]
    else:
      fileList.append(arg)

  #Read in the vsini table
  vsini_data = ascii.read(vsini_file)[10:]

  if len(fileList) == 0:
    fileList = [f for f in os.listdir("./") if f.endswith("telluric_corrected.fits")]
  for fname in fileList:
    orders = HelperFunctions.ReadFits(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
    
    
    #Find the vsini of this star
    header = fits.getheader(fname)
    starname = header["object"]
    found = False
    for data in vsini_data:
      if data[0] == starname:
        vsini = float(data[1])
        found = True
    if not found:
      outfile = open("Warnings.log", "a")
      outfile.write("Cannot find %s in the vsini data: %s\n" %(starname, vsini_file))
      outfile.close()
      warnings.warn("Cannot find %s in the vsini data: %s" %(starname, vsini_file))
    print starname, vsini
    
    #Begin looping over the orders
    column_list = []
    header_list = []
    for i, order in enumerate(orders):
      print "Smoothing order %i/%i" %(i+1, len(orders))
      #Fix errors
      order.err[order.err > 1e8] = numpy.sqrt(order.y[order.err > 1e8])

      #Linearize
      xgrid = numpy.linspace(order.x[0], order.x[-1], order.x.size)
      order = FittingUtilities.RebinData(order, xgrid)
      
      dx = order.x[1] - order.x[0]
      smooth_factor = 0.8
      theta = roundodd(vsini/3e5 * order.x.mean()/dx * smooth_factor)
      denoised = SmoothData(order, 
                            windowsize=theta, 
                            smoothorder=3, 
                            lowreject=3, 
                            highreject=3,
                            expand=10, 
                            numiters=10)
      #denoised, theta = GPSmooth(order.copy())
      #denoised, theta = CrossValidation(order.copy(), 5, 2, 2, 10)
      #denoised, theta = OptimalSmooth(order.copy())
      #denoised.y *= order.cont/order.cont.mean()
      print "Window size = %.4f nm" %theta


      column = {"wavelength": denoised.x,
                "flux": order.y / denoised.y,
                "continuum": denoised.cont,
                "error": denoised.err}
      header_list.append((("Smoother", theta, "Smoothing Parameter"),))
      column_list.append(column)
      if plot:
        plt.figure(1)
        plt.plot(order.x, order.y/order.y.mean())
        plt.plot(denoised.x, denoised.y/denoised.y.mean())
        plt.title(starname)
        plt.figure(2)
        plt.plot(order.x, order.y/denoised.y)
        plt.title(starname)
        #plt.plot(order.x, (order.y-denoised.y)/numpy.median(order.y))
        #plt.show()
    if plot:
      plt.show()
    outfilename = "%s_smoothed.fits" %(fname.split(".fits")[0])
    print "Outputting to %s" %outfilename
    HelperFunctions.OutputFitsFileExtensions(column_list, fname, outfilename, mode='new', headers_info=header_list)
