import numpy
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as interp
from scipy.optimize import leastsq, brute
from scipy import mat
from scipy.linalg import svd, diagsvd
import sys
import os
import FitsUtils
from collections import defaultdict
import DataStructures
from astropy import units, constants
import MakeModel
import RotBroad
import FittingUtilities
import FindContinuum

"""
  Function to create a model spectrum from a list
    of line wavelengths and strengths.
  fullmodel:   xypoint structure with the full model.
                 This code will determine the important
                 lines from that
  xgrid:       numpy array of the x points in the data
  vsini:       In km/s
  epsilon:     Linear limb darkening parameter
  resolution:  Detector resolution ( lambda/(delta lambda) )
"""
def FitData(fullmodel, data, vsini, epsilon, resolution, threshold=1.0, vsys=0.0):
  #First, find the lines in the relevant part of the spectrum
  left = numpy.searchsorted(fullmodel.x, data.x[0])
  right = numpy.searchsorted(fullmodel.x, data.x[-1])
  segment = fullmodel[left:right]
  segment.cont = FittingUtilities.Continuum(segment.x, segment.y, fitorder=15, lowreject=1.5, highreject=10, numiter=10)
  lines, strengths = FindLines(segment, threshold=threshold, numlines=20)
  """
  plt.plot(segment.x, segment.y/segment.cont)
  for i in range(len(lines)):
    x = lines[i]
    y = strengths[i]
    plt.plot([x,x], [y-0.1, y-0.05], 'r-')
  plt.show()
  """
  #Prepare parameter list for fit
  pars = []
  const_pars = [lines.size,]
  for i in range(lines.size):
    const_pars.append(lines[i])
    pars.append(strengths[i])
  pars.append(vsini)
  pars.append(epsilon)
  pars.append(vsys)
  const_pars.append(resolution)
  fitpars, success = leastsq(ErrorFunction, pars, args=(const_pars, data))
  #model = GenerateModel(lines, strengths, data.x, vsini, epsilon, resolution, vsys)
  chisq = numpy.sum(ErrorFunction(fitpars, const_pars, data))/float(data.size() - len(pars))
  model = ErrorFunction(fitpars, const_pars, data, return_model=True)
  print "X^2 = %g" %chisq
  plt.plot(model.x, model.y, 'r-')
  
  
def ErrorFunction(pars, const_pars, data, return_model=False):
  #Unpack parameters
  numlines = const_pars[0]
  lines = numpy.zeros(numlines)
  strengths = numpy.zeros(numlines)
  for i in range(numlines):
    lines[i] = const_pars[i+1]
    strengths[i] = pars[i]
  resolution = const_pars[numlines+1]
  vsini = pars[numlines]
  epsilon = pars[numlines+1]
  vsys = pars[numlines+2]
  model = GenerateModel(lines, strengths, data.x, vsini, epsilon, resolution, vsys)
  
  if not return_model:
    return (data.y - model.y*data.cont)**2 / data.err**2
  else:
    return model
  
  
  
def GenerateModel(lines, strengths, xgrid, vsini, epsilon, resolution, vsys):
  #Make spectrum of delta functions with same xgrid as given
  strengths = 1.0 - strengths
  model = DataStructures.xypoint(x=xgrid, y=numpy.zeros(xgrid.size))
  delta_v = (xgrid[1] - xgrid[0])/xgrid[xgrid.size/2] * 3e5
  factor = 10./delta_v
  for i, line in enumerate(lines):
    idx = numpy.searchsorted(model.x, line*(1+vsys/constants.c.cgs.value))
    model.y[idx] = -strengths[i]*factor
  model.y += 1.0
  model = RotBroad.Broaden2(model.copy(), vsini*units.km.to(units.cm), linear=True, epsilon=epsilon)
  model = MakeModel.ReduceResolution(model, resolution)
  return model
  
"""  
  Function to find the strongest spectral lines
"""
def FindLines(model, threshold=1.0, numlines=10):
  slope = [(model.y[i]-model.y[i-1])/(model.x[i]-model.x[i-1]) for i in range(1, model.x.size)]
  slope = numpy.append(numpy.array(0.0), numpy.array(slope))
  linepoints = numpy.where(model.y/model.cont < threshold)[0]
  points = []
  lines = []
  strengths = []
  for line in linepoints:
    if len(points) == 0 or int(line) - 1 == points[-1]:
      points.append(int(line))
    else:
      if len(points) > 1:
        asign = numpy.sign(slope[points[0]:points[-1]])
        signchange = ((numpy.roll(asign, 1) - asign) < 0).astype(int)
        minima = numpy.where(signchange > 0.5)[0] + points[0]
        for minimum in minima:
          if model.y[minimum - 1] < model.y[minimum]:
            minimum -= 1
          lines.append(model.x[minimum])
          yval = model.y[minimum] / model.cont[minimum]
          strengths.append(yval)
        
      else:
        lines.append(model.x[points[0]])
        yval = model.y[points[0]] / model.cont[points[0]]
        strengths.append(yval)
          
      points = [int(line),]
  
  #Sort by line strength
  lines = numpy.array(lines)
  strengths = numpy.array(strengths)
  sortindices = [i[0] for i in sorted(enumerate(strengths), key=lambda x:x[1])]
  lines = lines[sortindices]
  strengths = strengths[sortindices]
  return lines[:numlines], strengths[:numlines]
  
  
  
if __name__ == "__main__":
  model_dir = "%s/School/Research/Models/Sorted/Stellar/Vband/" %(os.environ["HOME"])
  modelfile = "%slte86-4.00+0.0-alpha0.KURUCZ_MODELS.dat.sorted" %model_dir
  threshold = 0.90
  print "Reading file %s" %modelfile
  x, y = numpy.loadtxt(modelfile, usecols=(0,1), unpack=True)
  model = DataStructures.xypoint(x=x*units.angstrom.to(units.nm)/1.00026, 
                                 y=10**y)
  #model.cont = FittingUtilities.Continuum(model.x, model.y, fitorder=21, lowreject=1.5, highreject=20)
  
  for fname in sys.argv[1:]:
    orders = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", errors="error")
    for order in orders:
      #Linearize
      datafcn = interp(order.x, order.y)
      x = numpy.linspace(order.x[0], order.x[-1], order.size())
      order = DataStructures.xypoint(x=x, y=datafcn(x))
      order.cont = FittingUtilities.Continuum(order.x, order.y)
      plt.plot(order.x, order.y/order.cont, 'k-')
      FitData(model, order, 140, 0.6, 60000)
    plt.show()
      
      
      
      
      
      
      
      
      
      
      
      
      
