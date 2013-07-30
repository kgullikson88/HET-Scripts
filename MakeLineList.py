import numpy
import matplotlib.pyplot as plt
import FitsUtils
import pyfits
import Units
import os
import FittingUtilities
from astropy import units, constants
import DataStructures
import sys


bclength = 1000  #Boxcar smoothing length


def main1():
  model_dir = "%s/School/Research/Models/Sorted/Stellar/Vband/" %(os.environ["HOME"])
  modelfile = "%slte86-4.00+0.0-alpha0.KURUCZ_MODELS.dat.sorted" %model_dir
  threshold = 0.90
  print "Reading stellar model"
  x,y = numpy.loadtxt(modelfile, usecols=(0,1), unpack=True)
  x *= units.angstrom.to(units.nm)
  y = 10**y
  model = DataStructures.xypoint(x=x, y=y)
  model.cont = FittingUtilities.Continuum(model.x, model.y, fitorder=21, lowreject=1.5, highreject=20)
  plt.plot(model.x, model.y)
  plt.plot(model.x, model.cont)
  plt.show()


  print "Finding lines"
  linepoints = numpy.where(model.y/model.cont < threshold)[0]
  points = []
  lines = []
  strengths = []
  for line in linepoints:
    #print len(points)
    if len(points) == 0 or int(line) - 1 == points[-1]:
      points.append(int(line))
    else:
      index = int(numpy.median(points) + 0.5)
      if len(points) > 1:
        minindex = model.y[points[0]:points[-1]].argmin() + points[0]
      else:
        minindex = points[0]
      lines.append(model.x[minindex])
      yval = model.y[minindex] / model.cont[minindex]
      strengths.append(yval)
      points = [int(line),]

  """
  #Make sure there are no points too close to each other
  tol = 0.05
  lines = sorted(lines)
  for i in range(len(lines) - 2, 0, -1):
    if numpy.abs(lines[i] - lines[i-1]) < tol:
      del lines[i]
      del lines[i-1]
    elif numpy.abs(lines[i] - lines[i+1]) < tol:
      del lines[i+1]
      del lines[i]
    else:
      index = numpy.searchsorted(x,lines[i]) - 1
      yval = trans[index]
      plt.plot((lines[i], lines[i]), (yval-0.05, yval-0.1), 'r-')
  """
  plt.plot(model.x, model.y/model.cont, 'k-')
  for line in lines:
    idx = numpy.searchsorted(model.x, line)
    x = model.x[idx]
    y = model.y[idx] / model.cont[idx]
    plt.plot([x, x], [y-0.05, y-0.1], 'r-')
  plt.show()
  numpy.savetxt("Linelist3.dat", numpy.transpose((lines, strengths)), fmt="%.8f")



def main2(modelfile, threshold=0.8, reverse=False, xfactor=1.0, yfactor=1.0, unlog=False, fitcontinuum=False, fitorder=3):
  print "Reading model"
  x,y = numpy.loadtxt(modelfile, usecols=(0,1), unpack=True)
  if reverse:
    x = x[::-1]*xfactor
    y = y[::-1]*yfactor
  if unlog:
    y = 10**y
  model = DataStructures.xypoint(x=x, y=y)
  if fitcontinuum:
    model.cont = FittingUtilities.Continuum(model.x, model.y, fitorder=fitorder, lowreject=1.5, highreject=20)
  

  print "Finding lines"
  slope = [(model.y[i]-model.y[i-1])/(model.x[i]-model.x[i-1]) for i in range(1, model.x.size)]
  slope = numpy.append(numpy.array(0.0), numpy.array(slope))
  linepoints = numpy.where(model.y/model.cont < threshold)[0]
  points = []
  lines = []
  strengths = []
  for line in linepoints:
    #print len(points)
    if len(points) == 0 or int(line) - 1 == points[-1]:
      points.append(int(line))
    else:
      print "\n", model.y[points[0]:points[-1]] / model.cont[points[0]:points[-1]]
      print slope[points[0]:points[-1]]
      if len(points) > 1:
        asign = numpy.sign(slope[points[0]:points[-1]])
        signchange = ((numpy.roll(asign, 1) - asign) < 0).astype(int)
        minima = numpy.where(signchange > 0.5)[0] + points[0]
        print signchange
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

  #Make sure no lines are really close to other ones
  for i in range(len(lines)-2, -1, -1):
    line = lines[i]
    after = lines[i+1]
    if numpy.abs(line - after) < 0.03:
      #Find the stronger line
      if strengths[i] > strengths[i+1]:
        lines.pop(i+1)
        strengths.pop(i+1)
      else:
        lines.pop(i)
        strengths.pop(i)
        
    

      
  plt.plot(model.x, model.y/model.cont, 'k-')
  for line in lines:
    idx = numpy.searchsorted(model.x, line)
    x = model.x[idx]
    y = model.y[idx] / model.cont[idx]
    plt.plot([x, x], [y-0.05, y-0.1], 'r-')
  plt.show()
  numpy.savetxt("Linelist.dat", numpy.transpose((lines, strengths)), fmt="%.8f")

  

if __name__ == "__main__":
  model_dir = "%s/School/Research/Models/Sorted/Stellar/Vband/" %(os.environ["HOME"])
  modelfile = "%slte86-4.00-0.5-alpha0.KURUCZ_MODELS.dat.sorted" %model_dir
  if len(sys.argv)>1:
    modelfile = sys.argv[1]
  main2(modelfile, threshold=0.95, reverse=True)

