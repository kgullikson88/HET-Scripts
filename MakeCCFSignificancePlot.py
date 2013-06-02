import numpy
import matplotlib.pyplot as plt
import sys
import os
from collections import defaultdict
import MakeModel
import DataStructures
import FindContinuum
import FitsUtils
#import Units
from astropy import units, constants
import RotBroad


Corr_dir = "Cross_Correlations/"
model_dir = "%s/School/Research/Models/Sorted/Stellar/Vband/" %(os.environ["HOME"])
modelfiles = os.listdir(model_dir)

if __name__ == "__main__":
  basename = sys.argv[1]
  velocity = float(sys.argv[2])  #Turn this into an optional argument later...
  scale = 0.1

  if basename.endswith(".fits"):
    basename = basename[:-5]
  
  files = os.listdir(Corr_dir)
  fileList = defaultdict(list)
  for fname in files:
    if fname.startswith(basename) and fname.endswith("-0.0"):
      vsini = float(fname.split("kps")[0].split(".")[-1])
      fileList[vsini].append(fname)
  for val in sorted(fileList.keys()):
    files = fileList[val]
    Temperatures = []
    Significances = []
    logg = []
    metal = []
    for fname in files:
      vels, corr = numpy.loadtxt(Corr_dir + fname, unpack=True)
      fit = numpy.poly1d( numpy.polyfit(vels, corr, 2) )
      std = numpy.std(corr - fit(vels))
      index = numpy.searchsorted(vels, velocity)
      index = min(index, vels.size - 1)
      print fname
      Temperatures.append(float(fname.split("_")[-1].split("K")[0]))
      logg.append(float(fname.split("K")[-1][:4]))
      metal.append(float(fname[-4:]))
      Significances.append((corr[index] - fit(vels[index]))/std)
    Temperatures = numpy.array(Temperatures)
    Significances = numpy.array(Significances)
    logg = numpy.array(logg)
    metal = numpy.array(metal)
    sorter = sorted(range(len(Temperatures)),key=lambda x:Temperatures[x])
    T = Temperatures[sorter]
    S = Significances[sorter]
    G = logg[sorter]
    Z = metal[sorter]
    plt.plot(T, S, '-o', color='black')
    plt.xlabel("Secondary Temperature (K)")
    plt.ylabel("CCF Value at v = %g" %velocity)
    plt.title("vsini = %g" %val)

    #Find the maximum value
    maxindex = numpy.argmax(S)
    Tmax = Temperatures[maxindex]
    loggmax = G[maxindex]
    metalmax = Z[maxindex]
    #Tmax = 4000
    #loggmax = 4.0
    #metalmax=0.5
    print "T = %g\nlogg = %g\nZ = %g" %(Tmax, loggmax, metalmax)

    #Find model file with best fit
    modelstart = "lte%.2i-%.2f%+.1f" %(Tmax/100.0, loggmax, metalmax)
    modelfile = [f for f in modelfiles if f.startswith(modelstart)][0]
    #modelfile = "lte40-4.0+0.5.Cond.PHOENIX2004.tab.7.sorted"
    x,y = numpy.loadtxt("%s%s" %(model_dir, modelfile), usecols=(0,1), unpack=True)
    x *= units.angstrom.to(units.nm)
    #x /= 1.00026
    x *= 1 - velocity*units.km.to(units.cm) / constants.c.cgs.value
    y = 10**y
    model = DataStructures.xypoint(x=x, y=y)
    print "vsini = ", val

    #Find original filename from the current directory
    filename = "%s.fits" %(basename)
    if filename in os.listdir("./"):
      orders = FitsUtils.MakeXYpoints(filename, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
      for order in orders:
        #left = numpy.searchsorted(model.x, order.x[0])
        #right = numpy.searchsorted(model.x, order.x[-1])
        #segment = DataStructures.xypoint(x=model.x[left:right], y=model.y[left:right])
        segment = MakeModel.RebinData(model, order.x)
        segment.cont = FindContinuum.Continuum(segment.x, segment.y)
        segment = RotBroad.Broaden(segment, val*units.km.to(units.cm))
        segment = MakeModel.ReduceResolution(segment, 50000)
        order.cont = FindContinuum.Continuum(order.x, order.y)

        segment.y /= segment.cont
        segment.y = (segment.y - 1.0)*scale + 1.0
        
        plt.figure(2)
        plt.plot(order.x, order.y/order.cont, 'k-')
        plt.plot(segment.x, segment.y-0.05, 'g-')
      plt.show()
      
    else:
      print "Original file not found. Searched for %s in current directory!" %filename
        
    plt.show()
    
