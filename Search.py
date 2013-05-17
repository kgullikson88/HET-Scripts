import Correlate
import FitsUtils
import numpy
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import os
import sys
import DataStructures
import FindContinuum
import matplotlib.pyplot as plt
import Units

homedir = os.environ["HOME"]
modeldir = homedir + "/School/Research/Models/Sorted/Stellar/Vband/"

#Define regions contaminated by telluric residuals or other defects. We will not use those regions in the cross-correlation
badregions = [[588.8, 589.9],
              [627.1, 635.4]]

#Set up model list
model_list = [ modeldir + "lte30-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte32-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte34-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte35-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte36-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte37-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte38-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte39-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte30-4.00+0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte32-4.00+0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte34-4.00+0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte35-4.00+0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte36-4.00+0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte37-4.00+0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte38-4.00+0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte39-4.00+0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte30-4.00-0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte32-3.50-0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte34-3.50-0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte35-4.00-0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte36-3.50-0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte37-4.00-0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte38-4.00-0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte39-4.00-0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte40-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte42-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte44-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte46-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte48-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte50-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte51-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte52-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte53-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte54-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte55-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte56-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte57-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte58-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte59-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte60-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte61-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte62-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte63-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte64-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte65-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte66-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte67-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte68-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte69-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte69-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte70-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte70-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte72-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte74-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte74-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte76-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte78-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted"]
star_list = []
temp_list = []
gravity_list = []
metal_list = []
for fname in model_list:
  temp = int(fname.split("lte")[-1][:2])*100
  gravity = float(fname.split("lte")[-1][3:7])
  metallicity = float(fname.split("lte")[-1][7:11])
  star_list.append(str(temp))
  temp_list.append(temp)
  gravity_list.append(gravity)
  metal_list.append(metallicity)

if __name__ == "__main__":
  #Parse command line arguments:
  fileList = []
  extensions=True
  tellurics=False
  trimsize = 50
  for arg in sys.argv[1:]:
    if "-e" in arg:
      extensions=False
    if "-t" in arg:
      tellurics=True  #telluric lines modeled but not removed
    else:
      fileList.append(arg)

  for fname in fileList:
    if extensions:
      orders = FitsUtils.MakeXYpoints(fname, extensions=extensions, x="wavelength", y="flux", errors="error")
      if tellurics:
        model_orders = FitsUtils.MakeXYpoints(fname, extensions=extensions, x="wavelength", y="model")
        for i, order in enumerate(orders):
          orders[i].cont = FindContinuum.Continuum(order.x, order.y, lowreject=2, highreject=2)
          plt.plot(order.x, order.y/orders[i].cont, 'b-')
          orders[i].y /= model_orders[i].y
          plt.plot(order.x, orders[i].y/orders[i].cont-0.5, 'r-')
          plt.plot(order.x, model_orders[i].y + 0.5, 'g-')
        plt.show()
    else:
      orders = FitsUtils.MakeXYpoints(fname, errors=2)
    numremoved = 0
    for i, order in enumerate(orders):
      order.x = order.x[trimsize:-trimsize]
      order.y = order.y[trimsize:-trimsize]
      order.cont = order.cont[trimsize:-trimsize]
      order.err = order.err[trimsize:-trimsize]
      if order.x.size > 10:
        orders[i] = order.copy()
      else:
        orders[i].pop(i-numremoved)
        numremoved += 1
    

    xspacing = 1e9
    for i, order in enumerate(orders):
      orders[i].cont = FindContinuum.Continuum(order.x, order.y, lowreject=2, highreject=2)
      spacing = (order.x[-1] - order.x[0]) / float(order.size() )
      if spacing < xspacing:
        xspacing = spacing
      
    data = DataStructures.CombineXYpoints(orders, xspacing=xspacing)

    #Remove bad regions from the data
    for region in badregions:
      left = numpy.searchsorted(data.x, region[0])
      right = numpy.searchsorted(data.x, region[1])
      data.y[left:right] = data.cont[left:right]
    #plt.plot(data.x, data.y/data.cont)
    #plt.show()

    #Do the cross-correlation
    for vsini in [10, 20, 30, 40]:
      Correlate.PyCorr(data, combine=False, resolution=30000, outdir="Cross_correlations/%s" %(fname.split(".fits")[0]), models=model_list, stars=star_list, temps=temp_list, gravities=gravity_list, metallicities=metal_list, vsini=vsini*Units.cm/Units.km, pause=10)


