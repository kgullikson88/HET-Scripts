import numpy
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as interp
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


if __name__ == "__main__":
  #Parse command line arguments
  fileList = []
  vsini = 100.0
  Tmin, Tmax, Tstep = 8000, 8800, 200
  Zmin, Zmax, Zstep = -0.5, 0.5, 0.5
  loggmin, loggmax, loggstep = 4.0, 4.0, 1.0
  model_dir = "%s/School/Research/Models/Sorted/Stellar/Vband/" %(os.environ["HOME"])
  for arg in sys.argv[1:]:
    if "vsini" in arg:
      vsini = float(arg.split("=")[-1])
    elif "Tmin" in arg:
      Tmin = float(arg.split("=")[-1])
    elif "Tmax" in arg:
      Tmax = float(arg.split("=")[-1])
    elif "Tstep" in arg:
      Tstep = float(arg.split("=")[-1])
    elif "Zmin" in arg:
      Zmin = float(arg.split("=")[-1])
    elif "Zmax" in arg:
      Zmax = float(arg.split("=")[-1])
    elif "Zstep" in arg:
      Zstep = float(arg.split("=")[-1])
    elif "loggmin" in arg:
      loggmin = float(arg.split("=")[-1])
    elif "loggmax" in arg:
      loggmax = float(arg.split("=")[-1])
    elif "loggstep" in arg:
      loggstep = float(arg.split("=")[-1])
    elif "modeldir" in arg:
      model_dir = arg.split("=")[-1]
    else:
      fileList.append(arg)

  if not model_dir.endswith("/"):
    model_dir = model_dir + "/"

  #Read in all of the necessary model files
  allmodels = os.listdir(model_dir)
  file_dict = defaultdict( lambda : defaultdict( lambda: defaultdict( str) ) )
  for model in allmodels:
    if "BSTAR_MODELS" in model:
      T = float(model[3:6])*100
      logg = float(model[7:11])
      Z = float(model[11:15])
    elif "PHOENIX-ACES" in model:
      T = float(model[3:5])*100
      logg = float(model[6:10])
      Z = float(model[10:14])
    elif "PHOENIX2004" in model:
      T = float(model[3:5])*100
      logg = float(model[6:9])
      Z = float(model[9:13])
    elif "KURUCZ" in model:
      T = float(model[3:5])*100
      logg = float(model[6:10])
      Z = float(model[10:14])
    else:
      continue

    #Only save the filenames if in the correct T, logg, and Z range
    if (T >= Tmin and T <= Tmax and
        logg >= loggmin and logg <= loggmax and
        Z >= Zmin and Z <= Zmax):

      if file_dict[T][logg][Z] == "":
        file_dict[T][logg][Z] = model
      elif "PHOENIX-ACES" in model and "PHOENIX2004" in file_dict[T][logg][Z]:
        #Prefer PHOENIX_ACES
        file_dict[T][logg][Z] = model
      else:
        print "Two models with the same T, logg, and Z!"
        print "(1):", file_dict[T][logg][Z]
        print "(2):", model
        inp = raw_input("Which one do you want to use? ")
        if inp == "2":
          file_dict[T][logg][Z] = model

  #Now, actually read in the models we saved and store as xypoints
  #model_dict = defaultdict( lambda : defaultdict( lambda: defaultdict( DataStructures.xypoint ) ) )
  model_dict = defaultdict( lambda : defaultdict( lambda: defaultdict( interp ) ) )
  for T in file_dict:
    for logg in file_dict[T]:
      for Z in file_dict[T][logg]:
        print "Reading file %s" %file_dict[T][logg][Z]
        x, y = numpy.loadtxt(model_dir + file_dict[T][logg][Z], usecols=(0,1), unpack=True)
        #model_dict[T][logg][Z] = DataStructures.xypoint(x=x*units.angstrom.to(units.nm)/1.00026,
        #                                                y=10**y)
        model = DataStructures.xypoint(x=x*units.angstrom.to(units.nm)/1.00026, y=10**y)
        model = RotBroad.Broaden(model, vsini*units.km.to(units.cm) )
        model_dict[T][logg][Z] = interp(model.x, model.y)

          
  #Done reading in the models. Now, loop over the actual data and try to fit
  for fname in fileList:
    orders = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", errors="error")

    #Loop over the models
    best_chisq = 9e99
    best_T = 0
    best_logg = 0
    best_Z = 0
    best_rv = 0
    for T in model_dict:
      for logg in model_dict[T]:
        for Z in model_dict[T][logg]:
          #First, find the rv from all the orders
          rv = []
          for order in orders:
            order.cont = FittingUtilities.Continuum(order.x, order.y)
            model = DataStructures.xypoint(x=order.x, y=model_dict[T][logg][Z](order.x) )
            model.cont = FittingUtilities.Continuum(model.x, model.y, lowreject=1.5, highreject=10)
            model = MakeModel.ReduceResolution(model, 60000)
            #model = RotBroad.Broaden(model, vsini*units.km.to(units.cm))
            offset= FittingUtilities.CCImprove(order, model, be_safe=False)
            rv.append(-offset/order.x.mean() * constants.c.cgs.value)

          #Apply the median rv to all, and determine X^2
          rv = numpy.median(rv)
          chisq = 0.0
          norm = 0.0
          for order in orders:
            order.cont = FittingUtilities.Continuum(order.x, order.y)
            model = DataStructures.xypoint(x=order.x,
                                           y=model_dict[T][logg][Z](order.x*(1+rv/constants.c.cgs.value)) )
            model.cont = FittingUtilities.Continuum(model.x, model.y, lowreject=1.5, highreject=10)
            model = MakeModel.ReduceResolution(model, 60000)
            chisq += numpy.sum( (order.y - model.y/model.cont*order.cont)**2 / order.err**2 )
            norm += order.size()
          chisq /= float(norm)
          print T, logg, Z, rv, chisq
          if chisq < best_chisq:
            best_chisq = chisq
            best_T = T
            best_logg = logg
            best_Z = Z
            best_rv = rv

    print "Best fit values:"
    print "T: %g\nlog(g): %g\n[Fe/H]: %g " %(best_T, best_logg, best_Z)

    #Subtract best model
    model_fcn = model_dict[best_T][best_logg][best_Z]
    for order in orders:
      order.cont = FittingUtilities.Continuum(order.x, order.y)
      model = DataStructures.xypoint(x=order.x,
                                     y=model_fcn(order.x*(1+best_rv/constants.c.cgs.value)) )
      model.cont = FittingUtilities.Continuum(model.x, model.y, lowreject=1.5, highreject=10)
      model = MakeModel.ReduceResolution(model, 60000)
      plt.figure(1)
      plt.plot(order.x, order.y/order.cont, 'k-')
      plt.plot(model.x, model.y/model.cont, 'r-')
      order.y -= model.y/model.cont*order.cont
      plt.figure(2)
      plt.plot(order.x, order.y/order.cont)
    plt.show()
          
    
