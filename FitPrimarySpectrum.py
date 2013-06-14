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



def BroadeningErrorFunction(pars, data, unbroadened):
  vsini, beta = pars[0]*units.km.to(units.cm), pars[1]
  #print vsini, beta
  #data, unbroadened = args[0], args[1]
  model = RotBroad.Broaden(unbroadened, vsini, beta=beta)
  model = MakeModel.RebinData(model, data.x)
  return numpy.sum( FittingUtilities.GeneralLSErrorFunction(data, model) )



#Fits the broadening profile using singular value decomposition
#oversampling is the oversampling factor to use before doing the SVD
#m is the size of the broadening function, in oversampled units
#dimension is the number of eigenvalues to keep in the broadening function. (Keeping too many starts fitting noise)
def Broaden(data, model, oversampling = 5, m = 201, dimension = 15):
  n = data.x.size*oversampling
  
  #n must be even, and m must be odd!
  if n%2 != 0:
    n += 1
  if m%2 == 0:
    m += 1
  
  #resample data
  Spectrum = interp(data.x, data.y/data.cont)
  Model = interp(model.x, model.y)
  xnew = numpy.linspace(data.x[0], data.x[-1], n)
  ynew = Spectrum(xnew)
  model_new = Model(xnew)

  #Make 'design matrix'
  design = numpy.zeros((n-m,m))
  for j in range(m):
    for i in range(m/2,n-m/2-1):
      design[i-m/2,j] = model_new[i-j+m/2]
  design = mat(design)
    
  #Do Singular Value Decomposition
  try:
    U,W,V_t = svd(design, full_matrices=False)
  except numpy.linalg.linalg.LinAlgError:
    outfilename = "SVD_Error.log"
    outfile = open(outfilename, "a")
    numpy.savetxt(outfile, numpy.transpose((data.x, data.y, data.cont)))
    outfile.write("\n\n\n\n\n")
    numpy.savetxt(outfile, numpy.transpose((model.x, model.y, model.cont)))
    outfile.write("\n\n\n\n\n")
    outfile.close()
    sys.exit("SVD did not converge! Outputting data to %s" %outfilename)
      
  #Invert matrices:
  #   U, V are orthonormal, so inversion is just their transposes
  #   W is a diagonal matrix, so its inverse is 1/W
  W1 = 1.0/W
  U_t = numpy.transpose(U)
  V = numpy.transpose(V_t)
  
  #Remove the smaller values of W
  W1[dimension:] = 0
  W2 = diagsvd(W1,m,m)
    
  #Solve for the broadening function
  spec = numpy.transpose(mat(ynew[m/2:n-m/2-1]))
  temp = numpy.dot(U_t, spec)
  temp = numpy.dot(W2,temp)
  Broadening = numpy.dot(V,temp)

  #Make Broadening function a 1d array
  spacing = xnew[2] - xnew[1]
  xnew = numpy.arange(model.x[0], model.x[-1], spacing)
  model_new = Model(xnew)
  Broadening = numpy.array(Broadening)[...,0]

  
  #If we get here, the broadening function looks okay.
  #Convolve the model with the broadening function
  model = DataStructures.xypoint(x=xnew)
  Broadened = interp(xnew, numpy.convolve(model_new,Broadening, mode="same") )
  model.y = Broadened(model.x)
  
  return MakeModel.RebinData(model, data.x)


def main4():
  linelist = "../Scripts/LineList.dat"
  lines, strengths = numpy.loadtxt(linelist, unpack=True)
  strengths = 1.0 - strengths

  fname = "HIP_70384.fits"
  orders = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", errors="error")

  for i, order in enumerate(orders):
    DATA = interp(order.x, order.y)
    order.x, xspacing = numpy.linspace(order.x[0], order.x[-1], order.x.size, retstep=True)
    order.y = DATA(order.x)
    order.cont = FittingUtilities.Continuum(order.x, order.y)
    
    left = numpy.searchsorted(lines, order.x[0])
    right = numpy.searchsorted(lines, order.x[-1])
    print right - left + 1
    unbroadened = order.copy()
    unbroadened.y = numpy.zeros(unbroadened.x.size)
    unbroadened.cont = numpy.ones(unbroadened.x.size)
    deltav = xspacing / numpy.median(order.x) * 3e5
    print deltav
    factor = 10./deltav
    for j, line in enumerate(lines[left:right]):
      x = line
      y = strengths[j+left]
      idx = numpy.searchsorted(unbroadened.x, line)
      unbroadened.y[idx] = -y*factor
    unbroadened.y += 1.0
    #model = Broaden(order, unbroadened, m=401, dimension=20)
    model2 = RotBroad.Broaden2(unbroadened.copy(), 100*units.km.to(units.cm), linear=True)
    model3 = RotBroad.Broaden2(unbroadened.copy(), 140*units.km.to(units.cm), linear=True)
    model2 = MakeModel.ReduceResolution(model2, 60000)
    model3 = MakeModel.ReduceResolution(model3, 60000)

    unbroadened.y = (unbroadened.y - 1.0)/factor + 1.0
    #plt.plot(model.x, model.y)
    plt.plot(order.x, order.y/order.cont)
    plt.plot(unbroadened.x, unbroadened.y)
    plt.plot(model2.x, model2.y)
    plt.plot(model3.x, model3.y)
    plt.show()



def main3():
  model_dir = "%s/School/Research/Models/Sorted/Stellar/Vband/" %(os.environ["HOME"])
  modelfile = "%slte86-4.00-0.5-alpha0.KURUCZ_MODELS.dat.sorted" %model_dir
  vsini = 150.0
  beta = 1.0
  x, y = numpy.loadtxt(modelfile, usecols=(0,1), unpack=True)
  MODEL = interp(x*units.angstrom.to(units.nm)/1.00026, 10**y)
  xlin = numpy.linspace(x[0], x[-1], x.size)
  #model = DataStructures.xypoint(x=xlin, y=MODEL(xlin))
  #model = DataStructures.xypoint(x=x*units.angstrom.to(units.nm)/1.00026, y=10**y)
  #model.cont = FittingUtilities.Continuum(model.x, model.y, fitorder=15, lowreject=1.5, highreject=10)
  #model2 = RotBroad.Broaden(model, vsini*units.km.to(units.cm))

  fname = "HIP_70384.fits"
  orders = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", errors="error")

  for i, order in enumerate(orders):
    DATA = interp(order.x, order.y)
    xlin = numpy.linspace(order.x[0], order.x[-1], order.x.size)
    order = DataStructures.xypoint(x=xlin, y=DATA(xlin))
    order.cont = FittingUtilities.Continuum(order.x, order.y)
    extended = numpy.append( numpy.append(order.y[::-1]/order.cont[::-1], order.y/order.cont), order.y[::-1]/order.cont[::-1])
    plt.plot(order.x, order.y)
    plt.plot(order.x, order.cont)
    plt.show()
    plt.plot(extended)
    plt.show()

    unbroadened = DataStructures.xypoint(x=xlin, y=MODEL(xlin))
    unbroadened.cont = FittingUtilities.Continuum(unbroadened.x, unbroadened.y, fitorder=4, lowreject=1.5, highreject=10)
    plt.plot(unbroadened.x, unbroadened.y)
    plt.plot(unbroadened.x, unbroadened.cont)
    plt.show()

    #Fit broadening
    #left = numpy.searchsorted(model.x, 2*order.x[0] - order.x[-1])
    #right = numpy.searchsorted(model.x, 2*order.x[-1] - order.x[0])
    #unbroadened = model[left:right]
    size = unbroadened.size()
    
    #model2 = Broaden(order, unbroadened, m=401, dimension=20)
    #model2 = MakeModel.RebinData(model2, order.x)
    ycorr = numpy.correlate(extended-1.0, unbroadened.y/unbroadened.cont-1.0, mode='same')[size:-size]
    #ycorr -= ycorr.min()
    plt.plot(ycorr)
    plt.show()

    model2 = order.copy()
    model2.y = numpy.correlate(extended, ycorr/ycorr.sum(), mode='same')[size:-size]
    
    model2.cont = FittingUtilities.Continuum(model2.x, model2.y, lowreject=1.5, highreject=10)
    model2.y = (model2.y/model2.cont-1)*10.0 + model2.cont

    plt.figure(1)
    plt.plot(order.x, order.y/order.cont, 'k-')
    plt.plot(unbroadened.x, unbroadened.y/unbroadened.cont, 'g-')
    plt.plot(model2.x, model2.y/model2.cont, 'r-')
    plt.figure(3)
    plt.plot(order.x, order.y / (model2.y*order.cont))
    order.y -= model2.y/model2.cont*order.cont
    plt.figure(2)
    plt.plot(order.x, order.y/order.cont)
    plt.show()
    


def main2():
  model_dir = "%s/School/Research/Models/Sorted/Stellar/Vband/" %(os.environ["HOME"])
  modelfile = "%slte86-4.00+0.5-alpha0.KURUCZ_MODELS.dat.sorted" %model_dir
  vsini = 150.0
  beta = 1.0
  x, y = numpy.loadtxt(modelfile, usecols=(0,1), unpack=True)
  model = DataStructures.xypoint(x=x*units.angstrom.to(units.nm)/1.00026, y=10**y)
  model.cont = FittingUtilities.Continuum(model.x, model.y, fitorder=15, lowreject=1.5, highreject=10)
  #model2 = RotBroad.Broaden(model, vsini*units.km.to(units.cm))

  fname = "HIP_70384.fits"
  orders = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", errors="error")

  for i, order in enumerate(orders):
    order.cont = FittingUtilities.Continuum(order.x, order.y)

    #Fit broadening
    left = numpy.searchsorted(model.x, 2*order.x[0] - order.x[-1])
    right = numpy.searchsorted(model.x, 2*order.x[-1] - order.x[0])
    unbroadened = model[left:right]
    pars = [vsini, beta]
    grid = [[120, 200, 10],
            [0.2, 1.8, 0.1]]
    #fitpars, fitval = leastsq(BroadeningErrorFunction, pars, args=(order, unbroadened))
    #fitvsini, fitbeta = fitpars[0], fitpars[1]
    fitvsini, fitbeta = brute(BroadeningErrorFunction, grid, args=(order, unbroadened), finish=None )
    print "Fitted values for order %i:\n\tvsini = %g\n\tbeta = %g" %(i+1, fitvsini, fitbeta)
    model2 = RotBroad.Broaden(unbroadened, fitvsini*units.km.to(units.cm), beta=fitbeta)
    model2 = MakeModel.RebinData(model2, order.x)
    model2.cont = FittingUtilities.Continuum(model2.x, model2.y, lowreject=1.5, highreject=10)

    plt.figure(1)
    plt.plot(order.x, order.y/order.cont, 'k-')
    plt.plot(model2.x, model2.y/model2.cont, 'r-')
    plt.figure(3)
    plt.plot(order.x, order.y / (model2.y*order.cont))
    order.y -= model2.y/model2.cont*order.cont
    plt.figure(2)
    plt.plot(order.x, order.y/order.cont)
  plt.show()
    
    

def main1():
  #Parse command line arguments
  fileList = []
  vsini = 100.0
  resolution = 60000
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
    elif "resolution" in arg:
      resolution = float(arg.split("=")[-1])
    else:
      fileList.append(arg)

  if not model_dir.endswith("/"):
    model_dir = model_dir + "/"

  reduce_resolution = True
  v_res = 3e5/resolution
  if v_res < 0.1*vsini:
    reduce_resolution = False
    print "Will not reduce detector resolution: %g\t%g" %(v_res, vsini)

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
      elif "KURUCZ" in model:
        #Prefer the KURUCZ models that I make over everything else
        file_dict[T][logg][Z] = model
      elif "KURUCZ" in file_dict[T][logg][Z]:
        continue
      elif "PHOENIX-ACES" in model and "PHOENIX2004" in file_dict[T][logg][Z]:
        #Prefer PHOENIX_ACES over PHOENIX2004 (ACES was made in 2009)
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
        model.cont = FittingUtilities.Continuum(model.x, model.y, fitorder=15, lowreject=1.5, highreject=10)
        model = RotBroad.Broaden(model, vsini*units.km.to(units.cm),  )
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
            if reduce_resolution:
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
            if reduce_resolution:
              model.y /= model.cont
              model = MakeModel.ReduceResolution(model, 60000)
              model.y *= model.cont
            model.y *= model.cont
            chisq += numpy.sum( (order.y - model.y/model.cont*order.cont)**2 / order.err**2 )
            norm += order.size()
            #plt.plot(order.x, order.y/order.cont, 'k-')
            #plt.plot(model.x, model.y/model.cont, 'r-')
          #plt.show()
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
      if reduce_resolution:
        model = MakeModel.ReduceResolution(model, 60000)
      plt.figure(1)
      plt.plot(order.x, order.y/order.cont, 'k-')
      plt.plot(model.x, model.y/model.cont, 'r-')
      plt.figure(3)
      plt.plot(order.x, order.y / (model.y*order.cont))
      order.y -= model.y/model.cont*order.cont
      plt.figure(2)
      plt.plot(order.x, order.y/order.cont)
    plt.show()
          
    
if __name__ == "__main__":
  main4()
