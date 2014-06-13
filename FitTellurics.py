import numpy
import sys
import os
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from astropy import units, constants
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import TelluricFitter
import DataStructures
import Units
from astropy import units, constants
import HelperFunctions
import FittingUtilities
import GetAtmosphere


homedir = os.environ["HOME"]
weather_file = homedir + "/School/Research/Useful_Datafiles/Weather.dat"

badregions = [[588.98, 589.037],   #Na D line 1
              [589.567, 589.632],  #Na D line 2
              [627.4, 629.0],  #O2 band
              [686.4, 690.7]]  #O2 band
              
              
namedict = {"pressure": ["PRESFIT", "PRESVAL", "Pressure"],
                  "temperature": ["TEMPFIT", "TEMPVAL", "Temperature"],
                  "angle": ["ZD_FIT", "ZD_VAL", "Zenith Distance"],
                  "resolution": ["RESFIT", "RESVAL", "Detector Resolution"],
                  "h2o": ["H2OFIT", "H2OVAL", "H2O abundance"],
                  "co2": ["CO2FIT", "CO2VAL", "CO2 abundance"],
                  "o3": ["O3FIT", "O3VAL", "O3 abundance"],
                  "n2o": ["N2OFIT", "N2OVAL", "N2O abundance"],
                  "co": ["COFIT", "COVAL", "CO abundance"],
                  "ch4": ["CH4FIT", "CH4VAL", "CH4 abundance"],
                  "o2": ["O2FIT", "O2VAL", "O2 abundance"],
                  "no": ["NOFIT", "NOVAL", "NO abundance"],
                  "so2": ["SO2FIT", "SO2VAL", "SO2 abundance"],
                  "no2": ["NO2FIT", "NO2VAL", "NO2 abundance"],
                  "nh3": ["NH3FIT", "NH3VAL", "NH3 abundance"],
                  "hno3": ["HNO3FIT", "HNO3VAL", "HNO3 abundance"]}  







def FindOrderNums(orders, wavelengths):
  """
    Given a list of xypoint orders and
    another list of wavelengths, this
    finds the order numbers with the
    requested wavelengths
  """
  nums = []
  for wave in wavelengths:
    for i, order in enumerate(orders):
      if order.x[0] < wave and order.x[-1] > wave:
        nums.append(i)
        break
  return nums




if __name__ == "__main__":
  #Initialize fitter
  fitter = TelluricFitter.TelluricFitter()
  fitter.SetObservatory("McDonald")
 
  fileList = []
  start = 0
  end = 999
  makenew = True
  edit_atmosphere=False
  humidity_low = 1.0
  humidity_high = 99.0
  for arg in sys.argv[1:]:
    if "-atmos" in arg:
      edit_atmosphere = True
    elif "-hlow" in arg:
      humidity_low = float(arg.split("=")[1])
    elif "-hhigh" in arg:
      humidity_high = float(arg.split("=")[1])
    else:
      fileList.append(arg)


  #START LOOPING OVER INPUT FILES
  for fname in fileList:
    #Make sure this file is an object file
    header = pyfits.getheader(fname)
    if header['imagetyp'].strip() != 'object' or "solar" in header['object'].lower():
      print "Skipping file %s, with imagetype = %s and object = %s" %(fname, header['imagetyp'], header['object'])
      continue

    logfile = open("fitlog_%s.txt" %(fname.split(".fits")[0]), "a")
    logfile.write("Fitting file %s\n" %(fname))
    name = fname.split(".fits")[0]
    outfilename = "Corrected_%s.fits" %name

    #Read file
    orders = HelperFunctions.ReadFits(fname, errors="error", extensions=True, x="wavelength", y="flux")

    #Get the observation time
    date = header["DATE-OBS"]
    time = header["UT"]
    t_seg = time.split(":")
    time = 3600*float(t_seg[0]) + 60*float(t_seg[1]) + float(t_seg[2])
    

    #Read in weather information (archived data is downloaded from weather.as.utexas.edu)
    infile = open(weather_file)
    lines = infile.readlines()
    infile.close()
    times = []
    RH = []
    P = []
    T = []
    idx = 0
    bestindex = 0
    difference = 9e9
    for line in lines[1:]:
      segments = line.split()
      if date in segments[0]:
        segments = segments[1].split(",")
        t = segments[0]
        t_seg = t.split(":")
        weather_time = 3600*float(t_seg[0]) + 60*float(t_seg[1]) + float(t_seg[2])
        if numpy.abs(time - weather_time) < difference:
          difference = numpy.abs(time - weather_time)
          bestindex = idx
        times.append(segments[0])
        T.append(float(segments[3]))
        RH.append(float(segments[4]))
        P.append(float(segments[5]))
        idx += 1
    
    
    angle = float(header["ZD"])
    resolution = 60000.0
    humidity = max(5, RH[bestindex])
    T_fahrenheit = T[bestindex]
    pressure = P[bestindex]*Units.hPa/Units.inch_Hg
    temperature = (T_fahrenheit - 32.0)*5.0/9.0 + 273.15

    if edit_atmosphere:
      filenames = [f for f in os.listdir("./") if "GDAS" in f]      
      height, Pres, Temp, h2o = GetAtmosphere.GetProfile(filenames, header['date-obs'], header['ut'])

      fitter.EditAtmosphereProfile("Temperature", height, Temp)
      fitter.EditAtmosphereProfile("Pressure", height, Pres)
      fitter.EditAtmosphereProfile("H2O", height, h2o)
      
    #Adjust fitter values
    fitter.AdjustValue({"angle": angle,
                        "pressure": pressure,
                        "resolution": resolution,
                        "temperature": temperature,
			                  "o2": 2.12e5})
    fitter.FitVariable({"h2o": humidity})
    #                    "temperature": temperature})
    fitter.SetBounds({"h2o": [humidity_low, humidity_high],
                      "temperature": [temperature-10, temperature+10],
                      "o2": [5e4, 1e6],
                      "resolution": [50000, 90000]})
    
    #Ignore the interstellar sodium D lines and parts of the O2 bands
    fitter.IgnoreRegions(badregions)
    
    # Determine the H2O abundance
    resolution = []
    h2o = []
    T = []
    o2 = []
    waveshifts = []
    wave0 = []
    chisquared = []
    fitter.DisplayVariables()
    #for i in [30, 35, 39, 41]:
    for i in FindOrderNums(orders, [650, 700, 717, 727]):
      print "\n***************************\nFitting order %i: " %(i)
      order = orders[i]
      fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                          "waveend": order.x[-1] + 20.0})
      order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=3, lowreject=1.5, highreject=10)
      primary = DataStructures.xypoint(x=order.x, y=numpy.ones(order.x.size))
      primary, model, R = fitter.Fit(data=order.copy(), 
                                     resolution_fit_mode="gauss", 
                                     fit_source=True, 
                                     return_resolution=True,
                                     adjust_wave="model",
                                     wavelength_fit_order=4)
      resolution.append(R)
      waveshifts.append(fitter.shift)
      wave0.append(fitter.data.x.mean())
      h2o.append(fitter.GetValue("h2o"))
      T.append(fitter.GetValue("temperature"))
      
      chisquared.append((1.0-min(model.y))/fitter.chisq_vals[-1])

    # Determine the average humidity (weight by chi-squared)
    humidity = numpy.sum(numpy.array(h2o)*numpy.array(chisquared)) / numpy.sum(chisquared)
    temperature = numpy.sum(numpy.array(T)*numpy.array(chisquared)) / numpy.sum(chisquared)
    logfile.write("Humidity/Temperature values and their chi-squared values:\n")
    for h, t, c in zip(h2o, T, chisquared):
      logfile.write("%g\t%g\t%g\n" %(h, t, 1.0/c))
    logfile.write("Best humidity = %.4f\n" %humidity)
    logfile.write("Best temperature = %.4f\n" %temperature)
    fitter.AdjustValue({"h2o": humidity,
                        "temperature": temperature})
    
    # Now, determine the O2 abundance
    fitter.FitVariable({"o2": 2.12e5})
    #for i in [33, 38]:
    for i in FindOrderNums(orders, [630, 695]):
      order = orders[i]
      fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                          "waveend": order.x[-1] + 20.0})
      order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=3, lowreject=1.5, highreject=10)
      primary = DataStructures.xypoint(x=order.x, y=numpy.ones(order.x.size))
      primary, model, R = fitter.Fit(data=order.copy(), 
                                     resolution_fit_mode="gauss", 
                                     fit_source=True,
                                     return_resolution=True,
                                     adjust_wave="model",
                                     wavelength_fit_order=4)
      resolution.append(R)
      waveshifts.append(fitter.shift)
      wave0.append(fitter.data.x.mean())
      idx = fitter.parnames.index("o2")
      o2.append(fitter.const_pars[idx])
      chisquared.append((1.0-min(model.y))/fitter.chisq_vals[-1])

    # Determine the average of the other parameter values
    chi2 = numpy.array(chisquared)
    o2 = numpy.array(o2)
    resolution = numpy.array(resolution)
    waveshifts = numpy.array(waveshifts)
    wave0 = numpy.array(wave0)
    velshifts = waveshifts/wave0 * constants.c.cgs.value*units.cm.to(units.km)
    vel = numpy.sum(velshifts*chi2) / numpy.sum(chi2)
    logfile.write("resolution, velocity shifts and their chi-squared\n")
    for R, v, c in zip(resolution, velshifts, chi2):
      logfile.write("%g\t%g\t%g\n" %(R, v, 1.0/c))
    logfile.write("O2 abundance and their chi-squared:\n")
    for o, c in zip(o2, chi2[-2:]):
      logfile.write("%g\t%g\n" %(o, 1.0/c))
    o2 = numpy.sum(o2*chi2[-2:])/numpy.sum(chi2[-2:])
    resolution = numpy.sum(resolution*chi2)/numpy.sum(chi2)
    logfile.write("Best o2 = %.4f ppmv\n" %o2)
    logfile.write("Best resolution = %.5f\n" %resolution)
    logfile.write("Best velocity shift = %.4f km/s" %vel)
    """
    
    o2 = 224773
    humidity = 42.8873
    resolution = 88472.7266
    vel = 4.8
    """

    # Finally, apply these parameters to all orders in the data
    for i, order in enumerate(orders):
      if (order.x[0] < 470 and order.x[-1] > 470) or max(order.y) < 0.01:
        model = order.copy()
        model.y = numpy.ones(order.size())
        data = order.copy()
        data.cont = numpy.ones(data.size())
      else:
        print "\n\nGenerating model for order %i of %i\n" %(i, len(orders))
        fitter.AdjustValue({"wavestart": order.x[0] - 5.0,
                            "waveend": order.x[-1] + 5.0,
                            "o2": o2,
                            "h2o": humidity,
                            "resolution": resolution})
        fitpars = [fitter.const_pars[j] for j in range(len(fitter.parnames)) if fitter.fitting[j] ]
        order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=3, lowreject=1.5, highreject=10)
        fitter.ImportData(order)
        fitter.resolution_fit_mode = "gauss"
        #wave0 = order.x.mean()
        #fitter.shift = vel/(constants.c.cgs.value*units.cm.to(units.km)) * wave0
        print "fitter.shift = ", fitter.shift
        primary, model = fitter.GenerateModel(fitpars, 
                                              separate_primary=True, 
                                              return_resolution=False)

        data = fitter.data
        if min(model.y) > 0.98:
          #The wavelength calibration might be off
          wave0 = order.x.mean()
          fitter.shift = vel/(constants.c.cgs.value*units.cm.to(units.km)) * wave0
          model = fitter.GenerateModel(fitpars, separate_primary=False, nofit=True)
          model.x /= (1.0 + vel/(constants.c.cgs.value*units.cm.to(units.km)))
          model = FittingUtilities.RebinData(model, order.x)
          data = order.copy()
          data.cont = FittingUtilities.Continuum(data.x, data.y, fitorder=3, lowreject=2, highreject=5)


      # Set up data structures for OutputFitsFile
      columns = {"wavelength": data.x,
                 "flux": data.y,
                 "continuum": data.cont,
                 "error": data.err,
                 "model": model.y,
                 "primary": primary.y}
      
      header_info = []
      numpars = len(fitter.const_pars)
      for j in range(numpars):
        try:
          parname = fitter.parnames[j]
          parval = fitter.const_pars[j]
          fitting = fitter.fitting[j]
          header_info.append([namedict[parname][0], fitting, namedict[parname][2] ])
          header_info.append([namedict[parname][1], parval, namedict[parname][2] ])
        except KeyError:
          print "Not saving the following info: %s" %(fitter.parnames[j])
      
      
      if (i == 0 and makenew) or not exists:
        HelperFunctions.OutputFitsFileExtensions(columns, fname, outfilename, headers_info=[header_info,], mode="new")
        exists = True
      else:
        HelperFunctions.OutputFitsFileExtensions(columns, outfilename, outfilename, headers_info=[header_info,], mode="append")
      
      
    logfile.close()
