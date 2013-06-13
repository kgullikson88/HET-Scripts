import numpy
import scipy.signal as sig
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys
import FitsUtils
import DataStructures
from astropy import units, constants



def MakeFilter(x, x1, x2, passmode="high"):
  retarray = numpy.zeros(x.size)
  a = numpy.pi/(x2-x1)
  b = x1*numpy.pi/(x2-x1)
  if passmode == "high":
    retarray[x>=x2] = 1.0
  elif passmode == "low":
    retarray[x<=x1] = 1.0
  retarray[numpy.logical_and(x > x1, x < x2)] = 0.5*(1+numpy.cos(a*x[numpy.logical_and(x > x1, x < x2)]+b))

  return retarray


def Filter(data, vel, delta = 0.005):
  #Figure out cutoff frequency from the velocity.
  featuresize = 2*data.x.mean()*vel/constants.c.cgs.value    #vel MUST be given in units of cm
  dlam = data.x[1] - data.x[0]   #data.x MUST have constant x-spacing
  Npix = featuresize / dlam
  lowfreq = 1.0/Npix
  highfreq = lowfreq + delta
  print "\nNumber of pixels in feature: %i\nCutoff frequency = %g - %g" %(Npix, lowfreq, highfreq)

  #Extend the y axis to remove edge effects, and fourier transform
  ylong = numpy.append( numpy.append(data.y[::-1], data.y), data.y[::-1])
  fft = numpy.fft.fft(ylong)
  freq = numpy.fft.fftfreq(ylong.size)

  #Low-pass filter the fft.
  filt = MakeFilter(numpy.abs(freq), lowfreq, lowfreq+0.005, passmode="low")
  fft *= filt

  ysmooth = numpy.fft.ifft(fft)
  sign = numpy.sign(ysmooth.real)
  smoothed = ( sign * numpy.sqrt(ysmooth.real**2 + ysmooth.imag**2) )[data.size():-data.size()]
  data.y = smoothed
  return data
  

def main1():
  fileList = []
  vsini = 40.0
  for arg in sys.argv[1:]:
    if 'vsini' in arg:
      vsini = float(arg.split("=")[-1])
    else:
      fileList.append(arg)
      
  for fname in fileList:
    fig = plt.figure()
    plotgrid = gridspec.GridSpec(3,1)
    mainaxis = plt.subplot(plotgrid[0:2])
    reducedaxis = plt.subplot(plotgrid[2], sharex = mainaxis, sharey = mainaxis)

    column_dicts = []
    orders = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", errors="error")
    for order in orders:
      flux = interp(order.x, order.y)
      x = numpy.linspace(order.x[0], order.x[-1], order.x.size)
      data = DataStructures.xypoint(x=x, y=flux(x))

      smoothed = Filter(data.copy(), vsini*units.km.to(units.cm))
      norm = smoothed.y.max()
      smoothed.y /= norm

      reduced = data.y/smoothed.y

      #Make columns for output
      columns = {"wavelength": data.x,
                 "flux": reduced,
                 "continuum": data.cont,
                 "error": data.err}
      column_dicts.append(columns)
      
      reducedaxis.plot(data.x,reduced)
      #smoothed = y/y2 * y.mean()
      mainaxis.plot(data.x,data.y, 'k-')
      mainaxis.plot(data.x,smoothed.y*norm, 'r-')
    plt.show()

    #Output
    outfilename = "%s_noprimary%i.fits" %(fname.split(".fits")[0], vsini)
    print "Outputting to %s" %outfilename
    FitsUtils.OutputFitsFileExtensions(column_dicts, fname, outfilename, mode="new")



if __name__ == "__main__":
  main1()
