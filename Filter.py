import scipy.signal as sig
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import DataStructures
from astropy import units, constants
import FittingUtilities

import FitsUtils
import HelperFunctions


def LowPass():
    fileList = []
    vsini = 40.0
    for arg in sys.argv[1:]:
        if 'vsini' in arg:
            vsini = float(arg.split("=")[-1])
        else:
            fileList.append(arg)

    for fname in fileList:
        column_list = []
        fig = plt.figure(1)
        plotgrid = gridspec.GridSpec(3, 1)
        mainaxis = plt.subplot(plotgrid[0:2])
        reducedaxis = plt.subplot(plotgrid[2], sharex=mainaxis)
        orders = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", errors="error",
                                        cont="continuum")
        for order in orders:
            # Linearize
            datafcn = interp(order.x, order.y, k=1)
            errorfcn = interp(order.x, order.err, k=1)
            linear = DataStructures.xypoint(order.x.size)
            linear.x = np.linspace(order.x[0], order.x[-1], linear.size())
            linear.y = datafcn(linear.x)
            linear.err = errorfcn(linear.x)
            smoothed = HelperFunctions.LowPassFilter(linear, 50 * units.km.to(units.cm))
            smoothed /= smoothed.mean()
            mainaxis.plot(linear.x, linear.y)
            mainaxis.plot(linear.x, smoothed, 'r-', linewidth=1)
            reducedaxis.plot(linear.x, linear.y / smoothed)
            columns = {"wavelength": linear.x,
                       "flux": linear.y,
                       "error": linear.err,
                       "continuum": FittingUtilities.Continuum(linear.x, linear.y, fitorder=3, lowreject=3,
                                                               highreject=3)}
            column_list.append(columns)
        outfilename = "%s_filtered.fits" % (fname.split(".fits")[0])
        print "Outputting to %s" % outfilename
        HelperFunctions.OutputFitsFileExtensions(column_list, fname, outfilename, mode='new')
        plt.show()


def HighPass():
    fileList = []
    vsini = 40.0
    for arg in sys.argv[1:]:
        if 'vsini' in arg:
            vsini = float(arg.split("=")[-1])
        else:
            fileList.append(arg)

    for fname in fileList:
        column_list = []
        fig = plt.figure(1)
        plotgrid = gridspec.GridSpec(3, 1)
        mainaxis = plt.subplot(plotgrid[0:2])
        reducedaxis = plt.subplot(plotgrid[2], sharex=mainaxis)
        orders = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", errors="error",
                                        cont="continuum")
        for order in orders:
            # Linearize
            datafcn = interp(order.x, order.y, k=1)
            errorfcn = interp(order.x, order.err, k=1)
            linear = DataStructures.xypoint(order.x.size)
            linear.x = np.linspace(order.x[0], order.x[-1], linear.size())
            linear.y = datafcn(linear.x)
            linear.err = errorfcn(linear.x)
            linear.cont = FittingUtilities.Continuum(linear.x, linear.y)
            smoothed = HelperFunctions.HighPassFilter(linear, vsini * units.km.to(units.cm))
            mean = np.mean(smoothed)
            std = np.std(smoothed)
            badindices = np.where(np.abs((smoothed - mean) / std > 3.0))[0]
            plt.figure(2)
            plt.plot(linear.x, (smoothed - mean) / std)
            plt.figure(3)
            plt.plot(linear.x, linear.y - smoothed)
            plt.figure(1)
            smoothed[badindices] = 0.0
            smoothed += np.median(linear.cont)
            smoothed /= np.median(linear.cont)
            #linear.y[badindices] = smoothed[badindices]
            mainaxis.plot(linear.x, linear.y / linear.cont, 'k-')
            mainaxis.plot(linear.x, smoothed, 'r-', linewidth=1)
            reducedaxis.plot(linear.x, smoothed)
            columns = {"wavelength": linear.x,
                       "flux": smoothed,
                       "error": linear.err,
                       "continuum": FittingUtilities.Continuum(linear.x, linear.y, fitorder=3, lowreject=3,
                                                               highreject=3)}
            column_list.append(columns)
        outfilename = "%s_filtered.fits" % (fname.split(".fits")[0])
        print "Outputting to %s" % outfilename
        plt.show()
        HelperFunctions.OutputFitsFileExtensions(column_list, fname, outfilename, mode='new')


if __name__ == "__main__":
    HighPass()
