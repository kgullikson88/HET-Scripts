import sys
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from scipy.optimize import leastsq
import os
import warnings

from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
import DataStructures
import FittingUtilities
import numpy as np
import MakeModel
import TelluricFitter

import HelperFunctions


tf = TelluricFitter.TelluricFitter()


def ReadCorrectedFile(fname, yaxis="model"):
    orders = []
    headers = []
    hdulist = pyfits.open(fname)
    numorders = len(hdulist)
    for i in range(1, numorders):
        order = hdulist[i].data
        xypt = DataStructures.xypoint(x=order.field("wavelength"),
                                      y=order.field(yaxis),
                                      cont=order.field("continuum"),
                                      err=order.field("error"))

        orders.append(xypt)
        headers.append(hdulist[i].header)
    return orders, headers


# def FixWavelength(data, model, fitorder=3):
#  modelfcn = tf.FitWavelength(data, model, fitorder=fitorder)
#  return modelfcn

def FitGaussian(data):
    A = 1.0 - min(data.y)
    sigma = 0.004
    w0 = data.x[data.y == min(data.y)]
    Baseline = 1.0
    pars = [A, w0, sigma, Baseline]
    gauss = lambda x, pars: pars[3] - pars[0] * np.exp(-(x - pars[1]) ** 2 / (2 * pars[2] ** 2))
    errfcn = lambda pars, d: gauss(d.x, pars) - d.y
    pars, success = leastsq(errfcn, pars, args=(data), diag=1.0 / np.array(pars), epsfcn=1e-10)
    #if A > 0.7:
    #  data.output("Testdata.dat")
    #  sys.exit()
    return pars, success


def WavelengthErrorFunction(pars, wavelengths, offsets, weights, maxdiff):
    fcn = np.poly1d(pars)
    prediction = fcn(wavelengths)
    penalty = np.sum(np.abs(prediction[np.abs(prediction) > maxdiff]))
    return (offsets - prediction) ** 2 + penalty


def FitPoly(wavelengths, offsets, weights, maxdiff=0.5, fitorder=3):
    pars = np.zeros(fitorder + 1)
    pars, success = leastsq(WavelengthErrorFunction, pars, args=(wavelengths, offsets, weights, maxdiff))
    return np.poly1d(pars)


def FixWavelength(data, model, fitorder=3, tol=20, numiters=5):
    #Find the lines from the input model
    linelist = FittingUtilities.FindLines(model, debug=False, tol=0.985, linespacing=0.05)
    linelist = model.x[linelist]

    gauss = lambda x, pars: pars[3] - pars[0] * np.exp(-(x - pars[1]) ** 2 / (2 * pars[2] ** 2))

    model_lines = []
    dx = []
    #Loop over lines
    for line in linelist:
        idx = np.searchsorted(model.x, line)
        if idx < tol or model.size() - idx < tol:
            continue
        pars, model_success = FitGaussian(model[idx - tol:idx + tol])
        if model_success < 5 and pars[0] > 0 and pars[0] < 1:
            model_lines.append(pars[1])
        else:
            continue

        idx = np.searchsorted(data.x, line)
        pars, data_success = FitGaussian(data[idx - tol:idx + tol])
        if data_success < 5 and pars[0] > 0 and pars[0] < 1:
            dx.append(pars[1] - model_lines[-1])
        else:
            model_lines.pop()

    model_lines = np.array(model_lines)
    dx = np.array(dx)

    mean, std = np.mean(dx), np.std(dx)
    badindices = np.where(np.abs(dx) > 0.01)[0]
    model_lines = np.delete(model_lines, badindices)False
    dx = np.delete(dx, badindices)

    print "Found %i lines" % len(model_lines)
    if len(model_lines) < 3 * fitorder:
        return lambda x: x, 0.0


    #Iteratively fit
    done = False
    mean = data.x.mean()
    iternum = 0
    mean = 0
    fit = lambda x: 0.0  #Default value
    while not done and len(model_lines) >= fitorder and iternum < numiters:
        iternum += 1
        done = True
        fit = np.poly1d(np.polyfit(model_lines - mean, dx, fitorder))
        residuals = fit(model_lines - mean) - dx
        std = np.std(residuals)
        badindices = np.where(np.abs(residuals) > 3 * std)[0]
        if badindices.size > 0 and model_lines.size - badindices.size > 2 * fitorder:
            done = False
            model_lines = np.delete(model_lines, badindices)
            dx = np.delete(dx, badindices)
    return lambda x: x + fit(x - mean), 0


def Correct_HET(original, corrected, offset=None, plot=True, get_primary=False):
    #Read in the data and model
    original_orders = HelperFunctions.ReadFits(original, extensions=True, x="wavelength", y="flux", errors="error",
                                               cont="continuum")
    corrected_orders, corrected_headers = ReadCorrectedFile(corrected)
    primary_orders, primary_headers = ReadCorrectedFile(corrected, yaxis="primary")
    print len(original_orders), len(corrected_orders)
    if offset == None:
        offset = len(original_orders) - len(corrected_orders)
    offset = 0

    if plot:
        fig = plt.figure(1)
        ax = fig.add_subplot(111)

    for i in range(offset, len(original_orders)):
        print "Order %i/%i" % (i + 1, len(original_orders))
        data = original_orders[i]
        data.cont = FittingUtilities.Continuum(data.x, data.y)
        try:
            model = corrected_orders[i - offset]
            #primary = primary_orders[i-offset]
            header = corrected_headers[i - offset]
            if i == 0:
                print "Order = %i\nHumidity: %g\nO2 concentration: %g\n" % (i, header['h2oval'], header['o2val'])
        except IndexError:
            model = DataStructures.xypoint(x=data.x, y=np.ones(data.x.size))
            #primary = DataStructures.xypoint(x=data.x, y=np.ones(data.x.size))
            print "Warning!!! Telluric Model not found for order %i" % i
        primary = data.copy()
        primary.y = FittingUtilities.Iterative_SV(data.y / model.y, 51, 3)

        if model.size() < data.size():
            left = np.searchsorted(data.x, model.x[0])
            right = np.searchsorted(data.x, model.x[-1])
            if right < data.size():
                right += 1
            data = data[left:right]
        elif model.size() > data.size():
            sys.exit("Error! Model size (%i) is larger than data size (%i)" % (model.size(), data.size()))

        #if np.sum((model.x-data.x)**2) > 1e-8:
        #  model = FittingUtilities.RebinData(model, data.x)

        #Adjust the wavelength calibration if necessary
        if min(model.y[1:-1]) < 0.98:
            data2 = data.copy()
            if plot:
                ax.plot(model.x, model.y, 'k--')
            data2.y /= primary.y
            #data2.y /= FittingUtilities.Continuum(data2.x, data2.y, lowreject=1, fitorder=1)
            modelfcn, mean = FixWavelength(data2, model, fitorder=2)
            test = modelfcn(model.x - mean)
            xdiff = [test[j] - test[j - 1] for j in range(1, len(test) - 1)]
            if min(xdiff) > 0 and min(test) > 0 and abs(test[0] - data.x[0]) < 0.15 and abs(
                            test[-1] - data.x[-1]) < 0.15:
                model.x = test.copy()
            model = FittingUtilities.RebinData(model, data.x)
            #if shift < 0:
            #  data = data[:shift]
            #  model = model[-shift:]
            #elif shift > 0:
            #  data = data[shift:]
            #  model = model[:-shift]

        if plot:
            ax.plot(data.x, data.y / data.cont)
            #ax.plot(primary.x, primary.y)
            #      if min(model.y) >= 0.99:
            #        plt.plot(data.x, data.y/data.cont)
            #      else:
            #        plt.plot(data2.x, data2.y)
            ax.plot(model.x, model.y)

        data.y[data.y / data.cont < 1e-5] = 1e-5 * data.cont[data.y / data.cont < 1e-5]
        badindices = np.where(np.logical_or(data.y <= 0, model.y < 0.05))[0]
        model.y[badindices] = data.y[badindices] / data.cont[badindices]

        data.y /= model.y
        data.err /= model.y
        if get_primary:
            data.y /= primary.y
        original_orders[i] = data.copy()
    if plot:
        plt.show()
        #sys.exit()
    return original_orders



def Correct(original, corrected, offset=None, get_primary=False, interpolate=True, adjust=True, plot=False):
    # Read in the data and model
    original_orders = HelperFunctions.ReadFits(original, extensions=True, x="wavelength", y="flux", errors="error",
                                               cont="continuum")
    corrected_orders, corrected_headers = ReadCorrectedFile(corrected)
    test_orders, header = ReadCorrectedFile(corrected, yaxis="flux")

    if plot:
        for order, model in zip(test_orders, corrected_orders):
            plt.plot(order.x, order.y / order.cont)
            plt.plot(model.x, model.y)
        plt.title("Correction in corrected file only")
        plt.show()

    if get_primary:
        primary_orders = ReadCorrectedFile(corrected, yaxis="primary")[0]
    if offset == None:
        offset = len(original_orders) - len(corrected_orders)
    print "Offset = ", offset
    for i in range(len(original_orders) - offset):
        data = original_orders[i]
        data.cont = FittingUtilities.Continuum(data.x, data.y)
        try:
            model = corrected_orders[i]
            header = corrected_headers[i]
            if get_primary:
                primary = primary_orders[i]
            if i == 0:
                print "Order = %i\nHumidity: %g\nO2 concentration: %g\n" % (i, header['h2oval'], header['o2val'])
        except IndexError:
            model = DataStructures.xypoint(x=data.x, y=np.ones(data.x.size))
            print "Warning!!! Telluric Model not found for order %i" % i

        if plot:
            plt.figure(1)
            plt.plot(data.x, data.y / data.cont)
            plt.plot(model.x, model.y)
            if get_primary:
                plt.plot(primary.x, primary.y)

        if model.size() < data.size():
            left = np.searchsorted(data.x, model.x[0])
            right = np.searchsorted(data.x, model.x[-1])
            if right < data.size():
                right += 1
            data = data[left:right]
        elif model.size() > data.size() and not interpolate:
            sys.exit("Error! Model size (%i) is larger than data size (%i)" % (model.size(), data.size()))

        if interpolate:
            fcn = spline(model.x, model.y, k=1)
            model = data.copy()
            model.y = fcn(data.x)
            if get_primary:
                fcn = spline(primary.x, primary.y, k=1)
                primary = data.copy()
                primary.y = fcn(primary.x)

        data.y[data.y / data.cont < 1e-5] = 1e-5 * data.cont[data.y / data.cont < 1e-5]
        badindices = np.where(np.logical_or(data.y <= 0, model.y < 0.05))[0]
        model.y[badindices] = data.y[badindices] / data.cont[badindices]
        model.y[model.y < 1e-5] = 1e-5

        if get_primary:
            data.y /= primary.y

        if adjust:
            model.cont = np.ones(model.size())
            lines = FittingUtilities.FindLines(model, tol=0.95).astype(int)
            if len(lines) > 5:
                scale = np.median(np.log(data.y[lines] / data.cont[lines]) / np.log(model.y[lines]))
            else:
                scale = 1.0
            print i, scale
            model.y = model.y ** (1.0/scale)

        #plt.plot(data.x, data.y / model.y)
        data.y /= model.y
        data.err /= model.y
        if get_primary:
            data.y *= primary.y
        original_orders[i] = data.copy()
    if plot:
        plt.show()
    return original_orders



def main1():
    primary = True
    plot = True
    if len(sys.argv) > 2:
        original = sys.argv[1]
        corrected = sys.argv[2]
        if len(sys.argv) > 3 and "prim" in sys.argv[3]:
            primary = True

        outfilename = "%s_telluric_corrected.fits" % (original.split(".fits")[0])
        print "Outputting to %s" % outfilename

        corrected_orders = Correct(original, corrected, offset=None, get_primary=primary, plot=plot, adjust=True)

        column_list = []
        if plot:
            plt.figure(2)
        for i, data in enumerate(corrected_orders):
            if plot:
                plt.plot(data.x, data.y / data.cont)
                #plt.plot(data.x, data.cont)
            #Set up data structures for OutputFitsFile
            columns = {"wavelength": data.x,
                       "flux": data.y,
                       "continuum": data.cont,
                       "error": data.err}
            column_list.append(columns)
        if plot:
            plt.title("Corrected data")
            plt.show()
        HelperFunctions.OutputFitsFileExtensions(column_list, original, outfilename, mode="new")

    else:
        allfiles = os.listdir("./")
        corrected_files = [f for f in allfiles if "Corrected_HRS" in f and f.endswith("-0.fits")]
        #original_files = [f for f in allfiles if any(f in cf for cf in corrected_files)]

        for fname in corrected_files:
            try:
                original = fname.split("Corrected_")[1]
                header = pyfits.getheader(original)
                if header['IMAGETYP'].lower() != "object":
                    warnings.warn("Skipping file %s with image type %s" % (original, header['IMAGETYP']))
                corrected = Correct(original, fname, offset=None, plot=plot, get_primary=primary)
            except IOError:
                warnings.warn("No original file matching the Corrected file %s" % fname)
                continue

            outfilename = "%s_telluric_corrected.fits" % (original.split(".fits")[0])
            print "Outputting to %s" % outfilename

            column_list = []
            for i, data in enumerate(corrected):
                if plot:
                    plt.plot(data.x, data.y / data.cont)
                #Set up data structures for OutputFitsFile
                columns = {"wavelength": data.x,
                           "flux": data.y,
                           "continuum": data.cont,
                           "error": data.err}
                column_list.append(columns)
            HelperFunctions.OutputFitsFileExtensions(column_list, original, outfilename, mode="new")

            if plot:
                plt.title(original)
                plt.xlabel("Wavelength (nm)")
                plt.ylabel("Flux")
                plt.show()


if __name__ == "__main__":
    main1()
