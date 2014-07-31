"""
  This script loops through all directories and finds the co-added spectra files.
  It then pulls the x-axis of each order, and plots the x-axis in a different 
  figure for each order, and with different lines for each file. It is used
  as a check to see how stable the CHIRON instrument wavelength calibration is,
  and if I can safely use the faster search method on all data at once.
"""

import os

import matplotlib.pyplot as plt
import astropy.time
from astropy.io import fits

import HelperFunctions


dirs = ['20130326',
        '20130407',
        '20130409',
        '20130412',
        '20130414',
        '20130415',
        '20130416',
        '20130421',
        '20130422',
        '20130423',
        '20130425',
        '20130426',
        '20130427',
        '20130802',
        '20130803',
        '20130806',
        '20130808',
        '20130812',
        '20130813',
        '20130814',
        '20130816',
        '20130818',
        '20130819']

if __name__ == "__main__":
    orders = []
    # dirs = [d for d in os.listdir("./") if d.startswith("2014") and len(d) == 8]
    for d in dirs:
        object_files = [f for f in os.listdir(d) if f.startswith("HIP") and f.endswith(".fits")]

        print d
        for f in object_files:
            header = fits.getheader("%s/%s" % (d, f))
            if header['IMAGETYP'].strip().lower() != "object":
                continue
            data = HelperFunctions.ReadExtensionFits("%s/%s" % (d, f))
            print "\t", f, len(data)
            if len(orders) < 1:
                for order in data:
                    orders.append({d: order.x})
            else:
                for i, order in enumerate(data):
                    orders[i][d] = order.x

    pixel = 1000
    for i, order in enumerate(orders):
        for date in sorted(order.keys()):
            date2 = "%s-%s-%s" % (date[:4], date[4:6], date[6:])
            jd = astropy.time.Time(date2, scale='utc', format='iso').jd
            dx = order[date][pixel + 1] - order[date][pixel]
            plt.plot(jd, dx, 'ro')
            #plt.plot(jd, order[date][pixel], 'ro')
        plt.xlabel("Julian Date")
        plt.ylabel("Delta - Wavelength at pixel %i (nm)" % (pixel))
        plt.title("Order %i/%i" % (i + 1, len(orders)))
        #plt.ylabel("Wavelength at pixel %i (nm)" %(pixel))
        ax = plt.gca()
        ax.ticklabel_format(style='sci', useOffset=False)
        plt.show()
