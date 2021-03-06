import sys

from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
import HelperFunctions
import FittingUtilities


if __name__ == "__main__":
    fileList = []
    tellurics = False
    normalize = False
    byorder = False  # Plots one order at a time
    pixelscale = False
    oneplot = False
    for arg in sys.argv[1:]:
        if "tellcorr" in arg:
            tellurics = True
        elif "-norm" in arg:
            normalize = True
        elif "-order" in arg:
            byorder = True
        elif "-pix" in arg:
            pixelscale = True
            # byorder = True
        elif "-one" in arg:
            oneplot = True
        else:
            fileList.append(arg)

    linestyles = ['k-', 'r-', 'b-', 'g-']

    for fnum, fname in enumerate(fileList):
        ls = linestyles[fnum % len(linestyles)]
        orders = HelperFunctions.ReadExtensionFits(fname)
        print fname, len(orders)
        if not oneplot:
            plt.figure(fnum)
            plt.title(fname)
        if tellurics:
            model = HelperFunctions.ReadFits(fname, extensions=True, x="wavelength", y="model")
        for i, order in enumerate(orders):

            # order.cont = FindContinuum.Continuum(order.x, order.y, lowreject=3, highreject=3)
            if pixelscale:
                order.x = np.arange(order.size())
            if tellurics:
                plt.plot(order.x, order.y / order.cont, 'k-')
                plt.plot(order.x, model[i].y, 'r-')
            else:
                if normalize:
                    plt.plot(order.x, order.y / order.cont, 'k-')
                    plt.text(order.x.mean(), 1.1, str(i + 1))
                else:
                    plt.plot(order.x, order.y, ls)
                    #plt.plot(order.x, order.cont)
            plt.xlabel("Wavelength (nm)")
            plt.ylabel("Flux")
            plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
            if byorder:
                plt.title("Order %i" % i)
                plt.show()
    if not byorder:
        plt.show()
