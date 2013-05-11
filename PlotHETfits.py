import pyfits
import matplotlib.pyplot as plt
import sys
import FitsUtils


if __name__ == "__main__":
  fileList = []
  for arg in sys.argv[1:]:
    fileList.append(arg)

  for fname in fileList:
    orders = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
    for order in orders:
      print order.x.mean()
      plt.plot(order.x, order.y/order.cont)
    plt.show()
