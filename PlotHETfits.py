import pyfits
import matplotlib.pyplot as plt
import sys
import FitsUtils


if __name__ == "__main__":
  fileList = []
  tellurics = False
  for arg in sys.argv[1:]:
    if "tell" in arg:
      tellurics = True
    else:
      fileList.append(arg)

  for fname in fileList:
    orders = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
    if tellurics:
      model = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="model")
    for i, order in enumerate(orders):
      if tellurics:
        plt.plot(order.x, order.y/order.cont/model[i].y)
      else:
        plt.plot(order.x, order.y/order.cont)
      
    plt.show()
