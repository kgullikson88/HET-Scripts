import pyfits
import matplotlib.pyplot as plt
import sys
import FitsUtils
import FindContinuum


if __name__ == "__main__":
  fileList = []
  tellurics = False
  for arg in sys.argv[1:]:
    if "tellcorr" in arg:
      tellurics = True
    else:
      fileList.append(arg)

  for fname in fileList:
    orders = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
    print fname, len(orders)
    if tellurics:
      model = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="model")
    for i, order in enumerate(orders):
      order.cont = FindContinuum.Continuum(order.x, order.y)
      if tellurics:
        plt.plot(order.x, order.y/order.cont, 'k-')
        plt.plot(order.x, model[i].y, 'r-')
      else:
        plt.plot(order.x, order.y/order.cont)
      
    plt.show()
