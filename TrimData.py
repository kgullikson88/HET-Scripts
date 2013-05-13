import FitsUtils
import numpy
import sys

trimming = {
           21: [0, 576.3],
           22: [582.6, 9e9],
           23: [583.45, 9e9],
           37: [0, 680.15]}
           

if __name__ == "__main__":
  for fname in sys.argv[1:]:
    if "-" in fname:
      num = int(fname.split("-")[-1].split(".fits")[0])
      outfilename = "%s-%i.fits" %(fname.split("-")[0], num+1)
    else:
      outfilename = "%s-0.fits" %(fname.split(".fits")[0])
    orders = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", errors="error", cont="continuum")
    for i, order in enumerate(orders):
      if i in trimming.keys():
        left = numpy.searchsorted(order.x, trimming[i][0])
        right = numpy.searchsorted(order.x, trimming[i][1])
        order.x = order.x[left:right]
        order.y = order.y[left:right]
        order.cont = order.cont[left:right]
        order.err = order.err[left:right]
        orders[i] = order.copy()
      
      columns = {"wavelength": order.x,
	         "flux": order.y,
                 "continuum": order.cont,
                 "error": order.err}
      
      if i == 0:
        FitsUtils.OutputFitsFileExtensions(columns, fname, outfilename, mode="new")
      else:
        FitsUtils.OutputFitsFileExtensions(columns, outfilename, outfilename, mode="append")
        
