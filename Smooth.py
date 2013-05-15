import FitsUtils
import numpy
import matplotlib.pyplot as plt
import sys


if __name__ == "__main__":
  window_len = 100
  numstd = 2.0
  numiters = 10
  for fname in sys.argv[1:]:
    orders = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")

    for order in orders:
      s = numpy.r_[order.y[window_len/2:0:-1], order.y, order.y[-1:-window_len/2:-1]]
      done = False
      i = 0
      while not done and i < numiters:
        done = True
        i += 1
        #w1 = numpy.ones(window_len)
        w2 = numpy.hanning(window_len)
        #y1 = numpy.convolve(w1/w1.sum(), s, mode='valid')
        y2 = numpy.convolve(w2/w2.sum(), s, mode='valid')

        reduced = order.y/y2
        sigma = numpy.std(reduced)
        mean = numpy.mean(reduced)
        badindices = numpy.where((reduced - mean)/sigma < -numstd)[0]
        #print badindices
        #plt.plot((reduced - mean)/sigma)
        #plt.plot(numpy.ones(reduced.size)*-numstd)
        #plt.show()
        #inp = raw_input("done")
        if badindices.size > 0:
          done = False
          order.y[badindices] = y2[badindices]
        
      #print order.x.size, y.size
      plt.plot(order.x, order.y, 'k-')
      #plt.plot(order.x, y1, 'r-')
      plt.plot(order.x, y2, 'g-')
      plt.show()

