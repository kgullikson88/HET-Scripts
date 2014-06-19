import FitsUtils
import numpy
import sys
import matplotlib.pyplot as plt
import HelperFunctions

trimming = {
           21: [0, 576.3],
           22: [582.6, 9e9],
           23: [583.45, 9e9],
           37: [0, 680.15]}


class Trimmer:
  def __init__(self, data=None):
    if data != None:
      self.data = data.copy()
    self.clicks = []
    logfile = open("trimlog.dat", "w")
    logfile.close()
    
  def InputData(self, data):
    self.data = data.copy()
    
  def Plot(self):
    self.fig =  plt.figure(1, figsize=(11,10))
    cid = self.fig.canvas.mpl_connect('key_press_event', self.keypress)
    plt.plot(self.data.x, self.data.y)
    plt.plot(self.data.x, self.data.cont)
    plt.show()
    return self.data.copy()

  def keypress(self, event):
    if event.key == "r":
      print "Set to remove points. Click on the bounds"
      self.clipmode = "remove"
      self.clickid = self.fig.canvas.mpl_connect('button_press_event', self.mouseclick)
    elif event.key == "i":
      print "Set to interpolate between points. Click the bounds"
      self.clipmode = "interpolate"
      self.clickid = self.fig.canvas.mpl_connect('button_press_event', self.mouseclick)
    

  def mouseclick(self, event):
    self.clicks.append(event.xdata)

    if len(self.clicks) == 2:
      left = max(0, numpy.searchsorted(self.data.x, min(self.clicks)))
      right = min(self.data.size()-1, numpy.searchsorted(self.data.x, max(self.clicks)))

      
      logfile = open("trimlog.dat", "a")
      if self.clipmode == "remove":
        logfile.write("Removing:\t%.3f  to %.3f\n" %(min(self.clicks), max(self.clicks)))
        self.data.x = numpy.delete(self.data.x, numpy.arange(left,right+1))
        self.data.y = numpy.delete(self.data.y, numpy.arange(left,right+1))
        self.data.cont = numpy.delete(self.data.cont, numpy.arange(left,right+1))
        self.data.err = numpy.delete(self.data.err, numpy.arange(left,right+1))
      elif self.clipmode == "interpolate":
        logfile.write("Interpolating:\t%.3f  to %.3f\n" %(min(self.clicks), max(self.clicks)))
        x1, x2 = self.data.x[left], self.data.x[right]
        y1, y2 = self.data.y[left], self.data.y[right]
        m = (y2 - y1) / (x2 - x1)
        self.data.y[left:right] = m*(self.data.x[left:right] - x1) + y1
        self.data.cont[left:right] = m*(self.data.x[left:right] - x1) + y1
        
      self.fig.clf()
      cid = self.fig.canvas.mpl_connect('key_press_event', self.keypress)
      plt.plot(self.data.x, self.data.y)
      plt.plot(self.data.x, self.data.cont)
      plt.draw()
      
      self.fig.canvas.mpl_disconnect(self.clickid)
      self.clicks = []
      logfile.close()
           

def main1():
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
        HelperFunctions.OutputFitsFileExtensions(columns, fname, outfilename, mode="new")
      else:
        HelperFunctions.OutputFitsFileExtensions(columns, outfilename, outfilename, mode="append")
        



def Interactive():
  trim = Trimmer()
  for fname in sys.argv[1:]:
    if "-" in fname:
      num = int(fname.split("-")[-1].split(".fits")[0])
      outfilename = "%s-%i.fits" %(fname.split("-")[0], num+1)
    else:
      outfilename = "%s-0.fits" %(fname.split(".fits")[0])
    orders = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", errors="error", cont="continuum")
    #trim = Trimmer(orders[0])
    
    logfile = open("trimlog.dat", "a")
    logfile.write("\n\n\n******************************************\n")
    logfile.write("\nTrimming file %s\n\n" %(fname))
    logfile.write("******************************************\n")
    logfile.close()
    
    for i, order in enumerate(orders):
      logfile = open("trimlog.dat", "a")
      logfile.write("********   Order %i  ******************\n" %(i+1))
      logfile.close()
      trim.InputData(order)
      order = trim.Plot()

      columns = {"wavelength": order.x,
                 "flux": order.y,
                 "continuum": order.cont,
                 "error": order.err}
      
      if i == 0:
        HelperFunctions.OutputFitsFileExtensions(columns, fname, outfilename, mode="new")
      else:
        HelperFunctions.OutputFitsFileExtensions(columns, outfilename, outfilename, mode="append")



def Automatic():
  badorders = [16, 39, 42, 45, 46, 67, 69]
  badwaves = [[461.2, 465.8],
              [590.7, 593],
              [703.4, 707.7],
              ]
  trimsize = 10
  for fname in sys.argv[1:]:
    orders = HelperFunctions.ReadExtensionFits(fname)
    
    # Remove the bad orders
    for n in badorders[::-1]:
      orders.pop(n)

    #Trim the bad parts of the spectrum
    column_list = []
    for i, order in enumerate(orders):
      order = order[trimsize:-trimsize]
      for badsection in badwaves:
        left = numpy.searchsorted(order.x, badsection[0])
        right = numpy.searchsorted(order.x, badsection[1])
        order.x = numpy.delete(order.x, numpy.arange(left, right))
        order.y = numpy.delete(order.y, numpy.arange(left, right))
        order.cont = numpy.delete(order.cont, numpy.arange(left, right))
        order.err = numpy.delete(order.err, numpy.arange(left, right))
      
      columns = {"wavelength": order.x,
                 "flux": order.y,
                 "continuum": order.cont,
                 "error": order.err}
      column_list.append(columns)


    #Output:
    if "-" in fname:
      n = int(fname.split("-")[-1].split(".fits")[0])
      n = n + 1
      base = fname.split("-")[0]
    else:
      n = 0
      base = fname.split(".fits")[0]
    outfilename = "%s-%i.fits" %(base, n)
    print "Outputting to %s" %outfilename

    HelperFunctions.OutputFitsFileExtensions(column_list, fname, outfilename, mode="new")



if __name__ == "__main__":
  #Interactive()
  Automatic()
