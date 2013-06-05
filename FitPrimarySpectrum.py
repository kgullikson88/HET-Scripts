import numpy
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import sys
import os
import FitsUtils


if __name__ == "__main__":
  #Parse command line arguments
  fileList = []
  Tmin, Tmax, Tstep = 8000, 8800, 200
  Zmin, Zmax, Zstep = -0.5, 0.5, 0.5
  loggmin, loggmax, loggstep = 4.0, 4.0, 1.0
  for arg in sys.argv[1:]:
    if "Tmin" in arg:
      Tmin = float(arg.split("=")[-1])
    elif "Tmax" in arg:
      Tmax = float(arg.split("=")[-1])
    elif "Tstep" in arg:
      Tstep = float(arg.split("=")[-1])
    elif "Zmin" in arg:
      Zmin = float(arg.split("=")[-1])
    elif "Zmax" in arg:
      Zmax = float(arg.split("=")[-1])
    elif "Zstep" in arg:
      Zstep = float(arg.split("=")[-1])
    elif "loggmin" in arg:
      loggmin = float(arg.split("=")[-1])
    elif "loggmax" in arg:
      loggmax = float(arg.split("=")[-1])
    elif "loggstep" in arg:
      loggstep = float(arg.split("=")[-1])
    else:
      fileList.append(arg)

  for fname in fileList:
    orders = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", errors="error")
    
