import pyfits
import FitsUtils
import sys
import numpy


if __name__ == "__main__":
  fileList = []
  outfileset = False
  for arg in sys.argv[1:]:
    if '-out' in arg:
      outfilename = arg.split("=")[-1]
      outfileset = True
    else:
      fileList.append(arg)

  if not outfileset:
    outfilename = raw_input("Enter name of the output file: ")


  f1, f2 = fileList[0], fileList[1]
  hdulist1 = pyfits.open(f1)
  hdulist2 = pyfits.open(f2)

  #Need to make sure there are not duplicate orders
  orders1 = []
  orders2 = []

  for i in range(1, len(hdulist1)):
    orders1.append(hdulist1[i].data.field("wavelength").mean())
    print orders1[-1]
    
  print "\n\n"
  for i in range(1, len(hdulist2)):
    orders2.append(hdulist2[i].data.field("wavelength").mean())
    print orders2[-1]

  if numpy.mean(orders1) > numpy.mean(orders2):
    temp = orders2
    orders2 = orders1
    orders1 = temp

    temp = hdulist2
    hdulist2 = hdulist1
    hdulist1 = temp

  print numpy.mean(orders1), numpy.mean(orders2)
  done = False
  while not done:
    if orders1[-1] >= orders2[0]:
      print "duplicate order found!"
      orders1.pop()
      hdulist1.pop()
    else:
      done = True
  
  #Now, append all orders of hdulist2 onto hdulist1, and save as outfilename
  for order in hdulist2[1:]:
    hdulist1.append(order)

  hdulist1.writeto(outfilename, clobber=True)


  hdulist1.close()
  hdulist2.close()
