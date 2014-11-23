import sys

import GenericSearch



#Define regions contaminated by telluric residuals or other defects. We will not use those regions in the cross-correlation
badregions = [[585.4, 595],
              [627.5, 632],
              [686, 706],
              [759, 9e9],
              #[655, 657],   #H alpha
              #[485, 487],   #H beta
              #[433, 435],   #H gamma
              #[409, 411],   #H delta
              #[396, 398],   #H epsilon
              #[388, 390],   #H zeta
              ]


if "darwin" in sys.platform:
    modeldir = "/Volumes/DATADRIVE/Stellar_Models/Sorted/Stellar/Vband/"
elif "linux" in sys.platform:
    modeldir = "/media/FreeAgent_Drive/SyntheticSpectra/Sorted/Stellar/Vband/"
else:
    modeldir = raw_input("sys.platform not recognized. Please enter model directory below: ")
    if not modeldir.endswith("/"):
        modeldir = modeldir + "/"

if __name__ == "__main__":
    # Parse command line arguments:
    fileList = []
    extensions = True
    tellurics = False
    trimsize = 1
    for arg in sys.argv[1:]:
        if "-e" in arg:
            extensions = False
        if "-t" in arg:
            tellurics = True  # telluric lines modeled but not removed
        else:
            fileList.append(arg)

    GenericSearch.CompanionSearch(fileList,
                                  extensions=extensions,
                                  resolution=60000.0,
                                  trimsize=trimsize,
                                  vsini_values=[1.0, 10.0, 20.0, 30.0, 40.0],
                                  observatory="McDonald",
                                  vbary_correct=True,
                                  debug=False,
                                  badregions=badregions,
                                  modeldir=modeldir)


