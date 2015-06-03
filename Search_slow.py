import sys

import GenericSearch
import StarData

# Define regions contaminated by telluric residuals or other defects. We will not use those regions in the cross-correlation
badregions = [[567.5, 575.5],
              [588.5, 598.5],
              [627, 632],
              [647, 655],
              [686, 706],
              [716, 734],
              [759, 9e9],
              # [655, 657],  # H alpha
              # [485, 487],  #H beta
              # [433, 435],  #H gamma
              # [409, 411],  #H delta
              #[396, 398],  #H epsilon
              #[388, 390],  #H zeta
]
interp_regions = []
trimsize = 10

if "darwin" in sys.platform:
    modeldir = "/Volumes/DATADRIVE/Stellar_Models/PHOENIX/Stellar/Vband/"
    hdf5_filename = '/Volumes/DATADRIVE/PhoenixGrid/HRS_Grid.hdf5'
elif "linux" in sys.platform:
    modeldir = "/media/FreeAgent_Drive/SyntheticSpectra/Sorted/Stellar/Vband/"
    hdf5_filename = '/media/ExtraSpace/PhoenixGrid/HRS_Grid.hdf5'
else:
    modeldir = raw_input("sys.platform not recognized. Please enter model directory below: ")
    if not modeldir.endswith("/"):
        modeldir = modeldir + "/"

if __name__ == '__main__':
    # Parse command line arguments:
    fileList = []
    for arg in sys.argv[1:]:
        if 1:
            fileList.append(arg)

    # Get the primary star vsini values
    prim_vsini = StarData.get_vsini(fileList)

    GenericSearch.slow_companion_search(fileList, prim_vsini,
                                        hdf5_file=hdf5_filename,
                                        extensions=True,
                                        resolution=None,
                                        trimsize=trimsize,
                                        modeldir=modeldir,
                                        badregions=badregions,
                                        metal_values=(0.0, -0.5, 0.5),
                                        vsini_values=(1, 5.0, 10.0, 20.0, 30.0),
                                        Tvalues=range(3000, 9000, 100),
                                        observatory='McDonald',
                                        debug=False,
                                        vbary_correct=True,
                                        addmode='simple',
                                        output_mode='hdf5')

