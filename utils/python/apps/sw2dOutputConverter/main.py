import os
import pathlib
import sys

import createCDF
import createTxt
import readBinResult

# Read binary result
bin_result = readBinResult.BinResult()
bin_result.readFile(sys.argv[1] + "/output/ResultHydro.bin")

# Convert result to the required format
format = sys.argv[2].lower()
if (format == "netcdf") :
    createCDF.create_netCDF(bin_result, sys.argv[1])
elif (format == "text") :
    createTxt.create_text(bin_result, sys.argv[1])
else :
    print("Unknown format : " + format)
    print("Allowed format are : netcdf, text")
