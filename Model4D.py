#!/usr/bin/env python3

# Model4D.py
# Python base class to read some model output for 4-D visualization.
#
# Author: Carl Drews, Atmospheric Chemistry Observations & Modeling
# Created: March 2020
# Copyright 2020 by the University Corporation for Atmospheric Research



import sys



# index constants for accessing coordinate pairs
LAT = 0
LON = 1



# Base class atmospheric model.
class AtmosModel:

   def __init__(self):
      self.BASE_DIRECTORY = "/model-output/"
      return

   # Display progress message on the console.
   # endString = set this to '' for no return
   def progress(self, message, endString='\n'):
      if (True):   # disable here in production
         print(message, end=endString)
         sys.stdout.flush()



   # return friendly name or short acronym for this model
   def getModelName(self):
      return("Atmos Model")

   def setBaseDirectory(self, newDir):
      self.BASE_DIRECTORY = newDir
      return

   # return friendly name for chemical species
   def getChemName(self, species):
      if (species.lower() == "o3"):
         return("Ozone")
      if (species.lower() == "co_fire"):
         return("CO fires")
      if (species.lower() == "extinctdn"):
         return("Aerosol extinction")
      if (species.lower() == "co"):
         return("Carbon monoxide")
      if (species.lower() == "co_bdry"):
         return("CO boundary")

      return("ChemName???")



   # return possible sub-directories and filename below BASE_DIRECTORY
   def getFilename(self, year, month, day, hour,
      minute=0):
      return("AtmosModel-{:4d}{:02d}{:02d}-{:02}000.nc"
         .format(year, month, day, hour))



   # return list of chemical thresholds for isosurfaces
   # species = name of chemical compound as model variable
   # return list of [[threshLow, threshMiddle, threshHigh], units]
   def getChemThresholds(self, species):
      thresholds = [12, 34, 56]
      units = "ppbv"

      if (species.upper() == "O3"):
         thresholds = [67, 70, 75]
      if (species.lower() == "co_fire"):
         thresholds = [200, 500, 1000]
      if (species.lower() == "co"):
         thresholds = [500, 700, 900]

      return([thresholds, units])



   # return the time step for this model in hours
   def getHourStride(self):
      return(12.0)



   # Read subset of model data over a rectangular portion of earth.
   # species = name of chemical or aerosol
   # modelFilename = path and name of some model output (NetCDF)
   # frameTime = UTC date and time to retrieve
   # latBox, lonBox = lat-lon bounds to retrieve
   # yStride = cell spacing to sample in latitude direction
   # xStride = cell spacing to sample in longitude direction
   # zStride = cell spacing to sample in vertical levels
   # return tuple of [values, lats, lons, heights, units, terrain]
   #		heights = meters above sea level
   def readModelBox(self, species, modelFilename, frameTime,
      latBox, lonBox, yStride, xStride, zStride):
      return([None, None, None, None, None, None])



   # Convert model units for display.
   # varValues = chemical concentrations, scaled on output
   # varUnits = unit string read from model output file
   # return newUnits for human display
   def convertDataUnits(self, varValues, varUnits):
      newVarUnits = varUnits

      if (varUnits.lower() == "mol/mol"):
         varValues *= 1e+9
         newVarUnits = "ppbv"

      if (varUnits.lower() == "ppmv"):
         varValues *= 1e+3
         newVarUnits = "ppbv"

      return(newVarUnits)



   # Find the highest surface chemical concentration over a time span.
   # species = looking at this chemical
   # startUTC, endUTC = looking across this time span
   # hourStep = number of hours between time steps
   # return tuple of [[latitude, longitude], max value found, at date time, units]
   def findMaxChemValue(self, species, startUTC, endUTC, hourStep):
      latLon = [40.0, -105.25]
      maxFound = 73.456

      timeSpan = endUTC - startUTC
      timeSpan /= 2
      maxDateTime = startUTC + timeSpan

      units = "ppbv"

      return([latLon, maxFound, maxDateTime, units])

