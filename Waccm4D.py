#!/usr/bin/env python3

# Waccm4D.py
# Python routines to read WACCM model output for 4-D visualization.
#
# Author: Carl Drews, Atmospheric Chemistry Observations & Modeling
# Created: March 2020
# Copyright 2020 by the University Corporation for Atmospheric Research



import Model4D

import datetime
import netCDF4
import numpy



class WaccmModel(Model4D.AtmosModel):
   
   def __init__(self):
      self.BASE_DIRECTORY = "/waccm-output/"
      return



   # return friendly name or short acronym for this model
   def getModelName(self):
      return("WACCM")

   def getFilename(self, year, month, day, hour):
      return("f.e22.beta02.FWSD.f09_f09_mg17.cesm2_2_beta02.forecast.002.cam.h3.{:4d}-{:02d}-{:02d}-00000.nc"
         .format(year, month, day))



   # Read subset of WRF-Chem data over a rectangular portion of earth.
   # species = name of chemical or aerosol
   # waccmFilename = path and name of WACCM model output (NetCDF)
   # frameTime = UTC date and time to retrieve
   # latBox, lonBox = lat-lon bounds to retrieve
   # yStride = cell spacing to sample in latitude direction
   # xStride = cell spacing to sample in longitude direction
   # zStride = cell spacing to sample in vertical levels
   # return tuple of [values, lats, lons, heights, units, terrain]
   #            heights = meters above sea level
   def readModelBox(self, species, waccmFilename, frameTime,
      latBox, lonBox, yStride, xStride, zStride):
   
      # open NetCDF file
      fp = netCDF4.Dataset(waccmFilename, 'r')
   
      # retrieve latitudes and longitudes
      latitude = fp.variables["lat"][:]
      longitude = fp.variables["lon"][:]
   
      # longitude range for WACCM is 0..360 degrees
      #self.progress("latitudes = {}".format(latitude))
      #self.progress("longitudes = {}".format(longitude))
   
      # find the indexes for the bounding box
      latBoxIndex = [0, len(latitude)-1]
      lonBoxIndex = [0, len(longitude)-1]
   
      # latitude
      for li in range(len(latitude)-1, -1, -1):
         if (latitude[li] <= latBox[0]):
            latBoxIndex[0] = li
            self.progress("Actual lat0 = {}".format(latitude[li]))
            break
      for li in range(latBoxIndex[0], len(latitude), yStride):
         if (latitude[li] >= latBox[1]):
            latBoxIndex[1] = li
            self.progress("Actual lat1 = {}".format(latitude[li]))
            break
      self.progress("latBoxIndex = {}".format(latBoxIndex))

      # longitude
      for li in range(len(longitude)-1, -1, -1):
         if (longitude[li] <= lonBox[0]):
            lonBoxIndex[0] = li
            self.progress("Actual lon0 = {}".format(longitude[li]))
            break
      for li in range(lonBoxIndex[0], len(longitude), xStride):
         if (longitude[li] >= lonBox[1]):
            lonBoxIndex[1] = li
            self.progress("Actual lon1 = {}".format(longitude[li]))
            break
      self.progress("lonBoxIndex = {}".format(lonBoxIndex))

      # retrieve the actual latitudes and longitudes
      lats = latitude[latBoxIndex[0]:latBoxIndex[1]+1:yStride]
      lons = longitude[lonBoxIndex[0]:lonBoxIndex[1]+1:xStride]
      #self.progress("Actual latitudes = {}".format(lats))
      #self.progress("Actual longitudes = {}".format(lons))
   
      # find the requested frame time
      dateIndex = 0
      date = fp.variables["date"][:]
      datesec = fp.variables["datesec"][:]
      for di in range(len(date)-1, -1, -1):
         dateStr = "{}".format(date[di])
   
         daySeconds = datesec[di]
         hour = int(daySeconds / 3600)
         daySeconds -= hour * 3600
         minute = int(daySeconds / 60)
         daySeconds -= minute * 60
   
         myDateTime = datetime.datetime(int(dateStr[0:4]), int(dateStr[4:6]),
            int(dateStr[6:8]), hour=hour, minute=minute, second=daySeconds)
         if (myDateTime <= frameTime):
            dateIndex = di
            self.progress("Actual date = {}".format(myDateTime))
            break
   
      # WACCM geopotential height is listed from top to bottom of atmosphere
      geoHeight = fp.variables["Z3"][dateIndex, :,
         latBoxIndex[0]:latBoxIndex[1]+1:yStride,
         lonBoxIndex[0]:lonBoxIndex[1]+1:xStride]
      geoHeight = geoHeight[::-1, :]       # reverse along vertical axis
   
      # altitude stride always begins at the surface,
      heights = geoHeight[::zStride, :]
      #self.progress("heights = {}".format(heights[:,0,0]))
      terrain = heights[0, :, :]
      #self.progress("terrain = {}".format(terrain[0:5, 0:5]))
   
      # finally retrieve the chemical values
      values = fp.variables[species][dateIndex, :,
         latBoxIndex[0]:latBoxIndex[1]+1:yStride,
         lonBoxIndex[0]:lonBoxIndex[1]+1:xStride]
      values = values[::-1, :]     # reverse along vertical axis
      values = values[::zStride, :]
   
      rawUnits = getattr(fp.variables[species], "units")
      #self.progress("{} raw values = {} {}".format(species, values[:,0,0], rawUnits))
   
      units = self.convertDataUnits(values, rawUnits)
      self.progress("units = {}".format(units))
      #self.progress("{} converted values = {} {}".format(species, values[:,0,0], units))
   
      # close NetCDF file
      fp.close()

      # expand the coordinates into two dimensions
      #self.progress("lats = {} {}".format(lats.shape, lats))
      #self.progress("lons = {} {}".format(lons.shape, lons))
      tempLats = numpy.zeros([lats.shape[0], lons.shape[0]])
      for loni in range(0, lons.shape[0]):
         tempLats[:, loni] = lats[:]
      #self.progress("tempLats = {} {}".format(tempLats.shape, tempLats))

      tempLons = numpy.zeros([lats.shape[0], lons.shape[0]])
      for lati in range(0, lats.shape[0]):
         tempLons[lati, :] = lons[:]
      #self.progress("tempLons = {} {}".format(tempLons.shape, tempLons))

      lats = tempLats
      lons = tempLons

      return([values, lats, lons, heights, units, terrain])
   
