#!/usr/bin/env python3

# WrfChem4D.py
# Python routines to read WRF-Chem model output for 4-D visualization.
#
# Author: Carl Drews, Atmospheric Chemistry Observations & Modeling
# Created: March 2020
# Copyright 2020 by the University Corporation for Atmospheric Research



import Model4D

import datetime
import netCDF4
import wrf
import numpy
import math



class WrfChemModel(Model4D.AtmosModel):
   
   def __init__(self):
      self.BASE_DIRECTORY = "/WRF-Chem/"
      return

   # return friendly name or short acronym for this model
   def getModelName(self):
      return("WRF-Chem")

   # return the time step for this model in hours
   def getHourStride(self):
      return(1.0)



   def getFilename(self, year, month, day, hour,
      minute=0):
      subDirectory = "{:4d}{:02d}{:02d}/wrf/".format(year, month, day)
      nameOnly = ("wrfout_hourly_d01_{:4d}-{:02d}-{:02d}_{:02d}:00:00"
         .format(year, month, day, hour))

      # apply suffix for regional extractions
      if (False):
         suffix = "_SouthCA.nc"		# Los Angeles NASA/JPL
         nameOnly += suffix

      # if requested date is a forecast, look in today's sub-directory instead
      requestedDate = datetime.datetime(year=year, month=month, day=day, hour=hour)
      todayNow = datetime.datetime.utcnow()
      #todayNow += datetime.timedelta(hours=-24)	# bogus - re-run today after 6:00 PM
      if (requestedDate > todayNow):
         #self.progress("\tForecast requested.")
         subDirectory = ("{:4d}{:02d}{:02d}/wrf/"
            .format(todayNow.year, todayNow.month, todayNow.day))

      return(subDirectory + nameOnly)



   # Search *inward* toward one corner of the bounding box.
   # Call this routine four times to get all four corners,
   # then probably take the maximum of those limits.
   # latBounds = requested south and north latitude extremes
   # lonBounds = requested west and east longitude extremes
   # lats, lons = coordinates of all grid cells, indexed by [lati, loni]
   # latStride, lonStride = number of cells to jump between samples
   # latSearchDirection, lonSearchDirection = look in this direction (-1 or 1)
   # return array of [latIndex, latIndex]
   def findOneBoxCorner(self, latBounds, lonBounds, lats, lons,
      latStride, lonStride,
      latSearchDirection, lonSearchDirection):

      # default to the southwest
      latBoxIndex = [0]
      lonBoxIndex = [0]

      # start searching at one corner of the entire model domain
      testLatIndex = 0
      searchForLat = latBounds[0]
      if (latSearchDirection < 0):
         testLatIndex = lats.shape[0] - 1
         searchForLat = latBounds[1]

      testLonIndex = 0
      searchForLon = lonBounds[0]
      if (lonSearchDirection < 0):
         testLonIndex = lons.shape[1] - 1
         searchForLon = lonBounds[1]

      # search inward toward the model cell that contains the corner
      while (True):
         # calculate the test distance
         testDistance = math.hypot(searchForLat - lats[testLatIndex, testLonIndex],
            searchForLon - lons[testLatIndex, testLonIndex])
         #self.progress("testDistance1 = {} at {}, {}"
         #   .format(testDistance, testLatIndex, testLonIndex))

         # look in the requested direction
         latDistance = math.hypot(
            searchForLat - lats[testLatIndex + latSearchDirection * latStride, testLonIndex],
            searchForLon - lons[testLatIndex + latSearchDirection * latStride, testLonIndex])
         lonDistance = math.hypot(
            searchForLat - lats[testLatIndex, testLonIndex + lonSearchDirection * lonStride],
            searchForLon - lons[testLatIndex, testLonIndex + lonSearchDirection * lonStride])

         # have we reached the closest point we can get?
         moveDistance = min(latDistance, lonDistance)
         if (moveDistance > testDistance):
            break

         # step toward the biggest improvement
         if (latDistance < lonDistance):
            testLatIndex += latSearchDirection * latStride
         else:
            testLonIndex += lonSearchDirection * lonStride

      #self.progress("testDistance2 = {} at {}, {}"
      #   .format(testDistance, testLatIndex, testLonIndex))

      # back off if we overshot,
      # using search direction to reverse the comparison operator
      if (lats[testLatIndex, testLonIndex] * latSearchDirection
         > searchForLat * latSearchDirection):
         testLatIndex -= latSearchDirection * latStride

         # stay inside the model domain
         if (testLatIndex >= lats.shape[0]):
            testLatIndex = lats.shape[0] - 1
         if (testLatIndex < 0):
            testLatIndex = 0

      if (lons[testLatIndex, testLonIndex] * lonSearchDirection
         > searchForLon * lonSearchDirection):
         testLonIndex -= lonSearchDirection * lonStride

         # stay inside the model domain
         if (testLonIndex >= lons.shape[1]):
            testLonIndex = lats.shape[1] - 1
         if (testLonIndex < 0):
            testLonIndex = 0

      #self.progress("testDistance3 = {} at {}, {}"
      #   .format(testDistance, testLatIndex, testLonIndex))

      # record this corner
      latBoxIndex = testLatIndex
      lonBoxIndex = testLonIndex
      #self.progress("Actual corner = {} North and {} East."
      #   .format(lats[latBoxIndex, lonBoxIndex], lons[latBoxIndex, lonBoxIndex]))

      return([latBoxIndex, lonBoxIndex])



   # Locate and return the smallest indexes of a box that spans
   # all four corners of the requested lat-lon bounds.
   # Because of Lambert Conformal rotation, this box
   # may not be oriented strictly north-south and east-west.
   # The algorithm employs a search for any map projection,
   # and is necessary for any projection not on a lat-lon grid.
   # latBounds = requested south and north latitude extremes
   # lonBounds = requested west and east longitude extremes
   # lats, lons = coordinates of all grid cells, indexed by [lati, loni]
   # latitudeStride, longitudeStride = number of cells to jump between samples
   # return tuple of [[lat0index, lat1index], [lon0index, lon1index]]
   def findBoxCorners(self, latBounds, lonBounds, lats, lons,
      latitudeStride, longitudeStride):

      # look for all four corners
      southWest = self.findOneBoxCorner(latBounds, lonBounds, lats, lons,
         latitudeStride, longitudeStride, 1, 1)
      #self.progress("southWest = {}".format(southWest))

      northWest = self.findOneBoxCorner(latBounds, lonBounds, lats, lons,
         latitudeStride, longitudeStride, -1, 1)
      #self.progress("northWest = {}".format(northWest))

      northEast = self.findOneBoxCorner(latBounds, lonBounds, lats, lons,
         latitudeStride, longitudeStride, -1, -1)
      #self.progress("northEast = {}".format(northEast))

      southEast = self.findOneBoxCorner(latBounds, lonBounds, lats, lons,
         latitudeStride, longitudeStride, 1, -1)
      #self.progress("southEast = {}".format(southEast))

      # determine the outermost corners
      lat0Index = min(southWest[0], southEast[0])
      lat1Index = max(northWest[0], northEast[0])
      lon0Index = min(southWest[1], northWest[1])
      lon1Index = max(southEast[1], northEast[1])
      indexBounds = [[lat0Index, lat1Index], [lon0Index, lon1Index]]
      #self.progress("Final index bounds to read are lat {} and lon {}"
      #   .format(indexBounds[0], indexBounds[1]))

      # display the actual lat-lon corners
      #self.progress("Final box corners are = \n\tSW ({},{}) \n\tNW ({},{}) \n\tNE ({},{}) \n\tSE ({},{})."
      #   .format(lats[lat0Index, lon0Index], lons[lat0Index, lon0Index],
      #      lats[lat1Index, lon0Index], lons[lat1Index, lon0Index],
      #      lats[lat1Index, lon1Index], lons[lat1Index, lon1Index],
      #      lats[lat0Index, lon1Index], lons[lat0Index, lon1Index]))

      return(indexBounds)



   # Retrieve units separately because wrf module doesn't get them.
   # species = name of chemical species
   # return units attribute of that variable (not converted)
   def getUnits(self, wrfFilename, species):
      fp = netCDF4.Dataset(wrfFilename, 'r', format='nc')
      rawUnits = getattr(fp.variables[species], "units", None)
      fp.close()

      return(rawUnits)



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

      # open NetCDF file for reading
      wrfData = netCDF4.Dataset(waccmFilename)

      # find the requested frame time
      dateIndex = 0
      dateTimes = wrf.getvar(wrfData, "Times", squeeze=False).data		# squeeze doesn't seem to work here
      #self.progress("dateTimes = {} {}".format(dateTimes, type(dateTimes)))
      if (len(dateTimes.shape) == 0):
         # force single date-time into an array
         dateTimesArray = numpy.array([dateTimes])
      else:
         dateTimesArray = dateTimes
      dateTimes = None

      # convert numpy datetime64 to datetime.datetime
      unixEpoch = numpy.datetime64(0, 's')
      oneSecond = numpy.timedelta64(1, 's')
      for di in range(len(dateTimesArray)-1, -1, -1):
         secondsSinceEpoch = (dateTimesArray[di] - unixEpoch) / oneSecond
         myDateTime = datetime.datetime.utcfromtimestamp(secondsSinceEpoch)
         if (myDateTime <= frameTime):
            dateIndex = di
            break

      self.progress("Actual date = {} at file index {}"
         .format(dateTimesArray[dateIndex], dateIndex))

      # retrieve latitudes and longitudes
      lats = wrf.getvar(wrfData, "lat", timeidx=dateIndex).data
      lons = wrf.getvar(wrfData, "lon", timeidx=dateIndex).data
      #self.progress("{} shape1 = {}".format("latitude", lats.shape))
      #self.progress("{} shape2 = {}".format("longitude", lons.shape))

      # find the indexes for the bounding box
      corners = self.findBoxCorners(latBox, lonBox, lats, lons, yStride, xStride)
      latBoxIndex = corners[0]
      lonBoxIndex = corners[1]

      # retrieve the model terrain height
      terrain = wrf.getvar(wrfData, "ter", timeidx=dateIndex, units="m").data
      #self.progress("terrain shape1 = {}".format(terrain.shape))
      terrain = terrain[
         latBoxIndex[0]:latBoxIndex[1]+1:yStride,
         lonBoxIndex[0]:lonBoxIndex[1]+1:xStride]
      #self.progress("terrain shape2 = {}".format(terrain.shape))
      #self.progress("terrain = {}".format(terrain[0:5, 0:5]))

      # retrieve model height for mass grid
      heights = wrf.getvar(wrfData, "height", units="m").data
      #self.progress("heights1 shape = {}".format(heights.shape))
      #self.progress("heights1 = {}".format(heights[:, 0, 0]))	# water point in the Pacific Ocean

      heights = heights[::zStride,
         latBoxIndex[0]:latBoxIndex[1]+1:yStride,
         lonBoxIndex[0]:lonBoxIndex[1]+1:xStride]
      #self.progress("heights2 shape = {}".format(heights.shape))
      #self.progress("heights2 = {}".format(heights[0, 0:5, 0:5]))

      # finally retrieve the chemical values
      values = wrf.getvar(wrfData, species, timeidx=dateIndex).data
      #self.progress("{} shape3 = {}".format(species, values.shape))
      values = values[::zStride,
         latBoxIndex[0]:latBoxIndex[1]+1:yStride,
         lonBoxIndex[0]:lonBoxIndex[1]+1:xStride]
      #self.progress("{} shape4 = {}".format(species, values.shape))

      # close NetCDF file
      wrfData.close()

      # get the units separately
      rawUnits = self.getUnits(waccmFilename, species)
      #self.progress("{} raw values = {} {}".format(species, values[:,0,0], rawUnits))

      units = self.convertDataUnits(values, rawUnits)
      #self.progress("{} converted values = {} {}".format(species, values[:,0,0], units))
      #self.progress("{} shape5 = {}".format(species, values.shape))
      #self.progress("values = {} {}".format(values[0, 0:5, 0:5], units))
 
      # subset the coordinates  
      lats = lats[
         latBoxIndex[0]:latBoxIndex[1]+1:yStride,
         lonBoxIndex[0]:lonBoxIndex[1]+1:xStride]
      lons = lons[
         latBoxIndex[0]:latBoxIndex[1]+1:yStride,
         lonBoxIndex[0]:lonBoxIndex[1]+1:xStride]

      #self.progress("{} shape6 = {}".format(species, values.shape))
      return([values, lats, lons, heights, units, terrain])
 


   # Find the highest surface chemical concentration over a time span.
   # species = looking at this chemical
   # startUTC, endUTC = looking across this time span
   # hourStep = number of hours between time steps
   # return tuple of [[latitude, longitude], max value found, at date time, units]
   def findMaxChemValue(self, species, startUTC, endUTC, hourStep):
      self.progress("findMaxChemValue():")

      # initialize the return values
      latLon = [-999.9, -999.9]
      maxFound = -999.9
      maxDateTime = datetime.datetime(1, 1, 1)

      maxFile = None
      units = None
      terrainMask = None

      # step through the date-time interval at even hours
      frameTime = datetime.datetime(startUTC.year,
         startUTC.month, startUTC.day, startUTC.hour)
      timeStep = datetime.timedelta(hours=hourStep)
      endTime = datetime.datetime(endUTC.year,
         endUTC.month, endUTC.day, endUTC.hour)

      while (frameTime <= endTime):
         self.progress("\tFrame time {}".format(frameTime))

         # formulate the filename
         frameFile = self.BASE_DIRECTORY + self.getFilename(
            frameTime.year, frameTime.month, frameTime.day, frameTime.hour)
         self.progress("\tFrame file {}".format(frameFile))

         # open NetCDF file for reading
         wrfData = netCDF4.Dataset(frameFile)

         dateIndex = 0
         lats = wrf.getvar(wrfData, "lat", timeidx=dateIndex).data
         lons = wrf.getvar(wrfData, "lon", timeidx=dateIndex).data

         if (terrainMask is None):
            # mask  out the oceans and seas
            terrainHeight = wrf.getvar(wrfData, "HGT", timeidx=dateIndex).data
            terrainMask = numpy.zeros(terrainHeight.shape, dtype=bool)
            for y in range(terrainMask.shape[0]):
               for x in range(terrainMask.shape[1]):
                  if (terrainHeight[y, x] <= 0.01):	# 1 centimeter above mean sea level
                     terrainMask[y, x] = True		# True means ignore

         # retrieve the numeric values
         values = wrf.getvar(wrfData, species, timeidx=dateIndex).data
         levelIndex = 0
         values = values[levelIndex, :, :]

         # close NetCDF file
         wrfData.close()

         if (units is None):
            units = self.getUnits(frameFile, species)

         # only examine chemical concentrations over land
         values = numpy.ma.masked_array(values, mask=terrainMask)

         # Use Python to extract max value over entire array.
         frameMax = values.max()
         #self.progress("\tframeMax = {}".format(frameMax))
         if (frameMax > maxFound):
            maxFound = frameMax
            maxDateTime = frameTime

            # now locate the lat-lon of the individual cell
            # This case occurs for about 1/10 of the frames,
            # so re-locating the new lat-lon is not inefficient.
            cellMax = -999.9
            for y in range(values.shape[0]):
               for x in range(values.shape[1]):
                  cellValue = values[y, x]
                  if (cellValue > cellMax):
                     cellMax = cellValue
                     latLon[Model4D.LAT]= lats[y, x]
                     latLon[Model4D.LON]= lons[y, x]

                     #self.progress("\tterrain height = {} mask = {} value = {}"
                     #   .format(terrainHeight[y, x], terrainMask[y, x], cellMax))
                     
            self.progress("\tNew max {} {} found at {} on {}"
               .format(maxFound, units, latLon, frameTime))

         frameTime += timeStep

      self.progress("\tMaximum value {} {} occurred at {}"
         .format(maxFound, units, maxDateTime))

      return([latLon, maxFound, maxDateTime, units])

