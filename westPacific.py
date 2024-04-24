#!/usr/bin/env python3

# westPacific.py
# Python script to create 4-D Google Earth animation of air parcels
# and their trajectories over the Western Pacific. This script
# demonstrates research by Ren Smith, ACOM/NCAR.
#
# Author: Carl Drews, Atmospheric Chemistry Observations & Modeling
# Created: April 2021
# Copyright 2021 by the University Corporation for Atmospheric Research



import sys
import datetime
import xml.etree.ElementTree as ET
import os
import math

import utilsLite

import acomKml
import Model4D
import Waccm4D
import WrfChem4D

import colorbar
import subprocess
import copy

import MarchingCubes
import utilsCollada
import shapes3D

import netCDF4
import contextlib
import pandas
import shutil



# Google Earth -> Tools -> Options -> 3D View -> Terrain -> Elevation Exaggeration
terrainExaggeration = 3.0	# set within Google Earth
heightExaggeration = 100	# for this demonstration
heightExaggeration = 30		# for ACCLIP



# Display progress message on the console.
# endString = set this to '' for no return
def progress(message, endString='\n'):
   if (True):   # disable here in production
      print(message, end=endString)
      sys.stdout.flush()



# Make sure that longitude is in range(-180, 180) for Google Earth.
def wrapLongitude(rawLongitude):
   if (rawLongitude <= 180.0):
      return(rawLongitude)

   return(rawLongitude - 360.0)



# return the average of a pair of numbers
def average(numPair):
   return((numPair[0] + numPair[1]) / 2.0)



# Create placemarks within one time-step of the animation.
# startTime, endTime = time bounds of the frame
# northFolder,... = KML folders for each side of the plot
#                 each containing many time steps
# latBounds, lonBounds = extent of the domain, including plot margins
# altBounds = extent of plot above sea level (km), including margins
# imageIndex = the images are organized by integer
# stageDirectory = where to write the .DAE Collada files created
# northLatBounds northLonBounds = extent of the north wall
# westLatBounds, westLonBounds = extent of the west wall
# return [wallCount, overlayCount] created here
def createOneFrame(startTime, endTime,
   northFolder, westFolder, surfaceFolder,
   latBounds, lonBounds, altBounds, imageIndex, stageDirectory,
   northLatBounds, northLonBounds, westLatBounds, westLonBounds):

   wallCount = 0
   overlayCount = 0

   imageFile = "images/ST_XY/{}.png".format(imageIndex)
   progress("imageFile = {}".format(imageFile))

   # create surface overlay
   overlayKML = acomKml.createImageOverlay(imageFile,
      latBounds, lonBounds, -1.0, startTime, endTime)
   surfaceFolder.append(overlayKML)
   overlayCount += 1

   # create northern wall
   imageFile = "images/ST_XZ/{}.png".format(imageIndex)
   daeFilename = ("NorthWall-{:%Y%m%d-%H%M}.dae"
      .format(startTime))
   utilsCollada.writeImageWallDAE(imageFile,
      northLatBounds, northLonBounds,
      [altBounds[0] * heightExaggeration, altBounds[1] * heightExaggeration],
      stageDirectory + daeFilename)
   wallKML = acomKml.createImageWall(daeFilename,
      average(northLatBounds), average(northLonBounds),
      0.0, startTime, endTime)
   northFolder.append(wallKML)
   wallCount += 1

   # create western wall
   imageFile = "images/ST_YZ/{}.png".format(imageIndex)
   daeFilename = ("WestWall-{:%Y%m%d-%H%M}.dae"
      .format(startTime))
   utilsCollada.writeImageWallDAE(imageFile,
      westLatBounds, westLonBounds,
      [altBounds[0] * heightExaggeration, altBounds[1] * heightExaggeration],
      stageDirectory + daeFilename)
   wallKML = acomKml.createImageWall(daeFilename,
      average(westLatBounds), average(westLonBounds),
      0.0, startTime, endTime)
   westFolder.append(wallKML)
   wallCount += 1

   return([wallCount, overlayCount])



# Create initial view of domain.
# latBox, lonBox = domain span
# return KML element for document
def viewCenter(latBox, lonBox):
   centerLat = (latBox[0] + latBox[1]) / 2.0
   centerLon = (lonBox[0] + lonBox[1]) / 2.0

   latRange = (latBox[1] - latBox[0]) / 2.0
   latRange *= 111.32 * 1000.0		# degrees to meters

   lookAt = acomKml.kmlElement("LookAt")
   lookAt.append(acomKml.kmlElement("latitude", "{}".format(centerLat)))
   lookAt.append(acomKml.kmlElement("longitude", "{}".format(centerLon)))
   lookAt.append(acomKml.kmlElement("altitude", "{}".format(0.0)))
   lookAt.append(acomKml.kmlElement("range", "{}".format(latRange)))
   lookAt.append(acomKml.kmlElement("heading", "{}".format(0.0)))
   lookAt.append(acomKml.kmlElement("tilt", "{}".format(70.0)))
   lookAt.append(acomKml.kmlElement("altitudeMode", "absolute"))

   return(lookAt)



# Read parcel tracks from NetCDF file.
# return list of parcel [indexs, times, latitudes, longitudes, altitudes (km)]
def readTracks(filepath):
   progress("Reading parcel tracks from {}".format(filepath))

   # open the trajectory file
   fp = netCDF4.Dataset(filepath, 'r')

   # retrieve the variables therein
   indexes = fp.variables["Parcel"][:]
   times = fp.variables["Time"][:]
   latitudes = fp.variables["Latitude"][:]
   longitudes = fp.variables["Longitude"][:]
   altitudes = fp.variables["Altitude"][:]

   # close the trajectory file
   fp.close()

   # process the variables a bit
   indexes = indexes.astype(int)

   # convert date-time stamps to Python datetime
   dateTimes = []
   for time in times:
      #progress("datetime = {:.0f} ".format(time), endString='')
      pass

   # TODO: The date-time stamps have to be strings, not floats.
   # Carl Drews - May 9, 2021
   startDate = datetime.datetime(2014, 1, 25, 12)
   endDate = datetime.datetime(2014, 1, 29, 7)
   hourStride = 3
   timeDelta = datetime.timedelta(days=0, seconds = hourStride * 3600)
   frame = endDate
   while (frame >= startDate):
      dateTimes.append(frame)
      frame -= timeDelta

   return([indexes, dateTimes, latitudes, longitudes, altitudes])



# Read parcel trajectories from NetCDF file.
# return list of parcel [indexs, times, latitudes, longitudes, altitudes (km)]
def readTrajectories(filepath):
   progress("Reading parcel trajectories from {}".format(filepath))

   # open the trajectory file
   fp = netCDF4.Dataset(filepath, 'r')

   # retrieve the variables therein
   numParticles = fp.dimensions["Particle"].size
   progress("numParticles = {}".format(numParticles))
   numTimes = fp.dimensions["Time"].size
   progress("numTimes = {}".format(numTimes))

   julianDays = fp.variables["Julian_day"][:]
   seconds = fp.variables["Seconds"][:]

   latitudes = fp.variables["Latitude"][:].transpose()
   longitudes = fp.variables["Longitude"][:].transpose()
   altitudes = fp.variables["Altitude"][:].transpose()

   # close the trajectory file
   fp.close()

   # process the variables a bit
   indexes = range(numParticles)

   # adjust trajectory times according to the ACCLIP flight
   timeAdjust = -datetime.timedelta(hours=12)

   # convert date-time stamps to Python datetime
   dateTimes = []
   secondsPerDay = 24.0 * 3600
   for day,sec in zip(julianDays, seconds):
      daySec = day + sec / secondsPerDay
      frame = pandas.to_datetime(daySec, origin="julian", unit="D")
      frame += timeAdjust
      dateTimes.append(frame)

   # convert pressure in hPa to height in kilometers
   # https://en.wikipedia.org/wiki/Scale_height
   R = 8.31446261815324 	# gas constant kg⋅m2⋅s−2⋅K−1⋅mol−1
   T = 250			# mean temperature in degrees K for Earth
   M = 0.029			# molar mass kg/mol for Earth
   g = 9.80665			# acceleration of gravity m/s2
   P0 = 1010			# surface pressure hPa
   scaleHeight = (R * T) / (M  * g) / 1000
   progress("scaleHeight = {} km".format(scaleHeight))

   for pi in range(numParticles):
      for ti in range(numTimes):
         alt = altitudes[pi][ti]
         altitudes[pi][ti] = -scaleHeight * math.log(alt / P0)

   # add small adjustment to match airplane altitude
   altitudes *= 1.02937

   return([indexes, dateTimes, latitudes, longitudes, altitudes])



# Create flying tour of the bounding box.
# myFolder = place the tour within here
# tourTitle = short descriptive name of the tour, like Tour01
# startWhen, endWhen = date-time of flying time
# latBox, lonBox = box around which to fly, probably contains slabs
# runMode = "manual" or "auto"
def createFlight(myFolder, tourTitle, startWhen, endWhen, latBox, lonBox,
   runMode):

   progress("createFlight {} from {} to {}"
      .format(tourTitle, startWhen, endWhen))

   metersPerKm = 1000.0	# kilometers to meters
   kmPerDegree = 111.32

   # calculate center of box, maybe use for LookAt target
   latCenter = (latBox[0] + latBox[1]) / 2.0
   lonCenter = (lonBox[0] + lonBox[1]) / 2.0
   lookAtPoint = [latCenter, lonCenter]

   # calculate range from which to view central point
   latRange = abs(lookAtPoint[Model4D.LAT] - latBox[0])
   lonRange = abs(lookAtPoint[Model4D.LON] - lonBox[0])
   flightRange = math.hypot(latRange, lonRange) * kmPerDegree
   progress("flightRange = {} km".format(flightRange))

   # time spent flying each leg
   angleStart = 45.0
   numLegs = 12
   legTime = (endWhen - startWhen) / numLegs
   progress("legTime = {}".format(legTime))
   legAngle = 360.0 / numLegs

   tiltAngle = [75.0]
   if (runMode == "auto"):
      # fly in a circle around the domain
      lookAtPoint = [latCenter, lonCenter]
      latRange = abs(lookAtPoint[Model4D.LAT] - latBox[0])
      lonRange = abs(lookAtPoint[Model4D.LON] - lonBox[0])
      flightRange = math.hypot(latRange, lonRange) * kmPerDegree
      progress("flightRange 2 = {} km".format(flightRange))
      angleStart = 45.0
      tiltAngle = [60.0, 90.0, 75.0]

   flightRange *= metersPerKm

   # set up the flying tour
   tour = acomKml.kmlElement("gx:Tour")
   myFolder.append(tour)

   tourName = acomKml.kmlElement("name",
      "{0}-{1.year}{1.month:02}{1.day:02}to{2.year}{2.month:02}{2.day:02}"
      .format(tourTitle, startWhen, endWhen))
   tour.append(tourName)

   playlist = acomKml.kmlElement("gx:Playlist")
   tour.append(playlist)

   # begin flight at southwest corner of bounding box
   flyTo = acomKml.kmlElement("gx:FlyTo")
   playlist.append(flyTo)
   flyTo.append(acomKml.kmlElement("gx:duration", "5.0"))
   flyTo.append(acomKml.kmlElement("gx:flyToMode", "bounce"))

   lookAt = acomKml.kmlElement("LookAt")
   lookAt.append(acomKml.kmlElement("latitude", "{}".format(lookAtPoint[Model4D.LAT])))
   lookAt.append(acomKml.kmlElement("longitude", "{}".format(lookAtPoint[Model4D.LON])))
   lookAt.append(acomKml.kmlElement("altitude", "{}".format(0.0)))
   lookAt.append(acomKml.kmlElement("range", "{}".format(flightRange)))
   lookAt.append(acomKml.kmlElement("heading", "{}".format(angleStart)))
   lookAt.append(acomKml.kmlElement("tilt", "{}".format(tiltAngle[0])))
   lookAt.append(acomKml.kmlElement("altitudeMode", "absolute"))
   takeoffLanding = None

   flyTo.append(lookAt)

   when = acomKml.kmlElement("gx:TimeStamp")
   lookAt.append(when)
   when.append(acomKml.kmlElement("when",
      startWhen.strftime("%Y-%m-%dT%H:%M:%SZ")))

   # fly the legs of the tour
   direction = 1	# clockwise
   if (startWhen.day % 2 == 0):
      direction = -1	# counter-clockwise

   for legi in range(1, numLegs+1):
      legDone = startWhen + legi * legTime
      progress("legDone = {}".format(legDone))
      legHeading = (angleStart + legAngle * legi * direction) % 360
      legTilt = tiltAngle[0]

      # fly to the next waypoint
      # tour displays a single moving time frame
      playlist.append(acomKml.createFlyToLookAt(
         duration=4, arriveWhen=legDone,
         latitude=lookAtPoint[Model4D.LAT], longitude=lookAtPoint[Model4D.LON],
         altitude=0, flightRange=flightRange,
         heading=legHeading, tilt=legTilt))

      if (tourTitle == "Tour02"):
         # tour changes when parcels start ascending
         riseDate = datetime.datetime(2022, 8, 18, 15)
         if (legDone + legTime > riseDate):
            # tilt view back to show higher parcels
            playlist.append(acomKml.createFlyToLookAt(
               duration=4, arriveWhen=datetime.datetime(2022, 8, 18, 21),
               latitude=35.5, longitude=121.5,
               altitude=320000, flightRange=1250801,
               heading=15, tilt=80))

            # fly to the sample rendezvous
            sampleDate = datetime.datetime(2022, 8, 19, 2, 40)
            playlist.append(acomKml.createFlyToLookAt(
               duration=4, arriveWhen=sampleDate,
               latitude=34.328, longitude=124.365,
               altitude=450000, flightRange=125080,
               heading=30, tilt=60))

            # wait here while samples taken
            progress("Waiting at {}".format(sampleDate))
            waitHere = acomKml.kmlElement("gx:Wait")
            waitHere.append(acomKml.kmlElement("gx:duration", "4.0"))	# seconds in real time
            playlist.append(waitHere)

            # follow the airplane back to base
            playlist.append(acomKml.createFlyToLookAt(
               duration=8, arriveWhen=datetime.datetime(2022, 8, 19, 4),
               latitude=37.1, longitude=127.1,
               altitude=0, flightRange=1250801,
               heading=0, tilt=40))

            playlist.append(acomKml.createFlyToLookAt(
               duration=8, arriveWhen=datetime.datetime(2022, 8, 19, 5, 40),
               latitude=37.1, longitude=127.1,
               altitude=0, flightRange=1250801,
               heading=0, tilt=40))

            playlist.append(acomKml.createFlyToLookAt(
               duration=4, arriveWhen=datetime.datetime(2022, 8, 19, 6, 8),
               latitude=37.1, longitude=127.1,
               altitude=0, flightRange=125080,
               heading=0, tilt=40))

            break

   return



# Create hovering tour of the bounding box.
# We slowly move from lower left of south wall to upper right.
# myFolder = place the tour within here
# startWhen, endWhen = date-time of flying time
# latBox, lonBox, altBox = box around which to fly
# latCenter, lonCenter = center of the action
def createHover01(myFolder, startWhen, endWhen,
   latBox, lonBox, altBox,
   latCenter, lonCenter):

   metersPerKm = 1000.0	# kilometers to meters
   kmPerDegree = 111.32

   # calculate center of box, maybe use for LookAt target
   altCenter = (altBox[0] + altBox[1]) / 2.0 * heightExaggeration
   lookAtPoint = [latCenter, lonCenter, altCenter]

   # calculate range from which to view central point
   latRange = abs(lookAtPoint[Model4D.LAT] - latBox[0])
   lonRange = abs(lookAtPoint[Model4D.LON] - lonBox[0])
   flightRange = math.hypot(latRange, lonRange) * kmPerDegree
   progress("flightRange = {} km".format(flightRange))
   flightRange *= metersPerKm

   # set up the flying tour
   tour = acomKml.kmlElement("gx:Tour")
   myFolder.append(tour)

   tourName = acomKml.kmlElement("name",
      "Hover01-{0.year}{0.month:02}{0.day:02}to{1.year}{1.month:02}{1.day:02}"
      .format(startWhen, endWhen))
   tour.append(tourName)

   playlist = acomKml.kmlElement("gx:Playlist")
   tour.append(playlist)

   # begin flight at southwest lower corner of bounding box
   flyTo = acomKml.kmlElement("gx:FlyTo")
   playlist.append(flyTo)
   flyTo.append(acomKml.kmlElement("gx:duration", "5.0"))
   flyTo.append(acomKml.kmlElement("gx:flyToMode", "bounce"))

   lookAt = acomKml.populateLookAt(
      lookAtPoint[Model4D.LAT], lookAtPoint[Model4D.LON], lookAtPoint[2],
      flightRange, 45, 120)	# view up into the sky

   when = acomKml.kmlElement("gx:TimeStamp")
   when.append(acomKml.kmlElement("when",
      startWhen.strftime("%Y-%m-%dT%H:%M:%SZ")))

   lookAt.append(when)
   flyTo.append(lookAt)

   # end flight at southeast upper corner of bounding box
   flyTo = acomKml.kmlElement("gx:FlyTo")
   playlist.append(flyTo)
   flyTo.append(acomKml.kmlElement("gx:duration", "20.0"))
   flyTo.append(acomKml.kmlElement("gx:flyToMode", "smooth"))

   lookAt = acomKml.populateLookAt(
      lookAtPoint[Model4D.LAT], lookAtPoint[Model4D.LON], lookAtPoint[2],
      flightRange, -45, 50)	# look down toward earth

   when = acomKml.kmlElement("gx:TimeStamp")
   when.append(acomKml.kmlElement("when",
      endWhen.strftime("%Y-%m-%dT%H:%M:%SZ")))

   lookAt.append(when)
   flyTo.append(lookAt)

   return



# Create second hovering tour of the bounding box.
# We stay at the upper southeast corner, perhaps with minimumal motion.
# myFolder = place the tour within here
# startWhen, endWhen = date-time of flying time
# latBox, lonBox, altBox = box around which to fly
# latCenter, lonCenter = center of the action
def createHover02(myFolder, startWhen, endWhen,
   latBox, lonBox, altBox,
   latCenter, lonCenter):

   metersPerKm = 1000.0	# kilometers to meters
   kmPerDegree = 111.32

   # calculate center of box, maybe use for LookAt target
   altCenter = (altBox[0] + altBox[1]) / 2.0 * heightExaggeration
   lookAtPoint = [latCenter, lonCenter, altCenter]

   # calculate range from which to view central point
   latRange = abs(lookAtPoint[Model4D.LAT] - latBox[0])
   lonRange = abs(lookAtPoint[Model4D.LON] - lonBox[0])
   flightRange = math.hypot(latRange, lonRange) * kmPerDegree
   flightRange *= 1.2
   progress("flightRange 02 = {} km".format(flightRange))
   flightRange *= metersPerKm

   # set up the flying tour
   tour = acomKml.kmlElement("gx:Tour")
   myFolder.append(tour)

   tourName = acomKml.kmlElement("name",
      "Hover02-{0.year}{0.month:02}{0.day:02}to{1.year}{1.month:02}{1.day:02}"
      .format(startWhen, endWhen))
   tour.append(tourName)

   playlist = acomKml.kmlElement("gx:Playlist")
   tour.append(playlist)

   # hover at southeast upper corner of bounding box
   lookAngle = 45	# 180 = straight up, 0 = straight down
   flyTo = acomKml.kmlElement("gx:FlyTo")
   playlist.append(flyTo)
   flyTo.append(acomKml.kmlElement("gx:duration", "5.0"))
   #flyTo.append(acomKml.kmlElement("gx:flyToMode", "bounce"))
   flyTo.append(acomKml.kmlElement("gx:flyToMode", "smooth"))

   lookAt = acomKml.populateLookAt(
      lookAtPoint[Model4D.LAT], lookAtPoint[Model4D.LON], lookAtPoint[2],
      flightRange, -45, lookAngle)	# look down toward earth

   when = acomKml.kmlElement("gx:TimeStamp")
   when.append(acomKml.kmlElement("when",
      startWhen.strftime("%Y-%m-%dT%H:%M:%SZ")))

   lookAt.append(when)
   flyTo.append(lookAt)

   # end flight at same place
   flyTo = acomKml.kmlElement("gx:FlyTo")
   playlist.append(flyTo)
   flyTo.append(acomKml.kmlElement("gx:duration", "20.0"))
   flyTo.append(acomKml.kmlElement("gx:flyToMode", "smooth"))

   lookAt = acomKml.populateLookAt(
      lookAtPoint[Model4D.LAT], lookAtPoint[Model4D.LON], lookAtPoint[2],
      flightRange, -45, lookAngle)	# look down toward earth

   when = acomKml.kmlElement("gx:TimeStamp")
   lookAt.append(when)
   when.append(acomKml.kmlElement("when",
      endWhen.strftime("%Y-%m-%dT%H:%M:%SZ")))

   flyTo.append(lookAt)

   return



# Return a safe list of comma-separated command-line arguments.
# csvParams = with no spaces: 123.4,abc,today
# numeric = convert params to floating-point
# integer = convert params to integer
# return [123.4, "abc", "today"]
def safeParams(csvParams, numeric=False, integer=False):
   allParams = csvParams.split(",")
   params = []
   for param in allParams:
      if (integer):
         params.append(utilsLite.safeInt(utilsLite.sanitize(param)))
      elif (numeric):
         params.append(utilsLite.safeFloat(utilsLite.sanitize(param)))
      else:
         params.append(utilsLite.sanitize(param))

   return(params)



# Safely move to new directory and guarantee return.
@contextlib.contextmanager
def pushd(new_dir):
   previous_dir = os.getcwd()
   os.chdir(new_dir)
   try:
      yield
   finally:
      os.chdir(previous_dir)



STAGE_DIRECTORY = "RenSmith/stage/"

# Main program begins here.
def main():
   # set up the default region to plot

   # default to Western Pacific north of New Guinea
   south = -20.0
   north = 30.0
   west = 120.0
   east = 180.0
   # include the plot margins here
   latBounds = [south - 7.1, north + 7.8]	# degrees north
   lonBounds = [west - 14.2, east + 12.0]	# degrees east
   heightBounds = [0 - 2.530, 18 + 2.822]	# km above sea level

   # wall plots have different margins
   northLats = [north, north]			# degrees north
   northLons = [west - 9.8, east + 7.8]	# degrees east
   westLats = [south - 8.2, north + 6.4]
   westLons = [west, west]

   if (False):
      # Western Pacific
      startDate = datetime.datetime(2014, 1, 25, 12)
      endDate = datetime.datetime(2014, 1, 29, 3)
   if (True):
      # East Asia
      startDate = datetime.datetime(2022, 7, 20, 0)
      endDate = datetime.datetime(2022, 8, 19, 6, 8)
      latBounds = [29.5, 41.5]	# degrees north
      lonBounds = [112.0, 131]	# degrees east

   dataDir = "RenSmith/"

   # default to no airplane track
   airplane = None
   airplaneName = None
   airplaneTrackFile = None

   # retrieve the command-line arguments, if any
   for argPair in sys.argv:
      progress("argPair = " + argPair)
      # the arguments are: arg=value
      pairValue = argPair.split('=')
      if (len(pairValue) < 2):
         continue

      if (pairValue[0].lower() == "datadir"):
         dataDir = pairValue[1]

      if (pairValue[0].lower() == "airplane"):
         airplane = pairValue[1].split(",")
         airplaneName = airplane[0]
         airplaneTrackFile = airplane[1]

   # Display the command-line arguments just received.
   progress("dataDir = {}".format(dataDir))
   progress("Date range is {} to {}".format(startDate, endDate))

   if (airplane is not None):
      progress("airplane = {} flying track {}"
         .format(airplaneName, airplaneTrackFile))

   progress("terrainExaggeration = {}".format(terrainExaggeration))
   progress("heightExaggeration = {}".format(heightExaggeration))

   # get time step for this model
   hourStride = 3
   progress("hourStride = {}".format(hourStride))
   timeDelta = datetime.timedelta(days=0, seconds = hourStride * 3600)

   region = "WesternPacific"
   region = "EastAsia"
   # create reference name like this: o3-20210418to20210420
   kmlBasename = ("{}-{:04d}{:02d}{:02d}to{:04d}{:02d}{:02d}"
      .format(region,
      startDate.year, startDate.month, startDate.day,
      endDate.year, endDate.month, endDate.day))

   # create kml sub-directory in the staging area
   stageDir = STAGE_DIRECTORY + kmlBasename + "/"
   if (not os.path.isdir(stageDir)):
      os.mkdir(stageDir)

   # create kml structure for this chemical
   kmlRoot = ET.Element("kml")
   kmlRoot.set("xmlns", "http://www.opengis.net/kml/2.2")
   kmlRoot.set("xmlns:gx", "http://www.google.com/kml/ext/2.2")
   kmlTree = ET.ElementTree(element=kmlRoot)

   # document and its description
   kmlDoc = ET.SubElement(kmlRoot, "Document")

   kmlDescription = ET.Element("description")
   kmlDescription.text = (
      "{} animation for region {}.\n"
            .format("Parcel trajectory", region)
         + "Dates: {} to {} UTC.\n".format(startDate, endDate)
         + "Created by westPacific.py on {} local time.\n"
            .format(datetime.datetime.now())
         + "Authors: Carl Drews and Ren Smith, ACOM/NCAR/UCAR.\n"
      )
   kmlDoc.append(kmlDescription)

   # create a sensible initial view
   kmlInitialView = viewCenter(latBounds, lonBounds)

   # Look at the entire time span so that all images get loaded.
   extraHours = 3
   allTime = acomKml.timeSpan(startDate,
      endDate + datetime.timedelta(hours=extraHours),
      prefix="gx:")
   kmlInitialView.append(allTime)

   kmlDoc.append(kmlInitialView)

   # create folder for the plot walls
   northFolder = acomKml.folder("North wall")
   kmlDoc.append(northFolder)
   westFolder = acomKml.folder("West wall")
   kmlDoc.append(westFolder)

   # create folder for the surface layer
   bottomFolder = acomKml.folder("Surface overlay")
   kmlDoc.append(bottomFolder)

   # create folders for the parcel models
   parcelFolder = acomKml.folder("Parcel trajectories")
   kmlDoc.append(parcelFolder)

   parcelPathFolder = acomKml.folder("Parcel paths")
   kmlDoc.append(parcelPathFolder)
   visibility = acomKml.kmlElement("visibility", "0")
   parcelPathFolder.append(visibility)

   # calculate and manage the backward-running image indexe
   numFrames = int((endDate - startDate) / timeDelta) + 1	# this really works!
   progress("numFrames = {}".format(numFrames))
   plotIndex = numFrames - 1		# zero-based index

   # loop through dates and times
   frame = startDate
   walls = 0
   overlays = 0

   while (False and (frame <= endDate)):
      progress("\nFrame time: {}".format(frame))
      frameInfo = createOneFrame(frame, frame + timeDelta,
         northFolder, westFolder, bottomFolder,
         latBounds, lonBounds, heightBounds,
         plotIndex, stageDir,
         northLats, northLons, westLats, westLons)

      walls += frameInfo[0]
      overlays += frameInfo[1]

      # move on to the next time frame
      frame += timeDelta
      plotIndex -= 1

      #break		# bogus	- just create the first frame

   progress("Created {} walls and {} surface overlays."
      .format(walls, overlays))

   # create 3D model of a single orange parcel
   parcelFile = shapes3D.createAirParcel("Air Parcel", stageDir)
   progress("parcelFile = {}".format(parcelFile))

   # don't show the pushpin
   dontShowKML = acomKml.createDontShowStyle()
   kmlDoc.append(dontShowKML)
   dontShowKML = acomKml.createDontShowStyle(idStr="airplane", lineWidth=5)
   kmlDoc.append(dontShowKML)

   # read the trajectories of the parcels
   if (False):
      trackFile = dataDir + "RF07_ERA5_GoogleEarth.nc"
      progress("trackFile = {}".format(trackFile))
      trackInfo = readTracks(trackFile)

   if (True):
      trackFile = dataDir + "traj3d_GV_20220819T024000Z_20220720T024000Z.nc"
      progress("trackFile = {}".format(trackFile))
      trackInfo = readTrajectories(trackFile)

   trackIndexes = trackInfo[0]
   trackTimes = trackInfo[1]
   trackLats = trackInfo[2]
   trackLons = trackInfo[3]
   trackAlts = trackInfo[4]

   #progress("trackIndexes = {}".format(trackIndexes))
   #progress("trackTimes = {} {}".format(trackTimes, len(trackTimes)))
   #progress("trackAlts = {} {}".format(trackAlts[42], len(trackAlts)))

   # adjust the altitudes according to view settings
   if (True):
      trackAlts *= (1000.0 * heightExaggeration / terrainExaggeration)
   if (False):
      trackAlts *= 1000.0

   # Google Earth has to have time running forward
   trackTimes.reverse()

   # animate those parcels along their tracks
   for parcelIndex in trackIndexes:
      lats = trackLats[parcelIndex]
      lons = trackLons[parcelIndex]
      altitudes = trackAlts[parcelIndex]

      # reverse the time-based arrays
      lats = lats[::-1]
      lons = lons[::-1]
      altitudes = altitudes[::-1]

      parcelName = "Parcel {}".format(parcelIndex)
      parcelDescription = "Animated air parcel."

      # create one parcel track
      trackKML = acomKml.createModelMoving(
         parcelName, parcelDescription,
         parcelFile, trackTimes,
         lats, lons, altitudes)
      parcelFolder.append(trackKML)

      # create one fixed parcel path
      pathKML = acomKml.createFixedPath(
         parcelName, parcelDescription,
         lats, lons, altitudes, "fixedPath02")
      parcelPathFolder.append(pathKML)

   # create a folder for the tour around the smoke
   flightFolder = acomKml.folder("Tours")
   kmlDoc.append(flightFolder)

   # end flights after endDate so all frames are fully seen
   endDatePlus = endDate + datetime.timedelta(hours=extraHours)
   createFlight(flightFolder, "Tour01", startDate, endDatePlus,
      latBounds, lonBounds, "auto")
   createFlight(flightFolder, "Tour02",
      endDate + datetime.timedelta(days=-4), endDate,
      latBounds, lonBounds, "auto")
   createFlight(flightFolder, "Tour03",
      endDate + datetime.timedelta(days=-4), endDate,
      latBounds, lonBounds, "manual")

   # the hover tours focus on a central point
   focusLat = -5
   focusLon = 155
   if (True):
      # ACCLIP Yellow Sea
      focusLat = 37.34
      focusLon = 122.48

   createHover01(flightFolder, startDate, endDatePlus,
      latBounds, lonBounds, heightBounds, focusLat, focusLon)
   createHover02(flightFolder, startDate, endDatePlus,
      latBounds, lonBounds, heightBounds, focusLat, focusLon)

   # create research airplane flying along track
   if (airplane is not None):
      airplaneFolder = acomKml.folder("Research flights")
      kmlDoc.append(airplaneFolder)

      planeFile = shapes3D.createAirplane(airplaneName, stageDir,
         sizeExaggeration=800)
      progress("planeFile = {}".format(planeFile))

      # set up path thickness and color
      fixedPathKML = acomKml.createFixedPathStyle("fixedPath01", "ccffff00")	# cyan airplane
      kmlDoc.append(fixedPathKML)
      fixedPathKML = acomKml.createFixedPathStyle("fixedPath02", "4000ff00")	# green parcels
      kmlDoc.append(fixedPathKML)

      # Carl Drews and Ren Smith decided on a height adjustment
      # of +6% applied to the GV aircraft altitude. This adjustment
      # matches the airplane height to the air parcels sampled.
      # The source of the height error is probably a conversion
      # from pressure coordinates to meters above sea level.
      # September 9, 2021
      heightAdjust = 1.06
      if (True):
         # new correction moved into parcel trajectories for ACCLIP; May 4, 2023
         heightAdjust = 1.0

      # Create airplane track and fixed path,
      # with exaggerated height to match air parcels.
      airplaneKMLs = shapes3D.createAirplaneTrack(
         airplaneName, planeFile, airplaneTrackFile,
         heightAdjust * heightExaggeration / terrainExaggeration)
      airplaneFolder.append(airplaneKMLs[0])
      airplaneFolder.append(airplaneKMLs[1])

      # create customized hover tour
      if (True):
         altX = heightExaggeration / terrainExaggeration
         hoverStart = datetime.datetime(2022, 8, 17, 0, 0)
         hoverMiddle = datetime.datetime(2022, 8, 19, 2, 40)
         hoverEnd = datetime.datetime(2022, 8, 19, 6, 8)
         acomKml.createHoverTour(flightFolder,
            [hoverStart, hoverMiddle, hoverMiddle, hoverMiddle, hoverMiddle, hoverEnd],	# time bounds
            [focusLat, focusLat, 34.38, 34.38, 35.00, 34.65],			# latitude
            [focusLon, focusLon, 124.32, 124.32, 124.32, 126.47],		# longitude
            [12000 * altX, 12000 * altX, 450 * 1000, 450 * 1000, 0, 0],		# altitude
            [1500*1000, 1500*1000, 500*1000, 100*1000, 1300*1000, 1500*1000],	# range
            [75, 75, 30, 0, 0, 0],					# heading
            [87, 87, 60, 0, 0, 0],					# tilt
            [2, 25, 3, 3, 3, 8],					# seconds duration
            [0, 2, 0, 3, 0, 0])					# seconds wait after

   # create a folder for rulers that show altitude and other measurements
   if (True):
      imageFilename = "3000metersAltitude.png"
      rulerFolder = acomKml.folder("Rulers")
      kmlDoc.append(rulerFolder)
      acomKml.createRulers(rulerFolder, latBounds, lonBounds,
         imageFilename, terrainExaggeration)
      shutil.copy(imageFilename, stageDir)

   # write the formatted kml file
   acomKml.indent(kmlRoot)
   kmlName = kmlBasename + ".kml"
   progress("Writing to KML file {}{}".format(stageDir, kmlName))
   kmlTree.write(stageDir + kmlName)

   # move to the stage directory for relative file paths
   with pushd(stageDir):
      # link to the images because ../../ does not work within zip
      imageLink = "images"
      if (not os.path.exists(imageLink)):
         cmd = "ln --symbolic --verbose ../../ {}".format(imageLink)
         progress("{}".format(cmd))
         returnCode = subprocess.call(cmd, shell=True)

      # create KMZ compressed archive
      kmzName = kmlBasename + ".kmz"
      imagePattern = "images/ST_*/*.png"
      cmd = ("zip -FS {} {} *.dae *.png {}"
         .format(kmzName, kmlName, imagePattern))
      progress("{}".format(cmd))
      returnCode = subprocess.call(cmd, shell=True)

   # deposit KMZ in the data directory
   cmd = ("mv {}{} {}".format(stageDir, kmzName, dataDir))
   progress("{}".format(cmd))
   returnCode = subprocess.call(cmd, shell=True)

   return(0)	# no error



# call the main
progress("{}".format(__file__))
progress("Start time: {} Local".format(datetime.datetime.now()))
retValue = main()
progress("End time: {} Local".format(datetime.datetime.now()))

sys.exit(retValue)

