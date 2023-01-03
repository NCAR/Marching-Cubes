#!/usr/bin/env python3

# shapes3D.py
# Python module to build 3D shapes using module PyCollada.
#
# Author: Carl Drews, Atmospheric Chemistry Observations & Modeling
# Created: August 2021
# Copyright 2021 by the University Corporation for Atmospheric Research



import sys
import collada
import numpy
import math
import datetime
import os

import MarchingCubes
import utilsCollada
import acomKml

import pandas
import netCDF4



# Display progress message on the console.
# endString = set this to '' for no return
def progress(message, endString='\n'):
   if (True):   # disable here in production
      print(message, end=endString)
      sys.stdout.flush()



# Create 3D Collada model of an air parcel and save it to .dae file.
# parcelName = name for the model
#       parcelName becomes the filename, with spaces replaced
# return name of model file
def createAirParcel(parcelName, saveDir):
   filename = parcelName.replace(" ", "-") + ".dae"
   filepath = saveDir + filename
   progress("Writing air parcel file {}".format(filepath))

   # set up a pyramid on top of an inverted pyramid
   triangles = []
   xSpan = 10.0         # kilometers
   ySpan = 10.0
   zSpan = 1.0          # km

   if (True):
      # bogus - make parcels large enough to see
      xSpan *= 10
      ySpan *= 10
      zSpan *= 50

   # 3D shape origin will be at center of the solid body
   radiusX = xSpan / 2.0
   radiusY = ySpan / 2.0
   radiusZ = zSpan / 2.0

   # set up the pyramid vertices
   top = MarchingCubes.Vertex(0.0, 0.0, radiusZ)
   bottom = MarchingCubes.Vertex(0.0, 0.0, -radiusZ)

   southwest = MarchingCubes.Vertex(-radiusX, -radiusY, 0.0)
   northwest = MarchingCubes.Vertex(-radiusX, radiusY, 0.0)
   northeast = MarchingCubes.Vertex(radiusX, radiusY, 0.0)
   southeast = MarchingCubes.Vertex(radiusX, -radiusY, 0.0)
   
   # set up the normal vectors
   hypotenuse = math.hypot(xSpan, zSpan)
   xNormal = zSpan / hypotenuse
   zNormal = xSpan / hypotenuse
   #progress("Normal vector x,z = {},{}".format(xNormal, zNormal))
   southUp = MarchingCubes.Vertex(0, -xNormal, zNormal)
   eastUp = MarchingCubes.Vertex(xNormal, 0, zNormal)
   northUp = MarchingCubes.Vertex(0, xNormal, zNormal)
   westUp = MarchingCubes.Vertex(-xNormal, 0, zNormal)

   southDown = MarchingCubes.Vertex(0, -xNormal, -zNormal)
   eastDown = MarchingCubes.Vertex(xNormal, 0, -zNormal)
   northDown = MarchingCubes.Vertex(0, xNormal, -zNormal)
   westDown = MarchingCubes.Vertex(-xNormal, 0, -zNormal)

   # express the pyramids with 8 triangles
   triangles.append(MarchingCubes.Triangle(southwest, southeast, top)
      .setNormals(southUp, southUp, southUp))
   triangles.append(MarchingCubes.Triangle(southeast, northeast, top)
      .setNormals(eastUp, eastUp, eastUp))
   triangles.append(MarchingCubes.Triangle(northeast, northwest, top)
      .setNormals(northUp, northUp, northUp))
   triangles.append(MarchingCubes.Triangle(northwest, southwest, top)
      .setNormals(westUp, westUp, westUp))

   triangles.append(MarchingCubes.Triangle(southwest, southeast, bottom)
      .setNormals(southDown, southDown, southDown))
   triangles.append(MarchingCubes.Triangle(southeast, northeast, bottom)
      .setNormals(eastDown, eastDown, eastDown))
   triangles.append(MarchingCubes.Triangle(northeast, northwest, bottom)
      .setNormals(northDown, northDown, northDown))
   triangles.append(MarchingCubes.Triangle(northwest, southwest, bottom)
      .setNormals(westDown, westDown, westDown))

   orange = [1.0, 0.647, 0.0, 1.0]
   blue = [0.0, 0.0, 1.0, 1.0]
   daeColor = blue
   utilsCollada.writeDAEfile(triangles, daeColor, 0.0, filepath)

   return(filename)



# Create 3D Collada model of an airplane and save it to .dae file.
# planeName = name for the model
#       planeName becomes the filename, with spaces replaced
# sizeExaggeration = factor to make plane larger
# return name of model file
def createAirplane(planeName, saveDir,
   sizeExaggeration=200):
   filename = planeName.replace(" ", "-") + ".dae"
   filepath = saveDir + filename
   progress("Writing airplane file {}".format(filepath))

   # set up a very pointy pyramid facing positive Y
   triangles = []
   fuselageLength = 46 / 1000	# km long
   fuselageLength *= sizeExaggeration	# make it visible!

   # 3D shape origin will be at center of the solid body
   radiusY = fuselageLength / 2.0
   wingSpan = fuselageLength * 142.0 / 150.0
   fuselageRadius = fuselageLength * 0.05

   # set up the pyramid vertices
   front = MarchingCubes.Vertex(0.0, radiusY, 0.0)
   tail = MarchingCubes.Vertex(0.0, -radiusY, 0.0)

   tailLowerRight = MarchingCubes.Vertex(fuselageRadius, -radiusY, -fuselageRadius)
   tailUpperRight = MarchingCubes.Vertex(fuselageRadius, -radiusY, fuselageRadius)
   tailLowerLeft = MarchingCubes.Vertex(-fuselageRadius, -radiusY, -fuselageRadius)
   tailUpperLeft = MarchingCubes.Vertex(-fuselageRadius, -radiusY, fuselageRadius)

   # set up the normal vectors
   right = MarchingCubes.Vertex(1, 0, 0)
   left = MarchingCubes.Vertex(-1, 0, 0)
   up = MarchingCubes.Vertex(0, 0, 1)
   down = MarchingCubes.Vertex(0, 0, -1)

   # express the fuselage with 4 triangles
   # right
   triangles.append(MarchingCubes.Triangle(front, tailUpperRight, tailLowerRight)
      .setNormals(right, right, right))
   # top
   triangles.append(MarchingCubes.Triangle(front, tailUpperLeft, tailUpperRight)
      .setNormals(up, up, up))
   # left
   triangles.append(MarchingCubes.Triangle(front, tailLowerLeft, tailUpperLeft)
      .setNormals(left, left, left))
   # bottom
   triangles.append(MarchingCubes.Triangle(front, tailLowerRight, tailLowerLeft)
      .setNormals(down, down, down))

   # add the main wing
   wingRadius = wingSpan / 2
   wingWidth = wingRadius / 2
   wingFront = MarchingCubes.Vertex(0, 0, 0)
   rightWingtip = MarchingCubes.Vertex(wingRadius, -wingWidth, 0)
   leftWingtip = MarchingCubes.Vertex(-wingRadius, -wingWidth, 0)
   triangles.append(MarchingCubes.Triangle(wingFront, rightWingtip, leftWingtip)
      .setNormals(up, up, up))

   # add the horizontal stabilizer
   wingRadius /= 2
   wingWidth /= 2
   wingFront = MarchingCubes.Vertex(0, -radiusY + wingWidth, 0)
   rightWingtip = MarchingCubes.Vertex(wingRadius, -radiusY, 0)
   leftWingtip = MarchingCubes.Vertex(-wingRadius, -radiusY, 0)
   triangles.append(MarchingCubes.Triangle(wingFront, rightWingtip, leftWingtip)
      .setNormals(up, up, up))

   # Add two vertical stablizers facing left and right;
   # a single stabilizer will look strangely dark on the other side.
   topWingtip = MarchingCubes.Vertex(0, -radiusY, wingRadius)
   triangles.append(MarchingCubes.Triangle(wingFront, topWingtip, tail)
      .setNormals(right, right, right))

   offset = -wingWidth * 0.01
   wingFront = MarchingCubes.Vertex(offset, -radiusY + wingWidth, 0)
   tail = MarchingCubes.Vertex(offset, -radiusY, 0.0)
   topWingtip = MarchingCubes.Vertex(0, -radiusY, wingRadius)
   triangles.append(MarchingCubes.Triangle(wingFront, tail, topWingtip)
      .setNormals(left, left, left))

   # colors are RGBA (A = alpha opacity)
   orange = [1.0, 0.647, 0.0, 1.0]
   blue = [0.0, 0.0, 1.0, 1.0]
   white = [1.0, 1.0, 1.0, 1.0]
   cyan = [0.0, 1.0, 1.0, 1.0]
   yellow = [1.0, 1.0, 0.0, 1.0]

   # set airplane color here based on name
   daeColor = white
   if (planeName == "GV"):
      daeColor = cyan
   if (planeName == "WB57"):
      daeColor = yellow

   utilsCollada.writeDAEfile(triangles, daeColor, 0.0, filepath)

   return(filename)



# Read track and create as KML track along airplane flight path.
# planeModel = name of the .dae 3D model file of airplane
# trackFile = .ict file containing flight waypoints
# verticalExaggeration = should match terrain exaggeration in Google Earth
# return [track, path, tour] as KML elements
def createAirplaneTrack(planeName, planeModel, trackFile,
   verticalExaggeration=1.0):
   if (False):
      trackTimes = [datetime.datetime(2019,8,8, 0), datetime.datetime(2019,8,9, 4)]
      lats = [46.93, 48.93]
      lons = [-119.0, -116.0]
      altitudes = [900, 5000]	# meters above sea level

   progress("trackFile = {}".format(trackFile))

   # there are no standard variable names within NetCDF flight tracks
   latLonAlt = ["GGLAT", "GGLON", "GGALT"]	# NASA DC-8
   if (planeName == "GV"):
      latLonAlt = ["LATC", "LONC", "GGALT"]	# NSF/NCAR HIAPER Gulfstream 5
   if (planeName == "WB57"):
      latLonAlt = ["G_LAT_MMS", "G_LONG_MMS", "G_ALT_MMS"]	# NASA WB-57

   flightParts = None
   flightExtension = os.path.splitext(trackFile)[1].lower()
   if (flightExtension == ".ict"):
      # read the MetNav .ict file that recorded the plane flight
      flightParts = readMetNavFlight(trackFile, 1*60)	# 1-minute intervals
      #flightParts = readMetNavFlight(trackFile, 5*60)	# 5 minutes
   if (flightExtension == ".nc"):
      # read the NetCDF .nc file that recorded the plane flight
      flightParts = readNetCDFFlight(trackFile, 1*60, latLonAlt)	# 1-minute intervals

   trackTimes = flightParts[0]
   lats = flightParts[1]
   lons = flightParts[2]
   altitudes = numpy.array(flightParts[3], dtype="float")

   # Google Earth tracks and paths are already exaggerated.
   trackKML = acomKml.createModelMoving(
      planeName, "Flight of research aircraft.",
      planeModel, trackTimes,
      lats, lons, altitudes,
      holdFinal=False, idStr="airplane")

   pathKML = acomKml.createFixedPath(
      planeName + " flight path", "Flight path of research aircraft.",
      lats, lons, altitudes)

   # Apply terrain exaggeration to tour height.
   altitudes *= verticalExaggeration

   tourKML = acomKml.createFlightTour(
      planeName + " tour", "Tour along the flight path.",
      trackTimes, lats, lons, altitudes,
      duration=120)

   return([trackKML, pathKML, tourKML])



# Read MetNav .ict file and extract timestamp, latitude, longitude, altitude
# trackFileIct = full path and name of file
# timeStep = increment at which to extract (seconds)
# return [timestamps, lats, lons, heights]
def readMetNavFlight(trackFileIct, timeStep=5*60):
   progress("reading MetNav file {}".format(trackFileIct))
   navInfo = [None, None, None, None]

   # get information from header
   headerLines = 0
   startDate = datetime.datetime.utcnow()
   with open(trackFileIct, "r") as fp:
      # number of header lines
      oneLine = fp.readline()
      headerLines = int(oneLine.split(",")[0])
      progress("There are {} header lines.".format(headerLines))

      # skip 5 lines
      for li in range(0, 5):
         oneLine = fp.readline()

      # starting date
      oneLine = fp.readline()
      #progress("Dates = {}".format(oneLine))
      oneLineParts = oneLine.split(",")
      startDate = datetime.datetime(int(oneLineParts[0]),
         int(oneLineParts[1]), int(oneLineParts[2]))
      progress("Start date: {}".format(startDate))

   # re-open MetNav file as CSV and read specific columns
   csvContents = pandas.read_csv(trackFileIct, skiprows=headerLines-1)
   #progress("Latitude column = {}".format(csvContents["Latitude"]))

   # retrieve the raw columns for all time steps
   navTimesAll = csvContents["Time_Start"]
   navLatsAll = csvContents["Latitude"]
   navLonsAll = csvContents["Longitude"]
   navHeightsAll = csvContents["MSL_GPS_Altitude"]

   navTimes = []
   navLats = []
   navLons = []
   navHeights = []

   # begin with the starting location
   seconds = int(navTimesAll[0])
   startSeconds = seconds
   navTime = startDate + datetime.timedelta(seconds=seconds)
   navTimes.append(navTime)

   navLats.append(navLatsAll[0])
   navLons.append(navLonsAll[0])
   navHeights.append(navHeightsAll[0])

   # extract locations at even timeStep intervals
   stepIndex = int(startSeconds / timeStep)
   for navTime, navLat, navLon, navHeight in zip(
      navTimesAll, navLatsAll, navLonsAll, navHeightsAll):
      newStep = int(navTime / timeStep)
      if (newStep == stepIndex):
         continue

      # record this time step
      navDatetime = startDate + datetime.timedelta(seconds=navTime)
      navTimes.append(navDatetime)

      navLats.append(navLat)
      navLons.append(navLon)
      navHeights.append(navHeight)

      stepIndex = newStep

   # end with the final location
   seconds = int(navTimesAll[len(navTimesAll) - 1])
   navTime = startDate + datetime.timedelta(seconds=seconds)
   navTimes.append(navTime)

   navLats.append(navLatsAll[len(navLatsAll) - 1])
   navLons.append(navLonsAll[len(navLonsAll) - 1])
   navHeights.append(navHeightsAll[len(navHeightsAll) - 1])

   if (False):
      progress("navTimes: ")
      for navTime in navTimes:
         progress("\t{} ".format(navTime))

      progress("navLats = {}".format(navLats))
      progress("navHeights = {}".format(navHeights))

   navInfo = [navTimes, navLats, navLons, navHeights]
   return(navInfo)



# Read NetCDF .nc file and extract timestamp, latitude, longitude, altitude
# trackFileNc = full path and name of file
# timeStep = increment at which to extract (seconds)
# latLonAltitude = variables for [latitude, longitude, height]
# return [timestamps, lats, lons, heights]
def readNetCDFFlight(trackFileNc, timeStep=5*60,
   latLonAltitude=["GGLAT", "GGLON", "GGALT"]):
   progress("reading NetCDF file {}".format(trackFileNc))
   navInfo = [None, None, None, None]

   # open the NetCDF file for reading
   flightData = netCDF4.Dataset(trackFileNc)

   # date format here is 01/28/2014
   flightDateStr = flightData.getncattr("Flight_Date")
   #progress("flightDateStr = {}".format(flightDateStr))
   startDate = datetime.datetime.strptime(flightDateStr, "%m/%d/%Y")
   progress("Start date: {}".format(startDate))

   # retrieve the raw columns for all time steps
   navTimesAll = flightData.variables["Time_UTC"][:]
   navLatsAll = flightData.variables[latLonAltitude[0]][:]
   navLonsAll = flightData.variables[latLonAltitude[1]][:]
   navHeightsAll = flightData.variables[latLonAltitude[2]][:]

   flightData.close()

   progress("navTimesAll = {}".format(navTimesAll))
   progress("navLatsAll = {}".format(navLatsAll))
   progress("navLonsAll = {}".format(navLonsAll))
   progress("navHeightsAll = {}; {} to {}"
      .format(navHeightsAll, navHeightsAll.min(), navHeightsAll.max()))

   navTimes = []
   navLats = []
   navLons = []
   navHeights = []

   # begin with the starting location
   seconds = int(navTimesAll[0])
   startSeconds = seconds
   navTime = startDate + datetime.timedelta(seconds=seconds)
   navTimes.append(navTime)

   navLats.append(navLatsAll[0])
   navLons.append(navLonsAll[0])
   navHeights.append(navHeightsAll[0])

   # extract locations at even timeStep intervals
   stepIndex = int(startSeconds / timeStep)
   for navTime, navLat, navLon, navHeight in zip(
      navTimesAll, navLatsAll, navLonsAll, navHeightsAll):
      newStep = int(navTime / timeStep)
      if (newStep == stepIndex):
         continue

      # record this time step
      navDatetime = startDate + datetime.timedelta(seconds=int(navTime))
      navTimes.append(navDatetime)

      navLats.append(navLat)
      navLons.append(navLon)
      navHeights.append(navHeight)

      stepIndex = newStep

   # end with the final location
   seconds = int(navTimesAll[len(navTimesAll) - 1])
   navTime = startDate + datetime.timedelta(seconds=seconds)
   navTimes.append(navTime)

   navLats.append(navLatsAll[len(navLatsAll) - 1])
   navLons.append(navLonsAll[len(navLonsAll) - 1])
   navHeights.append(navHeightsAll[len(navHeightsAll) - 1])

   if (False):
      progress("navTimes: ")
      for navTime in navTimes:
         progress("\t{} ".format(navTime))

      progress("navLats = {}".format(navLats))
      progress("navHeights = {}".format(navHeights))

   navInfo = [navTimes, navLats, navLons, navHeights]
   return(navInfo)

