#!/usr/bin/env python3

# bushFires.py
# Python script to create 4-D Google Earth animation of Australian bush fires.
#
# Author: Carl Drews, Atmospheric Chemistry Observations & Modeling
# Created: January 2020
# Copyright 2020 by the University Corporation for Atmospheric Research



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
import CamChem4D

import colorbar
import subprocess
import copy

import MarchingCubes
import utilsCollada
import shapes3D

import shutil



# Google Earth -> Tools -> Options -> 3D View -> Terrain -> Elevation Exaggeration
terrainExaggeration = 3.0	# set within Google Earth
#terrainExaggeration = 10.0	# Google Earth maximum exaggeration is 3
heightExaggeration = 10.0	# set here for altitude of model slabs above terrain

metersPerKm = 1000.0	# kilometers to meters



# Display progress message on the console.
# endString = set this to '' for no return
def progress(message, endString='\n'):
   if (True):   # disable here in production
      print(message, end=endString)
      sys.stdout.flush()



# Determine if this slab should be rendered in 3-D.
# The inclusion algorithm may look at adjacent cells
# if we are developing a surface of constant concentration.
# varValues = 3-D grid of values at this frame time
# heightIndex = vertical location within grid
# latIndex, lonIndex = horizontal location within grid
# threshold = looking for concentrations above this value
# return True if this grid cell should be displayed
def includeSlab(varValues,
   heightIndex, latIndex, lonIndex,
   threshold):

   # count corners on both sides of threshold
   numLess = 0
   numGreater = 0

   # check values at the eight corners
   for z in range(2):
      for y in range(2):
         for x in range(2):
            cellValue = varValues[heightIndex + z, latIndex + y, lonIndex + x]
            if (cellValue < threshold):
               numLess += 1
            else:
               numGreater += 1

   if (False and numLess > 0 and numGreater > 0):
      progress("threshold value = {}".format(threshold))
      progress("corners are {}".format(varValues[heightIndex:heightIndex+2,
         latIndex:latIndex+2, lonIndex:lonIndex+2]))
      progress("numLess = {}   numGreater = {}".format(numLess, numGreater))

   return(numLess > 0 and numGreater > 0)



# index constants for the atmospheric layers
SURFACE = 0
TROPOSPHERE = 1
STRATOSPHERE = 2
MESOSPHERE = 3
THERMOSPHERE = 4

# Place new surface placemark into the correct layer of the atmosphere.
# surfaceMark = newly created KML placemark representing 3-D WACCM surface
# layerIndex = which layer; with 0=surface, 1=troposphere, etc.
# myLayerFolders = slabMark is appended to one of these...
# subFolderIndex = ...within this threshold/color sub-folder
def surfaceIntoLayer(surfaceMark, layerIndex, myLayerFolders,
   subFolderIndex):

   targetLayer = myLayerFolders[layerIndex]
   subFolders = targetLayer.findall("Folder")
   subFolders[subFolderIndex].append(surfaceMark)



# Calculate atmospheric layer based on height.
# zIndex = index into the WACCM height array
# height = bottom of layer in kilometers
# return one of the index constants defined above
def heightToLayer(zIndex, height):
   targetLayer = None

   if (zIndex == 0):
      targetLayer = SURFACE
   elif (height < 15.0):
      targetLayer = TROPOSPHERE
   elif (height < 50.0):
      targetLayer = STRATOSPHERE
   elif (height < 85.0):
      targetLayer = MESOSPHERE
   else:
      targetLayer = THERMOSPHERE

   return(targetLayer)



# Return short abbreviation for atmospheric layer.
def layerAbbrev(layerIndex):
   layerDict = {
      0: "surface",
      1: "trop",
      2: "strat",
      3: "meso",
      4: "thermo"
   }

   return(layerDict[layerIndex])



# Calculate index of value within threshold table.
# myValue = some floating-point value to look up
# myDivisions = minimum values for that color
# return index into the table
def valueToThresholdIndex(myValue, myDivisions):

   myIndex = -1
   for division in myDivisions:
      if (myValue >= division):
         myIndex += 1

   return(myIndex)



# Return the RGBA color and transparency of a chemical threshold.
# isoIndex = index into the chemical thresholds defined elsewhere,
#		probably in AtmosModel.getChemThresholds(self, species)
# return list of [[RBGA], transparent]
def isoColor(isoIndex):
   colors = [
      [0.0, 1.0, 0.0, 1.0],	# green
      [1.0, 1.0, 0.0, 1.0],	# yellow
      [1.0, 0.0, 0.0, 1.0]	# red
   ]

   # In theory we could specify these opacity values above,
   # in the Alpha channel. The Collada file format controls
   # the overall transparency separately, so we follow that
   # convention. Carl Drews - April 23, 2021
   transparency = [0.6, 0.4, 0.2]	# green, yellow, red

   return([colors[isoIndex], transparency[isoIndex]])



# Create safe units string for this chemical value.
# modelUnits = string from the original WACCM or WRF-Chem output
# chemValue = chemical value expressed in those units
def isoValueUnits(modelUnits, chemValue):
   safeUnits = modelUnits.replace("/", "_per_")		# safe for filenames
   isoValueStr = "{:.0f}".format(chemValue)
   if (chemValue < 1.0):
      isoValueStr = "{:.1e}".format(chemValue)

   return((safeUnits, isoValueStr))



# Create placemarks within one time-step of the animation.
# startTime, endTime = time bounds of the frame
# layerFolders = KML folders for each layer of the atmosphere,
#                 each containing many time steps
# cellCount = count of grid cells that cross the isosurface
# triangleCount = count of triangles created so far
# lineCount = count of multi-segment grid/contour lines created so far
# markCount = count of placemarks created so far
# myWaccmData = tuple of [values, lats, lons, heights]
#		values = chemical concentration in 3-D space
#		lats, lons = vector of horizontal grid locations
#		heights = vertical grid locations in 3-D space
# chemDivisions = minimum values for green, yellow, and red
# minPlotHeight = create 3D surfaces beginning at this height above sea level (meters)
# maxPlotHeight = create slabs up to this height above sea level (meters)
# stageDirectory = where to write the .DAE Collada files created
# species = looking at this chemical compound
# gridLines = create 3D grid/contour lines in a cage around the surface
# return updated [cellCount, triangleCount, lineCount, markCount]
def createOneFrame(startTime, endTime, layerFolders, cellCount,
   triangleCount, lineCount, markCount, myWaccmData, chemDivisions,
   minPlotHeight, maxPlotHeight,
   stageDirectory, species,
   gridLines=True):

   # collect the WACCM data
   chemValues = myWaccmData[0]
   if (False):
      # force a surface layer at 3.0 for testing of terrain following
      chemValues[0, :] = 2.5
      chemValues[1, :] = 3.5
   #progress("chemValues shape = {}".format(chemValues.shape))
   lats = myWaccmData[1]
   lons = myWaccmData[2]
   #progress("lats {}   lons {}".format(lats.shape, lons.shape))
   heights = myWaccmData[3]
   #progress("heights shape = {}".format(heights.shape))
   units = myWaccmData[4]
   ground = myWaccmData[5]

   progress("minPlotHeight = {} meters".format(minPlotHeight))
   progress("maxPlotHeight = {} meters".format(maxPlotHeight))

   # adjust heights (m) for Google Earth (km)
   meterToKm = 0.001
   heights *= meterToKm
   minPlotHeight *= meterToKm
   maxPlotHeight *= meterToKm

   # create a quick and easy way to access height for this z-level
   centerHeights = heights[:, int(heights.shape[1] / 2), int(heights.shape[2] / 2)]
   progress("centerHeights = {} km".format(centerHeights))

   # all 3D models have origin at center of lat-lon bounding box
   centerLat = (lats[0, 0] + lats[-1, -1]) / 2.0
   centerLon = (lons[0, 0] + lons[-1, -1]) / 2.0
   progress("3D model origin lat, lon = ({}, {})".format(centerLat, centerLon))
   centerCoords = MarchingCubes.Vertex(centerLon, centerLat, 0.0)

   if (True):
      # find the cell with the greatest value
      levIndex = 0
      latIndex = 0
      lonIndex = 0
      maxValue = chemValues[levIndex, latIndex, lonIndex]
      for z in range(0, chemValues.shape[0]):
         if (maxPlotHeight <= 0.0 and z > 0):
            # surface plot only
            continue

         # check if level within altitude range
         if (minPlotHeight >= 0.0 and centerHeights[z] < minPlotHeight):
            continue
         if (maxPlotHeight > 0.0 and centerHeights[z] > maxPlotHeight):
            continue

         for y in range(0, chemValues.shape[1]):
            for x in range(0, chemValues.shape[2]):
               oneValue = chemValues[z, y, x]
               if (oneValue > maxValue):
                  maxValue = oneValue
                  levIndex = z
                  latIndex = y
                  lonIndex = x

      progress("Maximum value {} {} occurs at [{}, {}, {}]".
         format(maxValue, units, levIndex, latIndex, lonIndex))

   # loop through the chemical thresholds for each surface
   for isoIndex, isoValue in enumerate(chemDivisions):
      #if (isoIndex != 1):
      #   continue		# bogus
      if (False and isoIndex == 0):
         # bogus - skip the lowest threshold
         continue

      # derive one isosurface, separated into layers of the atmosphere
      surfaceInfo = createOneSurface(chemValues, lats, lons, heights,
         units, ground, isoValue, minPlotHeight, maxPlotHeight,
         layerFolders, isoIndex, gridLines)
      cellCount += surfaceInfo[0]
      triangleCount += surfaceInfo[1]
      lineCount += surfaceInfo[2]
      layers = surfaceInfo[3]
      lineList = surfaceInfo[4]

      # write the DAE Collada files now
      yearMonthDay = startTime.strftime("%Y%m%d")
      hourMinute = startTime.strftime("%H%M")
      for atmosLayer, triList in enumerate(layers):
         if (len(triList) == 0):
            continue

         # Apply vertical exaggeration to match Google Earth terrain,
         # and convert absolute lat-lon to relative km.
         for triangle in triList:
            triangle.exaggerateHeight(terrainExaggeration)
            # Curve the triangle surfaces down from the origin.
            triangle.latLonToKm(centerCoords)

            # calculate normal vectors for each triangle
            # The normals have already been set to point outside the plume
            # by crawlOneEdgeCrossing() called within createOneSurface().
            # The normal direction was set by the winding direction
            # when we determined if the polygon should be reversed.
            triNormal = MarchingCubes.normalVector(triangle)
            triangle.setNormals(triNormal, triNormal, triNormal)
            #if (triNormal.z <= 0):
            #   progress("Warning: Suspicious triNormal vector points down.")
            #   progress("triangle in km = {}\n{}".format(triangle, triangle.toString()))

         #progress("First triangle in km = {}\n{}".format(triList[0], triList[0].toString()))

         # filename goes into placemark, so save it
         safeUnits, isoValueStr = isoValueUnits(units, isoValue)

         daeNameOnly = ("{}-{}-{}-{}{}-{}.dae"
            .format(species, yearMonthDay, hourMinute, isoValueStr, safeUnits,
            layerAbbrev(atmosLayer)))
         daePathName = stageDirectory + daeNameOnly
         daeColor = isoColor(isoIndex)
         utilsCollada.writeDAEfile(triList,
            daeColor[0], daeColor[1],
            daePathName)

         # create the corresponding KML placemark for that surface
         daeNameOnly = os.path.splitext(daeNameOnly)[0]
         surfaceMark = acomKml.renderSurfaceKML(daeNameOnly, markCount,
            centerLat, centerLon, 0.0,
            isoValue, startTime, endTime, chemDivisions)
         surfaceIntoLayer(surfaceMark, atmosLayer, layerFolders, isoIndex)
         markCount += 1

         # create the corresponding placemark for grid lines over that surface
         # The KML and DAE cases here are mutually exclusive because
         # the lines are modified during vertical exaggeration.
         if (gridLines and False):
            # Implement the wireframe as KML <LineString> elements.

            # Apply km-to-meters to match Google Earth meter units for KML,
            # and make adjustment for any exaggeration > Google Earth 3.
            gridExaggeration = metersPerKm * (terrainExaggeration / 3.0)

            # loop through the multi-segment lines
            for oneLine in lineList[atmosLayer]:
               for oneVertex in oneLine:
                  oneVertex.z *= gridExaggeration

            # create a lot of KML LineStrings
            gridNameOnly = daeNameOnly + "-grid"
            gridMark = acomKml.renderGridKML(gridNameOnly, markCount,
               lineList[atmosLayer], startTime, endTime)
            surfaceIntoLayer(gridMark, atmosLayer, layerFolders, isoIndex)
            markCount += 1

         if (gridLines and True):
            # Implement the wireframe as a 3D COLLADA model saved to a .DAE file;
            # the many grid lines will be expressed in a <lines> element.

            # loop through the multi-segment lines
            for oneLine in lineList[atmosLayer]:
               for oneVertex in oneLine:
                  # apply vertical exaggeration to the Z-coordinates
                  oneVertex.z *= terrainExaggeration

                  # Curve the vertex coordinates down from the origin.
                  oneVertex.latLonToKm(centerCoords)

            # filename goes into placemark, so save lines in COLLADA file
            safeUnits, isoValueStr = isoValueUnits(units, isoValue)

            daeNameOnly = ("{}-{}-{}-{}{}-{}-grid.dae"
               .format(species, yearMonthDay, hourMinute, isoValueStr, safeUnits,
               layerAbbrev(atmosLayer)))
            daePathName = stageDirectory + daeNameOnly
            daeColor = [[1, 1, 1, 1], 0.0]	# white, opaque
            utilsCollada.writeDAElines(lineList[atmosLayer],
               daeColor[0], daeColor[1], daePathName)

            # create the corresponding KML placemark for that wireframe
            daeNameOnly = os.path.splitext(daeNameOnly)[0]
            surfaceMark = acomKml.renderSurfaceKML(daeNameOnly, markCount,
               centerLat, centerLon, 0.0,
               isoValue, startTime, endTime, chemDivisions)
            surfaceIntoLayer(surfaceMark, atmosLayer, layerFolders, isoIndex)
            markCount += 1

   return([cellCount, triangleCount, lineCount, markCount])



# Create isosurface at one chemical threshold.
# chemValues = 3D grid of chemical concentrations
# lats, lons = coordinates associated with chemValues
# heights = z-values for the grid
# units = for the chemical species
# ground = terrain elevation
# isoValue = create 3D surface along this chemical value
# minPlotHeight = start plotting at this level
# maxPlotHeight = go up to this level and no higher
# layerFolders = KML folders for each layer of the atmosphere
# thresholdIndex = which isoValue within the atmospheric layer
# gridLines = create 3D grid/contour lines in a cage around the surface
# return array of [cells, triangles, lines, triLayers, lineLayers]
#	cells = number of cells that intersect the surface
#	triangles = how many triangles created here
#	lines = how many multi-segment grid lines created here
#	triLayers = list of triangles created for each atmospheric layer
#	lineLayers = list of grid lines created for each atmospheric layer
def createOneSurface(chemValues, lats, lons, heights, units, ground,
   isoValue, minPlotHeight, maxPlotHeight, layerFolders, thresholdIndex,
   gridLines=True):
   progress("Creating isoSurface at {} {}.".format(isoValue, units))

   # bogus - flatten the model layers for testing
   if (False):
      for hi in range(heights.shape[0]):
         heights[hi, :, :] = heights[hi, 0, 0]

   # create a quick and easy way to access height for this z-level
   centerHeights = heights[:, int(heights.shape[1] / 2), int(heights.shape[2] / 2)]
   #progress("centerHeights = {} km".format(centerHeights))

   cellCount = 0
   triangleCount = 0
   lineCount = 0

   # begin with a layer list of empty triangle lists
   trianglesInLayer = [ [] for _ in range(len(layerFolders)) ]
   linesInLayer = [ [] for _ in range(len(layerFolders)) ]

   for z in range(0, chemValues.shape[0] - 1):
      if (maxPlotHeight <= 0.0 and z > 0):
         # surface plot only
         continue

      # check if level within altitude range
      if (minPlotHeight >= 0.0 and centerHeights[z] < minPlotHeight):
         continue
      if (maxPlotHeight > 0.0 and centerHeights[z] > maxPlotHeight):
         continue

      for y in range(0, chemValues.shape[1] - 1):

         for x in range(0, chemValues.shape[2] - 1):
            # Does the isosurface pass through this cell?
            if (not includeSlab(chemValues, z, y, x, isoValue)):
               continue
            cellCount += 1

            # calculate the slab dimensions
            latSpan = lats[y + 1, x] - lats[y, x]
            lonSpan = lons[y, x + 1] - lons[y, x]

            # calculate slab height for this level
            heightSpan = (heights[z + 1, y, x] - heights[z, y, x])
            #progress("height = {} heightSpan = {}"
            #   .format(heights[z, y, x], heightSpan))
            #progress("heights = {}".format(heights[:, y, x]))

            # represent isosurface inside grid cell with triangles
            isoTriangles, isoLines = MarchingCubes.renderCube(
               lons, lats, heights, chemValues,
               x, y, z, isoValue, gridLines)
            triangleCount += len(isoTriangles)
            #progress("First triangle = {}\n{}".format(isoTriangles[0], isoTriangles[0].toString()))
            #progress("Second triangle = {}\n{}".format(isoTriangles[1], isoTriangles[1].toString()))

            # assign triangles to one of the atmospheric layers
            layerIndex = heightToLayer(z, heights[z, y, x])
            trianglesInLayer[layerIndex].extend(isoTriangles)

            # assign lines to one of the atmospheric layers
            lineCount += len(isoLines)
            linesInLayer[layerIndex].extend(isoLines)

            #break	# bogus
         #if (triangleCount > 0):
         #   break		# bogus
      #if (triangleCount > 0):
      #   break		# bogus

   # display what we found
   progress("For each atmospheric layer, we have:")
   layerIndex = 0
   for tris, lines in zip(trianglesInLayer, linesInLayer):
      if (len(tris) > 0):
         progress("\t{} has {} triangles and {} perimeter lines: "
            .format(layerIndex, len(tris), len(lines)))
      layerIndex += 1

   if (False):
      # create upper slab for adjusting height scaling factor
      if (levIndex < len(heights)-1):
         kmlMark = renderSlab(str(maxValue) + " upper", triangleCount,
            lats[latIndex, lonIndex], lons[latIndex, lonIndex],
            (heights[levIndex+1, latIndex, lonIndex] - ground[y, x]) * heightExaggeration,
            latSpan, lonSpan, heightSpan * heightExaggeration,
            chemValues[levIndex+1, latIndex, lonIndex],
            startTime, endTime, chemDivisions)
         triangleCount += 1
         layerFolders[0].append(kmlMark)

   progress("The isosurface passed through {} grid cells.".format(cellCount))

   return([cellCount, triangleCount, lineCount, trianglesInLayer, linesInLayer])



# Create a colorbar for the displayed chemical.
# barFolder = KML folder that wil contain the bar
def createColorbar(barFolder, chemicalName):
   screenOverlay = acomKml.kmlElement("ScreenOverlay")
   barFolder.append(screenOverlay)

   barName = acomKml.kmlElement("name", "Legend: {}".format(chemicalName))
   screenOverlay.append(barName)

   icon = acomKml.kmlElement("Icon")
   screenOverlay.append(icon)
   href = acomKml.kmlElement("href", "{}-colorbar.png".format(chemicalName))
   icon.append(href)

   overlayXY = acomKml.kmlElement("overlayXY")
   overlayXY.set("x", "1.0")
   overlayXY.set("y", "1.0")
   overlayXY.set("xunits", "fraction")
   overlayXY.set("yunits", "fraction")
   screenOverlay.append(overlayXY)

   screenXY = acomKml.kmlElement("screenXY")
   screenXY.set("x", "90")
   screenXY.set("y", "4")
   screenXY.set("xunits", "insetPixels")
   screenXY.set("yunits", "insetPixels")
   screenOverlay.append(screenXY)

   size = acomKml.kmlElement("size")
   size.set("x", "400")
   size.set("y", "0")
   size.set("xunits", "pixels")
   size.set("yunits", "pixels")
   screenOverlay.append(size)



kmPerDegree = 111.32	# kilometers per degree at equator

# Create initial view of domain.
# latBox, lonBox = domain span
# return KML element for document
def viewCenter(latBox, lonBox):
   centerLat = (latBox[0] + latBox[1]) / 2.0
   centerLon = (lonBox[0] + lonBox[1]) / 2.0
   centerLon = acomKml.wrapLongitude(centerLon)

   latRange = (latBox[1] - latBox[0]) / 2.0
   latRange *= kmPerDegree * 1000.0		# degrees to meters

   lookAt = acomKml.kmlElement("LookAt")
   lookAt.append(acomKml.kmlElement("latitude", "{}".format(centerLat)))
   lookAt.append(acomKml.kmlElement("longitude", "{}".format(centerLon)))
   lookAt.append(acomKml.kmlElement("altitude", "{}".format(0.0)))
   lookAt.append(acomKml.kmlElement("range", "{}".format(latRange)))
   lookAt.append(acomKml.kmlElement("heading", "{}".format(0.0)))
   lookAt.append(acomKml.kmlElement("tilt", "{}".format(70.0)))
   lookAt.append(acomKml.kmlElement("altitudeMode", "absolute"))

   return(lookAt)



# Create flying tour of the bounding box.
# myFolder = place the tour within here
# startWhen, endWhen = date-time of flying time
# latBox, lonBox = box around which to fly, probably contains slabs
# runMode = "manual" or "auto"
def createFlight(myFolder, startWhen, endWhen, latBox, lonBox,
   runMode):

   # calculate center of box, maybe use for LookAt target
   latCenter = (latBox[0] + latBox[1]) / 2.0
   lonCenter = (lonBox[0] + lonBox[1]) / 2.0
   lookAtPoint = [latCenter, lonCenter]
   if (False):
      lookAtPoint = [-35.30826, 149.12464]	# Canberra
   if (False):
      lookAtPoint = [40.08171, -104.81291]	# Fort Lupton
   if (False):
      lookAtPoint = [90.0, 0.0]			# Arctic North Pole
   if (False):
      lookAtPoint = [40.63, -105.82]		# Cameron Peak Fire 2020

   # calculate range from which to view central point
   latRange = abs(lookAtPoint[Model4D.LAT] - latBox[0])
   lonRange = abs(lookAtPoint[Model4D.LON] - lonBox[0])
   flightRange = math.hypot(latRange, lonRange) * kmPerDegree
   flightRange /= 2.0		# get closer to Canberra
   if (False):
      flightRange *= 1.8	# farther away from Mt. Elbert
   if (False):
      flightRange *= 1.9094	# Iceland to North Pole
   if (False):
      flightRange *= 0.5	# closer to Cameron Peak Fire
   progress("flightRange = {} km".format(flightRange))

   # time spent flying each leg
   angleStart = 45.0
   if (True):
      angleStart = -(-22.6055)	# Iceland (longitude)
   numLegs = 12
   legTime = (endWhen - startWhen) / numLegs
   progress("legTime = {}".format(legTime))
   legAngle = 360.0 / numLegs

   # set takeoff/landing and crusing altitudes
   tiltAngle = [88.0, 80.0]	# Canberra
   if (True):
      tiltAngle = [70.0, 60.0]	# Colorado
      tiltAngle = [80.0, 70.0]	# Front Range
   if (False):
      tiltAngle = [90.0, 60.0]	# Arctic
 
   if (runMode == "auto"):
      # fly in a circle around the domain
      lookAtPoint = [latCenter, lonCenter]
      latRange = abs(lookAtPoint[Model4D.LAT] - latBox[0])
      lonRange = abs(lookAtPoint[Model4D.LON] - lonBox[0])
      flightRange = math.hypot(latRange, lonRange) * kmPerDegree
      progress("flightRange 2 = {} km".format(flightRange))
      angleStart = 45.0
      tiltAngle = [80.0, 70.0]

   flightRange *= metersPerKm

   # set up the flying tour
   tour = acomKml.kmlElement("gx:Tour")
   myFolder.append(tour)

   tourName = acomKml.kmlElement("name",
      "Tour-{0.year}{0.month:02}{0.day:02}to{1.year}{1.month:02}{1.day:02}"
      .format(startWhen, endWhen))
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
   lookAt.append(acomKml.kmlElement("altitudeMode", "relative"))
   takeoffLanding = None

   keflavik = False
   if (keflavik and runMode != "auto"):
      # Arctic - launch from Keflavik airport
      takeoffLanding = acomKml.kmlElement("Camera")
      takeoffLanding.append(acomKml.kmlElement("latitude", "{}".format(63.9850)))
      takeoffLanding.append(acomKml.kmlElement("longitude", "{}".format(-angleStart)))
      takeoffLanding.append(acomKml.kmlElement("altitude", "{}".format(1000.0)))
      takeoffLanding.append(acomKml.kmlElement("altitudeMode", "relativeToGround"))
      takeoffLanding.append(acomKml.kmlElement("heading", "{}".format(70)))
      takeoffLanding.append(acomKml.kmlElement("tilt", "{}".format(90)))
      lookAt = copy.deepcopy(takeoffLanding)

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
      legHeading = (angleStart + legAngle * legi * direction) % 360
      legTilt = tiltAngle[1]
      if (legi == 1 or legi == numLegs - 1):
         # gradually increase/decrease altitude
         legTilt = (tiltAngle[0] + tiltAngle[1]) / 2
      if (legi == numLegs):
         # land
         legTilt = tiltAngle[0]

      # begin flight at southwest corner of bounding box
      flyTo = acomKml.kmlElement("gx:FlyTo")
      playlist.append(flyTo)
      flyTo.append(acomKml.kmlElement("gx:duration", "4.0"))	# seconds in real time
      flyTo.append(acomKml.kmlElement("gx:flyToMode", "smooth"))

      lookAt = acomKml.kmlElement("LookAt")
      lookAt.append(acomKml.kmlElement("latitude", "{}".format(lookAtPoint[Model4D.LAT])))
      lookAt.append(acomKml.kmlElement("longitude", "{}".format(lookAtPoint[Model4D.LON])))
      lookAt.append(acomKml.kmlElement("altitude", "{}".format(0.0)))
      lookAt.append(acomKml.kmlElement("range", "{}".format(flightRange)))
      lookAt.append(acomKml.kmlElement("heading", "{}".format(legHeading)))
      lookAt.append(acomKml.kmlElement("tilt", "{}".format(legTilt)))
      lookAt.append(acomKml.kmlElement("altitudeMode", "absolute"))

      if (legi == numLegs and keflavik and runMode != "auto"):
         # Arctic - land back at Keflavik airport
         lookAt = takeoffLanding

      flyTo.append(lookAt)

      when = acomKml.kmlElement("gx:TimeStamp")
      lookAt.append(when)
      when.append(acomKml.kmlElement("when",
         legDone.strftime("%Y-%m-%dT%H:%M:%SZ")))

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
      saniParam = utilsLite.sanitize(param)

      if (integer):
         params.append(utilsLite.safeInt(saniParam))
      elif (numeric):
         params.append(utilsLite.safeFloat(saniParam))
      else:
         params.append(saniParam)

   return(params)



STAGE_DIRECTORY = "stage/"

# Main program begins here.
def main():
   # set up the default model and species to plot
   model = "wrf-chem"
   species = ["o3"]
   minHeight = 0.0		# meters
   maxHeight = 5.0e3		# meters; use 1e6 for unlimited height, -1 for surface only
   runMode = "auto"
   userHourStride = None
   gridLines = True		# show contour lines in a cage around the isosurfaces

   # default to CONUS
   latBounds = [25.0, 50.0]	# degrees
   lonBounds = [-130.0, -60.0]	# degrees

   dates = None
   hours = [0, 0]
   dataDir = None

   # set up levels of the atmosphere
   spheres = ["Surface", "Troposphere", "Stratosphere", "Mesosphere", "Thermosphere"]

   # allow caller to specify chemical threshold values
   manualThresholds = None

   # default to no airplane research flights
   airplanes = []

   # retrieve the command-line arguments, if any
   for argPair in sys.argv:
      progress("argPair = " + argPair)
      # the arguments are: arg=value
      pairValue = argPair.split('=')
      if (len(pairValue) < 2):
         continue

      if (pairValue[0].lower() == "runmode"):
         runMode = utilsLite.sanitize(pairValue[1]).lower()
      if (pairValue[0].lower() == "model"):
         model = utilsLite.sanitize(pairValue[1]).lower()

      if (pairValue[0].lower() == "species"):
         species = safeParams(pairValue[1])

      if (pairValue[0].lower() == "minheight"):
         minHeight = utilsLite.safeFloat(pairValue[1])
      if (pairValue[0].lower() == "maxheight"):
         maxHeight = utilsLite.safeFloat(pairValue[1])

      if (pairValue[0].lower() == "latbounds"):
         latBounds = safeParams(pairValue[1], numeric=True)
      if (pairValue[0].lower() == "lonbounds"):
         lonBounds = safeParams(pairValue[1], numeric=True)

      if (pairValue[0].lower() == "dates"):
         dates = safeParams(pairValue[1], integer=True)
      if (pairValue[0].lower() == "hours"):
         hours = safeParams(pairValue[1], numeric=True)

      if (pairValue[0].lower() == "hourstride"):
         userHourStride = utilsLite.safeFloat(pairValue[1])

      if (pairValue[0].lower() == "datadir"):
         dataDir = pairValue[1]

      if (pairValue[0].lower() == "thresholds"):
         manualThresholds = safeParams(pairValue[1], numeric=True)

      if (pairValue[0].lower() == "airplane"):
         oneFlight = pairValue[1].split(",")
         airplanes.append(oneFlight)

      if (pairValue[0].lower() == "grid"):
         gridLines = (pairValue[1].lower() != "no")	# default is yes

   # Display the command-line arguments just received.
   progress("runMode = {}".format(runMode))
   progress("model = {}".format(model))
   progress("species = {}".format(species))
   progress("minHeight = {}".format(minHeight))
   progress("maxHeight = {}".format(maxHeight))

   progress("latBounds = {}".format(latBounds))
   progress("lonBounds = {}".format(lonBounds))

   progress("dates = {}".format(dates))
   progress("hours = {}".format(hours))

   progress("dataDir = {}".format(dataDir))

   # Process the command-line arguments.
   useModel = WrfChem4D.WrfChemModel()
   if (model == "waccm"):
      useModel = Waccm4D.WaccmModel()
   if (model == "cam-chem"):
      useModel = CamChem4D.CamChemModel()

   if (dataDir is not None):
      useModel.setBaseDirectory(dataDir)

   if (dates is None):
      # default time span is two days centered on right now
      startDate = datetime.datetime.utcnow()
      deltaHours = datetime.timedelta(hours=24)
      endDate = startDate + deltaHours
      startDate -= deltaHours
   else:
      # set up the time bounds
      minutes0 = int((hours[0] % 1) * 60 + 0.5)
      minutes1 = int((hours[1] % 1) * 60 + 0.5)

      dateStr = "{}:{}:{}".format(dates[0], int(hours[0]), minutes0)
      startDate = datetime.datetime.strptime(dateStr, "%Y%m%d:%H:%M")
      dateStr = "{}:{}:{}".format(dates[1], int(hours[1]), minutes1)
      endDate = datetime.datetime.strptime(dateStr, "%Y%m%d:%H:%M")

   progress("Date range is {} to {}".format(startDate, endDate))

   # get time step for this model
   hourStride = useModel.getHourStride()
   if (userHourStride is not None):
      hourStride = userHourStride
   progress("hourStride = {}".format(hourStride))

   # display cells in color above certain chemical concentrations
   chemThresholds = []
   chemUnits = []
   for chem in species:
      myThresholds = useModel.getChemThresholds(chem)
      if (manualThresholds is not None):
         chemThresholds.append(manualThresholds)
      else:
         chemThresholds.append(myThresholds[0])
      chemUnits.append(myThresholds[1])

   progress("Chemical thresholds:")
   for chem, thresh, units in zip(species, chemThresholds, chemUnits):
      progress("\t{} {} {}".format(chem, thresh, units))

   for flight in airplanes:
      progress("airplane = {} flying track {}"
         .format(flight[0], flight[1]))

   progress("terrainExaggeration = {}".format(terrainExaggeration))

   # set up how many model cells to skip between samples
   latStride = 1
   lonStride = 1
   verticalStride = 1

   maxChemValue = math.inf
   if (runMode == "auto"):
      # set up the daily automatic run
      progress("Run auto from {} to {}".format(startDate, endDate))

      maxChem = useModel.findMaxChemValue(
         species[0], startDate, endDate, hourStride)
      latMaxChem = maxChem[0][Model4D.LAT]
      lonMaxChem = maxChem[0][Model4D.LON]
      maxChemValue = maxChem[1]
      maxChemDate = maxChem[2]
      maxUnits = maxChem[3]
      progress("Maximum surface value of {} {} occurred at {} N, {} W on {}"
         .format(maxChemValue, maxUnits, latMaxChem, lonMaxChem, maxChemDate))

      # create domain around that site
      domainHalfSpan = 150.0		# kilometers
      kmToDegrees = 1.0/kmPerDegree	# convert km to degrees
      halfSpanDeg = domainHalfSpan * kmToDegrees

      latBounds = [latMaxChem - halfSpanDeg, latMaxChem + halfSpanDeg]
      lonBounds = [lonMaxChem - halfSpanDeg, lonMaxChem + halfSpanDeg]
      progress("Auto domain is lats {} and lons {}".format(latBounds, lonBounds))

      if (maxUnits == "ppmv"):
         maxChemValue *= 1e3
         maxUnits = "ppbv"

      progress("Maximum chemical value is {} {}".format(maxChemValue, maxUnits))

   timeDelta = datetime.timedelta(days=0, seconds = hourStride * 3600)

   chemIndex = -1
   for chemical in species:
      chemIndex += 1
      progress("Chemical: {}".format(chemical))

      # create reference name like this: o3-20210418to20210420
      kmlBasename = ("{}-{:04d}{:02d}{:02d}to{:04d}{:02d}{:02d}"
         .format(chemical,
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
         "{} animation for chemical {}.\n"
            .format(useModel.getModelName(), chemical)
         + "Dates: {} to {} UTC.\n".format(startDate, endDate)
         + "Created by bushFires.py on {} local time.\n"
            .format(datetime.datetime.now())
         + "Author: Carl Drews, ACOM/NCAR/UCAR.\n"
      )
      kmlDoc.append(kmlDescription)

      # create a sensible initial view
      kmlInitialView = viewCenter(latBounds, lonBounds)

      # look at the entire time span here
      allTime = acomKml.timeSpan(startDate,
         endDate + datetime.timedelta(hours=1),
         prefix="gx:")
      kmlInitialView.append(allTime)

      kmlDoc.append(kmlInitialView)

      # create a folder for the smoke / chemical concentration
      smokeFolder = acomKml.folder(chemical)
      kmlDoc.append(smokeFolder)

      # create folders for the atmospheric layers
      sphereFolders = []
      for sphere in spheres:
         sphereFolder = acomKml.folder(sphere)
         sphereFolders.append(sphereFolder)
         smokeFolder.append(sphereFolder)

         # create folders for thresholds
         for ti in range(len(chemThresholds[chemIndex])):
            upperLevel = maxChemValue
            if (ti < len(chemThresholds[chemIndex]) - 1):
               upperLevel = chemThresholds[chemIndex][ti + 1]

            folderName = ("Value range: {:.1f} < {} < {:.1f} ppbv".format(
               chemThresholds[chemIndex][ti], species[0], upperLevel))
            threshFolder = acomKml.folder(folderName)
            sphereFolder.append(threshFolder)

      # keep track of model complexity and algorithm
      cellCount = 0		# how many grid cells intersected any isosurface
      triangleCount = 0		# total number of 3D triangles used
      lineCount = 0		# total number of multi-segment lines constructed
      markCount = 0		# total number of KML placemarks created

      # loop through dates and times
      frame = startDate
      units = "ppbv"

      while (frame <= endDate):
         progress("\nFrame time: {}".format(frame))
         waccmData = None

         # locate and open WACCM model output
         filename = (useModel.BASE_DIRECTORY
            + useModel.getFilename(frame.year, frame.month, frame.day,
            frame.hour, frame.minute))
         if (not os.path.exists(filename)):
            progress("{} file {} does not exist."
               .format(useModel.getModelName(), filename))
            frame += timeDelta
            continue

         # Note: Depending on the time stride and configuration of the WACCM files,
         # we might repeatedly open and read from the same file. Here may be on
         # opportunity for optimization, depending on performance and flexibility
         # in adapting to WACCM ouput. Carl Drews - January 9, 2020

         # read only data for this bounding box
         progress("Reading {} from {} file {}"
            .format(chemical, useModel.getModelName(), filename))
         waccmData = useModel.readModelBox(chemical, filename, frame,
            latBounds, lonBounds, latStride, lonStride, verticalStride)
         units = waccmData[4]

         # create a single 3-D visual at this frame time
         counts = createOneFrame(frame, frame + timeDelta, sphereFolders,
            cellCount, triangleCount, lineCount, markCount, waccmData,
            chemThresholds[chemIndex], minHeight, maxHeight,
            stageDir, chemical, gridLines)
         cellCount = counts[0]
         triangleCount = counts[1]
         lineCount = counts[2]
         markCount = counts[3]
         progress(
            "cell count = {}   triangle count = {}   line count = {}   placemark count = {}"
            .format(cellCount, triangleCount, lineCount, markCount))

         # move on to the next time frame
         frame += timeDelta

         #if (triangleCount > 0):
         #   break		# bogus

      progress("")
      progress(
         "Examined {} grid cells; created {} triangles, {} lines, and {} placemarks."
         .format(cellCount, triangleCount, lineCount, markCount))
      kmlDescription.text += ("Contains {} triangles, {} lines, and {} placemarks.\n"
         .format(triangleCount, lineCount, markCount))
      kmlDescription.text += "Set terrain exaggeration = 3 in Google Earth.\n"
      kmlDescription.text += ("Nominal minHeight = {} maxHeight = {} meters.\n"
         .format(minHeight, maxHeight))
      kmlDescription.text += ("Thresholds for {} are {} {}."
         .format(chemical, chemThresholds[chemIndex], units))

      # create a folder for the colorbar
      colorbarFolder = acomKml.folder("Colorbar")
      smokeFolder.append(colorbarFolder)
      createColorbar(colorbarFolder, chemical)

      # create a folder for the tours around the smoke
      flightFolder = acomKml.folder("Tours")
      kmlDoc.append(flightFolder)

      # end flights after endDate so all frames are fully seen
      endDatePlus = endDate + datetime.timedelta(hours=1)  # one extra hour
      createFlight(flightFolder, startDate, endDatePlus,
         latBounds, lonBounds, runMode)

      # create customized hover tour
      if (True):
         # Williams Flats Fire
         hoverStart = datetime.datetime(2019, 8, 8, 23, 50)
         hover2 = datetime.datetime(2019, 8, 9, 1, 00)
         hover3 = datetime.datetime(2019, 8, 9, 2, 00)
         hoverEnd = datetime.datetime(2019, 8, 9, 3, 10)
         acomKml.createHoverTour(flightFolder,
            [hoverStart, hover2, hover3, hoverEnd],	# time bounds
            [48.2, 48.2, 48.2, 48.2],			# latitude
            [-117.7, -118.0, -118.0, -117.7],		# longitude
            [3000, 2000, 1000, 1000],			# altitude
            [200000, 200000, 180000, 150000],		# range
            [100, 30, 340, 300],			# heading
            [70, 65, 65, 65],				# tilt
            [4, 7, 7, 7],				# duration (seconds)
            [0, 0, 0, 0])

      # create a folder for rulers that show altitude and other measurements
      if (True):
         imageFilename = "3000metersAltitude.png"
         rulerFolder = acomKml.folder("Rulers")
         kmlDoc.append(rulerFolder)
         actualLatBounds = [waccmData[1].min(), waccmData[1].max()]
         actualLonBounds = [waccmData[2].min(), waccmData[2].max()]
         progress("actualLatBounds = {}".format(actualLatBounds))
         progress("actualLonBounds = {}".format(actualLonBounds))
         acomKml.createRulers(rulerFolder, actualLatBounds, actualLonBounds,
            imageFilename, terrainExaggeration)
         shutil.copy(imageFilename, stageDir)

      # set up for research flight paths
      if (len(airplanes) > 0):
         airplaneFolder = acomKml.folder("Research flights")
         kmlDoc.append(airplaneFolder)

      # create research airplane flying along track
      firstFlight = True
      for flight in airplanes:
         airplaneName = flight[0]
         airplaneTrackFile = flight[1]

         planeFile = shapes3D.createAirplane(airplaneName, stageDir)
         progress("planeFile = {}".format(planeFile))

         # don't show the pushpins
         if (firstFlight):
            dontShowKML = acomKml.createDontShowStyle(idStr="airplane", lineWidth=5)
            kmlDoc.append(dontShowKML)
            fixedPathKML = acomKml.createFixedPathStyle("fixedPath", "ccffff00")
            kmlDoc.append(fixedPathKML)
            outlineKML = acomKml.createOutlineStyle()
            kmlDoc.append(outlineKML)

         # create airplane track and fixed path
         airplaneKMLs = shapes3D.createAirplaneTrack(
            airplaneName, planeFile, airplaneTrackFile,
            terrainExaggeration)
         airplaneFolder.append(airplaneKMLs[0])
         airplaneFolder.append(airplaneKMLs[1])
         airplaneFolder.append(airplaneKMLs[2])

         firstFlight = False

      # write the formatted kml file
      acomKml.indent(kmlRoot)
      kmlName = stageDir + kmlBasename + ".kml"
      progress("Writing to KML file {}".format(kmlName))
      kmlTree.write(kmlName)

      # draw colorbar for auto-scaled chemical thresholds
      colorbar.drawColorbar(chemical, useModel.getChemName(chemical),
         chemThresholds[chemIndex], units, 100, stageDir)         # WRF-Chem near surface

      # create KMZ compressed archive
      kmzName = kmlBasename + ".kmz"
      dailyDir = "daily/"
      cmd = ("zip --junk-paths -FS {}{}-{} {} {}*.dae {}*.png"
         .format(dailyDir, runMode.capitalize(), kmzName, kmlName,
         stageDir, stageDir))
      progress("{}".format(cmd))
      returnCode = subprocess.call(cmd, shell=True)

      # end of chemical species loop

   return(0)	# no error



# call the main
progress("{}".format(__file__))
progress("Start time: {} Local".format(datetime.datetime.now()))
retValue = main()
progress("End time: {} Local".format(datetime.datetime.now()))

sys.exit(retValue)

