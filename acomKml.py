#!/usr/bin/env python3

# acomKml.py
# Python utility routines to create KML files of WACCM output for animation in Google Earth.
#
# Author: Carl Drews, Atmospheric Chemistry Observations & Modeling
# Created: January 2020
# Copyright 2020 by the University Corporation for Atmospheric Research



import xml.etree.ElementTree as ET
import sys
import datetime
import math



# Display progress message on the console.
# endString = set this to '' for no return
def progress(message, endString='\n'):
   if (True):   # disable here in production
      print(message, end=endString)
      sys.stdout.flush()



# Create a single KML element with body text.
def kmlElement(tag, text=None):
   elem = ET.Element(tag)
   if (text is not None):
      elem.text = text
   return(elem)



# Create a KML folder for containing other KML elements.
def folder(folderName):
   myFolder = ET.Element("Folder")
   nameElement = kmlElement("name", folderName)
   myFolder.append(nameElement)
   return(myFolder)



# Create a Location KML tag.
def location(latitude, longitude, altitude):
   loc = ET.Element("Location")

   loc.append(kmlElement("longitude", str(longitude)))
   loc.append(kmlElement("latitude", str(latitude)))
   loc.append(kmlElement("altitude", str(altitude)))

   return(loc)



# Create an Orientation KML tag.
# heading = angle of rotation right (clockwise) from due North (degrees)
def orientation(heading=0):
   ori = ET.Element("Orientation")

   ori.append(kmlElement("heading", str(heading)))
   ori.append(kmlElement("tilt", str(0)))
   ori.append(kmlElement("roll", str(0)))

   return(ori)



# Create a Scale KML tag.
def scale(scaleX, scaleY, scaleZ):
   sca = ET.Element("Scale")

   sca.append(kmlElement("x", str(scaleX)))
   sca.append(kmlElement("y", str(scaleY)))
   sca.append(kmlElement("z", str(scaleZ)))

   return(sca)



# Create a 3-D Model KML tag for Google Earth.
# modelFile = Collada .dae file
# heading = degrees to the right of due North
def model3d(modelId, latitude, longitude, altitude,
   scaleX, scaleY, scaleZ, modelFile, heading=0):
   mod = ET.Element("Model")
   mod.set("id", modelId)

   mod.append(kmlElement("altitudeMode", "absolute"))
   mod.append(location(latitude, longitude, altitude))
   if (heading != 0):
      mod.append(orientation(heading))
   mod.append(scale(scaleX, scaleY, scaleZ))

   link = ET.Element("Link")
   link.append(kmlElement("href", modelFile))
   mod.append(link)

   return(mod)



# Create a TimeSpan KML element for placemarks and overlays . . .
# startUTC = first valid date and time
# endUTC = last valid date and time
# prefix = could be "gx:"
# return KML element for appending
def timeSpan(startUTC, endUTC, prefix=None):
   spanTag = "TimeSpan"
   if (prefix is not None):
      spanTag = prefix + spanTag

   myTimeSpan = ET.Element(spanTag)
   myTimeSpan.append(kmlElement("begin", startUTC.strftime("%Y-%m-%dT%H:%M:%SZ")))
   myTimeSpan.append(kmlElement("end", endUTC.strftime("%Y-%m-%dT%H:%M:%SZ")))

   return(myTimeSpan)



# Make sure that longitude is in range(-180, 180) for Google Earth.
def wrapLongitude(rawLongitude):
   if (rawLongitude <= 180.0):
      return(rawLongitude)

   return(rawLongitude - 360.0)



# Create KML Placemark to be placed within a KML Document.
# name = name of the placemark, or None if not supplied
# modelId = numerical ID of 3-D Model tag
# latitude, longitude = coordinates of placemark
# altitude = height of placemark in meters above surface
# scaleX, scaleY = horizontal scaling factors for stretching the slab
# scaleZ = vertical scaling factor
# modelFile = filename of Collada .dae file
# timeUTC = time that this placemark is valid
# timeEndUtc = time that placemark is no longer displayed
# heading = rotation right of due North
def placemark(name, modelId, latitude, longitude, altitude,
   scaleX, scaleY, scaleZ, modelFile,
   timeUTC=None, timeEndUTC=None, heading=0):
   mark = ET.Element("Placemark")

   if (name is not None):
      kmlName = ET.Element("name")
      kmlName.text = name
      mark.append(kmlName)

   #mark.append(camera(latitude - 0.01, longitude + 0.001, altitude - 200))
   mark.append(model3d("model_{:07d}".format(modelId),
      latitude, longitude, altitude, scaleX, scaleY, scaleZ,
      modelFile, heading))

   if ((timeUTC is not None) and (timeEndUTC is None)):
      # just one time supplied
      kmlTimeStamp = ET.Element("TimeStamp")
      kmlTimeStamp.append(kmlElement("when", timeUTC.strftime("%Y-%m-%dT%H:%M:%SZ")))
      mark.append(kmlTimeStamp)

   if ((timeUTC is not None) and (timeEndUTC is not None)):
      # begin and end times are supplied
      kmlTimeSpan = timeSpan(timeUTC, timeEndUTC)
      mark.append(kmlTimeSpan)

   return(mark)



# Display one WACCM surface as a 3-D set of triangles,
# using the Placemark tag of Google Earth.
# surfaceName = short identifying string
# placeCount = integer used to create unique ID
# latitude, longitude = horizontal location of 3D model
# height = altitude of model in meters above sea level
# value = chemical concentration of slab
# startTime, endTime = time window in which slab is visible
# chemColorValues = minimum concentrations for green, yellow, and red
# heading = rotation right of due North
# return the KML placemark created
def renderSurfaceKML(surfaceName, placeCount,
   latitude, longitude, height,
   value, startTime, endTime,
   chemColorValues, heading=0):

   myMark = placemark(surfaceName, placeCount,
      latitude, wrapLongitude(longitude), height,
      1.0, 1.0, 1.0,            # scale X Y Z
      "{}.dae".format(surfaceName),
      startTime, endTime, heading)

   return(myMark)



# Create KML LineString Placemark to be placed within a KML Document.
# These are 3D lines that outline an isosurface (contours).
# name = name of the placemark, or None if not supplied
# outlineId = numerical ID of placemark
# multiLines = list of 3D lines, each with multiple segments
# timeUTC = time that this placemark is valid
# timeEndUtc = time that placemark is no longer displayed
def placemarkLines(name, outlineId, multiLines = [],
   timeUTC=None, timeEndUTC=None):
   mark = ET.Element("Placemark")
   mark.set("id", "outline_{:07d}".format(outlineId))

   if (name is not None):
      kmlName = ET.Element("name")
      kmlName.text = name
      mark.append(kmlName)

   # set up for line styling
   styleKML = kmlElement("styleUrl", "#outline")
   mark.append(styleKML)

   if ((timeUTC is not None) and (timeEndUTC is None)):
      # just one time supplied
      kmlTimeStamp = ET.Element("TimeStamp")
      kmlTimeStamp.append(kmlElement("when", timeUTC.strftime("%Y-%m-%dT%H:%M:%SZ")))
      mark.append(kmlTimeStamp)

   if ((timeUTC is not None) and (timeEndUTC is not None)):
      # begin and end times are supplied
      kmlTimeSpan = timeSpan(timeUTC, timeEndUTC)
      mark.append(kmlTimeSpan)

   return(mark)



# Display one grid as a 3-D set of LineStrings,
# using the Placemark tag of Google Earth.
# surfaceName = short identifying string
# placeCount = integer used to create unique ID
# multiLines = list of 3D lines, each with multiple segments
# startTime, endTime = time window in which slab is visible
# return the KML placemark created
def renderGridKML(surfaceName, placeCount, multiLines,
   startTime, endTime):

   myMark = placemarkLines(surfaceName, placeCount,
      multiLines, startTime, endTime)

   # add MultiGeometry and the LineStrings
   multiGeo = kmlElement("MultiGeometry")
   myMark.append(multiGeo)

   for oneLine in multiLines:
      # set up one multi-segment LineString
      myLine = kmlElement("LineString")
      myLine.append(kmlElement("altitudeMode", "absolute"))

      # populate the LineString with coordinates
      allCoordinates = ""
      for oneVertex in oneLine:
         # no spaces between coordinates!
         allCoordinates += ("{},{},{}\n".format(oneVertex.x, oneVertex.y, oneVertex.z))
      coordKML = kmlElement("coordinates", allCoordinates)
      myLine.append(coordKML)

      multiGeo.append(myLine)

   return(myMark)



# Add newlines and indentation.
# https://stackoverflow.com/questions/3095434/inserting-newlines-in-xml-file-generated-via-xml-etree-elementtree-in-python
def indent(elem, level=0):
    i = "\n" + level*"   "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "   "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i



# Create horizontal image overlay for Google Earth.
# image = filename to display over the earth
# latBounds, lonBounds = image extent in degree coordinates
#	These bounds include margins beyond the plot outline.
# height = kilometers above sea level; use -1 for ground
# startTime, endTime = time span for this overlay, or None
# colorRGBAhex = used to set transparency
# return KML fragment for adding to folder/document
def createImageOverlay(image,
      latBounds, lonBounds, height=-1.0,
      startTime=None, endTime=None,
      colorRGBAhex="60ffffff"):

      # create ground overlay with name
      overlayKML = kmlElement("GroundOverlay")
      nameKML = kmlElement("name", "Ground from {} to {}"
         .format(startTime, endTime))
      overlayKML.append(nameKML)

      if ((startTime is not None) and (endTime is not None)):
         # create time span for valid time
         overlayTimeSpan = timeSpan(startTime, endTime)
         overlayKML.append(overlayTimeSpan)

      if (height < 0.0):
         # surface
         modeKML = kmlElement("altitudeMode", "clampToGround")
         overlayKML.append(modeKML)
      else:
         # raised above surface
         modeKML = kmlElement("altitudeMode", "absolute")
         overlayKML.append(modeKML)
         altitudeKML = kmlElement("altitude", "{}".format(height))
         overlayKML.append(altitudeKML)

      # set up bounds of the image
      south = latBounds[0]
      north = latBounds[1]
      west = lonBounds[0]
      east = lonBounds[1]

      # bounding box
      latLonKML = kmlElement("LatLonBox")
      southKML = kmlElement("south", "{}".format(south))
      northKML = kmlElement("north", "{}".format(north))
      westKML = kmlElement("west", "{}".format(west))
      eastKML = kmlElement("east", "{}".format(east))
      latLonKML.extend([southKML, northKML, westKML, eastKML])
      overlayKML.append(latLonKML)

      # the image itself
      iconKML = kmlElement("Icon")
      hrefKML = kmlElement("href", image)
      iconKML.append(hrefKML)
      overlayKML.append(iconKML)

      # transparency
      colorKML = kmlElement("color", colorRGBAhex)
      overlayKML.append(colorKML)

      return(overlayKML)



# Create vertical wall with image mapped onto it as texture.
# modelFile = filename of Collada .dae model to display on the wall
# latCenter, lonCenter = in degree coordinates
# altitude = bottom of wall in kilometers above sea level
# startTime, endTime = time span for this 3D model, or None
# return KML fragment for adding to folder/document
def createImageWall(modelFile,
   latCenter, lonCenter, altitude=0.0,
   startTime=None, endTime=None):

   mark = kmlElement("Placemark")
   nameKML = kmlElement("name", "Wall for {}".format(modelFile))
   mark.append(nameKML)

   if (False):
      # for testing
      styleKML = kmlElement("styleUrl", "#pushpin")
      mark.append(styleKML)
      pointKML = kmlElement("Point")
      coordsKML = kmlElement("coordinates",
         "{}, {}, {}".format(lonCenter, latCenter, altitude))
      pointKML.append(coordsKML)
      mark.append(pointKML)

   # set valid time window
   if ((startTime is not None) and (endTime is not None)):
      # create time span for valid time
      wallTimeSpan = timeSpan(startTime, endTime)
      mark.append(wallTimeSpan)

   # Model section
   modelKML = kmlElement("Model")
   mark.append(modelKML)

   # altitude mode
   altModeKML = kmlElement("altitudeMode", "absolute")
   modelKML.append(altModeKML)

   # location, including altitude
   locationKML = kmlElement("Location")
   latKML = kmlElement("latitude", "{}".format(latCenter))
   lonKML = kmlElement("longitude", "{}".format(lonCenter))
   altKML = kmlElement("altitude", "{}".format((altitude)* 1000.0))
   locationKML.extend([latKML, lonKML, altKML])
   modelKML.append(locationKML)

   # the 3D model itself
   linkKML = kmlElement("Link")
   hrefKML = kmlElement("href", modelFile)
   linkKML.append(hrefKML)
   modelKML.append(linkKML)

   return(mark)



# Create style to hide the yellow pushpins.
# idStr = default ID string; override for airplanes
# lineWidth = width of path, probably in points
def createDontShowStyle(idStr="dontShow", lineWidth=1):
   styleKML = kmlElement("Style")
   styleKML.set("id", idStr)

   iconStyle = kmlElement("IconStyle")
   iconStyle.set("id", "{}Style".format(idStr))
   icon = kmlElement("Icon")
   icon.append(kmlElement("href"))
   iconStyle.append(icon)

   labelStyle = kmlElement("LabelStyle")
   scale = kmlElement("scale", "0")
   labelStyle.append(scale)

   lineStyle = kmlElement("LineStyle")
   width = kmlElement("width", "{}".format(lineWidth))
   lineStyle.append(width)

   styleKML.append(iconStyle)
   styleKML.append(labelStyle)
   styleKML.append(lineStyle)
   return(styleKML)



# Create 3D model moving along a track.
# modelName = name for the enclosing placemark
# modelDescription = text description of placemark
# modelFile = Collada .dae file to animate along track
# times[] = array of datetime objects for the track
# lats[], lons[], altitudes[] = arrays of coordinates same size as times[]
# holdFinal = remain in final location for one more time frame
# return KML section for inclusion in some folder
def createModelMoving(modelName, modelDescription, modelFile,
   times, lats, lons, altitudes,
   holdFinal=True, idStr="dontShow"):

   # create the enclosing Placemark
   mark = kmlElement("Placemark")
   if (modelName is not None):
      kmlName = kmlElement("name", modelName)
      mark.append(kmlName)

   mark.append(kmlElement("description", modelDescription))

   # hide the pushpin
   styleKML = kmlElement("styleUrl", "#{}".format(idStr))
   mark.append(styleKML)

   # create the track
   trackKML = kmlElement("gx:Track")
   trackKML.append(kmlElement("altitudeMode", "absolute"))

   # loop through the locations
   for time, lat, lon, altitude in zip(times, lats, lons, altitudes):
      whenKML = kmlElement("when", time.strftime("%Y-%m-%dT%H:%M:%SZ"))
      trackKML.append(whenKML)
      coordKML = kmlElement("gx:coord", "{} {} {}".format(lon, lat, altitude))
      trackKML.append(coordKML)

   if (holdFinal):
      # create one more final location to match overlay time spans
      delta = times[-1] - times[-2]
      finalTime = times[-1] + delta
      whenKML = kmlElement("when", finalTime.strftime("%Y-%m-%dT%H:%M:%SZ"))
      trackKML.append(whenKML)
      coordKML = (kmlElement("gx:coord", "{} {} {}"
         .format(lons[-1], lats[-1], altitudes[-1])))
      trackKML.append(coordKML)

   # add reference to the 3D model file
   modelID = modelName.replace(" ", "_")
   modelKML = kmlElement("Model")
   modelKML.set("id", modelID)

   link = ET.Element("Link")
   link.append(kmlElement("href", modelFile))
   modelKML.append(link)

   trackKML.append(modelKML)
   mark.append(trackKML)

   return(mark)



# Create style for airplane flight path.
# styleID = ID for this style
# styleColor = ABGR value like "ccffff00"	# 80% opaque cyan
def createFixedPathStyle(styleID, styleColor):
   styleKML = kmlElement("Style")
   styleKML.set("id", styleID)

   iconStyle = kmlElement("IconStyle")
   iconStyle.set("id", "fixedPathStyle")
   icon = kmlElement("Icon")
   icon.append(kmlElement("href"))
   iconStyle.append(icon)

   labelStyle = kmlElement("LabelStyle")
   scale = kmlElement("scale", "0")
   labelStyle.append(scale)

   lineStyle = kmlElement("LineStyle")
   width = kmlElement("width", "2")
   lineStyle.append(width)
   color = kmlElement("color", styleColor)
   lineStyle.append(color)

   styleKML.append(iconStyle)
   styleKML.append(labelStyle)
   styleKML.append(lineStyle)
   return(styleKML)



# Create 3D model moving along a track.
# Create fixed path through 3D space.
# pathName = name for the enclosing placemark
# pathDescription = text description of placemark
# lats[], lons[], altitudes[] = arrays of coordinates same size as times[]
# styleID = line style for this path
# return KML section for inclusion in some folder
def createFixedPath(pathName, pathDescription,
   lats, lons, altitudes, styleID):

   # create the enclosing Placemark
   mark = kmlElement("Placemark")
   if (pathName is not None):
      kmlName = kmlElement("name", pathName)
      mark.append(kmlName)

   mark.append(kmlElement("description", pathDescription))

   # begin with fixed path turned off
   visibility = kmlElement("visibility", "0")
   mark.append(visibility)

   # set up for line styling
   styleKML = kmlElement("styleUrl", "#{}".format(styleID))
   mark.append(styleKML)

   # create the path as a LineString
   pathKML = kmlElement("LineString")
   pathKML.append(kmlElement("altitudeMode", "absolute"))

   # loop through the locations
   allCoordinates = ""
   for lat, lon, altitude in zip(lats, lons, altitudes):
      # no spaces between coordinates!
      allCoordinates += ("{},{},{}\n".format(lon, lat, altitude))
   coordKML = kmlElement("coordinates", allCoordinates)
   pathKML.append(coordKML)

   mark.append(pathKML)

   return(mark)



# Create and fill a KML LookAt tag.
# towardLat, towardLon, towardAlt = camera looks at this location
# fromRange, fromHeading, fromTilt = camera looks from this vantage point
# return KML tag suitable for tour
def populateLookAt(towardLat, towardLon, towardAlt,
   fromRange, fromHeading, fromTilt):

   lookAt = kmlElement("LookAt")
   lookAt.append(kmlElement("latitude", "{}".format(towardLat)))
   lookAt.append(kmlElement("longitude", "{}".format(towardLon)))
   lookAt.append(kmlElement("altitude", "{}".format(towardAlt)))
   lookAt.append(kmlElement("range", "{}".format(fromRange)))
   lookAt.append(kmlElement("heading", "{}".format(fromHeading)))
   lookAt.append(kmlElement("tilt", "{}".format(fromTilt)))     # view up into the sky
   lookAt.append(kmlElement("altitudeMode", "absolute"))

   return(lookAt)



# Create and fill a KML Camera tag.
# atLat, atLon, atAlt = camera is located at these coordinates
# heading = compass angle where 0 = due North, 90 = east
# return KML tag suitable for tour
def populateCamera(atLat, atLon, atAlt, heading):

   camera = kmlElement("Camera")
   camera.append(kmlElement("latitude", "{}".format(atLat)))
   camera.append(kmlElement("longitude", "{}".format(atLon)))
   camera.append(kmlElement("altitude", "{}".format(atAlt)))
   camera.append(kmlElement("altitudeMode", "absolute"))

   camera.append(kmlElement("heading", "{}".format(heading)))	# 0 is North
   camera.append(kmlElement("tilt", "{}".format(80.0)))		# 0 is straight down

   return(camera)



# Create hovering tour between several vantage points.
# We slowly move from start to end position.
# myFolder = place the tour within here
# whens[] = date-time stamps of tour waypoints
# lats[], lons[], alts[] = waypoint locations
# ranges[], headings[], tilts[] = wapoint views
# durations[] = how many seconds this segment should take
# waits[] = seconds to hold still after each segment
def createHoverTour(myFolder, whens,
   lats, lons, alts,
   ranges, headings, tilts,
   durations, waits):

   numWaypoints = len(whens)
   progress("Hover tour numWaypoints = {}".format(numWaypoints))

   # set up the flying tour
   tour = kmlElement("gx:Tour")
   myFolder.append(tour)

   progress("Hover from {} to {}".format(whens[0], whens[-1]))

   startStr = whens[0].strftime("%Y%m%d-%H%M")
   endStr = whens[-1].strftime("%Y%m%d-%H%M")
   tourNameStr = "Hover03-{}to{}".format(startStr, endStr)
   progress("tourName = {}".format(tourNameStr))
   tourName = kmlElement("name", tourNameStr)
   tour.append(tourName)

   # the PlayList is a series of FlyTo elements
   playlist = kmlElement("gx:Playlist")
   tour.append(playlist)

   firstWaypoint = True
   for when, lat, lon, alt, range, heading, tilt, duration, wait \
      in zip(whens, lats, lons, alts, ranges, headings, tilts, durations, waits):
      flyTo = kmlElement("gx:FlyTo")
      playlist.append(flyTo)

      flyTo.append(kmlElement("gx:duration", "{}".format(duration)))
      if (firstWaypoint):
         flyTo.append(kmlElement("gx:flyToMode", "bounce"))
         firstWaypoint = False
      else:
         flyTo.append(kmlElement("gx:flyToMode", "smooth"))

      lookAt = populateLookAt(
         lat, lon, alt, range, heading, tilt)

      timeStamp = kmlElement("gx:TimeStamp")
      timeStamp.append(kmlElement("when",
         when.strftime("%Y-%m-%dT%H:%M:%SZ")))
      lookAt.append(timeStamp)

      flyTo.append(lookAt)

      if (wait > 0):
         waitHere = kmlElement("gx:Wait")
         waitHere.append(kmlElement("gx:duration",
            "{}".format(wait)))		# seconds in real time
         playlist.append(waitHere)

   return



# Create timed tour along a flight path.
# The view will look out the cockpit front window.
# flightName = name for this tour
# flightDescription = text description of this tour
# times[] = date-time stamps of tour waypoints
# lats[], lons[], alts[] = waypoint locations
# duration = how many seconds this tour should take
# holdFinal = remain in final location for one more time frame
# idStr = default ID string; override for airplanes
# return KML section for inclusion in some folder
def createFlightTour(flightName, flightDescription,
   times, lats, lons, alts,
   duration=60, holdFinal=True, idStr="dontShow"):

   numWaypoints = len(times)
   progress("Flight tour numWaypoints = {} over {} seconds"
      .format(numWaypoints, duration))

   # set up the time intervals
   firstTimeStep = 5.0	# seconds in wall-clock time
   timeStep = (duration - firstTimeStep) / (numWaypoints - 1)
   progress("timeStep = {}, then {} seconds".format(firstTimeStep, timeStep))

   # set up the flying tour
   tour = kmlElement("gx:Tour")
   progress("Flight tour from {} to {}".format(times[0], times[-1]))

   tourName = kmlElement("name", flightName)
   tour.append(tourName)

   startStr = times[0].strftime("%Y%m%d-%H%M")
   endStr = times[-1].strftime("%Y%m%d-%H%M")
   tourDescStr = flightDescription + " {} to {} UTC".format(startStr, endStr)
   tourDesc = kmlElement("description", tourDescStr)
   tour.append(tourDesc)

   # the PlayList is a series of FlyTo elements
   playlist = kmlElement("gx:Playlist")
   tour.append(playlist)

   firstWaypoint = True
   altitudeDelta = 4000		# chase plane is slightly higher
   camHeading = 0.0
   for ti in range(len(times)):
      when = times[ti]
      lat = lats[ti]
      lon = lons[ti]
      alt = alts[ti]

      # look ahead to the next waypoint
      nextLat = None
      nextLon = None
      if (ti + 1 < len(times)):
         nextLat = lats[ti+1]
         nextLon = lons[ti+1]

      flyTo = kmlElement("gx:FlyTo")
      playlist.append(flyTo)

      # move smoothly to the first point
      if (firstWaypoint):
         flyTo.append(kmlElement("gx:duration", "{}".format(firstTimeStep)))
         flyTo.append(kmlElement("gx:flyToMode", "bounce"))
         firstWaypoint = False
      else:
         flyTo.append(kmlElement("gx:duration", "{}".format(timeStep)))
         flyTo.append(kmlElement("gx:flyToMode", "smooth"))

      # point camera toward the next waypoint
      if (nextLat is not None):
         if (lat != nextLat and lon != nextLon):
            camHeading = math.atan2(nextLon - lon, nextLat - lat) * 180.0 / math.pi
      camera = populateCamera(lat, lon, alt + altitudeDelta, camHeading)

      # chase plane flies slightly behind the research aircraft
      timeStamp = kmlElement("gx:TimeStamp")
      timeDelta = datetime.timedelta(seconds=3*60)
      whenCamera = when + timeDelta
      timeStamp.append(kmlElement("when",
         whenCamera.strftime("%Y-%m-%dT%H:%M:%SZ")))
      camera.append(timeStamp)

      flyTo.append(camera)

   return(tour)



# Create style for gridded outline of isosurface.
def createOutlineStyle():
   styleKML = kmlElement("Style")
   styleKML.set("id", "outline")

   lineStyle = kmlElement("LineStyle")
   width = kmlElement("width", "1")	# pixels
   lineStyle.append(width)
   color = kmlElement("color", "ffffffff")	# 100% opaque white
   lineStyle.append(color)

   styleKML.append(lineStyle)
   return(styleKML)



# Create folder containing various measuring sticks.
# myFolder = place the KML rulers within here
# latBox, lonBox = domain in which to measure
# imageFile = .PNG image file with grid and labels
# terrainExaggeration = multiplier for terrain height in Google Earth
def createRulers(myFolder, latBox, lonBox, imageFile,
   terrainExaggeration=3):
   overlay = kmlElement("GroundOverlay")
   myFolder.append(overlay)

   overlayAltitude = 3000       # meters above sea level
   imageFile = "3000metersAltitude.png"

   # add overlay name and useful instructions
   name = kmlElement("name",
      "Show altitude of 3000 meters above sea level.")
   overlay.append(name)

   descrip = kmlElement("description",
      "Height includes 3x terrain exaggeration."
      + " Use Properties -> Altitude -> slider to adjust altitude,"
      + " or enter number directly.")
   overlay.append(descrip)

   # begin with height overlay turned off
   visibility = kmlElement("visibility", "0")
   overlay.append(visibility)

   # initial altitude and mode
   altKml = kmlElement("altitude", "{}"
      .format(overlayAltitude * terrainExaggeration))
   overlay.append(altKml)

   modeKml = kmlElement("altitudeMode", "absolute")
   overlay.append(modeKml)

   # lat-lon bounding box
   latLonBox = kmlElement("LatLonBox")
   south = kmlElement("south", "{}".format(latBox[0]))
   north = kmlElement("north", "{}".format(latBox[1]))
   west = kmlElement("west", "{}".format(lonBox[0]))
   east = kmlElement("east", "{}".format(lonBox[1]))

   latLonBox.append(south)
   latLonBox.append(north)
   latLonBox.append(west)
   latLonBox.append(east)

   overlay.append(latLonBox)

   # the overlay image itself
   icon = kmlElement("Icon")
   overlay.append(icon)
   href = kmlElement("href", imageFile)
   icon.append(href)

   # transparency
   color = kmlElement("color", "60ffffff")
   overlay.append(color)

   return



# Create FlyTo section with camera LookAt sub-section.
# duration = seconds to take, in real time
# return constructed KML fragment
def createFlyToLookAt(duration=4, flyMode="smooth",
   arriveWhen=datetime.datetime.utcnow(),
   latitude=40.01, longitude=-105.3, altitude=5e5,
   flightRange=1.2e5, heading=30, tilt=75,
   altMode="absolute"):

   flyTo = kmlElement("gx:FlyTo")
   flyTo.append(kmlElement("gx:duration", "{}".format(duration)))
   flyTo.append(kmlElement("gx:flyToMode", flyMode))

   lookAt = kmlElement("LookAt")
   when = kmlElement("gx:TimeStamp")
   lookAt.append(when)
   when.append(kmlElement("when",
      arriveWhen.strftime("%Y-%m-%dT%H:%M:%SZ")))

   lookAt.append(kmlElement("latitude", "{}".format(latitude)))
   lookAt.append(kmlElement("longitude", "{}".format(longitude)))
   lookAt.append(kmlElement("altitude", "{}".format(altitude)))
   lookAt.append(kmlElement("range", "{}".format(flightRange)))
   lookAt.append(kmlElement("heading", "{}".format(heading)))
   lookAt.append(kmlElement("tilt", "{}".format(tilt)))
   lookAt.append(kmlElement("altitudeMode", altMode))

   flyTo.append(lookAt)
   return(flyTo)

