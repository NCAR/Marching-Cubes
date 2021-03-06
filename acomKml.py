#!/usr/bin/env python3

# acomKml.py
# Python utility routines to create KML files of WACCM output for animation in Google Earth.
#
# Author: Carl Drews, Atmospheric Chemistry Observations & Modeling
# Created: January 2020
# Copyright 2020 by the University Corporation for Atmospheric Research



import xml.etree.ElementTree as ET
import sys



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
def placemark(name, myId, latitude, longitude, altitude,
   scaleX, scaleY, scaleZ, modelFile,
   timeUTC=None, timeEndUTC=None, heading=0):
   mark = ET.Element("Placemark")

   if (name is not None):
      kmlName = ET.Element("name")
      kmlName.text = name
      mark.append(kmlName)

   #mark.append(camera(latitude - 0.01, longitude + 0.001, altitude - 200))
   mark.append(model3d("model_{:07d}".format(myId),
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
# id = default ID string; override for airplanes
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
def createFixedPathStyle():
   styleKML = kmlElement("Style")
   styleKML.set("id", "fixedPath")

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
   color = kmlElement("color", "ccffff00")	# 80% opaque cyan
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
# return KML section for inclusion in some folder
def createFixedPath(pathName, pathDescription,
   lats, lons, altitudes):

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
   styleKML = kmlElement("styleUrl", "#fixedPath")
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



# Create hovering tour between several vantage points.
# We slowly move from start to end position.
# myFolder = place the tour within here
# whens[] = date-time stamps of tour waypoints
# lats[], lons[], alts[] = waypoint locations
# ranges[], headings[], tilts[] = wapoint views
# duration = how many seconds this tour should take
def createHoverTour(myFolder, whens,
   lats, lons, alts,
   ranges, headings, tilts,
   duration):

   numWaypoints = len(whens)
   progress("numWaypoints = {}".format(numWaypoints))

   # set up the time intervals
   firstTimeStep = 2.0	# seconds in wall-clock time
   timeStep = (duration - firstTimeStep) / (numWaypoints - 1)
   progress("timeStep = {}, then {} seconds".format(firstTimeStep, timeStep))

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
   for when, lat, lon, alt, range, heading, tilt \
      in zip(whens, lats, lons, alts, ranges, headings, tilts):
      flyTo = kmlElement("gx:FlyTo")
      playlist.append(flyTo)

      if (firstWaypoint):
         flyTo.append(kmlElement("gx:duration", "{}".format(firstTimeStep)))
         flyTo.append(kmlElement("gx:flyToMode", "bounce"))
         firstWaypoint = False
      else:
         flyTo.append(kmlElement("gx:duration", "{}".format(timeStep)))
         flyTo.append(kmlElement("gx:flyToMode", "smooth"))

      lookAt = populateLookAt(
         lat, lon, alt, range, heading, tilt)

      timeStamp = kmlElement("gx:TimeStamp")
      timeStamp.append(kmlElement("when",
         when.strftime("%Y-%m-%dT%H:%M:%SZ")))
      lookAt.append(timeStamp)

      flyTo.append(lookAt)

   return

