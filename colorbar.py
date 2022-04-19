#!/usr/bin/env python3

# Python script to create a single colorbar for Google Earth.
#
# Author: Carl Drews, Atmospheric Chemistry Observations & Modeling
# Created: February 2020
# Copyright 2020 by the University Corporation for Atmospheric Research

# set up environment variables needed by the imports
import os
os.environ['MPLCONFIGDIR']="/tmp"
os.environ['NCARG']="/usr/local/ncarg"
os.environ['NCARG_ROOT']="/usr/local/ncarg"

import warnings
import matplotlib

import numpy
import datetime
import sys

# set up for PNG image output
matplotlib.use('Agg')
from matplotlib.backends.backend_agg import FigureCanvasAgg
import matplotlib.pyplot



# Display progress message on the console.
# endString = set this to '' for no return
def progress(message, endString='\n'):
   if (True):   # disable here in production
      print(message, end=endString)
      sys.stdout.flush()



# Create color map to use for CO fire emissions in Google Earth.
# The color levels are equally spaced.
def coFiresColorMapEqual():

   mapColors = [
      [0,1,0,1],        # green
      [1,1,0,1],        # yellow
      [1,0,0,1],        # red
   ]
   numColors = len(mapColors)

   myMap = matplotlib.colors.ListedColormap(
      mapColors[0:numColors-1], name='CO-Fires', N=2)

   # set extension colors
   myMap.set_over(mapColors[numColors-1])

   return(myMap)



# Create color map to use for CO fire emissions in Google Earth.
# The discrete levels are: 300, 600, 1000
# centerFraction = center divider falls here (normalized to 1.0)
# numSegments = divisions for center to land on tick mark
def coFiresColorMap(centerFraction, numSegments):

   # set up color scale with irregular spacing
   numPts = 3
   points = [0.0, centerFraction, 1.0]

   # set up color dictionary
   cdict = {'red': ((points[0], 0.0, 0.0),
                    (points[1], 0.0, 1.0),	# turn on red at this anchor
                    (points[2], 1.0, 1.0)),
            'green': ((points[0], 1.0, 1.0),    # green
                      (points[1], 1.0, 1.0),    # yellow
                      (points[2], 1.0, 1.0)),   # yellow
            'blue': ((points[0], 0.0, 0.0),
                     (points[1], 0.0, 0.0),
                     (points[2], 0.0, 0.0)),
            'alpha': ((points[0], 1.0, 1.0),
                      (points[1], 1.0, 1.0),
                      (points[2], 1.0, 1.0))
   }

   myMap = matplotlib.colors.LinearSegmentedColormap(
      'CO-Fires', cdict, numSegments)

   # set extension colors
   myMap.set_over([1,0,0,1])	# red

   return(myMap)



# Create one colorbar for the requested chemical.
# chemical = species name in WACCM
# chemName = more friendly name
# contours = divisions for the colors
# units = of the numerical values
# numColorDivs = divisions to separate color on tick mark
# colorDir = create image in thiis directory
def drawColorbar(chemical, chemName, contours, units, numColorDivs,
   colorDir="./"):
   numContours = len(contours)

   a = numpy.array([[contours[0], contours[numContours-1]]])

   matplotlib.pyplot.figure(figsize=(9, 1.5))
   middleFraction = float(contours[1] - contours[0]) / (contours[2] - contours[0])
   img = matplotlib.pyplot.imshow(a,
      cmap=coFiresColorMap(middleFraction, numColorDivs))
   matplotlib.pyplot.gca().set_visible(False)

   cax = matplotlib.pyplot.axes([0.1, 0.2, 0.8, 0.3])
   matplotlib.pyplot.colorbar(orientation="horizontal", cax=cax,
      ticks=contours, extend="max", extendfrac=0.3)

   matplotlib.pyplot.figtext(0.01, 0.38, chemName)
   matplotlib.pyplot.figtext(0.01, 0.24, units)

   filename = "{}{}-colorbar.png".format(colorDir, chemical)
   progress("Saving colorbar image as {}".format(filename))
   matplotlib.pyplot.savefig(filename,
      dpi=100, bbox_inches="tight", pad_inches=0.04)

   return(0)	# no error

