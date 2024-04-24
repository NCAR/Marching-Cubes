# Python script to create a single colorbar for Google Earth.
#
# Author: Carl Drews, Atmospheric Chemistry Observations & Modeling
# Created: February 2020
# Copyright 2020 by the University Corporation for Atmospheric Research

# set up environment variables needed by the imports
import os
#os.environ['MPLCONFIGDIR']="/tmp"
#os.environ['NCARG']="/usr/local/ncarg"
#os.environ['NCARG_ROOT']="/usr/local/ncarg"

#import warnings
#import matplotlib

#import numpy
import datetime
import sys

import colorbar

# set up for PNG image output
#matplotlib.use('Agg')
#from matplotlib.backends.backend_agg import FigureCanvasAgg
#import matplotlib.pyplot



# Display progress message on the console.
# endString = set this to '' for no return
def progress(message, endString='\n'):
   if (True):   # disable here in production
      print(message, end=endString)
      sys.stdout.flush()



# Main routine begins here.
def main():
   colorbar.drawColorbar("CO01-Canberra", "CO Fires", [300, 600, 1000], "ppbv", 7)
   colorbar.drawColorbar("O3-stratosphere", "Ozone", [8500, 8700, 9000], "ppbv", 5)	# WACCM stratosphere
   colorbar.drawColorbar("o3-surface", "Ozone", [60, 62, 64], "ppbv", 2)	# WRF-Chem surface
   colorbar.drawColorbar("o3-5000m", "Ozone", [70, 72, 75], "ppbv", 5)	# WRF-Chem near surface
   colorbar.drawColorbar("O3-Arctic", "Ozone", [6000, 6200, 6500], "ppbv", 5)	# WACCM Arctic
   return



# call the main
progress("Start time: {} Local\n".format(datetime.datetime.now()))
retValue = main()
progress("End time: {} Local\n".format(datetime.datetime.now()))

sys.exit(retValue)

