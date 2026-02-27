#!/usr/bin/env python3

# Era4D.py
# Python routines to read ERA5 model output for 4-D visualization.
#
# Author: Carl Drews, Atmospheric Chemistry Observations & Modeling
# Created: February 2026
# Copyright 2026 by the University Corporation for Atmospheric Research



import CamChem4D

import datetime
import netCDF4
import numpy



class EraModel(CamChem4D.CamChemModel):
   
   def __init__(self):
      self.BASE_DIRECTORY = "/home/drews/GoogleEarth/RenSmith/ACCLIP/ERA5/"
      return



   # return friendly name or short acronym for this model
   def getModelName(self):
      return("ERA5")

   # return the time step for this model in hours
   def getHourStride(self):
      return(3.0)

   def getFilename(self, year, month, day, hour,
      minute=0):
      hourFloat = hour + minute/60.0
      if (hourFloat < 0.5):
         # the zero hour (midnight) is in the previous day
         year, month, day = self.previousDay(year, month, day)

      return("FCnudged_f09.carmats16_cesm2.2.carma_trop_strat.2010_2021_cams5.3_so2_meic_withbb_RichConv.002.Aug19_E5avg.cam.h1.{:04d}-{:02d}-{:02d}-10800.nc"
         .format(year, month, day))

