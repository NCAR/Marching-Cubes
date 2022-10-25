#!/usr/bin/env python3

# CamChem4D.py
# Python routines to read CAM-Chem model output for 4-D visualization.
#
# Author: Carl Drews, Atmospheric Chemistry Observations & Modeling
# Created: July 2022
# Copyright 2022 by the University Corporation for Atmospheric Research



import Waccm4D

import datetime
import netCDF4
import numpy



class CamChemModel(Waccm4D.WaccmModel):
   
   def __init__(self):
      self.BASE_DIRECTORY = "/cam-chem-output/"
      return



   # return friendly name or short acronym for this model
   def getModelName(self):
      return("CAM-Chem")

   # return the time step for this model in hours
   def getHourStride(self):
      return(1.0)

   def getFilename(self, year, month, day, hour,
      minute=0):
      hourFloat = hour + minute/60.0
      seconds = 0
      if (hourFloat >= 18.5):
         seconds = 66600
      elif (hourFloat >= 12.5):
         seconds = 45000
      elif (hourFloat >= 6.5):
         seconds = 23400
      elif (hourFloat >= 0.5):
         seconds = 1800

      else:
         # the zero hour (midnight) is in the previous day
         fileDate = datetime.datetime(year, month, day)
         fileDate += datetime.timedelta(-1)
         year = fileDate.year
         month = fileDate.month
         day = fileDate.day
         seconds = 66600

      return("b.e22.beta02.BWHIST.f09_g17.honga_tonga.so2_h2o25-35km_triplearea180Tg_5days.cam.h2.{:04d}-{:02d}-{:02d}-{:05d}.nc"
         .format(year, month, day, seconds))

