#!/usr/bin/env python3

# utilsCollada.py
# Python module to provide Collada utilities with module PyCollada.
#
# Author: Carl Drews, Atmospheric Chemistry Observations & Modeling
# Created: April 2021
# Copyright 2021 by the University Corporation for Atmospheric Research



import sys
import collada
import numpy
import math



# Display progress message on the console.
# endString = set this to '' for no return
def progress(message, endString='\n'):
   if (True):   # disable here in production
      print(message, end=endString)
      sys.stdout.flush()



# Write DAE file containing triangles.
# No texture, just a color and opacity.
# triangles = list of MarchingCubes.triangle objects
# colorRBGA = array of normalized rgba values[0.0-1.0]
# transparency = float 0.0-1.0, from 100% opaque to entirely clear
# filename = path and name.dae of file
def writeDAEfile(triangles, colorRGBA, transparency, filename):
   numTriangles = len(triangles)
   progress("Writing DAE Collada file {} with {} triangles."
      .format(filename, numTriangles))

   if (False):
      fp = open(filename, 'w')
      fp.write("{} triangles\n".format(len(triangles)))

      for triangle in triangles:
         fp.write("{}\n{}".format(triangle, triangle.toString()))

      fp.close()

   # set up the triangular mesh
   mesh = collada.Collada()

   # set up the asset section
   axis = collada.asset.UP_AXIS.Z_UP
   myContributor = collada.asset.Contributor(author="Carl Drews")
   assetSection = collada.asset.Asset(
      title="ACOM 3D model with color and opacity.",
      contributors=[myContributor],
      unitname="kilometer", unitmeter=1000,
      upaxis=axis)
   mesh.assetInfo = assetSection
   #progress("mesh dir = {}".format(dir(mesh)))
   #progress("mesh vars = {}".format(vars(mesh)))
   #progress("assetInfo = {}".format(mesh.assetInfo))

   # get positions of the triangle vertices
   # There are 3 vertices per triangle * 3 dimensions per vertex.
   coordsPerTriangle = 3 * 3
   mesh1Position = numpy.zeros(numTriangles * coordsPerTriangle)
   for ti, triangle in enumerate(triangles):
      pi = ti * coordsPerTriangle
      #progress("triangle {}: {}".format(ti, triangle.toString()))
      mesh1Position[pi:pi+coordsPerTriangle] = triangle.getVertexCoordinates()
   #progress("mesh1Position = {}".format(mesh1Position))

   # get normal vectors of the triangle vertices
   mesh1Normal = numpy.zeros(numTriangles * coordsPerTriangle)
   for ti, triangle in enumerate(triangles):
      pi = ti * coordsPerTriangle
      mesh1Normal[pi:pi+coordsPerTriangle] = triangle.getNormalCoordinates()
   #progress("mesh1Normal = {}".format(mesh1Normal))

   # convert positions and normals into sources
   mesh1Position_src = collada.source.FloatSource("mesh1-geometry-position",
      mesh1Position, ('X', 'Y', 'Z'))
   mesh1Normal_src = collada.source.FloatSource("mesh1-geometry-normal",
      mesh1Normal, ('X', 'Y', 'Z'))

   # set up the geometry
   geom = collada.geometry.Geometry(mesh, "mesh1-geometry", "mesh1-geometry",
      [mesh1Position_src, mesh1Normal_src])

   # create input list of vertex positions
   input_list = collada.source.InputList()
   input_list.addInput(0, 'VERTEX', "#mesh1-geometry-position")
   # I am not sure if the following line is in the right place,
   # but it does link in the normal vectors for display as desired.
   # Carl Drews - May 10, 2021
   input_list.addInput(0, 'NORMAL', "#mesh1-geometry-normal")

   # set up indexes to the triangle vertices
   indexes1 = numpy.array(range(0, numTriangles * 3))	# 3 vertices per triangle
   #progress("indexes1 = {}".format(indexes1))

   # collect triangles into the geometry
   triSet1 = geom.createTriangleSet(indexes1, input_list, "material_0_1")
   geom.primitives.append(triSet1)
   mesh.geometries.append(geom)

   # set up the surface material
   effect1 = (collada.material.Effect("material_0_1-effect",
      [], "phong",
      emission=(0.0, 0.0, 0.0, 1),
      ambient=(0.05, 0.05, 0.05, 1),
      diffuse=(colorRGBA[0], colorRGBA[1], colorRGBA[2], colorRGBA[3]),
      specular=(0.0, 0.0, 0.0, 1.0),
      shininess=16.0,
      reflective=(0.0, 0.0, 0.0, 0.0),
      transparent=None,
      transparency=transparency,
      double_sided=True))

   material1 = collada.material.Material("material_0_1ID", "material_0_1", effect1)

   mesh.effects.append(effect1)
   mesh.materials.append(material1)

   # create material and geometry node
   matNode1 = collada.scene.MaterialNode("material_0_1", material1, inputs=[])
   geomNode = collada.scene.GeometryNode(geom, [matNode1])

   # go ahead - make a scene
   node = collada.scene.Node("Model", children=[geomNode])
   myScene = collada.scene.Scene("ACOM_Scene", [node])
   mesh.scenes.append(myScene)
   mesh.scene = myScene		# bogus - do we need both?

   # write the Collada structure to DAE file
   mesh.write(filename)

   return



# Write 3D Collada model of vertical wall with image textured onto it.
# myImage = image to map onto the wall
# latBounds, lonBounds = extent of possibly diagonal wall
# heightBounds = extent of image above sea level (km)
# filename = name of the .dae file
def writeImageWallDAE(myImage,
   latBounds, lonBounds, heightBounds, filename):
 
   progress("Writing image wall {}".format(filename))

   # wall is tangent to earth's surface at wall center
   centerLat = (latBounds[0] + latBounds[1]) / 2
   centerLon = (lonBounds[0] + lonBounds[1]) / 2
   centerAltitude = 0.0

   # set up Collada assets
   mesh = collada.Collada()
   axis = collada.asset.UP_AXIS.Z_UP
   assetSection = collada.asset.Asset(
      title="Image mapped onto vertical wall: {}".format(myImage),
      unitname="kilometer", unitmeter=1000,
      upaxis=axis)
   mesh.assetInfo = assetSection
   
   #progress("mesh dir = {}".format(dir(mesh)))
   #progress("mesh vars = {}".format(vars(mesh)))
   #progress("assetInfo = {}".format(mesh.assetInfo))
   
   # set up major components of the model
   image = collada.material.CImage("material_0_1_0-image", myImage)
   surface = collada.material.Surface("material_0_1_0-image-surface", image)
   sampler2d = collada.material.Sampler2D("material_0_1_0-image-sampler", surface)
   map1 = collada.material.Map(sampler2d, "UVSET0")
   
   # set up texture effect
   effect2 = collada.material.Effect("material_0_1_0-effect", [surface, sampler2d],
      "lambert", emission=(0.0, 0.0, 0.0, 1),
      ambient=(0.0, 0.0, 0.0, 1),  diffuse=map1, transparent=map1, transparency=0.0,
      double_sided=True)
   
   # texture material
   mat2 = collada.material.Material("material_0_1_0ID", "material_0_1_0", effect2)
   mesh.effects.append(effect2)
   mesh.materials.append(mat2)
   mesh.images.append(image)
   
   #red x-axis
   #green z-axis
   #blue y-axis

   # How far does this wall extend?
   latLonSpan = math.hypot(lonBounds[1] - lonBounds[0],
      latBounds[1] - latBounds[0])
   progress("writeImageWallDAE() latLonSpan = {}".format(latLonSpan))

   # calculate the increments for the linear segments
   degreeIncrement = 2.0	# degrees between segments
   #degreeIncrement = 50.0	# bogus, for testing bent vertical grid lines
   latSpan = latBounds[1] - latBounds[0]
   lonSpan = lonBounds[1] - lonBounds[0]
   latIncrement = degreeIncrement * latSpan / latLonSpan
   lonIncrement = degreeIncrement * lonSpan / latLonSpan

   # Create the wall by marching along the bottom, then points along the top.
   # It's a lot easier to generate uv coordinate at same time in same order.
   m1position = []
   topPosition = []
   degreePos = 0.0
   latPos = latBounds[0]
   lonPos = lonBounds[0]

   m1uv = []
   m1uvTop = []
   posCount = 0
   while (True):
      # supply x, y, z in geo-coordinates
      m1position.extend([lonPos, latPos, heightBounds[0]])	# bottom row
      topPosition.extend([lonPos, latPos, heightBounds[1]])	# top row
      m1uv.extend([degreePos / latLonSpan, 0])			# bottom
      m1uvTop.extend([degreePos / latLonSpan, 1])		# top

      # move on to the next segment
      degreePos += degreeIncrement
      latPos += latIncrement
      lonPos += lonIncrement
      posCount += 1

      # have we reached the end of the wall?
      if (degreePos >= latLonSpan):
         # end with the right edge of the wall
         m1position.extend([lonBounds[1], latBounds[1], heightBounds[0]])
         topPosition.extend([lonBounds[1], latBounds[1], heightBounds[1]])
         m1uv.extend([1, 0])
         m1uvTop.extend([1, 1])
         posCount += 1
         break

   progress("There are {} horizontal positions along the wall.".format(posCount))

   # append top row to bottom row
   m1position.extend(topPosition)
   m1uv.extend(m1uvTop)
   #progress("m1position = {}   length = {}".format(m1position, len(m1position)))
   #progress("m1uv = {}   length {}".format(m1uv, len(m1uv)))

   # convert geo-coordinates to cartesian
   dimCount = 3		# number of dimensions (x,y,z)
   for vi in range(int(len(m1position) / dimCount)):		# vertex index
      posIndex = vi * dimCount					# m1position index
      cartesian = sphericalToCartesian(
         [centerLat, centerLon, centerAltitude],
         [m1position[posIndex + 1], m1position[posIndex], m1position[posIndex + 2]])

      # update vertices to Cartesian coordinates
      m1position[posIndex] = cartesian[CART_X]
      m1position[posIndex + 1] = cartesian[CART_Y]
      m1position[posIndex + 2] = cartesian[CART_Z]

   # indexes referring into m1position(xyz) and m1uv(uv)
   # vertices per triangle * two triangles (lower and upper) * (m1pos + m1uv) = 12
   indexesPerPosition = 12
   indices2 = numpy.zeros((posCount - 1) * indexesPerPosition, dtype=int)
   for pi in range(0, posCount - 1):
      baseIndex = pi * indexesPerPosition
      # lower triangle winding counter-clockwise
      indices2[baseIndex : baseIndex+2] = pi
      indices2[baseIndex+2 : baseIndex+4] = pi + 1
      indices2[baseIndex+4 : baseIndex+6] = posCount + pi + 1

      # upper triangle winding counter-clockwise
      indices2[baseIndex+6 : baseIndex+8] = posCount + pi + 1
      indices2[baseIndex+8 : baseIndex+10] = posCount + pi
      indices2[baseIndex+10 : baseIndex+12] = pi

   #progress("indices2 = {}   length {}".format(indices2, len(indices2)))

   # calculate normals of the positions along the now-bent Cartesian wall
   # For now, just use the normal of the entire wall. Carl Drews - July 2, 2021
   incrementSpan = math.hypot(latIncrement, lonIncrement)
   normalX = latIncrement / incrementSpan
   normalY = -lonIncrement / incrementSpan
   normalZ = 0.0

   m1normal = numpy.zeros(posCount * dimCount * 2)	# bottom and top row
   m1normal[0::3] = normalX
   m1normal[1::3] = normalY
   m1normal[2::3] = normalZ
   #progress("m1normal = {}   length {}".format(m1normal, len(m1normal)))

   # encapsulate those coordinates
   m1position_src = collada.source.FloatSource("mesh1-geometry-position",
      numpy.array(m1position), ('X', 'Y', 'Z'))
   m1normal_src = collada.source.FloatSource("mesh1-geometry-normal",
      numpy.array(m1normal), ('X', 'Y', 'Z'))
   m1uv_src = collada.source.FloatSource("mesh1-geometry-uv",
      numpy.array(m1uv), ('S', 'T'))
   
   geom1 = collada.geometry.Geometry(mesh,"mesh1-geometry1",
      "mesh1-geometry1",[m1position_src, m1normal_src, m1uv_src])
   
   input_list1 = collada.source.InputList()
   input_list1.addInput(0, 'VERTEX', "#mesh1-geometry-position")
   input_list1.addInput(1, 'TEXCOORD', "#mesh1-geometry-uv", set="0")
   
   triset2 = geom1.createTriangleSet(indices2, input_list1, "material_0_1_0")
   geom1.primitives.append(triset2)
   mesh.geometries.append(geom1)
   
   matnode2 = collada.scene.MaterialNode("material_0_1_0", mat2, inputs=[])
   geomnode1 = collada.scene.GeometryNode(geom1, [matnode2])
   
   node = collada.scene.Node("Model", children=[geomnode1])
   myscene = collada.scene.Scene("ACOM-Scene", [node])
   mesh.scenes.append(myscene)
   mesh.scene = myscene
   
   mesh.write(filename)
   return



# Constants for accessing the 3D arrays.
LAT = 0
LON = 1
ALT = 2

CART_X = 0
CART_Y = 1
CART_Z = 2

# Convert geo-coordinates to 3D Cartesian coordinates
# within the 3D Collada model. Large 3D models have to
# account for the curvature of the earth and the poleward
# curve of latitude lines.
# origin[lat, lon, altitude] = point of tangent to the earth's surface
#	The origin is customarily the center of the 3D model.
# geoCoord[lat, lon, altitude] = a point somewhere in the 3D model
# approximate = use approximation if degree span < 1.0 degrees
# return cartesian[x, y, z] is that equivalent point in xyz space
#	The Cartesian coordinates are curved downward from the origin.
def sphericalToCartesian(origin, geoCoord, approximate=True):
   cartesian = [0.0, 0.0, 0.0]

   kmPerDegree = 111.32

   # How far does this request extend? Maybe we can approximate.
   latLonSpan = math.hypot(geoCoord[LON] - origin[LON],
      geoCoord[LAT] - origin[LAT])
   #progress("sphericalToCartesian0 origin = {} geoCoord = {}"
   #   .format(origin, geoCoord))
   #progress("latLonSpan = {}".format(latLonSpan))

   if (approximate and latLonSpan < 0.1):
      # The following approximation is valid for lat-lon domains
      # spanning less than 0.1 degrees (about 10 km).
      cartesian[CART_X] = ((geoCoord[LON] - origin[LON])
         * kmPerDegree * math.cos(math.radians(geoCoord[LAT])))
      cartesian[CART_Y] = ((geoCoord[LAT] - origin[LAT])
         * kmPerDegree)
      cartesian[CART_Z] = (geoCoord[ALT] - origin[ALT])

      return(cartesian)

   # Algorithm steps, as viewed by an airplane flying due north at the origin.
   # 1. Rotate the geo-spherical coordinates together so that the airplane is heading straight north along 90° west longitude, with the Prime Meridian on my right and the International Date Line on my left. We are going to make the positive x-axis extending to the right, and the positive y-axis extending forward. The z-axis shall be up.
   geoCoord[LON] += -90.0 - origin[LON]
   origin[LON] = -90.0
   #progress("sphericalToCartesian1 origin = {} geoCoord = {}"
   #   .format(origin, geoCoord))

   # 2. Convert the geo-spherical coordinates to 3D Cartesian.
   # https://math.libretexts.org/Bookshelves/Calculus/Book%3A_Calculus_(OpenStax)/12%3A_Vectors_in_Space/12.7%3A_Cylindrical_and_Spherical_Coordinates
   radiusEarth = 6371	# earth mean radius in kilometers
   theta = math.radians(origin[LON])
   phi = math.radians(90.0 - origin[LAT])
   originCartesian = [
      (radiusEarth + origin[ALT]) * math.sin(phi) * math.cos(theta),
      (radiusEarth + origin[ALT]) * math.sin(phi) * math.sin(theta),
      (radiusEarth + origin[ALT]) * math.cos(phi)]
   theta = math.radians(geoCoord[LON])
   phi = math.radians(90.0 - geoCoord[LAT])
   coordCartesian = [
      (radiusEarth + geoCoord[ALT]) * math.sin(phi) * math.cos(theta),
      (radiusEarth + geoCoord[ALT]) * math.sin(phi) * math.sin(theta),
      (radiusEarth + geoCoord[ALT]) * math.cos(phi)]
   #progress("sphericalToCartesian2 originCartesian = {}\n\tcoordCartesian = {}"
   #   .format(originCartesian, coordCartesian))

   # 3. Rotate the system around the x-axis, so that my airplane ends up on the North Pole. Now we have a Cartesian coordinate system matching the 3D Collada models.
   # Use a rotation matrix:
   # https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
   theta = math.radians(origin[LAT] - 90.0)
   #theta = math.radians(-90.0)		# bogus, for testing simple rotation
   rotX = numpy.array([
      [1,               0,                0],
      [0, math.cos(theta), -math.sin(theta)],
      [0, math.sin(theta),  math.cos(theta)]
   ])
   #progress("rotX = {}".format(rotX))

   originArray = numpy.array(originCartesian)
   originRotated = numpy.matmul(rotX, originArray)

   coordArray = numpy.array(coordCartesian)
   coordRotated = numpy.matmul(rotX, coordArray)
   #progress("sphericalToCartesian3 originRotated = {}\n\tcoordRotated = {}"
   #   .format(originRotated, coordRotated))

   # 4. Subtract the radius of the earth from all z-coordinates. This will move the Cartesian origin to the airplane itself, and all other coordinates will now be expressed relative to the airplane in xyz space (right-hand rule preserved).
   originRotated[CART_Z] -= radiusEarth
   coordRotated[CART_Z] -= radiusEarth
   #progress("sphericalToCartesian4 originRotated = {}\n\tcoordRotated = {}"
   #   .format(originRotated, coordRotated))

   return(coordRotated)

