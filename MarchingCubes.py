#!/usr/bin/env python3

# MarchingCubes.py
# Python module to implement Marching Cubes algorithm for 3D isosurface.
#
# Author: Carl Drews, Atmospheric Chemistry Observations & Modeling
# Created: April 2021
# Copyright 2021 by the University Corporation for Atmospheric Research



import sys
import math
import numpy

import utilsCollada		# for converting geo-coordinates
import copy



# Display progress message on the console.
# endString = set this to '' for no return
def progress(message, endString='\n'):
   if (True):   # disable here in production
      print(message, end=endString)
      sys.stdout.flush()



# A data point in 3D space.
class Vertex:
   # Create Vertex at a certain location.
   def __init__(self, xLocation=0.0, yLocation=0.0, zLocation=0.0):
      self.x = xLocation
      self.y = yLocation
      self.z = zLocation
      self.value = 0.0
      return

   def toString(self):
      return("x:{:06f} y:{:06f} z:{:06f} val:{:06f}"
         .format(self.x, self.y, self.z, self.value))

   # Return array of x,y,z coordinates.
   def getCoordinates(self):
      return([self.x, self.y, self.z])

   # accessors for the chemical value
   def setValue(self, chemValue):
      self.value = chemValue
      return(self)	# so we can chain constructor with setValue()
   def getValue(self):
      return(self.value)

   # Convert absolute lat-lon to km relative to 3D model origin.
   # Model origin is probably at the center of the WRF-Chem domain.
   # originVertex = center of the model in earth coordinates (lat, lon, alt)
   def latLonToKm(self, originVertex):
      kmVertex = utilsCollada.sphericalToCartesian(
         [originVertex.y, originVertex.x, originVertex.z],
         [self.y, self.x, self.z])

      self.x = kmVertex[utilsCollada.CART_X]
      self.y = kmVertex[utilsCollada.CART_Y]
      self.z = kmVertex[utilsCollada.CART_Z]

      return



# A grid cell in 3D space. The coordinates are integers.
class IntVertex(Vertex):
   # Create integer Vertex at a certain location (column, row, level).
   def __init__(self, column=0, row=0, level=0):
      super().__init__(int(column), int(row), int(level))
      return



# Triangular surface in 3D space.
class Triangle:

   def __init__(self, vertex1, vertex2, vertex3):
      self.vertices = [vertex1, vertex2, vertex3]

      # set up default normal vectors of zero length
      upNormal = Vertex(0.0, 0.0, 0.0)
      self.normals = [upNormal, upNormal, upNormal]
      return

   # specify the normal vectors in the same order as the vertices
   def setNormals(self, normal1, normal2, normal3):
      self.normals = [normal1, normal2, normal3]
      return(self)	# so we can chain constructor with setNormals()
 
   def toString(self):
      myString = ""
      for vertex in self.vertices:
         myString += "Vertex {}\n".format(vertex.toString())
      for vertex in self.normals:
         myString += "Normal {}\n".format(vertex.toString())
      return(myString)


   # Return array of vertex x,y,z coordinates (9 floats).
   def getVertexCoordinates(self):
      myVertices = []
      for vertex in self.vertices:
         myVertices.extend(vertex.getCoordinates())
      return(myVertices)


   # Return array of normal vector x,y,z coordinates (9 floats).
   def getNormalCoordinates(self):
      myNormals = []
      for normal in self.normals:
         myNormals.extend(normal.getCoordinates())
      return(myNormals)


   # Multiple the height to match Google Earth terrain exaggeration.
   def exaggerateHeight(self, altitudeFactor):
      for vertex in self.vertices:
         vertex.z *= altitudeFactor
      return


   # Convert absolute lat-lon to km relative to 3D model origin.
   # Model origin is probably at the center of the WRF-Chem domain.
   # originVertex = center of the model in earth coordinates (lat, lon, alt)
   def latLonToKm(self, originVertex):
      for vertex in self.vertices:
         kmVertex = utilsCollada.sphericalToCartesian(
            [originVertex.y, originVertex.x, originVertex.z],
            [vertex.y, vertex.x, vertex.z])

         vertex.x = kmVertex[utilsCollada.CART_X]
         vertex.y = kmVertex[utilsCollada.CART_Y]
         vertex.z = kmVertex[utilsCollada.CART_Z]

      return



# Some component within the marching cube (Corner, Face, or EdgeCrossing).
class Component:
   def __init__(self, column=0, row=0, level=0):
      self.cell = IntVertex(column, row, level)
      return

   def toString(self):
      return("cell {}".format(self.cell.toString()))

   # return list of cube coordinates ordered [column, row, level]
   def getCoords(self):
      return([self.cell.x, self.cell.y, self.cell.z])



# Corner component of the marching cube.
class Corner(Component):
   def __init__(self, column=0, row=0, level=0,
      xLocation=0.0, yLocation=0.0, zLocation=0.0):
      super().__init__(column, row, level)

      # set the location in 3D floating-point space
      self.location = Vertex(xLocation, yLocation, zLocation)

   def toString(self):
      return("Corner: {}\n\t{}"
         .format(super().toString(), self.location.toString()))

   # accessors for chemical value stored in the location
   def setValue(self, chemValue):
      self.location.setValue(chemValue)
      return(self)	# so we can chain constructor with setValue()
   def getValue(self):
      return(self.location.getValue())

   # set floating-point location for this corner in 3D space
   def setLocation(self, xLocation, yLocation, zLocation):
      self.location.x = xLocation
      self.location.y = yLocation
      self.location.z = zLocation
      return



# Face of one of the marching cubes.
# Each cube face will have a linear contour plot in two colors.
# Interior lines are line segments drawn across the face
# to separate different colors of the face's contour plot.
# The contour lines connect one edge crossing to another.
# The algorithm seeks to minimize the total length of interior lines.
class Face(Component):

   def __init__(self, column=0, row=0, level=0):
      super().__init__(column, row, level)

      # how many valid edge crossings for this face
      self.crossingCount = 0
      # total lenth of all interior lines
      self.lengthInteriorLines = 0.0
      return

   def toString(self):
      return("Face: {} edge crossings and {} interior lines.".
         format(self.crossingCount, self.lengthInteriorLines))



# An edge crossing occurs along the edges of a cube
# where two adjacent vertices are less than and greater
# than the desired isoValue. The edge crossing is located
# by linear interpolation between the two vertex values.
# The isosurface will pass through this edge crossing,
# and the crossing site will form one vertex of a triangle.
# If the isosurface does not cross this edge, then crossAtVertex = False.
class EdgeCrossing(Component):

   # construct new EdgeCrossing at the specified vertex
   def __init__(self, column=0, row=0, level=0,
      xLoc=0.0, yLoc=0.0, zLoc=0.0):
      super().__init__(column, row, level)

      self.surfaceCrossesHere = False			# cube corners bracket the isosurface
      self.crossAtVertex = Vertex(xLoc, yLoc, zLoc)	# location between the cube corners
      self.calculated = False				# crossAtVertex has been calculated
      self.neighbors = []			# 0 or 2 other edge crossings on adjacent faces
      self.visited = False			# has been examined during crawl
      self.outside = Vertex()		# orthogonal unit vector; points outside surface
					# toward region of lower chemical concentration

   def toString(self):
      return("EdgeCrossing {}\n\toutside {}"
         .format(self.crossAtVertex.toString(), self.outside.toString()))

   def setLocation(self, newVertex):
      self.crossAtVertex.x = newVertex.x
      self.crossAtVertex.y = newVertex.y
      self.crossAtVertex.z = newVertex.z

   def setOutsideVector(self, xDelta, yDelta, zDelta):
      self.outside.x = xDelta
      self.outside.y = yDelta
      self.outside.z = zDelta



# Calculate the hypotenuse distance between two vertices.
def vertexDistance(vertex1, vertex2):
   deltaX = vertex2.x - vertex1.x
   deltaY = vertex2.y - vertex1.y
   deltaZ = vertex2.z - vertex1.z

   return(math.hypot(math.hypot(deltaX, deltaY), deltaZ))



# Determine if triangles made from vertexList would
# point in roughly same direction as outside normal.
# We use cube coordinates here instead of geo-coordinates
# or Cartesian to remove: Lambert conformal grid rotation,
# variation on the y-axes from terrain following,
# and loss of accuracy due to very thin atmospheric layers.
# crossingList = three or more edge crossings
# return True if vertexList should be reversed
def needsReversal(crossingList):

   # calculate cross product normal of potential triangle
   futureTriangle = Triangle(crossingList[0].cell,
      crossingList[1].cell, crossingList[2].cell)
   #progress("futureTriangle = {}".format(futureTriangle.toString()))
   triDirection = normalVector(futureTriangle)
   #progress("triDirection = {}".format(triDirection.toString()))

   # compare that vector with the valid component of crossing
   refCrossing = crossingList[0]
   #progress("refCrossing.outside = {}".format(refCrossing.outside.toString()))
   sameSign = (refCrossing.outside.x * triDirection.x
      + refCrossing.outside.y * triDirection.y
      + refCrossing.outside.z * triDirection.z)
   #progress("sameSign = {}".format(sameSign))

   return(sameSign < 0.0)	# we want same sign/direction



# Divide polygon into multiple triangles,
# all anchored at the center point of the polygon.
# return list of triangles
def traversePolygon(vertexList):
   # find the average location in 3D space
   anchorVertex = Vertex()
   anchorCount = 0
   for vertex in vertexList:
      anchorVertex.x += vertex.x
      anchorVertex.y += vertex.y
      anchorVertex.z += vertex.z
      anchorCount += 1

   anchorVertex.x /= anchorCount
   anchorVertex.y /= anchorCount
   anchorVertex.z /= anchorCount
   #progress("anchorVertex = {}".format(anchorVertex.toString()))

   # walk around the polygon perimeter
   triList = []
   previousVertex = None
   for myVertex in vertexList:
      if (previousVertex is None):
         previousVertex = myVertex
         continue

      # make perimeter triangle from achor and two polygon vertices
      myTriangle = Triangle(copy.deepcopy(anchorVertex),
         copy.deepcopy(previousVertex), copy.deepcopy(myVertex))
      #progress("Perimeter triangle = {}".format(myTriangle.toString()))
      triList.append(myTriangle)
      previousVertex = myVertex

   # close out with the final connecting triangle
   myTriangle = Triangle(copy.deepcopy(anchorVertex),
      copy.deepcopy(previousVertex), copy.deepcopy(vertexList[0]))
   #progress("Perimeter triangle final = {}".format(myTriangle.toString()))
   triList.append(myTriangle)

   return(triList)


# Calculate the normal vector of a triangle.
# Method: Take cross product of the two adjacent sides.
# The order of the vertices used in the calculation
# will affect the direction of the normal direction
# (in or out of the face with respect to winding).
# Carl Drews - May 18, 2021
def normalVector(myTriangle):
   # extract vectors of the adjacent sides
   aSide = Vertex(myTriangle.vertices[1].x - myTriangle.vertices[0].x,
      myTriangle.vertices[1].y - myTriangle.vertices[0].y,
      myTriangle.vertices[1].z - myTriangle.vertices[0].z)
   bSide = Vertex(myTriangle.vertices[2].x - myTriangle.vertices[0].x,
      myTriangle.vertices[2].y - myTriangle.vertices[0].y,
      myTriangle.vertices[2].z - myTriangle.vertices[0].z)

   # take the cross product
   normalX = aSide.y * bSide.z - aSide.z * bSide.y
   normalY = aSide.z * bSide.x - aSide.x * bSide.z
   normalZ = aSide.x * bSide.y - aSide.y * bSide.x

   # convert the vector to unit length
   vectorLength = math.hypot(math.hypot(normalX, normalY), normalZ)
   normalX /= vectorLength
   normalY /= vectorLength
   normalZ /= vectorLength

   return(Vertex(normalX, normalY, normalZ))



# Data structure for one grid cell within the Marching Cubes.
# This class contains the components to generate the list
# of triangles that form the isosurface within the cube.
# The cube need not be of uniform edge length, or even
# strictly orthogonal. Think of it as a mostly-rectangular
# slab oriented nearly horizontally in space.
# The Cube object stores the slab components in a 3x3x3
# mixed-type list(s) of Corners, Faces, and EdgeCrossings.
# The components are referenced by their 3x3x3 coordinates:
#	myCube[z][y][x]
class Cube:

   def __init__(self):
      # build the 27 cells of our 3D grid
      self.structure = [
         # array is built from the bottom upwards
         [
            # bottom face
            [Corner(0,0,0),       EdgeCrossing(1,0,0), Corner(2,0,0)],
            [EdgeCrossing(0,1,0), Face(1,1,0),         EdgeCrossing(2,1,0)],
            [Corner(0,2,0),       EdgeCrossing(1,2,0), Corner(2,2,0)]
         ],

         # center slice
         [
            [EdgeCrossing(0,0,1), Face(1,0,1), EdgeCrossing(2,0,1)],
            [Face(0,1,1),         Face(1,1,1), Face(2,1,1)],	# unused face in the center of cube
            [EdgeCrossing(0,2,1), Face(1,2,1),  EdgeCrossing(2,2,1)]
         ],

         # top face
         [
            [Corner(0,0,2),       EdgeCrossing(1,0,2), Corner(2,0,2)],
            [EdgeCrossing(0,1,2), Face(1,1,2),         EdgeCrossing(2,1,2)],
            [Corner(0,2,2),       EdgeCrossing(1,2,2), Corner(2,2,2)]
         ]
      ]
      return


   # Display text version of cube for debugging.
   def toString(self):
      cubeStr = "Cube;, top to bottom, front to back, left to right:\n"
      # display from top to bottom
      for z in range(2, -1, -1):
         for y in range(3):
            #cubeStr += "{}\n".format(self.structure[z][y])
            for x in range(3):
               cubeStr += "{} ".format(self.structure[z][y][x].toString())
            cubeStr += "\n"

         cubeStr += "----------------------\n"
      return(cubeStr)


   # Assign XYZ-locations in 3D space for all the vertices.
   # vertexLocX = 2x2x2 numpy array of X coordinates
   # vertexLocY = 2x2x2 numpy array Y coordinates
   # vertexLocZ = 2x2x2 numpy array Z coordinates
   def setVerticesXYZ(self, vertexLocX, vertexLocY, vertexLocZ):
      for z in range(2):
         for y in range(2):
            for x in range(2):
               zIndex = z * 3	# skip past the edge crossing
               yIndex = y * 3
               xIndex = x * 3

               self.structure[z * 2][y * 2][x * 2].setLocation(
                  vertexLocX[z,y,x], vertexLocY[z,y,x], vertexLocZ[z,y,x])
      return


   # Assign chemical values for all the vertices.
   # chemValues = 2x2x2 numpy array of float values
   def setVertexValues(self, chemValues):
      for z in range(2):
         for y in range(2):
            for x in range(2):
               zIndex = z * 3	# skip past the edge crossing
               yIndex = y * 3
               xIndex = x * 3

               self.structure[z * 2][y * 2][x * 2].setValue(chemValues[z,y,x])
      return


   # Calculate one edge crossing between adjacent cube corners.
   # chemVertex1, chemVertex2 = cube corners on both ends of the edge
   # myCrossing = the edge crossing between chemical vertices
   # crossValue = looking for this isosurface chemical value
   def calculateOneCrossing(self, chemVertex1, chemVertex2,
      myCrossing, crossValue):

      chemValue1 = chemVertex1.getValue()
      chemValue2 = chemVertex2.getValue()
      #progress("chemValues = {} and {}".format(chemValue1, chemValue2))

      if ((chemValue1 < crossValue and chemValue2 < crossValue)
         or (chemValue1 >= crossValue and chemValue2 >= crossValue)):
         # there is no value crossing along this edge
         #progress("\tNo crossing.")
         myCrossing.surfaceCrossesHere = False
         myCrossing.calculated = True
         return

      #progress("\tchemVertices = {} and {}"
      #   .format(chemVertex1.toString(), chemVertex2.toString()))

      # perform linear interpolation according to chemical value
      fraction1 = (crossValue - chemValue1) / (chemValue2 - chemValue1)
      #progress("\tfraction1 = {}".format(fraction1))
      crossX = chemVertex1.location.x + fraction1 * (chemVertex2.location.x - chemVertex1.location.x)
      crossY = chemVertex1.location.y + fraction1 * (chemVertex2.location.y - chemVertex1.location.y)
      crossZ = chemVertex1.location.z + fraction1 * (chemVertex2.location.z - chemVertex1.location.z)

      myCrossing.setLocation(Vertex(crossX, crossY, crossZ).setValue(crossValue))
      myCrossing.surfaceCrossesHere = True
      myCrossing.calculated = True

      # calculate the normal vector pointing outside to low chem values
      coords1 = chemVertex1.getCoords()
      coords2 = chemVertex2.getCoords()
      crossingCoords = [(a + b)/2 for (a,b) in zip(coords1, coords2)]
      if (chemValue1 <= crossValue):
         deltaVector = [(b - a) for (a,b) in zip(crossingCoords, coords1)]
      else:
         deltaVector = [(b - a) for (a,b) in zip(crossingCoords, coords2)]
      myCrossing.setOutsideVector(deltaVector[0], deltaVector[1], deltaVector[2])
      #progress("\tmyCrossing = {}".format(myCrossing.toString()))

      return


   # Calculate locations for all 12 the possible edge crossings.
   # Mark them as calculated, and indicate if no actual crossing there.
   def calculateEdgeCrossings(self, isoValue):
      #progress("isoValue = {}".format(isoValue))

      # bottom face
      self.calculateOneCrossing(self.structure[0][0][0], self.structure[0][2][0],
         self.structure[0][1][0], isoValue)
      self.calculateOneCrossing(self.structure[0][2][0], self.structure[0][2][2],
         self.structure[0][2][1], isoValue)
      self.calculateOneCrossing(self.structure[0][2][2], self.structure[0][0][2],
         self.structure[0][1][2], isoValue)
      self.calculateOneCrossing(self.structure[0][0][2], self.structure[0][0][0],
         self.structure[0][0][1], isoValue)
      #progress("")

      # around the middle slice
      self.calculateOneCrossing(self.structure[0][0][0], self.structure[2][0][0],
         self.structure[1][0][0], isoValue)
      self.calculateOneCrossing(self.structure[0][2][0], self.structure[2][2][0],
         self.structure[1][2][0], isoValue)
      self.calculateOneCrossing(self.structure[0][2][2], self.structure[2][2][2],
         self.structure[1][2][2], isoValue)
      self.calculateOneCrossing(self.structure[0][0][2], self.structure[2][0][2],
         self.structure[1][0][2], isoValue)
      #progress("")


      # top face
      self.calculateOneCrossing(self.structure[2][0][0], self.structure[2][2][0],
         self.structure[2][1][0], isoValue)
      self.calculateOneCrossing(self.structure[2][2][0], self.structure[2][2][2],
         self.structure[2][2][1], isoValue)
      self.calculateOneCrossing(self.structure[2][2][2], self.structure[2][0][2],
         self.structure[2][1][2], isoValue)
      self.calculateOneCrossing(self.structure[2][0][2], self.structure[2][0][0],
         self.structure[2][0][1], isoValue)
      #progress("")

      return


   # Create connections between edge crossings around one face.
   # cubeFace = one (center) face in the cube
   # edges[4] = set of edge crossings around that face
   def resolveOneFace(self, cubeFace, edges):
      # create sub-list of valid crossings
      crossList = []
      for edge in edges:
         if (edge.surfaceCrossesHere):
            crossList.append(edge)

      cubeFace.crossingCount = len(crossList)
      if (cubeFace.crossingCount == 0):
         #progress("cubeFace = {}".format(cubeFace.toString()))
         return		# no connections to make

      if (cubeFace.crossingCount == 2):
         # connect the only two crossings with each other
         crossList[0].neighbors.append(crossList[1])
         crossList[1].neighbors.append(crossList[0])
         cubeFace.lengthInteriorLines = vertexDistance(
            crossList[0].crossAtVertex, crossList[1].crossAtVertex)
         #progress("cubeFace = {}".format(cubeFace.toString()))
         return

      if (cubeFace.crossingCount == 4):
	 # determine which of the four edge crossings should be connected
         #progress("Note: Cube face {} has {} edge crossings."
         #   .format(cubeFace, cubeFace.crossingCount))
         # measure the four interior lines
         lengths = numpy.zeros(cubeFace.crossingCount)
         for ci, cross in enumerate(crossList):
            if (ci < len(crossList) -1):
               nextCross = crossList[ci + 1]
            else:
               nextCross = crossList[0]
            #progress("cross = {}".format(cross.toString()))
            #progress("nextCross = {}".format(nextCross.toString()))
            lengths[ci] = vertexDistance(cross.crossAtVertex, nextCross.crossAtVertex)

         #progress("\tContour lengths = {}".format(lengths))

         # find shortest total length of interior lines
         dist01and23 = lengths[0] + lengths[2]
         dist12and30 = lengths[1] + lengths[3]
         #progress("\tdist01and23 = {} dist12and30 = {}".format(dist01and23, dist12and30))
         if (dist01and23 <= dist12and30):
            #progress("\tConnect 0 with 1, 2 with 3")
            crossList[0].neighbors.append(crossList[1])
            crossList[1].neighbors.append(crossList[0])
            crossList[2].neighbors.append(crossList[3])
            crossList[3].neighbors.append(crossList[2])
            cubeFace.lengthInteriorLines = dist01and23
         else:
            #progress("\tConnect 1 with 2, 3 with 0")
            crossList[1].neighbors.append(crossList[2])
            crossList[2].neighbors.append(crossList[1])
            crossList[3].neighbors.append(crossList[0])
            crossList[0].neighbors.append(crossList[3])
            cubeFace.lengthInteriorLines = dist12and30

         #progress("cubeFace = {}".format(cubeFace.toString()))
         return

      progress("Warning: Cube face {} has odd number {} edge crossings."
         .format(cubeFace, cubeFace.crossingCount))

      return


   # Connect together the calculated edge crossings on each face.
   # The only ambiguity arises when there are four crossings.
   # Resolve that ambiguity by minimizing total length of
   # isosurface lines drawn across that face. When this routine
   # is complete, all cube faces will be divided into regions
   # below and above the isovalue.
   def resolveFaceContours(self):
      # bottom face
      # The indexes here are: X, Y, Z
      self.resolveOneFace(self.structure[1][1][0],
         [self.structure[0][1][0], self.structure[1][2][0],
         self.structure[2][1][0], self.structure[1][0][0]])

      # clockwise around the middle
      self.resolveOneFace(self.structure[0][1][1],
         [self.structure[0][1][0], self.structure[0][2][1],
         self.structure[0][1][2], self.structure[0][0][1]])
      self.resolveOneFace(self.structure[1][2][1],
         [self.structure[1][2][0], self.structure[2][2][1],
         self.structure[1][2][2], self.structure[0][2][1]])
      self.resolveOneFace(self.structure[2][1][1],
         [self.structure[2][1][0], self.structure[2][0][1],
         self.structure[2][1][2], self.structure[2][2][1]])
      self.resolveOneFace(self.structure[1][0][1],
         [self.structure[1][0][0], self.structure[0][0][1],
         self.structure[1][0][2], self.structure[2][0][1]])

      # top face
      self.resolveOneFace(self.structure[1][1][2],
         [self.structure[0][1][2], self.structure[1][2][2],
         self.structure[2][1][2], self.structure[1][0][2]])

      return


   # Verify that all edge crossings have two neighbor crossings.
   def checkCrossingNeighbors(self):
      retValue = True
      for z in range(3):
         for y in range(3):
            for x in range(3):
               part = self.structure[z][y][x]
               if ("EdgeCrossing" not in str(type(part))):
                  continue
               if (not part.surfaceCrossesHere):
                  continue

               #progress("Checking edge crossing {}".format(part))
               if (len(part.neighbors) != 2):
                   #progress("Error: EdgeCrossing {} has {} neighbors."
                   #   .format(part, len(part.neighbors)))
                   retValue = False

      if (not retValue):
         progress("Warning: checkCrossingNeighbors() = {}".format(retValue))
      return(retValue)


   # Crawl from this edge crossing around the polygon
   # and back home. Split the polygon into triangles.
   # Calculate normal vectors for every new triangle.
   # Mark all traversed edge crossings as visited.
   # homeCrossing = begin crawl here
   # gridLines = create polygons around the edge of each interior surface
   # return list of Triangles around the crawl
   def crawlOneEdgeCrossing(self, homeCrossing, gridLines = True):
      triList = []
      lineList = []

      # step around the perimeter to make a polygon
      myCrossing = homeCrossing
      polygon = []
      referenceCrossings = []
      while (True):
         # record this crossing vertex into the polygon list
         polygon.append(myCrossing.crossAtVertex)
         myCrossing.visited = True

         # collect reference crossings
         referenceCrossings.append(myCrossing)

         # move around to the next neighbor crossing
         nextCrossing = myCrossing.neighbors[0]
         if (nextCrossing.visited):
            nextCrossing = myCrossing.neighbors[1]
         if (nextCrossing.visited):
            break

         myCrossing = nextCrossing

      #progress("polygon = {}".format(polygon))
      #for poly in polygon:
      #   progress("\tx={} y={} z={}".format(poly.x, poly.y, poly.z))

      if (len(polygon) > 6):
         # this case represents an extreme saddle
         progress("Note: polygon beginning at {} has {} points."
            .format(polygon[0].toString(), len(polygon)))

      # make sure normal vectors will point outside the plume
      if needsReversal(referenceCrossings):
         #progress("Perform reversal on polygon.")
         polygon.reverse()
      else:
         #progress("No need for reversal.")
         pass

      if (gridLines):
         # record this perimeter, making sure to return to first vertex
         # Note: This step will actually create twice as many grid
         # lines as needed for the wireframe. This is because each
         # line segment (of the perimeter) exactly matches another line
         # segment on the adjacent cube face. One could optimize the
         # line segments by only creating them on the three cube faces
         # closest to the origin (and the final outer faces of the domain).
         # Carl Drews - December 21, 2022
         perimeter = copy.deepcopy(polygon)
         perimeter.append(copy.deepcopy(polygon[0]))
         lineList.append(perimeter)

      if (len(polygon) == 3):
         # simple case - take those three vertices
         oneTriangle = Triangle(polygon[0], polygon[1], polygon[2])
         triList.append(oneTriangle)
         return(triList, lineList)

      # divide polygon into several triangles anchored at center point
      polyTriangles = traversePolygon(polygon)
      triList.extend(polyTriangles)
      return(triList, lineList)


   # Follow neighbors of all edge crossings. Build list of
   # triangles including their normal vectors.
   # gridLines = create polygons around the edge of each interior surface
   # return list of Triangles for this cube, and list of 3D polygons
   def crawlEdgeCrossings(self, gridLines=True):
      cubeTriangles = []
      cubeLines = []
      for z in range(3):
         for y in range(3):
            for x in range(3):
               part = self.structure[z][y][x]
               if ("EdgeCrossing" not in str(type(part))):
                  continue
               if (not part.surfaceCrossesHere):
                  continue
               if (part.visited):
                  continue

               partTriangles, partLines = self.crawlOneEdgeCrossing(part, gridLines)
               cubeTriangles.extend(partTriangles)
               cubeLines.extend(partLines)

      return(cubeTriangles, cubeLines)


   # Calculate and return isosurface within this cube.
   # surfaceValue = build surface along this chemical concentration
   # gridLines = create polygons around the edge of each interior surface
   # return list of triangles
   def calculateIsosurface(self, surfaceValue, gridLines=True):
      triList = []
      lineList = []

      if (False):
         # test with a simple slanted triangle
         v1 = Vertex(self.structure[0][0][0].x, self.structure[0][0][0].y,
            self.structure[0][0][0].z)
         v2 = Vertex(self.structure[2][2][0].x, self.structure[2][2][0].y,
            self.structure[2][2][0].z)
         v3 = Vertex(self.structure[2][0][2].x, self.structure[2][0][2].y,
            self.structure[2][0][2].z)

         tri = Triangle(v1, v2, v3)
         triList.append(tri)
         return(triList)

      # calculate the edge crossings
      self.calculateEdgeCrossings(surfaceValue)

      # resolve and connect contour lines on each face
      self.resolveFaceContours()
      self.checkCrossingNeighbors()

      # crawl the edge crossings to create triangles and contour lines
      triangles, lines = self.crawlEdgeCrossings(gridLines)
      triList.extend(triangles)
      lineList.extend(lines)

      return(triList, lineList)



# Calculate list of triangles representing isosurface for this grid cell.
# This grid cell has been determined to contain the isosurface.
# This routine is your main entry point into the isosurface calculation.
# xLocation[y,x], yLocation[y,x] = positions of all cells in 3D space
#	The data columns are straight up and down (no leaning).
# zLocation[z,y,x] = altitude of all cells in 3D space
#	The horizontal slices may vary in altitude and height.
# chemValues[z,y,x] = chemical values in the entire 3D space
# xIndex, yIndex, zIndex = lower southwest corner
# surfaceValue = create triangulated surface along this value
# gridLines = create polygons around the edge of each interior surface
#		These will become contour lines along each grid cell.
# return list of one or more triangles, and list of one or more perimeter lines
def renderCube(xLocation, yLocation, zLocation, chemValues,
   xIndex, yIndex, zIndex, surfaceValue,
   gridLines=True):

   # verify what we received
   #progress("xLocation shape = {}".format(xLocation.shape))
   #progress("yLocation shape = {}".format(yLocation.shape))
   #progress("zLocation shape = {}".format(zLocation.shape))
   #progress("zLocation = {}".format(zLocation[:, 10, 20]))

   #progress("chemValues shape = {}".format(chemValues.shape))

   # construct the cube
   gridCell = Cube()

   # assign all the vertex locations
   xTemp = numpy.zeros([2,2,2])
   xTemp[0,:,:] = xLocation[yIndex:yIndex+2, xIndex:xIndex+2]
   xTemp[1,:,:] = xTemp[0]	# data columns are strictly vertical
   #progress("xTemp = {}".format(xTemp))

   yTemp = numpy.zeros([2,2,2])
   yTemp[0,:,:] = yLocation[yIndex:yIndex+2, xIndex:xIndex+2]
   yTemp[1,:,:] = yTemp[0]
   #progress("yTemp = {}".format(yTemp))

   zTemp = numpy.zeros([2,2,2])
   zTemp[:,:,:] = zLocation[zIndex:zIndex+2, yIndex:yIndex+2, xIndex:xIndex+2]
   #progress("zTemp = {}".format(zTemp))

   if (False):
      # create small gaps between grid cells
      gapFractionHorz = 1.0 / 100.0
      gapFractionVert = 1.0 / 20.0
      xSpan = xTemp[0,0,1] - xTemp[0,0,0]
      ySpan = yTemp[0,1,0] - yTemp[0,0,0]
      zSpan = zTemp[1,0,0] - zTemp[0,0,0]

      xMargin = xSpan * gapFractionHorz
      yMargin = ySpan * gapFractionHorz
      zMargin = zSpan * gapFractionVert

      # create a slightly smaller cube in all dimensions
      xTemp[:,:,0] += xMargin
      xTemp[:,:,1] -= xMargin
      yTemp[:,0,:] += yMargin
      yTemp[:,1,:] -= yMargin
      zTemp[0,:,:] += zMargin
      zTemp[1,:,:] -= zMargin

   gridCell.setVerticesXYZ(xTemp, yTemp, zTemp)

   # assign all the vertex chemical values
   vTemp = numpy.zeros([2,2,2])
   vTemp[:,:,:] = chemValues[zIndex:zIndex+2, yIndex:yIndex+2, xIndex:xIndex+2]
   #vTemp[1, 0, 1] = 234.5678		# bogus - force four edge crossings, multiply by 10 to span entire face
   #vTemp[1, 1, 0] = 234.5678 * 10	# bogus
   #vTemp[0, 1, 1] = 189.5678		# bogus
   #progress("vTemp = {}".format(vTemp))
   gridCell.setVertexValues(vTemp)

   #progress("gridCell cube = {}".format(gridCell.toString()))

   # calculate and retrieve triangulated surface at isoValue
   triList, lineList = gridCell.calculateIsosurface(surfaceValue, gridLines)

   return(triList, lineList)

