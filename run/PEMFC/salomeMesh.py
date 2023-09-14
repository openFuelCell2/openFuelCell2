#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.7.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'./')

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
Vertex_1 = geompy.MakeVertex(0, -0.001, -0.00225)
Vertex_2 = geompy.MakeVertex(0, -0.001, -0.00125)
Vertex_3 = geompy.MakeVertex(0, -0.001, -0.00025)
Vertex_4 = geompy.MakeVertex(0, -0.001, -5e-05)
Vertex_5 = geompy.MakeVertex(0, -0.001, -3e-05)
Vertex_6 = geompy.MakeVertex(0, -0.001, -1e-05)
Vertex_7 = geompy.MakeVertex(0, -0.001, 1e-05)
Vertex_8 = geompy.MakeVertex(0, -0.001, 3e-05)
Vertex_9 = geompy.MakeVertex(0, -0.001, 5e-05)
Vertex_10 = geompy.MakeVertex(0, -0.001, 0.00025)
Vertex_11 = geompy.MakeVertex(0, -0.001, 0.00125)
Vertex_12 = geompy.MakeVertex(0, -0.001, 0.00225)
Line_1 = geompy.MakeLineTwoPnt(Vertex_1, Vertex_2)
Line_2 = geompy.MakeLineTwoPnt(Vertex_2, Vertex_3)
Line_3 = geompy.MakeLineTwoPnt(Vertex_3, Vertex_4)
Line_4 = geompy.MakeLineTwoPnt(Vertex_4, Vertex_5)
Line_5 = geompy.MakeLineTwoPnt(Vertex_5, Vertex_6)
Line_6 = geompy.MakeLineTwoPnt(Vertex_6, Vertex_7)
Line_7 = geompy.MakeLineTwoPnt(Vertex_7, Vertex_8)
Line_8 = geompy.MakeLineTwoPnt(Vertex_8, Vertex_9)
Line_9 = geompy.MakeLineTwoPnt(Vertex_9, Vertex_10)
Line_10 = geompy.MakeLineTwoPnt(Vertex_10, Vertex_11)
Line_11 = geompy.MakeLineTwoPnt(Vertex_11, Vertex_12)
Partition_1 = geompy.MakePartition([Line_1, Line_2, Line_3, Line_4, Line_5, Line_6, Line_7, Line_8, Line_9, Line_10, Line_11], [], [], [], geompy.ShapeType["EDGE"], 0, [], 0)
Extrusion_1 = geompy.MakePrismDXDYDZ(Partition_1, 0, 0.002, 0)
Vertex_13 = geompy.MakeVertexWithRef(Vertex_1, 0, 0.0005, 0)
Vertex_14 = geompy.MakeVertexWithRef(Vertex_1, 0, 0.0015, 0)
Extrusion_1_vertex_78 = geompy.GetSubShape(Extrusion_1, [78])
Vertex_15 = geompy.MakeVertexWithRef(Extrusion_1_vertex_78, 0, 0.0005, 0)
Vertex_16 = geompy.MakeVertexWithRef(Extrusion_1_vertex_78, 0, 0.0015, 0)
Line_12 = geompy.MakeLineTwoPnt(Vertex_15, Vertex_13)
Line_13 = geompy.MakeLineTwoPnt(Vertex_14, Vertex_16)
Partition_2 = geompy.MakePartition([Extrusion_1, Line_12, Line_13], [], [], [], geompy.ShapeType["FACE"], 0, [], 0)

# Geometry starts here
Geo = geompy.MakePrismDXDYDZ(Partition_2, 0.04, 0, 0)

# boundaries
interconnect0 = geompy.CreateGroup(Geo, geompy.ShapeType["FACE"])
geompy.UnionIDs(interconnect0, [632, 615, 649])
interconnect1 = geompy.CreateGroup(Geo, geompy.ShapeType["FACE"])
geompy.UnionIDs(interconnect1, [76, 52, 28])
interconnectSides = geompy.CreateGroup(Geo, geompy.ShapeType["FACE"])
geompy.UnionIDs(interconnectSides, [550, 608, 4, 86, 598, 69, 656, 134, 604, 106, 140, 34, 662, 570, 82, 58, 628, 645, 643, 660, 104, 626, 602, 568, 56, 138, 80, 32])
airInlet = geompy.CreateGroup(Geo, geompy.ShapeType["FACE"])
geompy.UnionIDs(airInlet, [585])
fuelInlet = geompy.CreateGroup(Geo, geompy.ShapeType["FACE"])
geompy.UnionIDs(fuelInlet, [121])
airOutlet = geompy.CreateGroup(Geo, geompy.ShapeType["FACE"])
geompy.UnionIDs(airOutlet, [587])
fuelOutlet = geompy.CreateGroup(Geo, geompy.ShapeType["FACE"])
geompy.UnionIDs(fuelOutlet, [123])
cathodeSides = geompy.CreateGroup(Geo, geompy.ShapeType["FACE"])
geompy.UnionIDs(cathodeSides, [540, 434, 482, 424, 492, 376, 428, 529, 544, 546, 396, 394, 488, 452, 527, 469, 430, 413, 454, 510, 411, 471, 486, 512])
anodeSides = geompy.CreateGroup(Geo, geompy.ShapeType["FACE"])
geompy.UnionIDs(anodeSides, [222, 181, 295, 239, 164, 196, 312, 314, 198, 220, 297, 254, 256, 162, 179, 237, 278, 280, 308, 202, 192, 250, 144, 260])
electrolyteSides = geompy.CreateGroup(Geo, geompy.ShapeType["FACE"])
geompy.UnionIDs(electrolyteSides, [355, 372, 338, 318, 366, 370, 353, 336])

# zones
interconnect = geompy.CreateGroup(Geo, geompy.ShapeType["SOLID"])
geompy.UnionIDs(interconnect, [630, 548, 36, 2, 647, 60, 606, 84, 125, 589])
cchannel = geompy.CreateGroup(Geo, geompy.ShapeType["SOLID"])
geompy.UnionIDs(cchannel, [572])
achannel = geompy.CreateGroup(Geo, geompy.ShapeType["SOLID"])
geompy.UnionIDs(achannel, [108])
anode = geompy.CreateGroup(Geo, geompy.ShapeType["SOLID"])
geompy.UnionIDs(anode, [142, 166, 183])
cathode = geompy.CreateGroup(Geo, geompy.ShapeType["SOLID"])
geompy.UnionIDs(cathode, [514, 490, 531])
cmpl = geompy.CreateGroup(Geo, geompy.ShapeType["SOLID"])
geompy.UnionIDs(cmpl, [432, 456, 473])
ccl = geompy.CreateGroup(Geo, geompy.ShapeType["SOLID"])
geompy.UnionIDs(ccl, [374, 415, 398])
electrolyte = geompy.CreateGroup(Geo, geompy.ShapeType["SOLID"])
geompy.UnionIDs(electrolyte, [340, 316, 357])
acl = geompy.CreateGroup(Geo, geompy.ShapeType["SOLID"])
geompy.UnionIDs(acl, [282, 258, 299])
ampl = geompy.CreateGroup(Geo, geompy.ShapeType["SOLID"])
geompy.UnionIDs(ampl, [224, 241, 200])

# edges
v1 = geompy.CreateGroup(Geo, geompy.ShapeType["EDGE"])
geompy.UnionIDs(v1, [31])
v2 = geompy.CreateGroup(Geo, geompy.ShapeType["EDGE"])
geompy.UnionIDs(v2, [55])
v3 = geompy.CreateGroup(Geo, geompy.ShapeType["EDGE"])
geompy.UnionIDs(v3, [79])
h1 = geompy.CreateGroup(Geo, geompy.ShapeType["EDGE"])
geompy.UnionIDs(h1, [614, 13])
h2 = geompy.CreateGroup(Geo, geompy.ShapeType["EDGE"])
geompy.UnionIDs(h2, [556, 92])
h3 = geompy.CreateGroup(Geo, geompy.ShapeType["EDGE"])
geompy.UnionIDs(h3, [498, 150])
h4 = geompy.CreateGroup(Geo, geompy.ShapeType["EDGE"])
geompy.UnionIDs(h4, [440, 208])
h5 = geompy.CreateGroup(Geo, geompy.ShapeType["EDGE"])
geompy.UnionIDs(h5, [266, 382])
h6 = geompy.CreateGroup(Geo, geompy.ShapeType["EDGE"])
geompy.UnionIDs(h6, [324])
l1 = geompy.CreateGroup(Geo, geompy.ShapeType["EDGE"])
geompy.UnionIDs(l1, [6])
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Geo, 'Geo' )
geompy.addToStudyInFather( Geo, interconnect0, 'interconnect0' )
geompy.addToStudyInFather( Geo, interconnect1, 'interconnect1' )
geompy.addToStudyInFather( Geo, interconnectSides, 'interconnectSides' )
geompy.addToStudyInFather( Geo, airInlet, 'airInlet' )
geompy.addToStudyInFather( Geo, fuelInlet, 'fuelInlet' )
geompy.addToStudyInFather( Geo, airOutlet, 'airOutlet' )
geompy.addToStudyInFather( Geo, fuelOutlet, 'fuelOutlet' )
geompy.addToStudyInFather( Geo, cathodeSides, 'cathodeSides' )
geompy.addToStudyInFather( Geo, anodeSides, 'anodeSides' )
geompy.addToStudyInFather( Geo, electrolyteSides, 'electrolyteSides' )
geompy.addToStudyInFather( Geo, interconnect, 'interconnect' )
geompy.addToStudyInFather( Geo, cchannel, 'cchannel' )
geompy.addToStudyInFather( Geo, achannel, 'achannel' )
geompy.addToStudyInFather( Geo, anode, 'anode' )
geompy.addToStudyInFather( Geo, cathode, 'cathode' )
geompy.addToStudyInFather( Geo, cmpl, 'cmpl' )
geompy.addToStudyInFather( Geo, ccl, 'ccl' )
geompy.addToStudyInFather( Geo, electrolyte, 'electrolyte' )
geompy.addToStudyInFather( Geo, acl, 'acl' )
geompy.addToStudyInFather( Geo, ampl, 'ampl' )
geompy.addToStudyInFather( Geo, v1, 'v1' )
geompy.addToStudyInFather( Geo, v2, 'v2' )
geompy.addToStudyInFather( Geo, v3, 'v3' )
geompy.addToStudyInFather( Geo, h1, 'h1' )
geompy.addToStudyInFather( Geo, h2, 'h2' )
geompy.addToStudyInFather( Geo, h3, 'h3' )
geompy.addToStudyInFather( Geo, h4, 'h4' )
geompy.addToStudyInFather( Geo, h5, 'h5' )
geompy.addToStudyInFather( Geo, h6, 'h6' )
geompy.addToStudyInFather( Geo, l1, 'l1' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

geo = smesh.Mesh(Geo)
Regular_1D = geo.Segment()
interconnect0 = geo.GroupOnGeom(interconnect0,'interconnect0',SMESH.FACE)
interconnect1 = geo.GroupOnGeom(interconnect1,'interconnect1',SMESH.FACE)
interconnectSides = geo.GroupOnGeom(interconnectSides,'interconnectSides',SMESH.FACE)
airInlet = geo.GroupOnGeom(airInlet,'airInlet',SMESH.FACE)
fuelInlet = geo.GroupOnGeom(fuelInlet,'fuelInlet',SMESH.FACE)
airOutlet = geo.GroupOnGeom(airOutlet,'airOutlet',SMESH.FACE)
fuelOutlet = geo.GroupOnGeom(fuelOutlet,'fuelOutlet',SMESH.FACE)
cathodeSides = geo.GroupOnGeom(cathodeSides,'cathodeSides',SMESH.FACE)
anodeSides = geo.GroupOnGeom(anodeSides,'anodeSides',SMESH.FACE)
electrolyteSides = geo.GroupOnGeom(electrolyteSides,'electrolyteSides',SMESH.FACE)
interconnect = geo.GroupOnGeom(interconnect,'interconnect',SMESH.VOLUME)
cchannel = geo.GroupOnGeom(cchannel,'cchannel',SMESH.VOLUME)
achannel = geo.GroupOnGeom(achannel,'achannel',SMESH.VOLUME)
anode = geo.GroupOnGeom(anode,'anode',SMESH.VOLUME)
cathode = geo.GroupOnGeom(cathode,'cathode',SMESH.VOLUME)
cmpl = geo.GroupOnGeom(cmpl,'cmpl',SMESH.VOLUME)
ccl = geo.GroupOnGeom(ccl,'ccl',SMESH.VOLUME)
electrolyte = geo.GroupOnGeom(electrolyte,'electrolyte',SMESH.VOLUME)
acl = geo.GroupOnGeom(acl,'acl',SMESH.VOLUME)
ampl = geo.GroupOnGeom(ampl,'ampl',SMESH.VOLUME)

#geo.GetMesh().RemoveSubMesh( smeshObj_1 ) ### smeshObj_1 has not been yet created
Quadrangle_2D = geo.Quadrangle(algo=smeshBuilder.QUADRANGLE)
geo_3D = geo.Hexahedron(algo=smeshBuilder.Hexa)
Regular_1D_2 = geo.Segment(geom=v1)
v1_2 = Regular_1D_2.NumberOfSegments(5)
Regular_1D_3 = geo.Segment(geom=v2)
v2_2 = Regular_1D_3.NumberOfSegments(10)
Regular_1D_4 = geo.Segment(geom=v3)
v3_2 = Regular_1D_4.NumberOfSegments(5)
Regular_1D_5 = geo.Segment(geom=h1)
h1_2 = Regular_1D_5.NumberOfSegments(5)
Regular_1D_6 = geo.Segment(geom=h2)
h2_2 = Regular_1D_6.NumberOfSegments(10)
Regular_1D_7 = geo.Segment(geom=h3)
h3_2 = Regular_1D_7.NumberOfSegments(5)
Regular_1D_8 = geo.Segment(geom=h4)
h4_2 = Regular_1D_8.NumberOfSegments(5)
Regular_1D_9 = geo.Segment(geom=h5)
h5_2 = Regular_1D_9.NumberOfSegments(5)
Regular_1D_10 = geo.Segment(geom=h6)
h6_2 = Regular_1D_10.NumberOfSegments(10)
Regular_1D_11 = geo.Segment(geom=l1)
l1_2 = Regular_1D_11.NumberOfSegments(40)

Propagation_of_Node = Regular_1D_2.PropagationOfDistribution()

#status = geo.AddHypothesis(Propagation_of_Node,v1)
status = geo.AddHypothesis(Propagation_of_Node,v2)
status = geo.AddHypothesis(Propagation_of_Node,v3)
status = geo.AddHypothesis(Propagation_of_Node,h1)
status = geo.AddHypothesis(Propagation_of_Node,h2)
status = geo.AddHypothesis(Propagation_of_Node,h3)
status = geo.AddHypothesis(Propagation_of_Node,h4)
status = geo.AddHypothesis(Propagation_of_Node,h5)
status = geo.AddHypothesis(Propagation_of_Node,h6)
status = geo.AddHypothesis(Propagation_of_Node,l1)

isDone = geo.Compute()

try:
  geo.ExportUNV( r'./salome.unv' )
  pass
except:
  print('ExportUNV() failed. Invalid file name?')

v1_3 = Regular_1D_2.GetSubMesh()
v2_3 = Regular_1D_3.GetSubMesh()
v3_3 = Regular_1D_4.GetSubMesh()
h1_3 = Regular_1D_5.GetSubMesh()
h2_3 = Regular_1D_6.GetSubMesh()
h3_3 = Regular_1D_7.GetSubMesh()
h4_3 = Regular_1D_8.GetSubMesh()
h5_3 = Regular_1D_9.GetSubMesh()
h6_3 = Regular_1D_10.GetSubMesh()
l1_3 = Regular_1D_11.GetSubMesh()


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
