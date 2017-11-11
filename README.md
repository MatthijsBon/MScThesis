# Underground Localization
This repository holds the source code for my thesis for the fulfillment of the degree of MSc Science in Geomatics for the built environment. The MSc Title was awarded on November 9th 2017.    

!! Note: Some data is removed from this repository for copyright purposes. Therefore, code in this repository will not execute properly. Transformer data and validation data should be retrieved somewhere else.

## INPUT DATA
There are three main sources for the input data:    
1. BAG buildings as polygons and address points
2. NWB (National Road Network) as polylines
3. Transformers as points

These three datasources have to be processed and cleaned before usage in the method. This consists of the following steps:
1. Find all address and VBO points within each building
2. Add the id's of all points within a building as attribute to that building
3. Find the most occurring address point for each building and select that as 'building point'

Optional:    
3a. Find the nearest neighbors for all buildings    
3b. Group neighboring buildings together  

4. Save building points 

## MANUAL SKELETON CONSTRUCTION
5. Connect buildings and transformers to the NWB using GRASS GIS and QGIS Networks plugin
5a. Method A: Connect to Closest Point on the street network 
5b. Method B: Connect to Closest Junction Vertex on the street network
5c. Method C: Iteratively Connect to the Closest Junction Vertex on the street network

A: Order Closest Junction
1. Split NWB in segments (for example 75m)
2. Extract nodes and save as endpoints shapefile
3. Run connect2ClosestPoint in python, for blds and trafo's
4. Merge WKT blds and trafo's with NWB
5. v.clean -snap t=0.100000
6. v.clean -break t=0.100000
7. Build graph (Network)
8. Join on location (Network + blds + trafos)
9. Add type attribute
10. Add length attribute to NWB

B: Order Closest Point
1. Connect NWB with trafos & blds
2. Build graph (Network)
3. Join on location (Network + blds + trafos)
4. Add type attribute
5. Check
6. Add length attribute to NWB

C: Order Iterative Closest Junction
1. Import maps in GRASS
2. v.net -connect NWB with trafos and then blds
3. Build graph (Network)
4. Join on location (Network + blds + trafos)
5. Add type attribute
6. Check
7. Add length attribute to NWB

## GRAPH ANALYSIS
6. Read in Graph in NetworkX format
7. Find all paths from all buildings to the closest transformer
8. Calculate edge betweenness for all edges in the network
9. Calculate Thickness & Current
10. Quantify cables


