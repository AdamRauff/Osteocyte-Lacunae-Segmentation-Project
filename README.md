# Osteocyte-Lacunae-Segmentation-Project

To run the program, simply run LacunarDistance.m
This script calls on all the necessary files.

Briefly, here is the pipeline of the program

1) GetVoxelLocations
  -set some flags
  -segment lacunae, and call on GUI for manual quality control
  -Calculate basic voxel statistics about each lacuna

2) AnalyzeLacunae.m
  -Calculate the principal moments of intertia of each lacuna

3) 
