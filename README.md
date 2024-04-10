# TTT_Ti64
Time-Temperature-Transformation model for Ti-6Al-4V alloys

This is a model for Ti6Al4V microstructure development.

This version implemented by: Bryan Webler, Carnegie Mellon University. Contact: webler@cmu.edu.

The initial developer was Andrew Huck during his time as PhD student Carnegie Mellon University. Contact: ahuck@alumni.cmu.edu

The files are:

1.   TTT_Ti6Al4V.ipynb - this is a Jupyter notebook containing the basic code to open a thermal history file and run the model.
2.   NGFunctions.py - this file contains the model physics, it must be in the same folder as the notebook
3.   SinglePointSTK.py - This file contains functions that execute the model, it must be in the same folder as the notebook
4.   sample_z_time_temp.csv - This is a sample input file, it contains time/temperature data for every position. This file is for one position coordinate (z).
5.   preprint_Location_dependent_phase_transformation_kinetics_during_laser_wire_deposition_additive_manufacturing_of_Ti6Al4V.pdf - This is a preprint of a manuscript that details the model phyiscs and comparison to experimental data.

The model requires a position/time/temperature history. This must be supplied as a *.csv file with position in meters, time in seconds, and temperature in Celsius. The model is currently configured for one position dimension.

The output from the notebook is a *.csv file with the following microstructure attributes for each position at the end of the provided time/temperature history:

1.   Total α phase fraction
2.   GB - grain boundary allotriomorph α fraction
3.   BW - basketweave α fraction
4.   COL - colony α fraction
5.   MASS - masssive α fraction
6.   MART - martensite fraction
7.   lath - α lath spacing
