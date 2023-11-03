volcalc - calculate the volume of a solid given two cross sectional boundary curves

> python volcalc_2.2.py filename.csv

This software computes volumes of a region in space defined by two boundary curves that
are cross sections of the region at right angles.  A number of assumptions are employed
regarding how these two curves are aligned, and how the shape deforms through the regions
"between" the curves.  Gaps in the data are filled with Bezier curves.  It also assumes
the input is in a csv file with a specified format.

filename.csv is a csv file in the following format.  The first row is a row of headers.
After that the data come in pairs of rows, one row for a curve in the (x,z)-plane, and one
for a curve in the (y,z)-plane.  This pair makes up the data for one 3-dimensional volume.
There may be any number of pairs of data rows in the file.  There must be a header called
"NumOfPoints" and that column provides the number of points in the curves.  The columns
following "NumOfPoints" come in groups of five.  There are as many groups of five as the
number of points.  We only use the data in the first two colums of each of these groups 
to define points in (x,z) or (y,z) space for the curves.  It is assumed the data are in mm.

The program writes a table of values to the file: filename.txt.  One row for
each 3D volume measurement (pair of curves), with the volume given in mL.

A full description of the volume calculations is given in:

Xiu Ting Yiew, Samantha Clarke, Allan R. Willms, Shane Bateman, 2019. "Feasibility of a
Novel 3-Dimensional Mathematical Algorithmic Computation of Feline Bladder Volumes Using
Point-of-Care Longitudinal and Transverse Cysto-Colic Ultrasonographic Images", Canadian
Journal of Veterinary Research, 83, pp. 298-312.

A brief description of the gap filling algorithm is given in:

Sabrina Ayoub, Xiu Ting Yiew, Allan R. Willms, 2023.  "Urinary Bladder Volume Estimation
Using 2D Linear Dimension Formula is More Accurate than 3D Bladder Circumference Tracing
in Client-Owned Cats", Canadian Journal of Veterinary Research, (submitted).
