# volcalc_2.2.py  
########################################################################
# Copyright 2021 Allan R. Willms
#
# volcalc.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see
# <http://www.gnu.org/licenses/>.
# 
# Allan R. Willms
# Dept. of Mathematics and Statistics
# University of Guelph
# Guelph, ON N1G 2W1
# Canada
#########################################################################
"""
volcalc_2.2 - calculate the volume of a solid given two cross sectional boundary curves

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
"""

import csv,sys,copy,argparse
import numpy as np
from scipy import linalg

NPTS = 200

def extract_loc(curve, fixed_coord, level, remove=False):
	"""
	Extract locations of (x,z) curve at the given level of the fixed coord (0 or 1).

	Remove curve points if "remove" flag is set and if point and both neighbours below level.
	Returns: sorted list of locations of the non-fixed coordinate
	"""
	fc = fixed_coord
	vc = (fc - 1)%2
	loc = []
	i = 0
	while i < len(curve):
		if ((curve[i][fc] - level)*(curve[i-1][fc] - level) <= 0) and (
				curve[i][fc] != level or curve[i-1][fc] != level):
			loc.append(curve[i-1][vc] + (curve[i][vc] - curve[i-1][vc])*(
					(level - curve[i-1][fc])/(curve[i][fc] - curve[i-1][fc])))
			n = len(curve)
			if remove and n > 3:
				if curve[(i-1)%n][fc] < level and curve[(i-2)%n][fc] < level and  \
						curve[(i-3)%n][fc] < level:
					del curve[(i-2)%n]
					if i-2 >= 0:  # deletion will effectively advance us so remove increment on i
						i -= 1
				elif curve[i][fc] < level and curve[(i+1)%n][fc] < level and \
						curve[(i+2)%n][fc] < level:
					del curve[(i+1)%n]
					if (i+1)%n < i:  # deletion will effectively advance us so remove increment on i
						i -= 1
		i += 1
	# sort locations in increasing order
	loc.sort()
	return loc

def split_xloc(xloc):
	"""
	Split a x-location set into distances below zero and distances above.

	Also, average near-neighbours to ensure there are no more than three distances on either side of
	zero.
	"""
	dist = [ [], [] ]
	for i in range(len(xloc)):
		if xloc[i] >= 0.0:
			dist[0].append(xloc[i])
		else:
			dist[1].append(abs(xloc[i]))
	for j in range(2):
		dist[j].sort()
		while len(dist[j]) > 3:
			d = dist[j][1] - dist[j][0]
			k = 1
			for i in range(2,len(dist[j])):
				d2 = dist[j][i] - dist[j][i-1]
				if d2 < d:
					k = i
					d = d2
			dist[j][k-1] = 0.5*sum(dist[j][k-1:k+1])
			del dist[j][k]
		dist[j].sort()
	return dist

def comp_area(xdist,ydist):
	"""
	Compute the area in the quadrant
	"""
	pi = np.pi
	# reorder so that xdist is the longest list
	nx = len(xdist)
	ny = len(ydist)
	if nx < ny:
		temp = xdist
		xdist = ydist
		ydist = temp
		temp = nx
		nx = ny
		ny = temp
	if nx == 0:
		area = 0.0
	elif ny == 0:
		if nx == 1:  # half rose petal
			area = pi*xdist[0]**2/16
		elif nx == 2:  # half an ellipse, height no more than width and below 45 deg line
			xm = 0.5*sum(xdist)
			a = abs(xm - xdist[0])
			b = min(a,np.sqrt(xm**2 - a**2))
			area = pi*a*b/2
		else: # nx == 3:  # half rose petal with half ellipse removed
			area = pi*xdist[2]**2/16
			xm = 0.5*sum(xdist[0:2])
			a = abs(xm - xdist[0])
			# approximate cos theta at point on rose when x=xm
			costheta = 1/np.sqrt(2) + (xm/xdist[2])*(1-1/np.sqrt(2))
			sintheta = np.sqrt(1-costheta**2)
			# make ellipse no higher that half height of rose curve at xm
			b = min(a,0.5*xdist[2]*(2*costheta**2 - 1)*sintheta)
			area -= pi*a*b/2
	elif ny == 1:
		if nx == 1:  # quarter ellipse
			area = pi*xdist[0]*ydist[0]/4
		elif nx == 2:  # quarter ellipse with half rose removed
			area = pi*xdist[1]*ydist[0]/4
			area -= pi*xdist[0]**2/16
		else: # nx == 3:  # quarter ellipse with half ellipse removed
			area = pi*xdist[2]*ydist[0]/4
			xm = 0.5*sum(xdist[0:2])
			a = abs(xm - xdist[0])
			# make ellipse stay below 45 degree line
			b = min(a,np.sqrt(xm**2 - a**2))
			area -= pi*a*b/2
	elif ny == 2:
		if nx == 2:  # quarter ellipse with small quarter ellipse removed
			area = pi*xdist[1]*ydist[1]/4
			area -= pi*xdist[0]*ydist[0]/4
		else: # nx == 3:  # quarter ellipse with half ellipse and half rose removed
			area = pi*xdist[2]*ydist[1]/4
			xm = 0.5*sum(xdist[0:2])
			a = abs(xm - xdist[0])
			# make ellipse stay below 45 degree line
			b = min(a,np.sqrt(xm**2 - a**2))
			area -= pi*a*b/2
			# now remove half rose on y-axis
			area -= pi*ydist[0]**2/16
	else: # ny == 3:
		# quarter ellipse with two half ellipses removed
		area = pi*xdist[2]*ydist[2]/4
		# remove half ellipse on x-axis
		xm = 0.5*sum(xdist[0:2])
		a = abs(xm - xdist[0])
		# make ellipse stay below 45 degree line
		b = min(a,np.sqrt(xm**2 - a**2))
		area -= pi*a*b/2
		# remove half ellipse on y-axis
		ym = 0.5*sum(ydist[0:2])
		b = abs(ym - ydist[0])
		# make ellipse stay above 45 degree line
		a = min(b,np.sqrt(ym**2 - b**2))
		area -= pi*a*b/2
	return area

def angle_calc(pt0,pt1):
	""" 
	Compute the angle with the positive horizontal axis made by the secant from pt0 to pt1, 
	pt0 and pt1 should be np.arrays of length 2.
	"""
	v = pt1 - pt0
	return np.arctan2(v[1],v[0])

def wrapping_slice(lst, *args):
	"""
	A slice that can wrap around the ends of a list.
	"""
	return [lst[i%len(lst)] for i in range(*args)]

def angle_from_fit_poly(points):
	"""
	Compute an angle between points by fitting a polynomial.

	Input is three or more points, we rotate around points[-1] so that points[0] is on the positive x-axis
	then do a least squares fit of a quadratic (if there are four or more points) or a line (if there are three 
	points), then rotate back.  Return value is the fitted angle from points[0] to points[end]
	"""
	n = len(points)
	theta = angle_calc(points[-1],points[0])
	if n == 3:
		order = 1
	else:
		order = 2
	# translate points[-1] to the origin
	for i in range(n-1):
		points[i] -= points[-1]
	points[-1] = np.zeros(2)
	# rotate points by -theta so points[0] is on positive x-axis
	ctheta = np.cos(theta)
	stheta = np.sin(theta)
	A = np.array([[ctheta, stheta],[-stheta, ctheta]])
	for i in range(n-1):
		points[i] = A@points[i]
	A = np.ones((n,order+1))
	b = np.ones(n)
	for i in range(n):
		for j in range(1,order+1):
			A[i,j] = A[i,j-1]*points[i][0]
		b[i] = points[i][1]
	(x,res,rank,s) = linalg.lstsq(A,b)
	angle = np.arctan(x[1])
	# rotate by theta back again and add pi so that angle is from points[0] to points[end]
	angle = (angle + theta + np.pi)%(2.0*np.pi)
	if angle > np.pi:
		angle -= 2.0*np.pi
	return angle

def distance_calc(curve):
	"""
	Calculate distances between points, sort them, and find the median distance.
	dist[j] is the distance from point j-1 to point j
	"""
	dist = []
	npts = len(curve)
	for j in range(npts):
		dist.append(linalg.norm((curve[j-1] - curve[j])))
	distsorted = [[i,dist[i]] for i in range(npts)]
	distsorted.sort(key=lambda d: d[1], reverse=True)
	mediandist = distsorted[round(npts/2)][1]
	return dist, distsorted, mediandist

def gap_edge_angle_calc(j,gap,curve,dist,mediandist,filldist):
	"""
	Determine tangent directions from the edges of the gap inward.

	Determine tangent directions from the edges of the gap between points j-1 and j in the
	inward (into the gap) direction.  Close points to the edge point are itself and the
	those that are sequentially less than the max of twice the filldistance away from the
	one closer to the end point, and a tenth of the gap distance.  If there are at least 4
	close points, least squares fit a quadratic to them and the end point, and use the
	angle of the quadratic at the end point.  If there are only 3 close points, then least
	squares fit a line to them.  If there are two close points use the angle between them.
	If there is one close point then do a weighted average of the two angles on either
	side of the edge point.  We use a maximum of maxclose points.
	"""
	maxclose = 5   # code below assumes this is at least 4
	useangle = []
	for direc in [1,-1]:   
		# direc means direction in curve list, so direc = 1 means we are computing the forward
		# tangent at the point before the gap, while direc = -1 means we are computing the
		# backward tangent after the gap.  The gap is from the point j-1 to point j.
		# In computing the tangent we look at at most maxclose points.
		if direc == 1:
			edge = j - 1
			distances = wrapping_slice(dist, edge - (maxclose-2), edge + 1)
		else:
			edge = j
			distances = wrapping_slice(dist, edge + maxclose-1, edge, -1)
		points = wrapping_slice(curve, edge - (maxclose-1)*direc, edge + direc, direc)
		angle = np.zeros(maxclose-1)
		closepoints = 1
		nearby = max(0.1*gap,2*filldist)
		for m in range(maxclose-2,-1,-1):
			angle[m] = angle_calc(points[m], points[m+1]) # this must be calculated for m=maxclose-2
			if distances[m] <= nearby:
				closepoints += 1
			else:
				break
		if closepoints > 2:
			useangle.append(angle_from_fit_poly(copy.deepcopy(points[maxclose-closepoints:maxclose])))
		elif closepoints == 2:  # only two points 
			useangle.append(angle[maxclose-2])
		else:   # only end point, use a weighted average of secants on either side
			npts = len(curve)
			if direc == 1:
				distances[0] = dist[(edge + direc)%npts]
			else:
				distances[0] = dist[edge%npts]
			npts = len(curve)
			angle[0] = angle_calc(curve[(edge)%npts], curve[(edge + direc)%npts])
			# angle[maxclose-2] is from edge-direc to edge
			if angle[maxclose-2] < 0.0 and angle[0] > angle[maxclose-2] + np.pi:  
				angle[0] -= 2.0*np.pi
			elif angle[maxclose-2] > 0.0 and angle[0] < angle[maxclose-2] - np.pi:
					angle[0] += 2.0*np.pi
			s = distances[maxclose-2] + distances[0]
			useangle.append((distances[0]*angle[maxclose-2] + distances[maxclose-2]*angle[0])/s)
	return useangle

def Bezier_filler(curve):
	"""
	Fill gaps in the curve with a Bezier curve.

	Fill gaps with points generated from a Bezier curve.  A gap is any distance greater
	than twice the median distance.  
	"""
	dist, distsorted, mediandist = distance_calc(curve)
	filldist = 2.0*mediandist
	# Fill in gaps that are larger than twice the median distance, starting with largest.
	k = 0
	while distsorted[k][1] > filldist:
		j = distsorted[k][0]
		useangle = gap_edge_angle_calc(j,distsorted[k][1],curve,dist,mediandist,filldist)
		# Find the intersection of the two lines coming from points j-1 and j in the specified
		# directions.
		A = np.array([[np.cos(useangle[0]),-np.cos(useangle[1])], \
				[np.sin(useangle[0]),-np.sin(useangle[1])]])
		s = linalg.solve(A,curve[j] - curve[j-1])
		intersectpt = curve[j-1] + A[:,0]*s[0]
		# The control points are forward (into gap) along the rays a distance equal to 
		# somewhere between 0.25*dist[j] and 0.75*dist[j], depending on the distances s from 
		# the edges to the intersection point, but never more than 0.8 of the way to
		# the intersection point.  If one or both of s is negative, then the
		# control points are 0.25*dist[j] forward along rays.
		if all(s > 0):
			rel = min(2.0,(s[0] + s[1])/dist[j])  # rel is between 1 and 2
			x = (0.75 - 0.5*(rel - 1.0))*dist[j]
			suse = [min(0.8*s[0],x),min(0.8*s[1],x)]
		else:
			suse = [0.25*dist[j],0.25*dist[j]]
		controlpt = [curve[j-1] + A[:,0]*suse[0], curve[j] - A[:,1]*suse[1]]
		t = 0
		ttry = min(0.5,1.0/np.floor(dist[j]/mediandist))
		bracket = [0,1]
		newpts = [curve[j-1]]   # put in a copy of the j-1 point
		newdist = []
		n = 0
		while linalg.norm(newpts[n] - curve[j]) > mediandist:
			pt = ((1-ttry)**3)*curve[j-1] + 3*(1-ttry)**2*ttry*controlpt[0] + \
					3*(1-ttry)*ttry**2*controlpt[1] + ttry**3*curve[j]
			tempnewdist = linalg.norm(pt - newpts[n])
			if tempnewdist < 0.75*mediandist:
				bracket[0] = ttry
				ttry = 0.5*(bracket[0] + bracket[1])
			elif tempnewdist > 1.25*mediandist:
				bracket[1] = ttry
				ttry = 0.5*(bracket[0] + bracket[1])
			else:
				newpts.append(pt)
				newdist.append(tempnewdist)
				n = n+1
				dt = ttry - t
				t = ttry
				bracket = [t,1]
				ttry = min(t + dt,0.5*(bracket[0] + bracket[1]))
		# insert the new points into the curve array and the dist array
		del newpts[0]    # remove the copy of point j-1
		for i in range(len(newpts)):
			curve.insert(j + i,newpts[i].copy())
			dist.insert(j + i,newdist[i])
		# update distance from last new point to what was point j
		i = j + len(newpts)
		dist[i] = linalg.norm(curve[i] - newpts[-1])
		i = k + 1
		while distsorted[i][1] > filldist:
			if distsorted[i][0] > j:
				distsorted[i][0] += len(newpts)
			i += 1
		k += 1
	return curve

def get_data(filename):
	"""
	data are read from a csv file in a particular format
	The first row has headers.
	The rows following the header row are in pairs, one row for either a transverse curve (TC)
	and the other for a sagital curve (SC).  These two rows make up a single patient.
	There is a header "NumOfPoints".
	The columns following the "NumOfPoints" column have the points in sets of 5 columns;
	we use just the first two in each set.

	curveset[set][type] is a list of 2x1 arrays (points).
	set is an integer counting the number of TC and SC curve sets.
	type is either 0 or 1 indicating a TC curve or an SC curve.
	"""
	with open(filename, newline='') as f:
		reader = csv.reader(f,dialect='excel')
		row = next(reader)
		NumOfPointscol = row.index('NumOfPoints')
		curveset = []
		csetind = -1 
		curveind = 1
		for row in reader:
			if curveind == 1:
				curveind = 0
				csetind += 1
				curveset.append([[],[]])
			else:
				curveind = 1
			npts = int(row[NumOfPointscol])
			# We assume that after the NumOfPoints column, the columns repeat the set mmX, mmY, mmZ, pxX,
			# pxY.  But we will only use mmX and mmY.
			k = NumOfPointscol + 1
			try:
				for j in range(npts):
					curveset[csetind][curveind].append(np.array(list(map(float,row[k:k+2]))))
					k += 5
			except:
				print('Warning: failed to read all points from csv file for set {0}, curve {1}'.format(
					csetind+1,curveind+1))
				sys.exit()
		f.close()
		if curveind == 0:  # there were an odd number of rows in the data file, ignore the
			#last one but give a warning.
			print("WARNING: The file ",filename," has an odd number of data rows; ignoring last one.")
			del(curveset[csetind])
	return curveset

# START of the main program:

# set up argument parser and help text
parser = argparse.ArgumentParser(description='Volume Calculation.')
parser.add_argument('input_file', help='Input CSV file')
args = parser.parse_args()
infile = args.input_file

# Read in the data from the csv file.
curveset = get_data(infile)

# Fill in gaps in the data using Bezier curves
for i in range(len(curveset)):
	for k in range(2):
		curveset[i][k] = Bezier_filler(curveset[i][k])

# The common axis of the two curves is depth into the body (front of abdomen to spine).
# We call this direction "z".
# The sagittal curve is in the (y,z) plane, where y is the direction parallel to the spine.
# The transverse curve is in the (x,z) plane, where x is the direction across the abdomen.
# We will fill in the surface by assuming quarter ellipses in each of the four (x,y) 
# quadrants, for each fixed z value.

outfile = open(infile[0:-3] + 'txt',mode='w',newline=None)
print('Subject  Volume (mL)',file=outfile)
print('-------- -----------',file=outfile)
for i in range(len(curveset)):
	# compute max and min of both curves together in both coordinates
	cmax = [curveset[i][0][0].copy(), curveset[i][1][0].copy()]
	cmin = [cmax[0].copy(), cmax[1].copy()]
	for k in range(2):
		for j in range(len(curveset[i][k])):
			for c in range(2):
				if curveset[i][k][j][c] > cmax[k][c]:
					cmax[k][c] = curveset[i][k][j][c]
				elif curveset[i][k][j][c] < cmin[k][c]:
					cmin[k][c] = curveset[i][k][j][c]
		# centre curve temporarily at midpoint
		cmid = 0.5*(cmax[k] + cmin[k])
		cmin[k] -= cmid
		cmax[k] -= cmid
		for j in range(len(curveset[i][k])):
			curveset[i][k][j] -= cmid
	# Adjust alignment by finding best place plus or minus 10% of smallest curve from midpoint
	start = -0.1*min(cmax[0][0] - cmin[0][0],cmax[1][0] - cmin[1][0])
	n = round(NPTS/5)
	dx = 2*abs(start)/(n - 1)
	# Extract z locations
	tranzloc = []
	sagzloc = []
	for j in range(n):
		xlev = start + j*dx
		tranzloc.append(extract_loc(curveset[i][0],0,xlev))
		sagzloc.append(extract_loc(curveset[i][1],0,xlev))
	# move both curves through the xlevels and compute the best alignment, penalize for moving in x
	# direction more than 0.5*start and for moving in z direction
	bestobjective = None
	for j1 in range(n):
		for j2 in range(n):
			zmean = [0.5*(tranzloc[j1][-1] + tranzloc[j1][0]), 0.5*(sagzloc[j2][-1] + sagzloc[j2][0])]
			objective = abs((tranzloc[j1][-1] - tranzloc[j1][0]) - (sagzloc[j2][-1] - sagzloc[j2][0])) \
					+ abs(zmean[0] - zmean[1]) \
					+ np.tan(np.pi*0.48*max(0.0,abs(start + j1*dx)/(0.5*abs(start)) - 1.0)) \
					+ np.tan(np.pi*0.48*max(0.0,abs(start + j2*dx)/(0.5*abs(start)) - 1.0))
			if bestobjective is None or bestobjective > objective:
				bestobjective = objective
				bestxcentre = [start + j1*dx, start + j2*dx]
				bestzcentre = [0.5*(zmean[0] - zmean[1])]
				bestzcentre.append(-bestzcentre[0])
	# centre curves at the above best point
	for k in range(2):
		for j in range(len(curveset[i][k])):
			curveset[i][k][j][0] -= bestxcentre[k]
			curveset[i][k][j][1] -= bestzcentre[k]
		cmin[k][1] -= bestzcentre[k]
		cmax[k][1] -= bestzcentre[k]
	# now loop through all z levels and extract x locations
	zmin = min(cmin[0][1],cmin[1][1])
	dz = (max(cmax[0][1],cmax[1][1]) - zmin)/(NPTS - 1)
	tranxloc = []
	sagxloc = []
	for j in range(NPTS):
		zlev = zmin + j*dz
		# for each curve find all places that cross this z level
		tranxloc.append(extract_loc(curveset[i][0],1,zlev,remove=True))
		sagxloc.append(extract_loc(curveset[i][1],1,zlev,remove=True))
	# Calculate the volume
	vol = 0.0
	for j in range(NPTS):
		zlev = zmin + j*dz
		# Divide distances into left and right of centre.
		split_trandist = split_xloc(tranxloc[j])
		split_sagdist = split_xloc(sagxloc[j])
		# for each quadrant, compute the area at this z level
		area = comp_area(split_trandist[0],split_sagdist[0]) + \
			comp_area(split_trandist[1],split_sagdist[0]) + \
			comp_area(split_trandist[1],split_sagdist[1]) + \
			comp_area(split_trandist[0],split_sagdist[1])
		# compute the volume of this slice
		if j == 0 or j == NPTS-1:
			vol += 0.5*area
		else:
			vol += area
	vol *= dz
	print("{0:<8} {1:>9,.2f}".format(i+1,float(vol)/1000),file=outfile)
outfile.close()
