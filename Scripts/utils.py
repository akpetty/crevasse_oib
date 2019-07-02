""" utils.py
	
	Common functions used by the various scripts within this repo 
	Model written by Alek Petty (June 2019)
	Contact me for questions (alek.a.petty@nasa.gov) or refer to the GitHub site

	Python dependencies:
		See below for the relevant module imports. Of note:
		matplotlib
		basemap

	Update history:
		06/2019: Version 1

"""

from osgeo import osr, gdal
import numpy as np
from glob import glob
from pylab import *
from scipy import ndimage
from matplotlib import rc
from scipy.interpolate import griddata as griddatascipy
import time
import h5py
from scipy.spatial import cKDTree as KDTree
import os
from matplotlib.path import Path
from scipy import stats
import matplotlib.patches as patches
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
from skimage.feature import peak_local_max
from skimage.morphology import watershed
from natsort import natsorted


def get_atm_from_dms(ATM_path, year, date, dms_time, res):
	"""
    Find ATM file path that started just before the DMS file time

    Args:
        ATM_path (str): ATM file path
        year (int): year of OIB campaign
        date (int): date of OIb campaign
        dms_time (string): DMS time string
        
	returns:
		the relevant ATM file
    """

    
	atm_files = glob(ATM_path+str(year)+'/'+date+'/*')
	
	for j in xrange(size(atm_files)):
		atm_time = float(atm_files[j][-17:-11])
		#print atm_time
		#print float(dms_time[0:-2])
		if (atm_time>float(dms_time[0:-2])):
			break

	return atm_files[j-1]

def get_pos_poly(xp, yp, idx1, idx2):
	"""
    Get polygon path"""

	gradx = xp[idx1+1]-xp[idx1-1]
	grady = yp[idx1+1]-yp[idx1-1]
	dx1 = 300*cos(arctan(grady/gradx)+(pi/2))
	dy1 = 300*sin(arctan(grady/gradx)+(pi/2))

	xp1 = xp[idx1]-dx1
	xp2 = xp[idx1]+dx1
	yp1 = yp[idx1]-dy1
	yp2 = yp[idx1]+dy1

	x1y1=[xp1, yp1]
	x2y2=[xp2, yp2]

	gradx = xp[idx2+1]-xp[idx2-1]
	grady = yp[idx2+1]-yp[idx2-1]
	dx2 = 300*cos(arctan(grady/gradx)+(pi/2))
	dy2 = 300*sin(arctan(grady/gradx)+(pi/2))

	xp3 = xp[idx2]+dx2
	xp4 = xp[idx2]-dx2
	yp3 = yp[idx2]+dy2
	yp4 = yp[idx2]-dy2

	x3y3=[xp3, yp3]
	x4y4=[xp4, yp4]

	hypot1 = hypot(xp2-xp1, yp2-yp1)
	hypot2 = hypot(xp3-xp2, yp3-yp2)
	hypot3 = hypot(xp4-xp3, yp4-yp3)
	hypot4 = hypot(xp1-xp4, yp1-yp4)

	return Path([x1y1, x2y2, x3y3, x4y4, x1y1], closed=True), [x1y1, x2y2, x3y3, x4y4, x1y1], [hypot1, hypot2, hypot3, hypot4]

def get_atm_files(atm_path_date, year):
	"""
    Get list of ATM files available in a given year

	returns:
		list of files in the given date
    """

	if (year<=2012):
		return glob(atm_path_date+'*.qi')
	else:
		return glob(atm_path_date+'*.h5')


def get_atmqih5(atm_file, year, utc_time=1, pole_str='GR'):

	"""
    Load the ATM file

    Args:
        ATM_path (str): ATM file path
        year (int): year of OIB campaign
        utc_time (int): calculate utc time from gps time (1=yes, 0=no)
        pole_str (string): campaign string 'GR'=Arctic/Greenland, 'AN'=Antarctic/Southern Ocean
        
	returns:
		ATM data
    """

	if (year>=2013):

		atm = h5py.File(atm_file, 'r')
		elevation = atm['elevation'][:]
		lon = atm['longitude'][:]
		lat = atm['latitude'][:]
		pitch = atm['instrument_parameters']['pitch'][:]
		roll = atm['instrument_parameters']['roll'][:]
		azi =atm['instrument_parameters']['azimuth'][:]
		#multiply by 1000 to put milliseconds in front of the decimal.
		gps_time=atm['instrument_parameters']['time_hhmmss'][:]*1000
		atm.close()
		if (utc_time==1):
			utc_time = gpshhmmss_to_utc_seconds(gps_time, year, pole_str)
			return lon, lat, elevation, azi, utc_time
		else:
			return lon, lat, elevation, azi

	else:
		fid = open(atm_file, 'rb')
		print pole_str, year
		if (pole_str=='GR'):
			if (year<=2010):
				#BIG ENDIAN
				dtypeT='>i4'
			else:
				dtypeT='<i4'
		elif (pole_str=='AN'):
			if (year==2009):
				#BIG ENDIAN
				dtypeT='>i4'
			else:
				dtypeT='<i4'
		#header = np.fromfile(fid, dtype='>u4', count=3888)
		numofwords = np.fromfile(fid, dtype=dtypeT, count=1)/4
		blankarray = np.fromfile(fid, dtype=dtypeT, count=numofwords[0]-1)
		initialword = np.fromfile(fid, dtype=dtypeT, count=1)
		skipBytes = np.fromfile(fid, dtype=dtypeT, count=1)
		print skipBytes[0]
		if (skipBytes[0]>20000.):
			if (year==2009):
				skipBytes=[2832]
			elif (year==2010):
				skipBytes=[2928]
			elif (year==2011):
				skipBytes=[3888]
			elif (year==2012):
				skipBytes=[4176]

		fid.seek(0)
		fid.seek(skipBytes[0], os.SEEK_SET)

		data = np.fromfile(fid, dtype=dtypeT)
		data = data.reshape(-1, 12)
		atm=np.zeros((data.shape))
		atm[:, 0] = data[:, 0]/1000.
		atm[:, 1] = data[:, 1]/1000000.
		atm[:, 2] = data[:, 2]/1000000.
		atm[:, 3] = data[:, 3]/1000.
		atm[:, 4] = data[:, 4]
		atm[:, 5] = data[:, 5]
		atm[:, 6] = data[:, 6]/1000.
		atm[:, 7] = data[:, 7]/1000.
		atm[:, 8] = data[:, 8]/1000.
		atm[:, 9] = data[:, 9]/10.
		atm[:, 10] = data[:, 10]
		atm[:, 11] = data[:, 11]

		lat = atm[:, 1]
		lon = atm[:, 2]
		elevation = atm[:, 3]
		pitch = atm[:, 7]
		roll = atm[:, 8]
		gps_time = atm[:, 11]
		azi = atm[:, 6]
		#pulse_s = data[:, 4]
		#ref_s = data[:, 5]
		#azi = data[:, 6]/1000.
		#pdop = data[:, 9]/10.
		#p_width = data[:, 10]

		if (utc_time==1):
			utc_time = gpshhmmss_to_utc_seconds(gps_time, year, pole_str)
			return lon, lat, elevation, azi, utc_time
		else:
			return lon, lat, elevation, azi

def gpshhmmss_to_utc_seconds(gps_timeT, year, pole_str='GR', nooffset=0):
	""" Convert GPS time to UTC time """

	# rough time of campaign
	if (pole_str=='GR'):
		month=4 #starts at 0, May
	elif (pole_str=='AN'):
		month=10

	gps_utc_time = get_gps_utc_offset(year, month)

	#2009 is out by 15 seconds for the Arctic data for some reason? Perhaps they changed the time twice??

	if (nooffset==1):
		gps_utc_time=0
	utc_timeT=[]

	for i in range(size(gps_timeT)):
		gps_time_str = "%09d" % int(gps_timeT[i])
		utc_timeT.append((int(gps_time_str[0:2])*60*60)+(int(gps_time_str[2:4])*60)+(int(gps_time_str[4:6]))+(int(gps_time_str[6:9]))/1000.-gps_utc_time)
		
	return array(utc_timeT)

def get_gps_utc_offset(year, month):
	""" Get the GPS/UTC time offset in second. Probably a much more elegant way of doing this or a library function already """

	if (year>=2017):
		gps_utc_time=18

	elif (year==2016):
		gps_utc_time=17

	elif ((year==2015) & (month>=6)):
		gps_utc_time=17

	elif ((year==2015)& (month<6)):
		gps_utc_time=16

	elif (year==2014):
		gps_utc_time=16
		
	elif (year==2013):
		gps_utc_time=16

	elif ((year>=2012)& (month>=6)):
		gps_utc_time=16

	elif ((year==2012)& (month<6)):
		gps_utc_time=15
	
	elif (year==2011):
		gps_utc_time=15

	elif (year==2010):
		gps_utc_time=15
		
	elif (year==2009):
		if (pole_str=='GR'):
			gps_utc_time=30
		else:
			gps_utc_time=15
	print gps_utc_time
	return gps_utc_time


def label_features(elevation2d_ridge_ma, xy_res, min_ridge_size, min_ridge_height):
	#GET BOOLEEN ARRAY WHERE TRUE IS VALID RIDGE DATA.
	num_points_req=min_ridge_size/(xy_res**2)
	elevation2d_ridge_mask = ~ma.getmask(elevation2d_ridge_ma)
	#SQURE OR DIAGONAL STRUCTIRE FOR LABELLING
	struct = ndimage.generate_binary_structure(2, 2)
	#scan across array and label connected componenets
	# returns a labelled array (integers for different labels) and the number of labels
	label_im, nb_labels = ndimage.label(elevation2d_ridge_mask, structure=struct)
	#find the unique label numbers and how many of each label are in the labelled array
	label_nums, label_num_counts = unique(label_im, return_counts=True)
	#make a note of the label numbers where there are less than a certain number of points in that label
	#labels_small = label_nums[where(label_num_counts<num_points_req)]

	label_imBIG = np.copy(label_im)
	mask_size = label_num_counts < num_points_req
	remove_pixel = mask_size[label_imBIG]

	label_imBIG[remove_pixel] = 0
	labels = np.unique(label_imBIG)
	label_imBIG = np.searchsorted(labels, label_imBIG)

	new_labels = np.zeros((label_im.shape))
	label_num_new=0
	for label_num in xrange(1, np.amax(label_imBIG)+1):
		elevation2d_ridge_maT=ma.masked_where(label_imBIG!=label_num, elevation2d_ridge_ma)
		imageT = elevation2d_ridge_maT.filled(0) - min_ridge_height
		local_maxi = peak_local_max(imageT, threshold_rel=0.25, min_distance=25/xy_res,exclude_border=False,indices=False)
		markers = ndimage.label(local_maxi)[0]

		if (np.amax(markers>0)):
			labelsT = watershed(-elevation2d_ridge_maT, markers)
			labelsT = ma.masked_where(ma.getmask(elevation2d_ridge_maT), labelsT)
			new_labels[where(labelsT!=0)]=labelsT[where(labelsT!=0)]+int(np.amax(new_labels))
		else:
			new_labels[where(label_imBIG==label_num)]=int(np.amax(new_labels)+1)

	new_labels=new_labels.astype('int')

	#return the new compressed, masked label array and the label numbers
	return new_labels

def grid_elevation(xS, yS,elevationS, xx2d, yy2d, kdtree=0):
	"""
    Interpolate elevation data onto a given grid 
	Mask where elevation is nan

    Args:
        xS (var): xpts of section
        yS (int): ypts of section
        elevationS (var): raw ATM elevation of section
        xx2d (var): 2d x grid (x coordinates)
        yy2d (var): 2d y grid (y coordinates)
	
	returns:
		gridded elevation

    """

	elevation2d = griddatascipy((xS, yS), elevationS, (xx2d, yy2d), method='nearest')
	if (kdtree==1):
		elevation2d = kdtree_clean(xx2d, yy2d, xS, yS, elevation2d)
	
	elevation2d_ma=ma.masked_where(isnan(elevation2d), elevation2d)

	return elevation2d_ma

def kdtree_clean(xx2d, yy2d, xS, yS, elevation2d):
	"""
    Remove bad data from the regridding using a kdtree
	
	dist is how far away the nearest neighbours are (default 4 m). 
    
    cleaned up gridded elevation data

    """

    # REMOVE DODGY ADDED DATA FROM THE REGRIDDING BASED ON KDTREE. 
	
	# need to decide on this threshold.
	# ONLY DO THIS FOR POINTS THAT HAVE ALREADY BEEN CLASSIFIED AS RIDGES
	grid_points = np.c_[xx2d.ravel(), yy2d.ravel()]
	tree = KDTree(np.c_[xS, yS])
	dist, _ = tree.query(grid_points, k=1)
	dist = dist.reshape(xx2d.shape)
	elevation2d_KD=ma.masked_where(dist > 4, elevation2d)
	return elevation2d_KD

def grid_atm(xS, yS, xy_res):
	xxS = np.arange(np.amin(xS),np.amax(xS), xy_res)
	yyS = np.arange(np.amin(yS),np.amax(yS), xy_res)
	
	xx2d, yy2d = meshgrid(xxS, yyS)

	return xx2d, yy2d
	
def calc_feature_stats(elevation2d_ridge_ma, num_ridges, label_im, xx2d, yy2d, level_elev, section_num, calc_orientation=0):
	"""
    Calculate ridge statistics

    Args:
        elevation2d_ridge_ma (var): masked elevation
        num_ridges (int): number of features
        label_im (var): labelled image
        xx2d (var): gridded x coordinates
        yy2d (var): gridded y coordinates
        level_elev (var): level ice elevation
        section_num (int): section number
        calc_orientation (int): If 1, calculate orientation of the feature
	
	returns:
		feature statistics

    """
	ridge_stats = ma.masked_all((num_ridges, 9))
	cov_matrix = ma.masked_all((num_ridges, 5))

	#make an empty gridded array to be filled with the ridge heights relative to level ice surface
	ridge_height_mesh = ma.masked_all((xx2d.shape))
	for i in xrange(1, num_ridges+1):
		#print i
		#get aray indices of each valid (big) label
		index_t = where(label_im==i)
		#height wrt to the lowest somehting percentile elevation
		ridge_height_mesh[index_t] = np.mean(elevation2d_ridge_ma[index_t])
		#mean x position  of ridge
		ridge_stats[i-1, 0] = mean(xx2d[index_t])
		#mean y position  of ridge
		ridge_stats[i-1, 1] = mean(yy2d[index_t])
		#mean x std of ridge points
		ridge_stats[i-1, 2] = std(xx2d[index_t])
		#mean y std of ridge points
		ridge_stats[i-1, 3] = std(yy2d[index_t])
		#mean height of ridge relative to level ice surface
		ridge_stats[i-1, 4] = np.mean(elevation2d_ridge_ma[index_t])
		#max (95th percentile) height of ridge relative to level ice surface
		ridge_stats[i-1, 5] =  np.percentile(elevation2d_ridge_ma[index_t], 95) 

		ridge_stats[i-1, 6] = np.amax(elevation2d_ridge_ma[index_t])
		#only want one coordinate size for number of points!
		ridge_stats[i-1, 7] = np.size(index_t[0])
		#section number ridge belongs to
		ridge_stats[i-1, 8] = section_num
		#CALCULATE ORIENTATION OF EACH RIDGE.
		if (calc_orientation==1):
			cov_matrix[i-1, 0:4]=np.cov(xx2d[index_t], yy2d[index_t]).flatten()
			cov_matrix[i-1, 4] = section_num

	return ridge_stats, ridge_height_mesh, cov_matrix, index_t

def calc_level_ice(elevationS, pint, pwidth, min_depth, lev_bounds=0):
	"""
    Calculate a level ice elevation (lowest turning point of the elevation distribution)

    Args:
        elevationS (var): raw ATM elevation data of a given section
        pint (int): percent of the distribution iterated through
        pwidth (int): percent of the distribution analyzed to find turning points
        min_depth (float): minimum depth added to the level ice surface to detect features
        lev_bounds (int): flag to return bounds of the level ice elevation
	
	returns:
		the level ice elevation and a few other things

    """
	difs = [np.percentile(elevationS, i+ pwidth)-np.percentile(elevationS, i) for i in range(0, int(100)-pwidth, pint)]
	difs2= diff(array(difs))
	difnegs = where(difs2<0)[0]
	if (size(difnegs)>0):
		min_index = difnegs[-1]+1
	else:
		min_index =where(difs==min(difs))[0][0]
	level_elev = np.percentile(elevationS, (pint*min_index)+(pwidth/2))
	level_elevl = np.percentile(elevationS, (pint*min_index))
	level_elevu = np.percentile(elevationS, (pint*min_index)+pwidth)
	print 'Level ice elevation:', level_elev
	thresh = level_elev+min_depth
	if (lev_bounds==1):
		return level_elev, level_elevl, level_elevu, min_index, thresh
	else:
		return level_elev, thresh, (pint*min_index)+(pwidth/2)

def getPlaneElev(elevation2d, xx2d, yy2d, order=2):
	"""
    Fit a plane (linear or quadratic) of the elevation data

    Args:
        elevation2d (var): 2d ATM elevation data of a given section
        xx2d (var): xpoints of the elevation grid
        yy2d (var): ypoints of the elevation grid
        order (int): 1=linear, 2=quadratic (default 2)
        
	returns:
		the elevation plane
    """

	egood=asarray(elevation2d[ma.nonzero(elevation2d)])
	xgood=xx2d[ma.nonzero(elevation2d)]
	ygood=yy2d[ma.nonzero(elevation2d)]
 
	if order == 1: # 1: linear, 2: quadratic
	    # best-fit linear plane
	    A = np.c_[xgood, ygood, np.ones(xgood.shape[0])]
	    C,_,_,_ = scipy.linalg.lstsq(A, egood)     # coefficients
	    
	    # evaluate it on grid
	    Z = C[0]*xx2d + C[1]*yy2d + C[2]
	    
	    # or expressed using matrix/vector product
	    #Z = np.dot(np.c_[XX, YY, np.ones(XX.shape)], C).reshape(X.shape)

	elif order == 2:
	    # best-fit quadratic curve
	    A = np.c_[np.ones(xgood.shape[0]), xgood, ygood, xgood*ygood, xgood**2, ygood**2]
	    C,_,_,_ = scipy.linalg.lstsq(A, egood)
	    
	    # evaluate it on a grid
	    Z = np.dot(np.c_[np.ones(xx2d.flatten().shape), xx2d.flatten(), yy2d.flatten(), xx2d.flatten()*yy2d.flatten(), xx2d.flatten()**2, yy2d.flatten()**2], C).reshape(yy2d.shape)

	elevation2dC = Z-elevation2d

	return elevation2dC

def get_pos_sections(posAV, m, along_track_res):
	""" Break up posAV flight location data into sections at the given resolution """

	latp = posAV[:, 1]
	lonp = posAV[:, 2]
	xp, yp = m(lonp, latp)
	dist = hypot(diff(xp), diff(yp))
	cum_dist = cumsum(dist)
	km_indices=[]
	for seg_length in xrange(along_track_res, int(cum_dist[-1]), along_track_res):
		km_indices.append(np.abs(cum_dist - seg_length).argmin())
	km_utc_times = posAV[km_indices, 0]

	return xp, yp, cum_dist, km_indices, km_utc_times

def posav_section_info(m, posAV_section, returnlonlat=0):
	""" Get data from a posAV section """
	mean_lat=mean(posAV_section[:, 1])
	mean_lon=mean(posAV_section[:, 2])
	mean_x, mean_y = m(mean_lon, mean_lat)
	mean_alt=mean(posAV_section[:, 3])
	mean_vel=mean(posAV_section[:, 4])
	mean_pitch=mean(posAV_section[:, 5])
	mean_roll=mean(posAV_section[:, 6]) 

	if (returnlonlat==1):
		return mean_x, mean_y, mean_lat, mean_lon, mean_alt, mean_pitch, mean_roll, mean_vel

	return mean_x, mean_y, mean_alt, mean_pitch, mean_roll, mean_vel

def get_atm_poly(xT, yT, elevationT, km_utc_times, utc_timeT, poly_path, i):
	""" Get ATM data iwithin a given polygon path"""
	pos_time1 = km_utc_times[i]-4
	pos_time2 = km_utc_times[i+1]+4
	xpos_rough = xT[where((utc_timeT>pos_time1)&(utc_timeT<pos_time2))]
	ypos_rough = yT[where((utc_timeT>pos_time1)&(utc_timeT<pos_time2))]
	elevation_rough = elevationT[where((utc_timeT>pos_time1)&(utc_timeT<pos_time2))]
	coords_rough = np.dstack((xpos_rough,ypos_rough))
	inpoly = poly_path.contains_points(coords_rough[0])

	xpos_km = xpos_rough[where(inpoly==True)]
	ypos_km = ypos_rough[where(inpoly==True)]
	elevation_km = elevation_rough[where(inpoly==True)]

	return xpos_km, ypos_km, elevation_km

def get_dms(image_path):
	""" Get RGB image from DMS system"""
	geo = gdal.Open(image_path) 
	band1 = geo.GetRasterBand(1)
	band2 = geo.GetRasterBand(2)
	band3 = geo.GetRasterBand(3)
	red = band1.ReadAsArray()
	green = band2.ReadAsArray()
	blue = band3.ReadAsArray()

	dms = (0.299*red + 0.587*green + 0.114*blue)
	dms = ma.masked_where(dms<1, dms)

	trans = geo.GetGeoTransform()
	width = geo.RasterXSize
	height = geo.RasterYSize

	x1 = np.linspace(trans[0], trans[0] + width*trans[1] + height*trans[2], width)
	y1 = np.linspace(trans[3], trans[3] + width*trans[4] + height*trans[5], height)
	xT, yT = meshgrid(x1, y1)
	return xT, yT, dms, geo

def reset_matplotlib():
    """
    Reset matplotlib to a common default.
    """
    
    # Set all default values.
    #mpl.rcdefaults()
    # Force agg backend.
    plt.switch_backend('agg')
    # These settings must be hardcoded for running the comparision tests and
    # are not necessarily the default values.
    rcParams['ytick.major.size'] = 2
    rcParams['axes.linewidth'] = .6
    rcParams['lines.linewidth'] = .6
    rcParams['patch.linewidth'] = .6
    rcParams['ytick.labelsize']=8
    rcParams['xtick.labelsize']=8
    rcParams['legend.fontsize']=9
    rcParams['font.size']=9
    rc('font',**{'family':'sans-serif','sans-serif':['Arial']})



