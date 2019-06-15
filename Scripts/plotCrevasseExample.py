""" plotCrevasseExample.py
	
	Script to plot an example of the crevasse detection scheme overlaid on a DMS image
	Model written by Alek Petty (June 2019)
	Contact me for questions (alek.a.petty@nasa.gov) or refer to the GitHub site

	Python dependencies:
		See below for the relevant module imports. Of note:
		

	Update history:
		06/2019: Version 1

"""

import matplotlib
matplotlib.use("AGG")

import numpy as np
from pylab import *
from scipy.io import netcdf
import numpy.ma as ma
import string

ut.reset_matplotlib()


#==================== Configuration ====================


# Projection of the DMS data
m=pyproj.Proj("+init=EPSG:3413")

datapath='.,/Data_output/'
rawdatapath = '../../../../DATA/ICEBRIDGE/'
ATM_path = rawdatapath+'ATM/ARCTIC/'
dms_path = rawdatapath+'DMS/GR/'
posAV_path =rawdatapath+'/POSAV/'
figpath='../Figures/'


pwidth=20
pint=5
min_depth = 0.5
xy_res=1
min_feature_size=100


date='20120417'
#dms_time='14530666'
dms_time='14530787'
# or...
#dms_image = 12
#date, dms_time = ro.pick_dms(dms_image)


#==================== Get DMS image ====================

year = int(date[0:4])
image_path = glob(dms_path+str(year)+'/'+date+'/*'+date+'_'+dms_time+'.tif')
xT, yT, dms, geo = ut.get_dms(image_path[0])
minX = xT[-1, 0]
minY = yT[-1, 0]
xT = xT - minX
yT = yT - minY

#==================== Get ATM elevation data ====================

print 'Finding atm...'
atm_file = ut.get_atm_from_dms(res=1)

lon, lat, elevation, azi = ut.get_atmqih5(atm_file, year, 1, utc_time=0)
xpts, ypts = m(lon, lat)
xpts = xpts - minX
ypts = ypts - minY

# Get elevation data within bounds of DMS image
dmsmask = where((xpts>np.amin(xT))&(xpts<np.amax(xT))&(ypts>np.amin(yT))&(ypts<np.amax(yT)))
lon=lon[dmsmask]
lat=lat[dmsmask]
xpts=xpts[dmsmask]
ypts=ypts[dmsmask]
elevation=elevation[dmsmask]

#elevation = -elevation
elevation=elevation-np.mean(elevation)

xx = np.arange(np.amin(xT),np.amax(xT), xy_res)
yy = np.arange(np.amin(yT),np.amax(yT), xy_res)
xx2d, yy2d = meshgrid(xx, yy)
elevation2d = ut.grid_elevation(xpts, ypts,elevation, xx2d, yy2d,kdtree=1)

# Elevation relative to a fitted plane (order 2=quadratic plane)
elevation2d_plane = ut.getPlaneElev(elevation2d, xx2d, yy2d, order=2)

# Find local level (modal) surface
level_elev, thresh, levpercent = ut.calc_level_ice(asarray(elevation2d_plane[ma.nonzero(elevation2d_plane)]), pint, pwidth, min_ridge_height)

# Elevation anomalies relative to a local level (modal) surface
elevation2d_anomalies=elevation2d_plane-level_elev
# Threhsold
thresh=thresh-level_elev

elevation2d_masked=ma.masked_where(elevation2d_anomalies<thresh, elevation2d_anomalies)
feature_area = ma.count(elevation2d_masked)

#================ Label the features ==================

labelled_image  = ut.label_features(elevation2d_masked, xy_res, min_ridge_size, min_ridge_height)
found_feature = 0
if (np.amax(labelled_image)>0):
	found_big_feature=1
	num_features = np.amax(labelled_image)
	#-------------- LABEL STATS ------------------
	ridge_stats, ridge_height_mesh, cov_matrix, index = ut.calc_feature_stats(elevation2d_masked, num_features, labelled_image, xx2d, yy2d, level_elev,0, calc_orientation=1)


#==================== Plotting ====================

minvalD = 0
maxvalD = 255
sizex = np.amax(xT) - np.amin(xT)
sizey = np.amax(yT) - np.amin(yT)
ratio = sizey/sizex
lowerp = 5
upperp = 99.5

minvalL = 0
maxvalL = np.round(np.percentile(ma.compressed(elevation2d_ridge_ma), 99), decimals=1)
minvalC = np.round(np.percentile(ma.compressed(elevation2dC), 1), decimals=1)
maxvalC = np.round(np.percentile(ma.compressed(elevation2dC), 99), decimals=1)
minvalE = np.round(np.percentile(elevation, 1), decimals=1)
maxvalE = np.round(np.percentile(elevation, 99), decimals=1)

textwidth=6

fig = figure(figsize=(textwidth,textwidth*ratio*0.5))

ax1 = subplot(131)
ax1.annotate('(a) Raw DMS + ATM' , xy=(0.03, 1.03), textcoords='axes fraction', color='k', horizontalalignment='middle', verticalalignment='middle')
res = 2
im1 = pcolormesh(xT[::res, ::res], yT[::res, ::res], dms[::res, ::res], vmin = minvalD, vmax = maxvalD, cmap = cm.gist_gray, rasterized=True)
im11 = scatter(xpts, ypts, c = elevation, vmin=minvalE, vmax=maxvalE, s=1, lw = 0, cmap = cm.RdYlBu_r, rasterized=True)

cax1 = fig.add_axes([0.25  , 0.37, 0.07, 0.03])
cbar = colorbar(im11,cax=cax1, orientation='horizontal', extend='both', use_gridspec=True)
cbar.set_label('H (m)', labelpad=1, rotation=0)
xticks1 = np.linspace(minvalE, maxvalE, 2)
cbar.set_ticks(xticks1)

ax2 = subplot(132)
ax2.annotate('(a) Quadratic plane - ATM' , xy=(0.03, 1.03), textcoords='axes fraction', color='k', horizontalalignment='middle', verticalalignment='middle')
res = 2
im2 = pcolormesh(xx2d, yy2d, elevation2dC, vmin = minvalC, vmax = maxvalC, cmap = cm.RdYlBu_r, rasterized=True)

cax2 = fig.add_axes([0.56  , 0.37, 0.08, 0.03])
cbar2 = colorbar(im2,cax=cax2, orientation='horizontal', extend='both', use_gridspec=True)
cbar2.set_label(r'H$_c$ (m)', labelpad=1, rotation=0)
xticks2 = np.linspace(minvalC, maxvalC, 2)
cbar2.set_ticks(xticks2)

ax3 = subplot(133)
ax3.annotate('(b) Features >'+str(min_ridge_height)+'m' , xy=(0.03, 1.03), textcoords='axes fraction', color='k', horizontalalignment='middle', verticalalignment='middle')
atm_2d_mean = '%.2f' %mean(ridge_stats[:, 6])

im3 = pcolormesh(xx2d, yy2d, elevation2d_ridge_ma, vmin = minvalL, vmax = maxvalL, cmap = cm.cubehelix, rasterized=True)
im31 = contour(xx2d, yy2d, label_im, np.arange(np.amax(label_im)+1), colors='k', linewidths=0.15)

cax3 = fig.add_axes([0.87 , 0.37, 0.08, 0.03])
cbar3 = colorbar(im3,cax=cax3, orientation='horizontal', extend='both', use_gridspec=True)
cbar3.set_label(r'$\Delta$H (m)', labelpad=1, rotation=0)
xticks3 = np.linspace(minvalL, maxvalL, 2)
cbar3.set_ticks(xticks3)

axesname=['ax1', 'ax2', 'ax3']
for ax in axesname:
	vars()[ax].set_xlim(np.amin(xT),np.amax(xT))
	vars()[ax].set_ylim(np.amin(yT),np.amax(yT))
	vars()[ax].set_xlabel('x (m)', labelpad=1)

ax2.set_yticklabels([])
ax3.set_yticklabels([])	

subplots_adjust(bottom=0.18, left=0.05, top = 0.9, right=0.98, hspace=0.22, wspace=0.05)

savefig(figpath+'ATMDMS_'+date+'_'+dms_time+str(int(min_ridge_height*100))+'m.png', dpi=300)





