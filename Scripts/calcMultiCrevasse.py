""" calcMultiCrevasse.py
	
	Script to run through the crevasse detection scheme
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

def plotCrevasseSection(xatm_sect, yatm_sect, xx2d, yy2d, elevation_sect, elevation2dC, elevation2d_ridge_ma, label_im, date, section):
	# Plot the given section

	#EXPRESS RAW ELEVATIONS RELATIVE TO THE MEAN
	elevation_sect=elevation_sect-mean(elevation_sect)

	sizex = np.amax(xatm_sect) - np.amin(xatm_sect)
	sizey = np.amax(yatm_sect) - np.amin(yatm_sect)
	ratio = sizex/sizey
	
	lowerp = 5
	upperp = 99.5

	
	minvalC = np.round(np.percentile(ma.compressed(elevation2dC), 1), decimals=1)
	maxvalC = np.round(np.percentile(ma.compressed(elevation2dC), 99), decimals=1)

	minvalE = np.round(np.percentile(elevation_sect, 1), decimals=1)
	maxvalE = np.round(np.percentile(elevation_sect, 99), decimals=1)

	textwidth=6

	fig = figure(figsize=(textwidth*ratio*0.5, textwidth))

	ax1 = subplot(311)
	ax1.annotate('(a) Raw ATM' , xy=(0.03, 1.03), textcoords='axes fraction', color='k', horizontalalignment='middle', verticalalignment='middle')
	
	lonStr='%.3f' %mean_lon
	latStr='%.3f' %mean_lat
	ax1.annotate('Lon:'+lonStr+'E'+' Lat:'+latStr+'N', xy=(0.7, 1.03), textcoords='axes fraction', color='k', horizontalalignment='middle', verticalalignment='middle')
	
	res = 2
	im11 = scatter(xatm_sect, yatm_sect, c = elevation_sect, vmin=minvalE, vmax=maxvalE, s=1, lw = 0, cmap = cm.RdYlBu_r, rasterized=True)

	cax1 = fig.add_axes([0.92  , 0.725, 0.01, 0.15])
	cbar = colorbar(im11,cax=cax1, orientation='vertical', extend='both', use_gridspec=True)
	cbar.set_label(r'H$_e$ (m)', labelpad=1, rotation=0)
	xticks1 = np.linspace(minvalE, maxvalE, 2)
	cbar.set_ticks(xticks1)

	ax2 = subplot(312)
	ax2.annotate('(b) Quadratic plane - gridded ATM' , xy=(0.03, 1.03), textcoords='axes fraction', color='k', horizontalalignment='middle', verticalalignment='middle')
	res = 2
	im2 = pcolormesh(xx2d, yy2d, elevation2dC, vmin = minvalC, vmax = maxvalC, cmap = cm.RdYlBu_r, rasterized=True)

	cax2 = fig.add_axes([0.92  , 0.425, 0.01, 0.15])
	cbar2 = colorbar(im2,cax=cax2, orientation='vertical', extend='both', use_gridspec=True)
	cbar2.set_label(r'H$_c$ (m)', labelpad=1, rotation=0)
	xticks2 = np.linspace(minvalC, maxvalC, 2)
	cbar2.set_ticks(xticks2)

	ax3 = subplot(313)

	ax3.annotate('(c) Features > '+str(min_ridge_height)+' m' , xy=(0.03, 1.03), textcoords='axes fraction', color='k', horizontalalignment='middle', verticalalignment='middle')
	
	if (found_big_ridge==1):
		minvalL = 0
		maxvalL = np.round(np.percentile(ma.compressed(elevation2d_ridge_ma), 99), decimals=1)

		atm_2d_mean = '%.2f' %mean(ridge_stats[:, 6])
		ax3.annotate('Mean D: '+atm_2d_mean+'m', xy=(0.98, 0.8), textcoords='axes fraction', color='k', horizontalalignment='right', verticalalignment='middle')


		im3 = pcolormesh(xx2d, yy2d, elevation2d_ridge_ma, vmin = minvalL, vmax = maxvalL, cmap = cm.viridis_r, rasterized=True)
		im31 = contour(xx2d, yy2d, label_im, np.arange(np.amax(label_im)+1), colors='k', linewidths=0.15)

		cax3 = fig.add_axes([0.92 , 0.15, 0.01, 0.15])
		cbar3 = colorbar(im3,cax=cax3, orientation='vertical', extend='both', use_gridspec=True)
		cbar3.set_label(r'D (m)', labelpad=1, rotation=0)
		xticks3 = np.linspace(minvalL, maxvalL, 2)
		cbar3.set_ticks(xticks3)

	axesname=['ax1', 'ax2', 'ax3']
	for ax in axesname:
		vars()[ax].set_xlim(np.amin(xatm_sect),np.amax(xatm_sect))
		vars()[ax].set_ylim(np.amin(yatm_sect),np.amax(yatm_sect))
		vars()[ax].set_ylabel('y (m)', labelpad=1)

	ax3.set_xlabel('x (m)', labelpad=1)
	ax1.set_xticklabels([])
	ax2.set_xticklabels([])	


	#subplots_adjust(bottom=0.05, left=0.05, top = 0.99, right=0.9, hspace=0.22, wspace=0.05)

	savefig(figpath+'Crevasse_'+str(date)+'_'+str(section)+'_'+str(along_track_res)+'m.png', bbox_inches='tight', dpi=300)
	close(fig)


datapath='.,/Output/'
rawdatapath = '../../../../DATA/ICEBRIDGE/'
ATM_path = rawdatapath+'ATM/ARCTIC/'
dms_path = rawdatapath+'DMS/GR/'
posAV_path =rawdatapath+'/POSAV/CREVASSE/'
figpath='../Figures/'

# Projection
m=pyproj.Proj("+init=EPSG:3413")


def main(atmFile, date)

	#==================== Parameter settings ====================
	pwidth=20
	pint=5
	min_depth = 0.5
	xy_res=2
	min_feature_size=100
	along_track_res=500
	pts_threshold=(along_track_res**2)/65
	num_points_req = min_feature_size/(xy_res**2)
	section_num=0

	#==================== Get PosAV (flight location) data ====================
	posAV = loadtxt(posAVPath+'/sbet_'+str(date)+'.out.txt', skiprows=1)
	xp, yp, dist, sect_idxs, sect_utc_times = ut.get_pos_sections(posAV, m, along_track_res)

	#==================== Get ATM elevation data ====================
	lonT, latT, elevationT, aziT, utc_timeT= ut.get_atmqih5(atmFile, date[0:4], utc_time=1)
	xT, yT = m(lonT, latT)

	# Get PosAV indices coinciding with start/end of ATM file.
	start_i = np.abs(sect_utc_times - utc_timeT[0]).argmin()
	end_i = np.abs(sect_utc_times - utc_timeT[-1]).argmin()
	print 'START/END:', start_i, end_i

	feature_statsALL=np.array([]).reshape(0,9)

	for i in xrange(start_i -1, end_i + 1):
		section_num+=1
		
		mean_x, mean_y, mean_lat, mean_lon, mean_alt, mean_pitch, mean_roll, mean_vel = ut.posav_section_info(m, posAV[sect_idxs[i]:sect_idxs[i+1]], returnlonlat=1	)
		print '    '
		print str(i)+'/'+str(end_i + 1)
		print 'Mean altitude:', mean_alt
		print 'Mean pitch:', mean_pitch
		print 'Mean roll:', mean_roll
		print 'Mean vel:', mean_vel
		
		if (abs(mean_alt-500)<2000) & (abs(mean_pitch)<50) & (abs(mean_roll)<50):
			
			poly_path, vertices, sides = ut.get_pos_poly(xp, yp, sect_idxs[i], sect_idxs[i+1])
			xptsA = [vertices[0][0], vertices[1][0], vertices[2][0], vertices[3][0]]
			yptsA = [vertices[0][1], vertices[1][1], vertices[2][1], vertices[3][1]]
			xminA= np.amin(xptsA)
			xmaxA= np.amax(xptsA)
			yminA= np.amin(yptsA)
			ymaxA= np.amax(yptsA)

			xatm_sect, yatm_sect, elevation_sect = ut.get_atm_poly(xT, yT, elevationT, sect_utc_times, utc_timeT, poly_path, i)
			xatm_sect = xatm_sect - np.amin(xatm_sect)
			yatm_sect = yatm_sect - np.amin(yatm_sect)

			#ro.plot_atm_poly(m, xatm_sect, yatm_sect, elevation_sect, poly_path, i, out_path, year)

			num_pts_section = size(xatm_sect)
			print 'Num pts in section:', size(xatm_sect)
			#if there are more than 15000 pts in the 1km grid (average of around 20000) then proceed
			if  (num_pts_section>pts_threshold):	
			
				# generate a 2d grid from the bounds of the ATM coordinates of this section 
				xx2d, yy2d = ut.grid_atm(xatm_sect, yatm_sect, xy_res)
				print 'Grid:', size(xx2d[0]), size(xx2d[1])

				# Elevation relative to a fitted plane (order 2=quadratic plane)
				elevation2d = ut.grid_elevation(xatm_sect, yatm_sect,elevation_sect, xx2d, yy2d,kdtree=1)

				# Elevation relative to a fitted plane (order 2=quadratic plane)
				elevation2d_plane = ut.getPlaneElev(elevation2d, xx2d, yy2d, order=2)

				# Find local level (modal) surface
				level_elev, thresh, levpercent = ut.calc_level_ice(asarray(elevation2d_plane[ma.nonzero(elevation2d_plane)]), pint, pwidth, min_ridge_height)
				
				# Elevation anomalies relative to a local level (modal) surface
				elevation2d_anomalies=elevation2d_plane-level_elev
				
				# Elevation threhsold
				thresh=thresh-level_elev

				elevation2d_masked=ma.masked_where(elevation2d_anomalies<thresh, elevation2d_anomalies)
				feature_area = ma.count(elevation2d_masked)


				if (feature_area>0):
					found_features=1
					
					labelled_image  = ut.label_features(elevation2d_masked, xy_res, min_ridge_size, min_ridge_height)

					# Big feature is one that has an area after labelling above the minimum area threshold
					found_big_feature=0
					if (np.amax(labelled_image)>0):
						found_big_feature=1
						num_features=np.amax(labelled_image)
						print 'Number of features detected:', num_features
						#-------------- LABEL STATS ------------------
						feature_stats, feature_mesh, cov_matrix, index = ut.calc_feature_stats(elevation2d_masked, num_features, labelled_image, xx2d, yy2d, level_elev,section_num, calc_orientation=1)

				plotCrevasseSection(xatm_sect, yatm_sect, xx2d, yy2d, elevation_sect, elevation2d_plane, elevation2d_masked, labelled_image, date, section_num)
			
			else:
				print 'No data - WHY?! --------------'
				print 'Num pts in section:', size(xatm_sect)
		
		#ASSIGN BULK STATISTICS AS WE HAVE NOT CARRIED OUT RIDGE CALCULATION AS PLANE IS DOING FUNNY THINGS
		#bulk_statsT = calc_bulk_stats()
		feature_statsAll = vstack([feature_statsAll, feature_stats])
		
			if not os.path.exists(outpath+str(year)):
				os.makedirs(outpath+str(year))

	feature_statsAll.dump(outpath+str(date)+'/feature_stats_'+str(int(along_track_res/1000))+'km_xyres'+str(xy_res)+'m_'+str(int(min_feature_height*100))+'cm_poly'+str(atm_path_date[-9:-1])+'_f'+str(atm_file).zfill(3)+'.txt')
	
	#CAN OUTPUT AS TEXT FILES INSTEAD - BIGGER BUT CAN OPEN RAW
	#savetxt(outpath+str(year)+'/ridge_stats_'+str(int(along_track_res/1000))+'km_xyres'+str(xy_res)+'m_'+str(int(min_ridge_height*100))+'cm_poly'+str(atm_path_date[-9:-1])+'_f'+str(atm_file)+'.txt', ridge_statsALL)
	#savetxt(outpath+str(year)+'/cov_matrix_'+str(int(along_track_res/1000))+'km_xyres'+str(xy_res)+'m_'+str(int(min_ridge_height*100))+'cm_poly'+str(atm_path_date[-9:-1])+'_f'+str(atm_file)+'.txt', covarALL)
	#


if __name__ == '__main__':

	# Just do this for a given data for now, but could loop over all years of data here..
	date=20130405
	atmPathDate = rawdatapath+'/ATM/ARCTIC/'+str(date[0:4])+'/'+str(date)+'/'
	atmFile = ut.get_atm_files(atmPathDate, date[0:4])[0]

	main(date, atmFile)







