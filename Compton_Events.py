#!/usr/bin/env python2.7
##	VIDUSHI, 23rd Feb 2018 
#######################################
## Banana pixels are sorted out and Veto tagged events not present in dblevt file
## e.x.: %run Compton_Events.py /data2/czti/special/GRB180416B/pol_analysis/AS1A04_199T01_9000002040_13792cztM0_level2_quad_clean_200_4.dblevt AS1cztbadpix201609##08v01.fits /data2/czti/special/GRB180416B/pol_analysis/"AS1CZT_##"GRB180416B"_200_4" --tmark 261612609.279
##=================================================================================================================
import numpy as np
import argparse
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
import math
from matplotlib.backends.backend_pdf import PdfPages
## Input Informations:---------------------------------------------------------------------------------------------
print ('\n \n Required Input files: \n (1)	*_quad_clean.dblevt \n (2)	Caldb badpix file, AS1cztbadpix20160908v01.fits \n')

## ----------------------------------------------------------------------------------------------------------------
parser	=	argparse.ArgumentParser()
parser.add_argument("inputfile_double_event",type=str,help='Enter the path of quad clean double evt file')
parser.add_argument("inputfile_banana_pix",type=str,help='Enter the path of CALDB banana pixel fits file')
parser.add_argument("outfile", nargs="?", type=str, help="Stem to be prepended to all output files")
parser.add_argument("--tmark",type=float,help='Trigger time for Veto lightcurve')
parser.add_argument("--tbin",type=float,help='Binning time for Veto lightcurve, default=1.0',default=1.0)
parser.add_argument("--emin",type=float,help='Lower energy limit, default=100',default=100.0)
parser.add_argument("--emax",type=float,help='Upper energy limit, default=400',default=400.0)
parser.add_argument("-p", "--plottype", help="Type of output file to produce (png, pdf)", type=str, default='pdf')
args	=	parser.parse_args()
args.tmark
tmin = int(args.tmark) - 300.0*args.tbin
tmax = int(args.tmark) + 400.0*args.tbin
print('trigger time = %.2f, tmin = %.2f, tmax = %.2f	\n' %(args.tmark, tmin, tmax))
tbins = np.arange(tmin, tmax + args.tbin, args.tbin)

hdu	=	fits.open(args.inputfile_double_event)
#badfile		=	fits.open(args.inputfile_banana_pix)

Q_hist_sum	=	[0.0]*int((tmax-tmin)/args.tbin)
Q_max	=	[0.0]*1

if args.plottype == 'pdf':
    plotfile = PdfPages("{stem}_compton_events.pdf".format(stem=args.outfile))


for quad in range(4):
# Open and read the Input file:----------------------------------------------------------------------------------
	print("		For Quadrant %d:" %(quad))
	Time	= hdu[quad+1].data['Time']
	detX	= hdu[quad+1].data['detx']
	detY	= hdu[quad+1].data['dety']
	DetId	= hdu[quad+1].data['DetID']
	pixId	= hdu[quad+1].data['pixID']
	energy	= hdu[quad+1].data['ENERGY']
	evtick	= hdu[quad+1].data['CZTNTICK']
	tmin = max([np.floor(min(Time)), tmin])
	tmax = min([np.floor(max(Time)), tmax])
	select_region = np.where( (Time >= tmin) & (Time <= tmax) )
	Vtable	=	Table([Time[select_region], detX[select_region], detY[select_region], DetId[select_region], pixId[select_region], energy[select_region], evtick[select_region]],	names=('Time', 'detX', 'detY', 'DetID','pixID','Energy', 'cztnTick'))
	print("Event length = %d" % len(Vtable))
	"""
	## banana pix correction: Whichever detID and pixel ID has flag 1 in caldb badpix fits file, remove the same detID and pixel ID from the input file
	badtable=	badfile[quad+1].data	
	badind	=	np.where(badtable['PIX_FLAG']==1)
	baddet	=	(badtable['DETID'])[badind[0]]
	badpix	=	(badtable['PIXID'])[badind[0]]
	badflag	=	(badtable['PIX_FLAG'])[badind[0]]
	print('Length of badpix with flag 1 = %d'%len(badind[0]))
	
	goodevt	=	[]
	badevt	=	[]
	evtdet	=	Vtable['DetID']
	evtpix	=	Vtable['pixID']
	for i in range (0,len(evtdet)):
		flag = 0
		for j in range (0,len(baddet)):
			if(	(evtdet[i] == baddet[j]) and (evtpix[i] == badpix[j])	):
				badevt.append([i,j])
				flag=1
				break
		if(flag==0):
			goodevt.append(i)
			
	badevt	=	np.array(badevt)
	goodevt	=	np.array(goodevt)
	#print("banana pixel length = %d" % len(badevt))
	print("cleaned event length = %d" % len(goodevt))
	cleaned_table0	=	Vtable[goodevt]

	ab	=	np.loadtxt("q"+str(quad)+"_bad_table.txt")
	#fig.savefig(str(GRBNAME)+"_compton_events.pdf", bbox_inches='tight')
	remove_detId	=	ab[:,0]
	remove_pixId	=	ab[:,1]
	evtdet_check	=	cleaned_table0['DetID']
	evtpix_check	=	cleaned_table0['pixID']
	
	evt	=	[]
	remove_evt	=	[]
	for i in range(0, len(evtdet_check)):
		flag = 0
		for j in range(0, len(remove_detId)):
			if(	(evtdet_check[i] == remove_detId[j])	and	(evtpix_check[i] == remove_pixId[j])	):
				remove_evt.append([i,j])
				flag = 1
				break
		if (flag == 0):
			evt.append(i)
	print("cleaned event length after remove less photons modules = %d" % len(evt))

	cleaned_table	=	cleaned_table0[evt]"""
	#cleaned_table	=	Vtable[goodevt]
	cleaned_table	=	Vtable
	####		===============================================================================================
	####	COMPTON CRITERIA	: 
	#### (1)	Finding the adjacent pixels with condition that distance between 2 pixel, 
	####	1.0 cm			<=			distance	<=	1.8 cm
	#### (2)	1	<	absorbed energy(h) / scattered energy(l)	<	6 
	####	Reason: Scattering pixel will have low energy, whereas absorbing pixel will absorb higher energy photon.
	Time	=	cleaned_table['Time']
	DetX	=	cleaned_table['detX']
	DetY	=	cleaned_table['detY']
	match_1st	=	[]
	match_2nd	=	[]
	dist	=	[]
	diff_x	=	[]
	diff_y	=	[]
	for i in range(0,	len(cleaned_table)-1):
		if (	(Time[i+1] - Time[i])	<= 0.00003	):
			dx	=	abs(int(DetX[i+1]) - int(DetX[i]))
			dy	=	abs(int(DetY[i+1]) - int(DetY[i]))
			diff_x.append(dx)
			diff_y.append(dy)
			dist.append( 	((dx)**2	+	(dy)**2)**(0.5)	)
			match_1st.append(i)
			match_2nd.append(i+1)
	print("Same CZT N TICK events (time difference between two events less than 30 micro sec), counted once = %d" %len(dist))	

	dist_min	=	1.0
	dist_max	=	1.8
	ratio_min	=	1.0
	ratio_max	=	6.0

	Event_table_1st	=	cleaned_table[match_1st]
	detId_1st	=	Event_table_1st['DetID']
	Event_table_2nd =	cleaned_table[match_2nd]
	detId_2nd	=	Event_table_2nd['DetID']
	Distance	=	np.array(dist)
	#print(Distance)

	## COMPTON EVENTS SORT OUT:
	adj_pix_stamp	=	[]
	remained_pix	=	[]
	for i in range(0, len(Event_table_1st)):
		#if ((dist_min <= Distance[i] <= dist_max) and (int(detId_1st[i]) == int(detId_2nd[i]))):
		if ((diff_x[i] > 0.0 or diff_y[i] > 0.0) and (diff_x[i] <= 1.0 and diff_y[i] <= 1.0) and (int(detId_1st[i]) == int(detId_2nd[i]))):
			adj_pix_stamp.append(i)
		else:
			remained_pix.append(i)
	print("Adjacent pixel events and same detId, 1 pair counted once = %d" %len(adj_pix_stamp))

	adj_pix_table1	=	Event_table_1st[adj_pix_stamp]
	Time_1	=	adj_pix_table1['Time']
	DetID_1	=	adj_pix_table1['DetID']
	PixID_1	=	adj_pix_table1['pixID']
	Energy_1	=	adj_pix_table1['Energy']
	DetX_1	=	adj_pix_table1['detX']
	DetY_1	=	adj_pix_table1['detY']
	dist_1_2	=	Distance[adj_pix_stamp]

	adj_pix_table2	=	Event_table_2nd[adj_pix_stamp]
	Time_2	=	adj_pix_table2['Time']
	DetID_2	=	adj_pix_table2['DetID']
	PixID_2	=	adj_pix_table2['pixID']
	Energy_2	=	adj_pix_table2['Energy']
	DetX_2	=	adj_pix_table2['detX']
	DetY_2	=	adj_pix_table2['detY']
	#print(len(dist_1_2), len(adj_pix_table1), len(adj_pix_table2))

	##### DECIDE WHICH PIX IS SCATTERED AND WHICH IS FOR ABSORBER IN EACH PAIR
	abs_1st_stamp	=	[]
	abs_2nd_stamp	=	[]
	for i in range(0, len(adj_pix_table1)):
		if (Energy_1[i]	>	Energy_2[i]):
			abs_1st_stamp.append(i)
		else:
			abs_2nd_stamp.append(i)	
	#print(len(abs_1st_stamp), len(abs_2nd_stamp))
	
	high_abs =	[]
	high_abs.extend(Energy_1[abs_1st_stamp])
	high_abs.extend(Energy_2[abs_2nd_stamp])
	high_abs	=	np.array(high_abs)

	low_scat	=	[]
	low_scat.extend(Energy_2[abs_1st_stamp])
	low_scat.extend(Energy_1[abs_2nd_stamp])
	low_scat	=	np.array(low_scat)

	total_energy	=	high_abs + low_scat
	ratio_h_by_l	=	high_abs/low_scat

	Time	=	[]
	Time.extend(Time_1[abs_1st_stamp])
	Time.extend(Time_2[abs_2nd_stamp])
	Time	=	np.array(Time)
	print("Time corresponding high absorber pixel, counted once = %d" % len(Time))

	evtdist	=	[]
	evtdist.extend(dist_1_2[abs_1st_stamp])
	evtdist.extend(dist_1_2[abs_2nd_stamp])
	evtdist	=	np.array(evtdist)
	#####################################
	h_detX	=	[]
	h_detX.extend(DetX_1[abs_1st_stamp])
	h_detX.extend(DetX_2[abs_2nd_stamp])
	h_detX	=	np.array(h_detX)
	h_detY	=	[]
	h_detY.extend(DetY_1[abs_1st_stamp])
	h_detY.extend(DetY_2[abs_2nd_stamp])
	h_detY	=	np.array(h_detY)
	######################################
	l_detX	=	[]
	l_detX.extend(DetX_2[abs_1st_stamp])
	l_detX.extend(DetX_1[abs_2nd_stamp])
	l_detX	=	np.array(l_detX)
	l_detY	=	[]
	l_detY.extend(DetY_2[abs_1st_stamp])
	l_detY.extend(DetY_1[abs_2nd_stamp])	
	l_detY	=	np.array(l_detY)
	#####################################

	comp_eve_stamp	=	[]
	for i in range(0, len(Time)):
		if(  (ratio_min < ratio_h_by_l[i] < ratio_max) and (args.emin <= total_energy[i] <= args.emax) and (low_scat[i] >= 0.0) ):
			comp_eve_stamp.append(i)

	print("Compton Events = %d" %len(comp_eve_stamp))
	##=============================================================================================================
	time	= 	Time[comp_eve_stamp]
	evtdist	=	evtdist[comp_eve_stamp]
	high_abs	=	high_abs[comp_eve_stamp]
	high_detX	=	h_detX[comp_eve_stamp]
	high_detY	=	h_detY[comp_eve_stamp]
	low_scat	=	low_scat[comp_eve_stamp]
	low_detX	=	l_detX[comp_eve_stamp]
	low_detY	=	l_detY[comp_eve_stamp]
	Compton_table	=	Table([time, evtdist, high_abs, high_detX, high_detY, low_scat, low_detX, low_detY],	names=('Time', 'Distance', 'HighE_or_abs', 'Abs_detX', 'Abs_detY', 'lowE_or_scat', 'Scat_detX', 'Scat_detY'))
	#print(Compton_table)
	Q1_hist, bin_edges = np.histogram(Compton_table['Time'], bins=tbins)
	Q_max	=	Q_max + max(Q1_hist)
	##=============================================================================================================
	Q_hist_sum = Q_hist_sum + Q1_hist
	center_time = (bin_edges[:-1] + bin_edges[1:])/2.0

##======================================================================================================================
Q_hist_sum_err	=	np.sqrt(Q_hist_sum)
#print(len(Q_hist_sum),len(Q_hist_sum_err))

all_Q_max	=	max(Q_hist_sum)
peak_stamp	=	np.where([Q_hist_sum == all_Q_max])[0]
peak_time	=	center_time[peak_stamp]
#print ("Peak compton count = %9.2f , Peak Time = %9.2f \n" %(all_Q_max, peak_time))

########################	SELECT PRE BACKGROUND INTERVAL		###################################
print("\n SELECT PRE GRB BACKGROUND INTERVAL *********************")
# Simple mouse click function to store coordinates
def onclick(event):
	global ix, iy
	ix, iy = event.xdata, event.ydata
	# assign global variable to access outside of function
	global coords
	coords.append((ix, iy))
	# Disconnect after 2 clicks
        if len(coords) == 2:
        	fig.canvas.mpl_disconnect(cid)
        	plt.close(1)
        return

fig = plt.figure(1)
ax = fig.add_subplot(111)
plt.step(center_time, Q_hist_sum)
ax.set_xlim(tmin, tmax)
coords = []
# Call click func
cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.show(1)
pre_bkg_min	=	coords[0][0]
pre_bkg_max	=	coords[1][0]
print("Pre-bkg start = %.2f , Pre-bkg stop = %.2f" %(pre_bkg_min, pre_bkg_max))

select_prebkg = np.where( (center_time >= pre_bkg_min) & (center_time <= pre_bkg_max) )
pre_grb	=	Q_hist_sum[select_prebkg]
pre_time	=	center_time[select_prebkg]
pre_time_interval	=	pre_time[-1] - pre_time[0]

pre_sum	=	0.0
for i in range(0, len(pre_grb)):
	pre_sum	=	pre_sum + pre_grb[i]	

pre_bkg_rate	=	pre_sum/float(pre_time_interval)
print("Pre-bkg total compton counts = %.2f, Pre-bkg time = %.2f , Pre-bkg rate (cts/s) = %.2f " %(pre_sum, pre_time_interval, pre_bkg_rate))

########################	SELECT GRB INTERVAL		#####################################
print("\n SELECT THE GRB INTERVAL ********************************")
def onclick(event):
	global ix1, iy1
	ix1, iy1 = event.xdata, event.ydata
	global coords_grb
	coords_grb.append((ix1, iy1))
	if len(coords_grb) == 2:
        	fig.canvas.mpl_disconnect(cid)
        	plt.close(1)
        return

fig = plt.figure(1)
ax = fig.add_subplot(111)
plt.step(center_time, Q_hist_sum)
ax.set_xlim(pre_bkg_max-10.0**args.tbin,tmax)
coords_grb = []
cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.show(1)
grb_start=coords_grb[0][0]
grb_stop=coords_grb[1][0]
print("GRB start time = %.2f , GRB stop time = %.2f" %(grb_start, grb_stop))

select_grb = np.where( (center_time >= grb_start) & (center_time <= grb_stop) )
hist_grb	=	Q_hist_sum[select_grb]
grb_time	=	center_time[select_grb]
grb_time_interval	=	grb_time[-1] - grb_time[0]
grb_sum	=	0.0
for i in range(0, len(hist_grb)):
	grb_sum	=	grb_sum + hist_grb[i]	

grb_rate	=	grb_sum/float(grb_time_interval)
print("Total compton events in GRB interval = %.2f, GRB Time = %.2f , GRB Rate (cts/s) = %.2f " %(grb_sum, grb_time_interval, grb_rate))

########################	SELECT POST GRB INTERVAL		#####################################
print("\n SELECT THE POST GRB BACKGROUND INTERVAL ****************")
def onclick(event):
	global ix2, iy2
	ix2, iy2 = event.xdata, event.ydata
	global coords2
	coords2.append((ix2, iy2))
	if len(coords2) == 2:
        	fig.canvas.mpl_disconnect(cid)
        	plt.close(1)
        return

fig = plt.figure(1)
ax = fig.add_subplot(111)
plt.step(center_time,Q_hist_sum)
ax.set_xlim(grb_stop-10.0**args.tbin,tmax)
coords2 = []
# Call click func
cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.show(1)
post_bkg_min=coords2[0][0]
post_bkg_max=coords2[1][0]
print("Post-bkg start = %.2f , Post-bkg stop = %.2f" %(post_bkg_min, post_bkg_max))

select_postgrb = np.where( (center_time >= post_bkg_min) & (center_time <= post_bkg_max) )
post_grb	=	Q_hist_sum[select_postgrb]
post_time	=	center_time[select_postgrb]
post_time_interval	=	post_time[-1] - post_time[0]
post_sum	=	0.0
for i in range(0, len(post_grb)):
	post_sum	=	post_sum + post_grb[i]	

post_bkg_rate	=	post_sum/float(post_time_interval)
print("Post-bkg total compton counts = %.2f, Post-bkg time = %.2f , Pre-bkg rate (cts/s) = %.2f " %(post_sum, post_time_interval, post_bkg_rate))

avg_bkg_rate	=	(	pre_bkg_rate + post_bkg_rate	)/2.0
print("Average background rate (cts/s) = %.2f" %avg_bkg_rate)
True_grb_counts	=	(	grb_rate - avg_bkg_rate	)*grb_time_interval
print("Background Corrected compton events rate (cts/s) = %.2f" %True_grb_counts)
bkg_contrast=float(grb_rate)/float(avg_bkg_rate)
print("Background contrast	= %.2f" % bkg_contrast)

outfile= args.outfile
#outfile.split()
GRBNAME=outfile.split('/')[4]
GRB_thresh=outfile.split('/')[6]
print(GRBNAME)

s = open("/data2/czti/special/"+str(GRBNAME)+"/"+str(GRBNAME)+".info", "a") #Ram args.outfile
s.write("\n"+str(GRB_thresh)+":"+"\n" +"Background Corrected compton events rate (cts/s) : "+ str(True_grb_counts) + "\n"+"Background contrast : " + str(bkg_contrast) +"\n" )
s.close() #Ram




fig = plt.figure()
plt.errorbar(center_time, Q_hist_sum , yerr = Q_hist_sum_err, color='b', fmt="*", ecolor='r', elinewidth=1.2, capsize=4, ms=11.0, label='Compton Events')
plt.step(center_time + 0.6*args.tbin, Q_hist_sum, 'b', lw = '1.75' )
plt.xlim((grb_start - 150.0*args.tbin), (grb_stop + 100.0*args.tbin))
plt.ylim(0, (all_Q_max + 0.25*all_Q_max))
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlabel('AstroSat Time (s)', color='red', fontweight= 'bold', fontsize=16)
plt.ylabel('Counts/s', color='red', fontweight= 'bold', fontsize=16)
plt.title('Compton Event Lightcurve', fontweight= 'bold', fontsize=16)
plt.axvline(grb_start, color='k', lw='1.25')
plt.axvline(grb_stop, color='k',  lw='1.25')
plt.axvline(args.tmark , color='r',ls='--' , lw='1.3')
plt.axvspan(grb_start, grb_stop, alpha=0.25, color='y')
plt.text((args.tmark -75.0*args.tbin),(all_Q_max+0.20*all_Q_max), 'Background Corrected Compton Events = %d' % True_grb_counts, color='red', fontsize=9, fontweight='bold')
plt.text((args.tmark -60.0*args.tbin),(all_Q_max+0.135*all_Q_max), 'S/N = %.2f' % bkg_contrast, color='red', fontsize=10, fontweight='bold')
#plt.legend(loc="best")
#plt.show()
#GRBNAME=args.GRBNAME
#fig.savefig(str(GRBNAME)+"_compton_events.pdf", bbox_inches='tight')
if args.plottype == 'pdf':

	plotfile.savefig()
else:
	plt.show()	
	plt.savefig(args.outfile + "_compton_events." + args.plottype)
plt.show()	
if args.plottype == 'pdf':
    plotfile.close()

hdu.close()
#badfile.close()
