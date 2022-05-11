##================================================================================================
#! /usr/bin/python2.7
## Calculate T90
## By Vidushi Sharma, 17 Dec 2017
##================================================================================================
## Step 1: Banana pixel sorted out and livetime corrected
## e.x.: %run Calculate_T90_CZTI.py GRB190928A/AS1CZT_GRB190928A_quad_clean.evt GRB190928A/AS1A06_002T02_9000003206_21631cztM0_level2_quad_badpix.fits GRB190928A/AS1A06_002T02_9000003206_21631cztM0_level2_quad_livetime.fits GRB190928A --tmark 307372337.0 
## e.x.: %run Calculate_T90_CZTI.py GRB181201A/AS1CZT_GRB181201A_quad_clean.evt GRB181201A/AS1A05_185T04_9000002548_17169cztM0_level2_quad_badpix.fits GRB181201A/AS1A05_185T04_9000002548_17169cztM0_level2_quad_livetime.fits GRB181201A --tmark 281327880.0 
##=================================================================================================

import numpy as np
import os,subprocess,argparse 
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from astropy.io import fits
import math
from scipy.ndimage import gaussian_filter1d
from matplotlib.backends.backend_pdf import PdfPages
import warnings
warnings.simplefilter('ignore', np.RankWarning)

# Input Informations:-----------------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("inputfile_single_event",type=str,help='Enter the path of quad clean evt file')
parser.add_argument("inputfile_banana_pix",type=str,help='Enter the path of CALDB banana pixel fits file')
parser.add_argument("inputfile_livetime",type=str,help='Enter the path of CALDB banana pixel fits file')
parser.add_argument("outfile", nargs="?", type=str, help="Stem to be prepended to all output files")
parser.add_argument("GRBNAME", type=str, help='enter the GRB NAME')
parser.add_argument("--tmark",type=float,help='Trigger time for lightcurve')
parser.add_argument("--tbin",type=float,help='Binning time for lightcurve, default=1.0',default=1.0)
parser.add_argument("--badpix",type=str,help='For badpix correction',default="0")
parser.add_argument("--livetime",type=str,help='For livetime correction on 1 s lightcurve',default="0")
parser.add_argument("-p", "--plottype", help="Type of output file to produce (png, pdf)", type=str, default='pdf')

args = parser.parse_args()
args.tmark
tmin = int(args.tmark) - 500.0*args.tbin
tmax = int(args.tmark) + 500.0*args.tbin
print ('\n \n Required Input files: \n (1)	*_quad_clean.dblevt \n (2)	Caldb badpix file, AS1cztbadpix20160908v01.fits \n (3)	Livetime fits file according tbin or type of GRB	\n')
print('bin time = %.2f, trigger time = %.2f, tmin = %.2f, tmax = %.2f	\n' %(args.tbin, args.tmark, tmin, tmax))
tbins = np.arange(tmin, np.float(tmax+0.0000001) , args.tbin)

evtfile	=	fits.open(args.inputfile_single_event)
badfile	=	fits.open(args.inputfile_banana_pix)
livefile	=	fits.open(args.inputfile_livetime)

#------------------------------- Badpix_correction --------------------------------
def goodpix_stamps(badpix_file, quadrant_no, event_table):
		badtable=	badpix_file[quadrant_no].data	
		badind	=	np.where(badtable['PIX_FLAG']==1)
		baddet	=	(badtable['DETID'])[badind[0]]
		badpix	=	(badtable['PIXID'])[badind[0]]
		badflag	=	(badtable['PIX_FLAG'])[badind[0]]
		print('Length of badpix with flag 1 = %d'%len(badind[0]))

		goodevt	=	[]
		badevt	=	[]
		for i in range (0,len(evtdet)):
			flag = 0
			for j in range (0,len(baddet)):
				if((evtdet[i]==baddet[j]) and (evtpix[i]==badpix[j])):
					badevt.append([i,j])
					flag=1
					break
			if(flag==0):
				goodevt.append(i)
			
		badevt	=	np.array(badevt)
		goodevt	=	np.array(goodevt)
		print("cleaned event length = %d, banana pixel length = %d" % (len(goodevt), len(badevt)))
		evtable	=	event_table[goodevt]
		return evtable
#------------------------------------------------------------------------------------

##========== LIVETIME correction, done for 1 s lightcurve: division of 1 s collected counts by fracexp (In livetime fits)
def livetime_corr_counts(livetime_file, quadrant_no):
	livetable	=	livefile[Q_n].data
	time	=	livetable['TIME']
	fracexp	=	livetable['FRACEXP']
	ltable	=	Table([time, fracexp], names=('Time', 'Fracexp'))

	lstart_match 	=	np.where(time >= tmin)[0]
	lstart	=	lstart_match[0]
	lstop_match 	=	np.where(time >= tmax)[0]
	lstop	=	lstop_match[0]

	evtlive	=	ltable[lstart:lstop]
	ltime	=	evtlive['Time']
	lfrac	=	evtlive['Fracexp']
	return lfrac
##==============================================================================================

####################### Lightcurve Counts: summing all quadrants ###############################
hist	=	[0.0]*int((tmax-tmin)/args.tbin)
for Q_n in range(1,5):
# Open and read the Input file: single event file
	print("		For Quadrant %d:" % Q_n)
	evtable	=	evtfile[Q_n].data
	evtime	=	evtable['Time']
	evtdet	=	evtable['DetID']
	evtpix	=	evtable['pixID']
	evtene	=	evtable['ENERGY']
	Vtable	=	Table([evtime,evtdet,evtpix,evtene], names=('Time','DetID','pixID','Energy'))
	tmin = max([np.floor(min(evtime)), tmin])
	tmax = min([np.floor(max(evtime)), tmax])
	start_match 	=	np.where(evtime >= tmin)[0]
	start	=	start_match[0]
	stop_match	=	np.where(evtime >= tmax)[0]
	stop	= stop_match[0]

	evtable	=	Vtable[start:stop]
	evtime	=	evtable['Time']
	evtdet	=	evtable['DetID']
	evtpix	=	evtable['pixID']
	energy	=	evtable['Energy']
	print("Event length = %d" % len(evtime))
	
	if(args.badpix == "1"):
		Event_table	=	goodpix_stamps(badfile, Q_n, evtable)
		print ('Badpixels are removed \n')
	else:
		Event_table	=	evtable
	
	Q1_hist, bin_edges = np.histogram(Event_table['Time'] , bins=tbins)
	center_time = (bin_edges[:-1] + bin_edges[1:])/2.0 
	Q1_hist	=	Q1_hist/float(args.tbin)
	
	if(args.livetime == "1"):
		Q1_hist	=	Q1_hist/livetime_corr_counts(livefile, Q_n)
		print ('Livetime correction done \n')
	else:
		Q1_hist	=	Q1_hist
	
	hist = hist + Q1_hist
	center_time = (bin_edges[:-1] + bin_edges[1:])/2.0 
############################################################################################

Counts	=	hist
Time	=	center_time
print(Counts[0],Time[0])	

############# SELECT PRE BACKGROUND + GRB + POST BACKGROUND INTERVALS #############
print("\n********* SELECT PRE GRB BACKGROUND START TO STOP INTERVAL ************")
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
plt.step(Time, Counts)
plt.xlabel('Time (AstroSat seconds)')
plt.ylabel('Count Rate')
plt.axvline([args.tmark], color='k', linestyle='dashed', label='Marked times',alpha=0.7)
ax.set_xlim(tmin, tmax)
coords = []
# Call click func
cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.show(1)
pre_bkg_min=coords[0][0]
pre_bkg_stop=coords[1][0]
pre_bkg_interval = pre_bkg_stop-pre_bkg_min	
print("Pre-bkg start = %.3f, Pre-bkg stop = %.3f, Pre-bkg Interval = %.3f" %(pre_bkg_min, pre_bkg_stop, pre_bkg_interval))

print("\n*************** SELECT THE GRB START TO STOP INTERVAL *****************")
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
plt.step(Time, Counts)
plt.xlabel('Time (AstroSat seconds)')
plt.ylabel('Count Rate')
plt.axvline([args.tmark], color='k', linestyle='dashed', label='Marked times',alpha=0.7)
ax.set_xlim(pre_bkg_min,tmax)
coords_grb = []
cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.show(1)
grb_start=coords_grb[0][0]
grb_stop=coords_grb[1][0]
grb_interval = 	grb_stop-grb_start
print("GRB start = %.3f , GRB stop = %.3f, GRB Interval = %.3f" %(grb_start,grb_stop, grb_interval))

## SELECT GRB INTERVAL
print("\n**************** SELECT THE POST BKG START TO STOP INTERVAL ****************")
def onclick(event):
	global ix2, iy2
	ix2, iy2 = event.xdata, event.ydata
	global coords_2
	coords_2.append((ix2, iy2))
	if len(coords_2) == 2:
        	fig.canvas.mpl_disconnect(cid)
        	plt.close(1)
        return

fig = plt.figure(1)
ax = fig.add_subplot(111)
plt.step(Time, Counts)
plt.xlabel('Time (AstroSat seconds)')
plt.ylabel('Count Rate')
plt.axvline([args.tmark], color='k', linestyle='dashed', label='Marked times',alpha=0.7)
ax.set_xlim(pre_bkg_min,tmax)
coords_2 = []
cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.show(1)
post_bkg_start=coords_2[0][0]
post_bkg_stop=coords_2[1][0]	
post_bkg_interval = post_bkg_stop-post_bkg_start
print("Post-bkg start = %.3f , Post-bkg stop = %.3f, Post-bkg Interval = %.3f" %(post_bkg_start, post_bkg_stop, post_bkg_interval))
##############################################################################################
def time_match(bin_time, match_time):
	s = np.where(bin_time >= match_time)[0]
	return s[0]	
t1	=	time_match(center_time, pre_bkg_min)
t2	=	time_match(center_time, pre_bkg_stop)

tgrb1	=	time_match(center_time, grb_start)
tgrb2	=	time_match(center_time, grb_stop)

t3	=	time_match(center_time, post_bkg_start)
t4	=	time_match(center_time, post_bkg_stop)

bkg_time=np.hstack((Time[t1:t2],Time[t3:t4])) 
bkg_counts=np.hstack((Counts[t1:t2],Counts[t3:t4])) 

grb_time=np.hstack((Time[tgrb1:tgrb2])) 
grb_counts=np.hstack((Counts[tgrb1:tgrb2]))

Time_sel = np.hstack((Time[t1:t4]))
Counts_sel = np.hstack((Counts[t1:t4]))
##############################################################################################

def t90_accum(Time, Counts, bkg_time, bkg_counts, grb_time, grb_counts, tgrb1, tgrb2):
	pfit = np.polyfit(bkg_time, bkg_counts, 2)  ## fitting with quadratic function to background
	bkg=np.polyval(pfit,Time)
	subs_data=Counts-bkg
	result=gaussian_filter1d(subs_data,sigma=0.001, axis=-1, order=0)
	
	mean_bkg = np.mean(bkg)
	peak_count=max(grb_counts)-mean_bkg;
	peak_time=grb_time[np.where(grb_counts==max(grb_counts))[0][0]]

	t_90_data=result[tgrb1:tgrb2]  
	t_90_accum=np.add.accumulate(t_90_data)
	t_90_time=Time[tgrb1:tgrb2] 
	max_acc=max(t_90_accum)*args.tbin

	fracout = np.array([max_acc*0.05, max_acc*0.95])
	tout = np.interp(fracout, t_90_accum, t_90_time)
	t90 = tout[1] - tout[0]

	return [mean_bkg, peak_count, peak_time, max_acc, t90, bkg, result, subs_data, t_90_accum, t_90_time, tout[0], tout[1]]

accurate_para = t90_accum(Time, Counts, bkg_time, bkg_counts, grb_time, grb_counts, tgrb1, tgrb2)	
print ("Mean background rate = %.3f" %accurate_para[0])
print ("Peak count rate above background = %.3f , at Time = %.3f s astrosat seconds" % (accurate_para[1], accurate_para[2]))
print ("Total accumulated counts = %.3f, for T90 interval = %.3f s\n" % (accurate_para[3], accurate_para[4]))

bkg = t90_accum(Time, Counts, bkg_time, bkg_counts, grb_time, grb_counts, tgrb1, tgrb2)[5]
result = t90_accum(Time, Counts, bkg_time, bkg_counts, grb_time, grb_counts, tgrb1, tgrb2)[6]
subs_data = t90_accum(Time, Counts, bkg_time, bkg_counts, grb_time, grb_counts, tgrb1, tgrb2)[7]
t_90_accum = t90_accum(Time, Counts, bkg_time, bkg_counts, grb_time, grb_counts, tgrb1, tgrb2)[8]
t_90_time = t90_accum(Time, Counts, bkg_time, bkg_counts, grb_time, grb_counts, tgrb1, tgrb2)[9]
t_90_start = t90_accum(Time, Counts, bkg_time, bkg_counts, grb_time, grb_counts, tgrb1, tgrb2)[10]
t_90_stop = t90_accum(Time, Counts, bkg_time, bkg_counts, grb_time, grb_counts, tgrb1, tgrb2)[11]


##============================== Simulations using Gaussian distributed Poisson errors on Counts==================================================

sim_runs = 10000

#Gaussian distributed Poisson errors
def fake_lc(counts):
    random_nums = np.random.randn(len(counts))
    random_signs = random_nums/abs(random_nums)
    counts_err = np.sqrt(counts)
    counts_err_fak = np.random.poisson(counts_err)
    counts_fak = counts + random_signs*counts_err_fak
    return counts_fak

bkg_matrix = [fake_lc(bkg_counts) for row in range(0, sim_runs)]
grb_matrix = [fake_lc(grb_counts) for row in range(0, sim_runs)]



sim_parameters = []
for i in range(0, sim_runs):
	bkg_sim_ct = np.array(bkg_matrix[i])
	grb_sim_ct = np.array(grb_matrix[i])
	para = t90_accum(Time, Counts, bkg_time, bkg_sim_ct, grb_time, grb_sim_ct, tgrb1, tgrb2)
	sim_parameters.append(para)

all_para_mean = []
all_para_std = []
for j in range(0, 5):
	same_element = []
	for i in range(0, sim_runs):
		same_element.append(sim_parameters[i][j])
	all_para_mean.append(np.mean(same_element))
	all_para_std.append(np.std(same_element))

print ("Parameter values & their uncertainties: Mean and Standard deviation for 5000 runs.")
print ("Bkg & grb counts are drawn randomly from a Poisson distribution obtained using true lightcurve.")
mean_bkg = all_para_mean[0]
mean_bkg_err = all_para_std[0]
print ("Mean background rate = %d $\pm$ %d" % (mean_bkg, mean_bkg_err))
peak_count = all_para_mean[1]
peak_count_err = all_para_std[1]
peak_time = all_para_mean[2]
print ("Peak count rate above background = %d $\pm$ %.1f, At Time = %.2f $\pm$ %.3f s AstroSat seconds" % (peak_count, peak_count_err, peak_time, all_para_std[2]))
max_acc = all_para_mean[3] 
max_acc_err = all_para_std[3]
t90_val = all_para_mean[4]
t90_val_err = all_para_std[4]
print ("Total accumulated counts = %d $\pm$ %d, With T90 interval = %.1f $\pm$ %.2f s" % (max_acc, max_acc_err, t90_val, t90_val_err))
print ("T90 Start = %.3f, T90 Stop = %.3f AstroSat seconds" % (t_90_start, t_90_stop))
##================================================================================================================================

def gen_plot(tmark, t1, t2, t3, t4, peak_count, peak_count_err, peak_time, max_acc, max_acc_err, t90_val, t90_val_err, mean_bkg, mean_bkg_err, Time, Counts,  bkg, result, subs_data, t_90_start, t_90_stop):
	plt.rc('text', usetex=True)
	if args.plottype == 'pdf':
		plotfile = PdfPages("{stem}_T90.pdf".format(stem=args.outfile))
 	# Set the font dictionaries (for plot title and axis titles)
	title_font = {'fontname':'Arial', 'size':'15', 'color':'black', 'weight':'normal',
              'verticalalignment':'bottom'} # Bottom vertical alignment for more space
	axis_font = {'fontname':'Arial', 'size':'12'}
	plt.figure(figsize=(8,7))
	plt.subplot(221)
	plt.title('Region Selection', **title_font)
	plt.step(bkg_time, bkg_counts, color='b')
	plt.step(grb_time, grb_counts, color='cyan', lw =1.5)
	plt.axvline(grb_time[0], color='darkcyan', lw =0.8, alpha=0.8)
	plt.axvline(grb_time[-1], color='darkcyan',  lw =0.8, alpha=0.8)
	plt.axvline(Time[t1], color='darkblue', ls='-.', alpha=0.8)
	plt.axvline(Time[t2], color='darkblue', ls='-.', alpha=0.8)
	plt.axvline(Time[t3], color='darkblue', ls='-.', alpha=0.8)
	plt.axvline(Time[t4], color='darkblue', ls='-.', alpha=0.8)
	plt.axvline([args.tmark], color='r', linestyle='dashed', label='Marked times', alpha=0.8, lw =2.0)
	plt.xlabel('AstroSat Time (s)', **axis_font)
	plt.ylabel('Count Rate', **axis_font)
	#---------------------------------------------------------------

	plt.subplot(222)
	plt.title('Cumulative Counts', **title_font)
	plt.tight_layout()
	plt.plot(t_90_time, t_90_accum, color = 'green', lw =1.75)
	plt.ylabel('Accumulated Count (Counts/sec)', **axis_font)
	plt.xlabel('AstroSat Time (s)', **axis_font)
	plt.axvline(t_90_start, color='darkgreen', ls='-.', lw = 1.75)
	plt.axvline(t_90_stop, color='darkgreen', ls='-.', alpha=1.75)
	plt.text(args.tmark, max(t_90_accum)-0.1*max(t_90_accum), r"\bf{$T_{90} = %.1f \pm %.2f$ s}" % (t90_val, t90_val_err), **axis_font)   
	plt.text(args.tmark, max(t_90_accum)-0.3*max(t_90_accum), r"\bf{$T_{90} \; Start \; Time= %.3f $ s}" % (t_90_start), **axis_font)  
	plt.axvspan(t_90_start, t_90_stop, alpha=0.1, color='y') 
	plt.axvline([args.tmark], color='r', linestyle='dashed', label='Marked times', alpha=0.8, lw =2.0)

	if args.plottype == 'pdf':
		plotfile.savefig()
	else:
		plt.show()	
		plt.savefig(args.outfile + "_T90." + args.plottype)
	#---------------------------------------------------------------

	plt.subplot(212)
	plt.ylabel('Count Rate', **axis_font)
	plt.xlabel('AstroSat Time (s)', **axis_font)
	plt.plot(Time, Counts)
	plt.xlim(Time[t1-3], Time[t4+3])
	plt.plot(Time, bkg, color='orange')
	plt.plot(Time, result, color='cyan')
	plt.plot(Time, subs_data, color='purple')
	plt.text(Time[t1+3], peak_count, r" Peak Count Rate = %d $\pm$ %.1f, at Astrosat time = %.3f s " % (peak_count, peak_count_err, peak_time), **axis_font)
	plt.text(Time[t1+3], peak_count-0.25*peak_count, r"Total Counts = %d $\pm$ %.1f " % (max_acc, max_acc_err), **axis_font)
	plt.text(Time[t1+3], peak_count-0.5*peak_count, r"Mean background Count = %d $\pm$ %.1f " % (mean_bkg, mean_bkg_err), **axis_font) 
	plt.grid()
	plt.axvline(t_90_start, color='darkgreen', ls='-.', lw = 1.75)
	plt.axvline(t_90_stop, color='darkgreen', ls='-.', alpha=1.75)
	plt.axvspan(t_90_start, t_90_stop, alpha=0.1, color='y')
	if args.tmark > 0: plt.axvline([args.tmark], color='r', linestyle='dashed', alpha=0.8, lw =2.0)
	plt.show()	
	if args.plottype == 'pdf':
    		plotfile.close()



gen_plot(args.tmark,  t1, t2, t3, t4, peak_count, peak_count_err, peak_time, max_acc, max_acc_err, t90_val, t90_val_err, mean_bkg, mean_bkg_err, Time,Counts, bkg, result, subs_data, t_90_start, t_90_stop)



