#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

import numpy
import matplotlib
from matplotlib import pyplot
from mpl_toolkits.axes_grid1 import make_axes_locatable
import csv
from MouseViewParser.readers import gorilla
import pandas as pd
#import pickle

# # # # #
# CONSTANTS

# BASIC CONTROLS
# Specify the file name.
#FILENAME = "mouseview_v3.csv"
# Overwrite the temporary data?
OVERWRITE_TMP = True

# EXPERIMENT SETTINGS
# Specify the name of the zone that starts the first trial.
TRIAL_START_ZONE = "Zone3"
# Specify the custom fields that should be loaded from the data.
CUSTOM_FIELDS = ["left_stim", "right_stim", "condition", "side"]

#at first I put them all in, but that won't work with later code. "Neutral" stimulus in pair is not listed in STIMULI
# Stimuli for each condition.
# STIMULI = { \
#     "1_ClothedWomanNudeWoman":  ["A1", "A11", "A2", "A22", "A3","A33","A4","A44"], \
#     "1_ClothedManNudeMan": ["B1", "B11", "B2", "B22", "B3","B33","B4","B44"], \
#     "2_NudeWomanNudeMan": ["D1", "D11", "D2", "D22", "D3","D33","D4","D44"], \
#     "2_ClothedWomanClothedMan": ["C1", "C11", "C2", "C22", "C3","C33","C4","C44"] \
#     }
STIMULI = { \
     "1_ClothedWomanNudeWoman": ["A1","A2","A3","A4"], \
     "1_ClothedManNudeMan": ["B1","B2", "B3","B4"], \
     "2_NudeWomanNudeMan": ["D1","D2","D3","D4"], \
     "2_ClothedWomanClothedMan": ["C1","C2","C3","C4"] \
     }

STIM_width = 300
STIM_height = 400
STIM_ratio = STIM_height/STIM_width
    # Conditions.
CONDITIONS = list(STIMULI.keys())
CONDITIONS.sort()
#CONDITIONS = ["1_ClothedWomanNudeWoman","1_ClothedManNudeMan","2_NudeWomanNudeMan","2_ClothedWomanClothedMan"]

CONDITIONS_2 = { \
    "Clothed-Nude": ["1_ClothedWomanNudeWoman","1_ClothedManNudeMan"], \
    "Pref-Nonpref": ["2_NudeWomanNudeMan","2_ClothedWomanClothedMan"] \
    }
            

# Areas of interest. This will be a little tricky. It will need to be redefined per condition. 
# Perhaps I can do figures, variables in a loop, with two conditions per cycle like before,
# and then just add conditional logic translating affective/neutral to current context
AOI = ["affective", "neutral", "other"]
# Number of trials per condition.
NTRIALS = 8
# Number of trials per stimulus (affective stimuli get repeated).
NTRIALSPERSTIM = 2
# Trial duration in milliseconds.
TRIAL_DURATION = 8000

# ANALYSIS SETTINGS
# Bin size for re-referencing samples to bins across the trial duration.
# This is in milliseconds.
BINWIDTH = 100.0 / 3.0

#Other factors
#GROUPS = ['fearful','nonfearful'] #between subjects
#EXPOSURE_CONDITIONS = ['FO','M','S'] #between subjects
#CELLS = ['fearfulFO','fearfulM','fearfulS','nonfearful']
#SESSIONS = ['A','B','C'] #within subjects

# FILES AND FOLDERS.
DIR = os.path.dirname(os.path.abspath(__file__))
DATAFILE = os.path.join(DIR, "data_exp_61655-v19_task-y41j_standard.csv")
TRIALDATAFOLDER = os.path.join(DIR, "uploads")
DATADIR = os.path.join(DIR, "data")
#FILEPATH = os.path.join(DATADIR, FILENAME)
DWELLPATH = os.path.join(DATADIR, "dwell_mouse_memmap.dat")
DWELLSHAPEPATH = os.path.join(DATADIR, "dwell_mouse_shape_memmap.dat")
OUTDIR = os.path.join(DIR, "output")
if not os.path.isdir(OUTDIR):
    os.mkdir(OUTDIR)

# # PLOTTING move this within loop?
# PLOTCOLMAPS = { \
#     "snake":  "Oranges", \
#     "pleasant": "Blues", \
#     }
# Font sizes for different elements.
FONTSIZE = { \
    "title":            32, \
    "axtitle":          24, \
    "legend":           14, \
    "bar":              14, \
    "label":            16, \
    "ticklabels":       12, \
    "annotation":       14, \
    }
# Set the y limits for various plots.
YLIM = {"dwell_p":  [-50, 50]}
# Set the locations of legends.
# LEGENDLOCS = { \
#     "snake":  "upper right", \
#     "pleasant": "lower right", \
#     }



    
# # # # #
# LOAD DATA

# Load the raw data.
if (not os.path.isfile(DWELLPATH)) or OVERWRITE_TMP:
    
    print("Loading raw data...")

    # Load the file.
     # Load the file.
    data = gorilla.read_file(DATAFILE, TRIALDATAFOLDER, \
    custom_fields=["condition", "left_stim", "right_stim"], \
    verbose=True)
    
    #get gendersex data and put in data dict
    df_gs = pd.read_csv("gendersex.csv")
    for ppname in data.keys():
        try:
            psex = df_gs.loc[df_gs['Participant Private ID'] == int(ppname), 'SEX'].iloc[0]
            data[ppname]["sex"] = psex
        except IndexError:
            data[ppname]["sex"] = "missing"
    #format time-lock time to start of trial (transform to range of 0-10000)
    for ppname in data.keys():
        for i in range(len(data[ppname]["trials"])):
            for ii in range(len(data[ppname]["trials"][i]["time"])):
                if ii ==0:
                    t_correct = data[ppname]["trials"][i]["time"][0]
                    data[ppname]["trials"][i]["time"][0] = 0
                else:
                    data[ppname]["trials"][i]["time"][ii] = data[ppname]["trials"][i]["time"][ii] - t_correct
        
    # Count the number of participants.
    participants = list(data.keys())
    participants.sort()
    n_participants = len(participants)
    
    print("\nLoaded {} participants".format(n_participants))
    
    # Empty matrices to start with.
    n_conditions = len(CONDITIONS)
    n_stimuli_per_condition = NTRIALS // NTRIALSPERSTIM
    n_aois = len(AOI) - 1
    n_bins = int(numpy.ceil(TRIAL_DURATION / BINWIDTH))
    shape = numpy.array([n_participants, n_conditions, \
        n_stimuli_per_condition, NTRIALSPERSTIM, n_aois, n_bins], \
        dtype=numpy.int32)
    dwell_shape = numpy.memmap(DWELLSHAPEPATH, dtype=numpy.int32, \
        mode="w+", shape=(len(shape)))
    dwell_shape[:] = shape[:]
    dwell_shape.flush()
    dwell = numpy.memmap(DWELLPATH, dtype=numpy.float32, mode="w+", \
        shape=tuple(shape))
    dwell[:] = 0.0
    
    # Compute the bin edges.
    bin_edges = numpy.arange(0, TRIAL_DURATION+BINWIDTH/2, BINWIDTH, \
        dtype=numpy.float32)
    
    # Run through all participants.
    print("\nLooping through participants...")
    for pi, participant in enumerate(participants):
        
        print("\tProcessing participant '{}' ({}/{})".format( \
            participant, pi+1, len(participants)))

        # Extract the participant's display resolution and view rect.
        dispsize = map(int, data[participant]["resolution"].split("x"))
        viewrect = map(int, data[participant]["viewport"].split("x"))
        
        # Keep track of the count for each affective stimulus.
        stim_count = {}
    
        # Run through all trials.
        for i in range(len(data[participant]["trials"])):
    
            # Start with some blank variables.
            condition = None
            stimuli = {"affective":None, "neutral":None}
            affective_stim = None
            aoi_rect = {"left_stim":None, "right_stim":None}
    
            # Reformat zone shape messages to match old format
            # Create dict with zone coordinates
            Zone_rect = {}
            coord_list = ['x','y','w','h']
            if len(data[participant]["trials"][i]["msg"]) >3: # >3 so it won't crash on empty trials without this information in msgs
                for zone in ['Zone1','Zone2']:
                    Zone_rect[zone] = []
                    if zone == 'Zone1':
                        Zonemsg = data[participant]["trials"][i]["msg"][6][1]
                    else:
                        Zonemsg = data[participant]["trials"][i]["msg"][7][1]
                    Zonemsg_split = Zonemsg.split(';')
                    for ji,coord in enumerate(coord_list):
                        if coord == 'x':
                              temp = Zonemsg_split[2]
                        elif coord == 'y':
                              temp = Zonemsg_split[3]
                        elif coord == 'w':
                              temp = Zonemsg_split[4]
                        elif coord == 'h':
                              temp = Zonemsg_split[5]
                       
                        temp_split = temp.split('=')
                        Zone_rect[zone].append(int(round(float(temp_split[1])))) 
                    
            #convert zone coordinates into image coordinates. Here, after consulting Gorilla documentation
            # I assume images are presented in center of zone, and are not stretched if zone bigger than image,
            #and if they don't fit in zone, are shrunk to width of zone w aspect ratio maintained 
            aoi_rect = {}
            data[participant]["trials"][i]["aoi_coord"] = []
            if len(data[participant]["trials"][i]["msg"]) >3:
                for zone in ['Zone1','Zone2']:
                    i_rect = []
                    z_rect = Zone_rect[zone]
                    zone_center = [(z_rect[0]+(z_rect[2]/2)),(z_rect[1]+(z_rect[3]/2))] #list with x,y for zone center
                    if z_rect[2] >= STIM_width :#presented in actual size
                        i_rect = [zone_center[0]-(STIM_width/2),(zone_center[1]-(STIM_height/2)),STIM_width,STIM_height] #x,y,w,h for actual size
                    else:#shrunk to width of zone w aspect ratio maintained
                        i_rect = [zone_center[0]-(z_rect[2]/2),zone_center[1]-((z_rect[2]*STIM_ratio)/2),z_rect[2],z_rect[2]*STIM_ratio]
                    for ir in range(len(i_rect)):
                        i_rect[ir] = int(round(i_rect[ir]))
                    if zone == "Zone1":
                        z_field = 'left_stim'
                    else:
                        z_field = 'right_stim'
                    aoi_rect[z_field] = i_rect
                    data[participant]["trials"][i][zone] = i_rect
            #this is not ideal. if trial is empty, use, aoi rect from first trial
            else: 
                
                aoi_rect['left_stim'] = data[participant]["trials"][0]['Zone1']
                aoi_rect['right_stim'] = data[participant]["trials"][0]['Zone2']
    
            # Run through all messages in this trial. 
            # Here is where we need to update for new conditions/stim names
            for j, msg in enumerate(data[participant]["trials"][i]["msg"]):
                    
                # Process the message that contains the condition info.
                if msg[0] == "condition":
                    condition = msg[1]
    
                # Process the messages related to the stimulus.
                # there is where we can address issue
                elif msg[0] in ["left_stim", "right_stim"]:
                    # Split the extension and the file name.
                    name, ext = os.path.splitext(msg[1])
                    #arbitrary assignment of stim in pair to affective/neutral
                    #for 1_ conditions, 1-digit = nude, 2-digit = clothed
                    #for 2_ conditions, 1-digit = man, 2-digit = woman
                    letter_code = name[0]
                    numeric_code = name[1:]
                    if len(numeric_code)==2:
                        stimuli["neutral"] = name
                    elif len(numeric_code)==1:
                        stimuli["affective"] = name
                        affective_stim = msg[0]
                
                # Process the AOI rect messages.
                elif msg[0] in ["zone1_shape", "zone2_shape"]:
                    # Extract the rect.
                    rect = msg[1].split("x")
                    rect = map(float, rect)
                    rect = map(round, rect)
                    rect = map(int, rect)
                    if rect[0] < viewrect[0]//2:
                        aoi_rect["left_stim"] = rect
                    else:
                        aoi_rect["right_stim"] = rect
            
            # Skip the practice trial (and trials in other conditions).
            if condition not in CONDITIONS:
                continue
            
            # Skip trials where the stimulus rects don't make sense with the
            # viewport, which should result in only one of the two AOI rects
            # being defined (as both will be lower than half the viewport).
            if (aoi_rect["left_stim"] is None) or \
                (aoi_rect["right_stim"] is None):
                print("\t\tAOI locations for participant " + \
                    "'{}' do not line up with viewrect {}".format(\
                    participant, viewrect))
                break
            
            # Add 1 to the stimulus count.
            if stimuli["affective"] not in stim_count.keys():
                stim_count[stimuli["affective"]] = 0
            else:
                stim_count[stimuli["affective"]] += 1
            #deals with counts over 4 preentations (>3)
            if stim_count[stimuli["affective"]] > 1:
                print("more than 4 presentations",participant)
                continue
            # Get the indices for the dwell matrix for the current condition
            # and stimulus.
            i_condition = CONDITIONS.index(condition)
            i_stimulus = STIMULI[condition].index(stimuli["affective"])

            # Go through all samples.
            for si in range(data[participant]["trials"][i]["time"].shape[0]):
                
                # Convenience renaming of the starting time for this sample.
                t0 = data[participant]["trials"][i]["time"][si]
                # Compute the end time of the current sample as the start
                # time of the next.
                if si < data[participant]["trials"][i]["time"].shape[0]-1:
                    t1 = data[participant]["trials"][i]["time"][si+1]
                # If this is the last sample, assume the median inter-sample
                # time.
                else:
                    t1 = t0 + numpy.nanmedian(numpy.diff( \
                        data[participant]["trials"][i]["time"]))
                # Convenience renaming of (x,y) coordinates.
                x = data[participant]["trials"][i]["x"][si]
                y = data[participant]["trials"][i]["y"][si]
                
                # Skip missing data.
                if numpy.isnan(x) or numpy.isnan(y):
                    continue
                # Skip timepoints beyond the edge.
                if t0 >= bin_edges[-1]:
                    continue

                # Check which area of interest this stimulus falls in.
                aoi = "other"
                for aoi_loc in aoi_rect.keys():
                    hor = (x > aoi_rect[aoi_loc][0]) and \
                        (x < aoi_rect[aoi_loc][0]+aoi_rect[aoi_loc][2])
                    ver = (y > aoi_rect[aoi_loc][1]) and \
                        (y < aoi_rect[aoi_loc][1]+aoi_rect[aoi_loc][3])
                    if hor and ver:
                        if aoi_loc == affective_stim:
                            aoi = "affective"
                        else:
                            aoi = "neutral"
                
                # Only record affective and neutral.
                if aoi == "other":
                    continue
                
                # Compute the first bin that this sample falls into. The last
                # bin in the range with smaller edges is the one that the
                # current sample fits within.
                si = numpy.where(bin_edges <= t0)[0][-1]
                
                # Compute the final bin that this sample falls into. The first
                # bin with a larger edge than the current sample will be the
                # bin after the current sample.
                if t1 < bin_edges[-1]:
                    ei = numpy.where(bin_edges > t1)[0][0] - 1
                else:
                    # Minus 2, as we need the last bin, which starts with the
                    # second-to-last bin edge.
                    ei = bin_edges.shape[0] - 2
                
                # If the sample falls within a bin.
                if si == ei:
                    dwell[pi, i_condition, i_stimulus, \
                        stim_count[stimuli["affective"]], AOI.index(aoi), si] \
                        += (t1 - t0) / BINWIDTH
                # If the sample falls in more than one bin.
                else:
                    # Compute the proportion of the first bin that is covered by
                    # the current sample.
                    p0 = (bin_edges[si+1]-t0) / BINWIDTH
                    # Compute the proportion of the last bin that is covered by
                    # the current sample.
                    p1 = (t1 - bin_edges[ei]) / BINWIDTH
                    
                    # Add the proportions to the dwell matrix.
                    dwell[pi, i_condition, i_stimulus, \
                        stim_count[stimuli["affective"]], AOI.index(aoi), si] \
                        += p0
                    dwell[pi, i_condition, i_stimulus, \
                        stim_count[stimuli["affective"]], AOI.index(aoi), ei] \
                        += p1
                    if ei - si > 1:
                        dwell[pi, i_condition, i_stimulus, \
                            stim_count[stimuli["affective"]], AOI.index(aoi), \
                            si+1:ei] = 1.0

# Load the data from the temporary file.
else:
    # Load the dwell data's shape from the dwell_shape file.
    dwell_shape = tuple(numpy.memmap(DWELLSHAPEPATH, dtype=numpy.int32, \
        mode="r"))
    # Load the dwell data from file.
    dwell = numpy.memmap(DWELLPATH, dtype=numpy.float32, mode="r", \
        shape=dwell_shape)
    # Recompute the bin edges.
    bin_edges = numpy.arange(0, TRIAL_DURATION+BINWIDTH/2, BINWIDTH, \
        dtype=numpy.float32)


# # # # #
# PLOT

# The dwell matrix has shape (n_participants, n_conditions, 
# n_stimuli_per_condition, n_stimulus_presentations, n_aois, n_bins)
# Plotting should happen as a time series (dwell[:,:,:,:,:,0:t). We're 
# interested in the difference between the neutral AOI  (dwell[:,:,:,:,1,:])
#  and the affective AOI (dwell[:,:,:,:,0,:]). This should be averaged across
# all stimuli (dwell[:,:,0:n,:,:,:]), and plotted separately for each stimulus
# presentation (dwell[:,:,:,0:NTRIALSPERSTIM,:,:]). The conditions 
# (dwell[:,i,:,:,:,:]) should be plotted in separate plots. The average dwell
# time difference and 95% confidence interval should be computed over all
# participants (dwell[0:n_participants,:,:,:,:,:]).

# Create a four-panel figure. (Top row for lines, bottom row for heatmaps;
# columns for different conditions.)
fig, axes = pyplot.subplots(nrows=2, ncols=len(CONDITIONS), \
    figsize=(8.0*len(CONDITIONS),9.0), dpi=300.0)
fig.subplots_adjust(left=0.05, bottom=0.06, right=0.95, top=0.95,
    wspace=0.15, hspace=0.2)

# Compute the bin centres to serve as x-values.
time_bin_centres = bin_edges[:-1] + numpy.diff(bin_edges)/2.0

# Loop through the conditions.
for ci, condition in enumerate(CONDITIONS):
    
    #This next section was previously hard-coded before loops
    cond_key1 = CONDITIONS[0]
    cond_key2 = CONDITIONS[1]
    cond_key3 = CONDITIONS[2]
    cond_key4 = CONDITIONS[3]
    
    # Add a legend.
    LEGENDLOCS = { \
    cond_key1: "upper right", \
    cond_key2: "lower right", \
    cond_key3: "upper right", \
    cond_key4: "lower right", \
    }
    
    # PLOTTING move this within loop?
    PLOTCOLMAPS = { \
    cond_key1:  "Oranges", \
    cond_key2: "Blues", \
    cond_key3:  "Greens", \
    cond_key4: "Purples", \
    }
        
    # Choose the top-row, and the column for this condition.
    ax = axes[0,ci]
    
    # Compute the difference between neutral and affective stimulus. This
    # quantifies approach (positive values) or avoidance (negative values).
    # Computed as affective minus neutral, only for the current condition.
    # d has shape (n_participants, n_stimuli, n_presentations, n_bins)
    d = dwell[:,ci,:,:,0,:] - dwell[:,ci,:,:,1,:]
    
    # Average over all different stimuli, and recode as percentages.
    # val has shape (n_participants, n_presentations, n_bins)
    val = 100 * numpy.nanmean(d, axis=1)
    
    # Compute the mean over all participants.
    m = numpy.nanmean(val, axis=0)
    # Compute the within-participant 95% confidence interval.
    nv = val - numpy.nanmean(val, axis=0) \
        + numpy.nanmean(numpy.nanmean(val, axis=0))
    sd = numpy.nanstd(nv, axis=0, ddof=1)
    sem = sd / numpy.sqrt(nv.shape[0])
    ci_95 = 1.96 * sem

    # Specify the colour map for the current condition.
    cmap = matplotlib.cm.get_cmap(PLOTCOLMAPS[condition])
    voffset = 3
    vmin = 0
    vmax = dwell.shape[3] + voffset
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    
    # Plot the separation line and background colours.
    ax.plot(time_bin_centres, numpy.zeros(time_bin_centres.shape), ':', lw=3, \
        color="black", alpha=0.5)
    ax.fill_between(time_bin_centres, \
        numpy.ones(time_bin_centres.shape)*YLIM["dwell_p"][0], \
        numpy.zeros(time_bin_centres.shape), color="black", alpha=0.05)
    annotate_x = time_bin_centres[0] + \
        0.02*(time_bin_centres[-1]-time_bin_centres[0])
    ax.annotate("Approach", (annotate_x, YLIM["dwell_p"][1]-5), \
        fontsize=FONTSIZE["annotation"])
    ax.annotate("Avoidance", (annotate_x, YLIM["dwell_p"][0]+3), \
        fontsize=FONTSIZE["annotation"])
    
    # LINES
    # Plot each stimulus presentation separately.
    for j in range(m.shape[0]):
        
        # Define the label.
        if j == 0:
            lbl = "First presentation"
        else:
            lbl = "Stimulus repeat {}".format(j)

        # PLOT THE LIIIIIIINE!
        ax.plot(time_bin_centres, m[j,:], "-", color=cmap(norm(j+voffset)), \
            lw=3, label=lbl)
        # Shade the confidence interval.
        ax.fill_between(time_bin_centres, m[j,:]-ci_95[j,:], \
            m[j,:]+ci_95[j,:], alpha=0.2, color=cmap(norm(j+voffset)))
        
    
    #legens
    if condition in LEGENDLOCS.keys():
        loc = LEGENDLOCS[condition]
    else:
        loc = "best"
    ax.legend(loc=loc, fontsize=FONTSIZE["legend"])

    # Finish the upper plot.
    ax.set_title(condition.capitalize(), fontsize=FONTSIZE["axtitle"])
    ax.set_xlabel("Time in trial (seconds)", fontsize=FONTSIZE["label"])
    ax.set_xlim([0, TRIAL_DURATION])
    ax.set_xticks(range(0, TRIAL_DURATION+1, 1000))
    ax.set_xticklabels(range(0, (TRIAL_DURATION//1000 + 1)), \
        fontsize=FONTSIZE["ticklabels"])
    if ci == 0:
        ax.set_ylabel(r"Dwell percentage difference", fontsize=FONTSIZE["label"])
    ax.set_ylim(YLIM["dwell_p"])
    ax.set_yticks(range(YLIM["dwell_p"][0], YLIM["dwell_p"][1]+1, 10))
    ax.set_yticklabels(range(YLIM["dwell_p"][0], YLIM["dwell_p"][1]+1, 10), \
        fontsize=FONTSIZE["ticklabels"])
    
    # HEATMAPS
    # Choose the heatmap row, and the column for the current condition.
    ax = axes[1,ci]

    # Create the colourmap for this condition.
    cmap = matplotlib.cm.get_cmap("coolwarm")
    vmin = YLIM["dwell_p"][0]
    vmax = YLIM["dwell_p"][1]
    vstep = YLIM["dwell_p"][1] // 2

    # Plot the heatmap.
    ax.imshow(m, cmap=cmap, vmin=vmin, vmax=vmax, interpolation="none", \
        aspect="auto", origin="upper")
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    divider = make_axes_locatable(ax)
    bax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = matplotlib.colorbar.ColorbarBase(bax, cmap=cmap, norm=norm, \
        ticks=range(vmin, vmax+1, vstep), orientation='vertical')
    if ci == axes.shape[1]-1:
        cbar.set_label(r"$\Delta$dwell percentage", \
            fontsize=FONTSIZE["bar"])

    #ax.set_title(CONDITIONS[condition], fontsize=FONTSIZE["axtitle"])
    ax.set_xlabel("Time in trial (seconds)", fontsize=FONTSIZE["label"])
    ax.set_xlim([-0.5, m.shape[1]-0.5])
    ax.set_xticks(numpy.linspace(-0.5, m.shape[1]+0.5, \
        num=TRIAL_DURATION//1000 + 1))
    ax.set_xticklabels(range(0, (TRIAL_DURATION//1000 + 1)), \
        fontsize=FONTSIZE["ticklabels"])
    if ci == 0:
        ax.set_ylabel("Presentation number", fontsize=FONTSIZE["label"])
    ax.set_ylim([NTRIALSPERSTIM-0.5, -0.5])
    ax.set_yticks(range(0, NTRIALSPERSTIM))
    ax.set_yticklabels(range(1, NTRIALSPERSTIM+1), \
        fontsize=FONTSIZE["ticklabels"])

# Save the figure.
fig.savefig(os.path.join(OUTDIR, "dwell_percentages.png"))
pyplot.close(fig)