# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 16:36 2018

@author: ejb

This program is intended to perform uncertainty analysis on externally 
post-processed results from SAPHIRE. 
It reads in exported cut sets from all seismic initiating events in the Vogtle 
PRA model, as well as an exported list of random samples of basic event 
probabilities, and quantifies end state frequencies for each sample in a 
way that is consistent with the original point estimate quantification.
It then generates an uncertainty distribution for each end state.

The quantification approach here is to identify basic events that are shared 
between multiple cut sets with the same end state, quantify them without the 
shared event, and then multiply the result by the shared event probability. 
SAPHIRE, on the other hand, would include the shared event in each cut set's 
probability before combining them with minimal cut set upper bound (MCUB).

MonteCarlo5 uses a set of 4575 random samples (out of 15000 requested, since 
SAPHIRE ended up with 10001 samples above 1.0 or below 0.0; not sure what 
happened to the other 424) that each include a sampled value for every basic 
event used in ANY end state with seismic cut sets.

"""

# Probability mismatch in cut set  21080 SAPHIRE: 5.208333333333334e-07 My calc: 0.00039965331559546006
# Bin 6 sample 2008 LCF

import os
import sys
import copy
import pickle

 
def calc_probs_dependent2(cut_sets, cutsetprobs, cutseteventsincl, basic_events, cut_set_indices):
    """Calculates combined probability for a list of cut sets.
 
    Uses minimal cut set upper bound calculation to combine probabilities
    from individual cut sets into a total probability.
    Unlike SAPHIRE, it takes into account dependence of cut sets, by
    removing common events and applying their probability to the group
    of all cut sets that contain them, instead of individually.
    This process occurs recursively, starting with the largest groups of cut 
    sets sharing a basic event and then operating on subgroups that share 
    another event. When the remaining group of cut sets are independent, 
    the MCUB calculation is performed normally (like SAPHIRE does).
    (Note that it also modifies the cut set data structure in the process.)
    
    I think this is O(n^2), but running time has not been a problem.
 
    Args:
        cut_sets: A list of cut set tuples (probability,
        list of basic event strings, sequence string, end state string)
        
        cutsetprobs: A corresponding list of just the probabilities
        
        cutseteventsincl: A list of lists of binary values, one list of binary
        values for each cut set in cut_sets. Each binary value corresponds to 
        a basic event in the cut set. Initially they are all True, but will 
        be changed to False instead of removing the event from the cut set list
     
        basic_events: A dictionary of basic events, with string as key and
        containing probability
     
        cut_set_indices: A list of indices in cut_sets that should be
        included in this probability calculation
     
    Returns:
        The combined probability for the given cut sets.
 
    """
 
    if not cut_set_indices: # func has been called on a group with no cut sets
        return 0.0
    
    # First go through and count how many times each basic event appears.
    event_count = {}
    event_prob = {}
    for cut_set_index in cut_set_indices:
        for i,event in enumerate(cut_sets[cut_set_index][1]):
            if not cutseteventsincl[cut_set_index][i]:
                continue
            if event[0] == "/":
                ev_prb = 1.0 - basic_events[event[1:]]
            else:
                ev_prb = basic_events[event]
            if event in event_count and ev_prb > 0.00:
                event_count[event] += 1
            else:
                event_count[event] = 1
                event_prob[event] = 0
            event_prob[event] += cut_sets[cut_set_index][0]
   
    # Find how often the most common event appears in this group of cut sets    
    if len(event_count) == 0:
        max_count = 0
    else:
        max_count = max(event_count.values())
 
    # If there's an event that appears in more than one cut set,
    # we need to treat those cut sets as dependent. 
    if max_count > 1:
        # Recursive case:
        # Make a list of only the cut sets that contain common_event.
        # Remove that event from each of those cut sets, and modify
        # their probabilities accordingly.
        # Run this function on the new list.
        # Return that result times the probability of the common event.

        # Find the most important common event (the one with the most 
        # probability in cut sets that include it).
        sorted_probs = [(x, event_prob[x]) for x in sorted(event_prob, 
                        key=event_prob.get, reverse=True)]
        max_prob_count = 0
        ind = -1
        while max_prob_count < 2:
            ind += 1
            max_prob_count = event_count[sorted_probs[ind][0]]
        most_important = sorted_probs[ind][0]
     
        # Alternative: Find the most common event.
#        most_common = [k for k, v in event_count.items() if v == max_count]
        # If two or more are equal, we'll just use the first one.
        # The other(s) will be dealt with on deeper recursions.
        # It would be faster to do them together but speed is not an issue.
#        common_event = most_common[0]

        common_event = most_important 
        
        # Find its probability.
        if common_event[0] == "/": # it's a success event
            # so set it to 1 - probability of failure
            common_prob = 1.0 - basic_events[common_event[1:]]
        else:
            common_prob = basic_events[common_event]
        
        # Make a list of which cut sets have that event, and modify
        # those cut sets' probabilities. And remove the event from them.
        # Changed to just mark it as inactive.
        common_cs = []
        other_cs = [] # list of which cut sets don't contain the common event
        for cut_set_index in cut_set_indices:
            if common_event in cut_sets[cut_set_index][1]:
#                if cut_set_index == 0:
#                    print(cutseteventsincl[cut_set_index], "old prob", cutsetprobs[cut_set_index])
#                    print("Removing",common_event, "with prob",common_prob)

                common_cs.append(cut_set_index)
                # Change its probability
                if common_prob != 0.0:
#                    cs = list(cut_sets[cut_set_index])
#                    cs[0] /= common_prob
#                    if cs[0] > 1.0: cs[0] = 1.0
                    cutsetprobs[cut_set_index] /= common_prob
                    if cutsetprobs[cut_set_index] > 1.0:
                        if cutsetprobs[cut_set_index] < 1.01:
                            cutsetprobs[cut_set_index] = 1.0
                        else:
                            print("Messed up event, cut set", cut_set_index, common_event, cutsetprobs[cut_set_index], common_prob, cut_sets[cut_set_index])
                            sys.exit(0)
                    # And remove the common event
#                    cs[1].remove(common_event)
#                    cut_sets[cut_set_index] = tuple(cs)
                    # (tuples are immutable so we create a replacement one)
                    i = cut_sets[cut_set_index][1].index(common_event)
                    cutseteventsincl[cut_set_index][i] = False
            else:
                other_cs.append(cut_set_index)
        
        # Find probability for the group containing the common event
        if common_prob == 0.0:
            batch_prob = 0.0
        else:
            batch_prob = calc_probs_dependent2(cut_sets, cutsetprobs, 
                                              cutseteventsincl, basic_events, 
                                              common_cs) * common_prob
        # Find probability for the group not containing the common event
        other_prob = calc_probs_dependent2(cut_sets, cutsetprobs, cutseteventsincl, basic_events, other_cs)
        
        # Return the MCUB combination of those two groups' probabilities
        return 1.0 - (1.0 - batch_prob) * (1.0 - other_prob)
       
     
    else:
        # Base case: Calculate MCUB in the normal way
        probcomp = 1.0 # complement of combined prob of cut sets
        for cut_set_index in cut_set_indices:
            probcomp *= 1.0 - cutsetprobs[cut_set_index]
        return 1.0 - probcomp
    
# end of calc_probs_dependent2




def cutsetcopy(cut_sets,indices):
    new_cut_sets = []
#    print(cut_sets[0])
    for i in indices:
        new_cut_sets.append(copy.deepcopy(cut_sets[i]))
#    print(new_cut_sets[0])
    return new_cut_sets




if __name__ == "__main__":
    """
    Note for reading this code:
        Most of the actual calculation happens in the for loop at the end,
        where calc_probs_dependent2 is called.
        Everything up to there is just setup (reading in files and creating 
        data structures for the basic events, cut sets, sequences, end states) 
        plus a little error checking.
    """
    
    start_sample = int(sys.argv[1])
    end_sample = int(sys.argv[2])
    
    rel_cats = ["BMT","CIF","CIF-SC","ECF","ICF-BURN","ICF-BURN-SC",
            "ISGTR","LCF","LCF-SC","NOCF","SGTR-C","SGTR-O",
            "SGTR-O-SC","V","V-F","V-F-SC"]
    bin_cdfs = [1.18e-5, 1.30e-6, 1.21e-6, 1.85e-6, 2.79e-6, 2.53e-6, 1.87e-6, 2.47e-7, 2.32e-9]
    bin_inits = ["1-IE-EQK-BIN-" + str(i) for i in range(0,9)]
    # Make dictionary of basic events with probabilities (and calculation type)
    # This info comes from MARD output, will be replaced later.
    with open("VOGT.BEI",'r') as infile:
        lines = infile.readlines()

    basic_events = {}
    for line in lines:
        pieces = line.split(",")
        if len(pieces) < 12 or line.startswith("*"): # skip irrelevant lines
            continue
        name = pieces[0].strip()
        FdT = pieces[1].strip() # I don't currently use this for anything.
        Prob = float(pieces[5].strip())
        CalcProb = float(pieces[12].strip())
        basic_events[name] = CalcProb
        
    basic_events_original = copy.deepcopy(basic_events)

    # list all the cut set files exported from SAPHIRE
    cdf_cutsetfilename = "1-CD-EQ_cut_sets.txt" # this one has seismic CDF cut sets for all bins together
    ft_cs_files = os.listdir("FTcutsets") # now all the cut set files for the bridge tree fault trees
    bin_cut_set_filenames = ["CutSetReport-1-EQK-BIN-" + str(bn) + ".txt" for bn in range(1,8)] # now the L2 cut sets for each bin
    
    
    ft_cut_sets = {}
    for ft_cs_filename in ft_cs_files:
        with open("FTcutsets/" + ft_cs_filename, "r") as cutsetfile:
            cutsetlines = cutsetfile.readlines()
        ft = ft_cs_filename[0:-4]
        ft_cut_sets[ft] = []
        for line in cutsetlines:
            pieces = line.split(",")
#            if not line[0] in ["1234567890"]: continue
            if not line[0] in "1234567890": continue
            prob = float(pieces[2])
            events = pieces[4:-1] + [pieces[-1].strip()]
            ft_cut_sets[ft].append((prob,events,"",ft))
#        print(ft, ft_cut_sets[ft])
#    print(ft_cut_sets)
    calc_ft_probs = {} # cut set probabilities that I will calculate here
    for ft in ft_cut_sets.keys():
        calc_ft_probs[ft] = []
        for cut_set in ft_cut_sets[ft]:
            prb = 1.0
            try:
                for event in cut_set[1]:
                    if event.startswith("/"):
                        prb *= (1 - basic_events[event[1:]])
                    else:
                        prb *= basic_events[event]
            except:
                print("Error in accessing basic event!")
                pass
            calc_ft_probs[ft].append(prb)
            # Check for match with SAPHIRE.
            # There's often a small difference, so report if > 1% and > 1e-4
            if abs(prb - cut_set[0]) > 1.e-4 and (prb == 0 or 
                  abs(prb - cut_set[0]) / prb > 1.e-2):
                print("Probability mismatch in cut set ", ft, len(calc_ft_probs[ft])-1, 
                      "SAPHIRE:",cut_set[0],"My calc:",prb)
#    sys.exit(0)
    cdf_cut_sets = {}
    with open(cdf_cutsetfilename, "r") as cutsetfile:
        cutsetlines = cutsetfile.readlines()
    for line in cutsetlines:
        pieces = line.split(",")
        if len(pieces) < 6: continue
        freq = float(pieces[2])
        initiator = pieces[4]
        init_freq = basic_events[initiator]
        prob = freq / init_freq
        events = pieces[5:-1]
        if initiator in cdf_cut_sets:
            index = len(cdf_cut_sets[initiator])
        else:
            index = 0
            cdf_cut_sets[initiator] = []
        cdf_cut_sets[initiator].append((prob,events,"","1-CD-EQ"))
    # Calculate cut set probabilities (imitating SAPHIRE)
    calc_cdf_probs = {}
    for key in cdf_cut_sets.keys():
        calc_cdf_probs[key] = [] # cut set probabilities that I will calculate here
    for initiator in cdf_cut_sets.keys():
        ccdp_comp = 1.0
        for cut_set in cdf_cut_sets[initiator]:
            prb = 1.0
            for event in cut_set[1]:
                if event.startswith("/"):
                    prb *= (1 - basic_events[event[1:]])
                else:
                    prb *= basic_events[event]
            calc_cdf_probs[initiator].append(prb)
            # Check for match with SAPHIRE.
            # There's often a small difference, so report if > 1% and > 1e-4
            if abs(prb - cut_set[0]) > 1.e-5 and (prb == 0 or 
                  abs(prb - cut_set[0]) / prb > 5.e-3):
                print("Probability mismatch in cut set ", initiator, len(calc_cdf_probs[initiator])-1, 
                      "SAPHIRE:",cut_set[0],"My calc:",prb)
            # update the CCDP
            ccdp_comp *= (1.0 - prb)
        binfreq = (1.0 - ccdp_comp) * basic_events[initiator]
        saphbinfreq = bin_cdfs[bin_inits.index(initiator)]
        if abs(binfreq - saphbinfreq) > 1.e-10 and (binfreq == 0 or abs(binfreq - saphbinfreq) / binfreq > 1.e-2):
            print("Probability mismatch in CDF for initiator", initiator, binfreq, saphbinfreq)
        else:
            print(initiator, binfreq, saphbinfreq, bin_inits.index(initiator))
    
    L2_cut_sets = [[]]
    L2_end_states = [{}] # each bin will have a dictionary of end states
    L2_calc_probs = [[]]
    # Read in each bin's L2 cut sets and store the list of them
    for cut_set_filename in bin_cut_set_filenames:
        with open(cut_set_filename, "r") as cutsetfile:
            cutsetlines = cutsetfile.readlines()
        cut_sets = [] # use a list here because we will reference cut sets by index
        end_states = {}

        # find initiator frequency for the first cut set
        # we'll check later to make sure they all have the same initiator
        initiator = cutsetlines[2].split(",")[4]
        init_freq = basic_events[initiator]
        binnum = int(cut_set_filename[-5])
#        if binnum > 1:
#            print("\n")
        print("Bin",binnum,"-",init_freq,"/ rcy - cut sets processed.", end=" ")
    
        # Make a list of each cut set's events, its probability and its end state.
        # Also for each sequence and end state, a list of indices of its cut sets.
        for line in cutsetlines:
            index = len(cut_sets)
            pieces = line.split(",")
            if len(pieces) < 6:
                continue  # I observe that any line with fewer than 5 commas is not
                            # a cut set
            
            # Read in cut set frequency and events list, calculate probability
            freq = float(pieces[2])
            prob = freq / init_freq
            events = pieces[5:-1]
            
            # Do some error checking
            if not pieces[4] == initiator:
                print("Error! File includes cut sets from multiple initiators:", 
                      initiator, "and", pieces[4])
    
            # Read in the name of the end state this cut set is associated with
            end_state = pieces[-1].split(">")[-1].strip() # get end state from file
            # (Show Origin must be checked when the cut sets are exported.)
            
            # Add this cut set's index to the list for its end state
            if end_state in end_states:
                end_states[end_state].append(index)
            else:
                end_states[end_state] = [index]
            # Add this cut set to the list of cut sets
            cut_sets.append((prob,events,"",end_state))
        
        L2_end_states.append(end_states)
        L2_cut_sets.append(cut_sets)
    
        # Calculate cut set probabilities (imitating SAPHIRE)
        # The main reason for this is to check that none of the basic event 
        # probabilities have been modified by flag sets; if they have, SAPHIRE's 
        # cut set probability will come out different than what we calculate here.
        calc_probs = [] # cut set probabilities that I will calculate here
        for cut_set in cut_sets:
            prb = 1.0
            for event in cut_set[1]:
                if event.startswith("/"):
                    prb *= (1 - basic_events[event[1:]])
                else:
                    prb *= basic_events[event]
            calc_probs.append(prb)
            # Check for match with SAPHIRE.
            # There's often a small difference, so report if > 1% and > 1e-4
            if abs(prb - cut_set[0]) > 1.e-4 and (prb == 0 or 
                  abs(prb - cut_set[0]) / prb > 1.e-2):
                print("Probability mismatch in cut set ", len(calc_probs)-1, 
                      "SAPHIRE:",cut_set[0],"My calc:",prb)
        L2_calc_probs.append(calc_probs)

    # Now we'll do the quantifications for the Monte Carlo samples
    max_samples = 4575
#    outfilename = "out.txt"
    samplefilename = "seismic_mc15000_all.csv"
    # highest level, it's a list of 5000 samples
    # each sample has a list of 7 bins, plus bin 0 is the total i guess
    # each bin has 16 RCs
    es_freqs = [[{},{},{},{},{},{},{},{}] for i in range(max_samples + 1)]
    # plus each sample has 6 1-BE fault tree failure events
    ft_probs = [{} for i in range(max_samples + 1)]
    # plus each sample has a CDF for each bin
    cd_freqs = [[0]*9 for i in range(max_samples + 1)]
    # Normalized end state frequencies for each end state in each sample
    normalized_es_freqs = [{} for i in range(max_samples + 1)]

    # Contains the end state frequencies SAPHIRE calculated
    # for each random sampling of possible basic event probabilities.
    saphire_sample_es_freqs = [{}]*(max_samples + 1) # a dict of end state frequencies for each sample
    point_es_freqs = [{}]*8 # point estimate bin frequencies

    print("\nCalculating point estimate end state frequencies for each bin...")
    for binnum in range(1,8):
        # Calculate a probability and frequency for each end state, taking into
        # account that the cut sets are not independent (SAPHIRE doesn't do this).
        end_state_probs = {}
        escutsetprobs = {}
        escutsetevents = {}
        end_states = L2_end_states[binnum]
        for end_state in end_states.items():
            escutsets = cutsetcopy(L2_cut_sets[binnum],end_state[1])
            escutsetprobs[end_state[0]] = [cut_set[0] for cut_set in escutsets]
            escutsetevents[end_state[0]] = [[True]*len(cut_set[1]) for cut_set in escutsets]
            esprob = calc_probs_dependent2(escutsets,
                                          copy.deepcopy(escutsetprobs[end_state[0]]),
                                          copy.deepcopy(escutsetevents[end_state[0]]),
                                          basic_events, 
                                          range(len(end_state[1])))
            point_es_freqs[binnum][end_state[0]] = esprob
            point_es_freqs[0][end_state[0]] += esprob*init_freq
        print("Bin",binnum,"point estimates calculated.")
    
    print("Processing Monte Carlo sample input.")
    # Now go through the sample file and run quantification on each end state's 
    # cut sets many times w/ the random basic event probabilities from SAPHIRE.
    with open(samplefilename,'r') as infile:
        lines = infile.readlines()
    # Update dictionary of basic events with probabilities from SAPHIRE uncertainty CSV output
    sampledata = False # indicates whether we are current in a sample data block
    samplenum = 0
    min_sample = 1 # don't change this except for debugging a specific sample
                   # it makes the printing at the end fail
#    end_state = "" # current end state we are running samples for
#    end_state_present = True # is the current end state present in the cut sets?
    used_basic_events = [] # basic events used in current end state's cut sets
    end_states_list = ["GROUP","1-REL-NOCF","1-REL-ECF","1-REL-ICF-BURN","1-REL-ICF-BURN-SC","1-REL-LCF","1-REL-LCF-SC","1-REL-BMT","1-REL-CIF","1-REL-CIF-SC","1-REL-SGTR-O","1-REL-SGTR-O-SC","1-REL-ISGTR"]
    for line in lines:
        if line.startswith("Uncertainty Method:MC"): # appears before each new end state
            sampledata = False # so the previous sample block is done
            samplenum = 0
            continue
        pieces = line.split(",")
        if line.startswith("1-REL-"): # SAPHIRE point estimates for each release category
            saphire_sample_es_freqs[0][pieces[0]] = float(pieces[1])
        if line.startswith(", GROUP,1-REL-"):
            # header line for the sample data
            # lists all the events used in the seismic L2 model
            used_basic_events = pieces[3:] # BEs used for all these end states
        if line.startswith("SAMPLEDATA"):
            sampledata = True # this line marks the start of the samples
#            print("\n")
            continue
        if sampledata:
            samplenum = int(pieces[0].strip())
            if samplenum == start_sample: print("Starting Monte Carlo sample calculations.")
        if sampledata and samplenum <= end_sample and samplenum >= start_sample:
            print("Sample",samplenum,end=" ")
            for i in range(len(end_states_list)):
                saphire_sample_es_freqs[samplenum][end_states_list[i]] = float(pieces[i+1])
                # SAPHIRE's quantification for this sample, which we will
                # improve on in our calculation below.
            # Modify basic event probabilities to the current sample values.
            for i in range(len(used_basic_events)):
                if len(used_basic_events[i])>0:
                    be = used_basic_events[i]
                    basic_events[be] = float(pieces[i+3])
                else:
                    break
            # Calculate probabilities for the 1-BE- basic events (fault trees)
            for ft in ft_cut_sets.keys():
                ft_comp = 1.0
                for cut_set in ft_cut_sets[ft]:
                    ft_comp *= (1.0 - prb)
                ft_prob = 1.0 - ft_comp
                ft_BE = "1-BE-" + ft[5:]
                if ft_BE in basic_events:
                    basic_events[ft_BE] = ft_prob
                else:
                    print("Problem! Basic event name for fault tree is wrong.")
                ft_probs[samplenum][ft_BE] = ft_prob
            print("FT probs done.",end=" ")
            
            norm_factors = [0]*8
            for binnum in range(1,8):
                print("Bin",binnum,end=": ")
                init_event = "1-IE-EQK-BIN-" + str(binnum)
                # Set the initiating event frequency for this sample
                init_freq = basic_events[init_event]
                # Calculate the CDF for just this bin, using this bin's CDF cut sets,
                # and with the basic events modified to this sample's values.
                bin_cdf_cut_sets = cutsetcopy(cdf_cut_sets[init_event], range(len(cdf_cut_sets[init_event])))
                cdfcsprobs = []
                ccdp_comp = 1.0
                for cut_set in bin_cdf_cut_sets:
                    prb = 1.0
                    for event in cut_set[1]:
                        if event.startswith("/"):
                            prb *= (1 - basic_events[event[1:]])
                        else:
                            prb *= basic_events[event]
                    cdfcsprobs.append(prb)
                    # update the CCDP
                    ccdp_comp *= (1.0 - prb)
                ccdp = 1.0 - ccdp_comp
                # The dependent version, which I won't use, in order to be
                # compatible with SAPHIRE
#                cutseteventsincluded = [[True]*len(cut_set[1]) for cut_set in bin_cdf_cut_sets]
#                ccdp = calc_probs_dependent2(bin_cdf_cut_sets,
#                                             cdfcsprobs,
#                                             cutseteventsincluded,
#                                             basic_events,
#                                             range(len(bin_cdf_cut_sets)))
                cd_freqs[samplenum][binnum] = ccdp * init_freq
                print("CDF,",end=" ")
                for end_state in L2_end_states[binnum].keys():
                    escutsets = cutsetcopy(L2_cut_sets[binnum], L2_end_states[binnum][end_state])
                    # Set the cut set probabilities for this sample
                    cutsetprobs = [] # cut set probabilities that I will calculate here
                    for cut_set in escutsets:
                        prb = 1.0
                        for event in cut_set[1]:
                            if event.startswith("/"):
                                prb *= (1 - basic_events[event[1:]])
                            else:
                                prb *= basic_events[event]
                        cutsetprobs.append(prb)
                    # No point checking for match with SAPHIRE, it won't match
                    cutseteventsincluded = [[True]*len(cut_set[1]) for cut_set in escutsets] 
                    # This line does the meat of the calculations:
                    esprob = calc_probs_dependent2(escutsets,
                                                  cutsetprobs,
                                                  cutseteventsincluded,
                                                  basic_events, 
                                                  range(len(L2_end_states[binnum][end_state])))
                    # Store bin-specific end state frequency for this sample
                    es_freqs[samplenum][binnum][end_state] = esprob*init_freq
                    # Add to total of all bins for this sample of this end state
                    if end_state in es_freqs[samplenum][0]:
                        es_freqs[samplenum][0][end_state] += esprob*init_freq
                    else:
                        es_freqs[samplenum][0][end_state] = esprob*init_freq
                print("RCs.",end=" ")
                # Total all the end state frequencies for this bin in this sample:
                es_freqs[samplenum][binnum]["All end states"] = sum(es_freqs[samplenum][binnum].values())
                # Comparing that value to the CDF we calculated for it gives the normalization factor
                norm_factors[binnum] = cd_freqs[samplenum][binnum] / es_freqs[samplenum][binnum]["All end states"]
            # normalize results (sample frequencies) by bin CDF
            for binnum in range(1,8):
                for end_state in es_freqs[samplenum][binnum].keys():
                    if end_state == "All end states": continue
                    if end_state in normalized_es_freqs[samplenum]:
                        normalized_es_freqs[samplenum][end_state] += es_freqs[samplenum][binnum][end_state] * norm_factors[binnum]
                    else:
                        normalized_es_freqs[samplenum][end_state] = es_freqs[samplenum][binnum][end_state] * norm_factors[binnum]
            print("Normalization done.")

#    print("\n")
#    for bin in range(1,8):
#        bin_total = [0]*max_samples
#        norm_factor = [0]*max_samples
#        for i in range(max_samples):
#            for es in my_sample_es_freqs[bin].items():
#                if len(es[1])>i:
#                    bin_total[i] += es[1][i]
#            norm_factor[i] = bin_cdfs[bin] / bin_total[i]
#            for es in my_sample_es_freqs[bin].items():
#                if es[0] not in normalized_es_freqs:
#                    normalized_es_freqs[es[0]] = [0]*max_samples
#                if my_sample_es_freqs[bin][es[0]]:
#                    if len(my_sample_es_freqs[bin][es[0]]) > i:
#                        normalized_es_freqs[es[0]][i] += norm_factor[i]*my_sample_es_freqs[bin][es[0]][i]
#                    else:
#                        normalized_es_freqs[es[0]][i] += 0
#                        print("Problem ",bin,es[0],i)
#        print("Bin",bin,"sample frequencies:",bin_total,"average:",sum(bin_total)/len(bin_total))
#    

        
#    # Print results
#    tot_freq[0] = sum(tot_freq[1:])
#    print("Point estimates - Total:",tot_freq[0])
#    for bin in range(1,8):
#        print("Bin",bin,":",tot_freq[bin])
#    for end_state in my_sample_es_freqs[0].items():
#        print(end_state[0])
#        for i in range(max_samples):
#            if i == len(saphire_sample_es_freqs[end_state[0]]):
#                break
#            saph = saphire_sample_es_freqs[end_state[0]][i]
#            me = my_sample_es_freqs[0][end_state[0]][i]
#            chg = 100*(me/saph-1.0)
#            if chg > 5 or chg < - 50:
#                print(saph, me, "     ", 100*(me/saph-1.0),"% change due to requant")
    
    # Pickle results for later use
    file_end = str(start_sample) + "-" + str(end_sample)
    pickle.dump(point_es_freqs, open("point_es_freqs" + file_end + ".p","wb"))
    pickle.dump(ft_probs, open("ft_probs" + file_end + ".p","wb"))
    pickle.dump(cd_freqs, open("cd_freqs" + file_end + ".p","wb"))
    pickle.dump(es_freqs, open("es_freqs" + file_end + ".p","wb"))
    pickle.dump(normalized_es_freqs, open("normalized_es_freqs" + file_end + ".p","wb"))
    print("Wrote pickle files.")
    
#    # write results to file
#    with open(outfilename,"a") as outfile:
#        outfile.write(samplefilename)
#        outfile.write("Point estimate end state frequencies,")
#        for end_state in my_sample_es_freqs[0].items():
#            outfile.write("\n" + end_state[0] + "\n")
#            for i in range(max_samples):
#                if i == len(saphire_sample_es_freqs[end_state[0]]):
#                    break
#                saph = saphire_sample_es_freqs[end_state[0]][i]
#                me = my_sample_es_freqs[0][end_state[0]][i]
#                try:
#                    norm = normalized_es_freqs[es[0]][i]
#                except:
#                    norm = 0
#                outfile.write(str(i+1) + " " + str(saph) + " " + str(me) + " " + str(norm) + "\n")
#                for bin in range(1,8):
#                    if my_sample_es_freqs[bin][end_state[0]]:
#                        outfile.write(str(my_sample_es_freqs[bin][end_state[0]][i]) + " ")
#                    else:
#                        outfile.write("0.00 ")
#                outfile.write("\n")
    
# =============================================================================
#     for end_state in sorted(end_state_probs, key=end_state_probs.__getitem__, 
#                             reverse=True):
#         print(end_state, "Frequency:", 
#               "{0:.4e}".format(end_state_probs[end_state] * init_freq))
#     
#     if outfilename:
#         with open(outfilename,"a") as outfile:
#             outfile.write("\n" + samplefilename)
#             outfile.write("\nEnd State Frequencies\n")
#             es_order = ["BMT","CIF","CIF-SC","ECF","ICF-BURN","ICF-BURN-SC",
#                         "ISGTR","LCF","LCF-SC","NOCF","SGTR-C","SGTR-O",
#                         "SGTR-O-SC","V","V-F","V-F-SC"]
#             for es in es_order:
#                 esname = "1-REL-" + es
#                 if esname in end_state_probs:
#                     outfile.write(esname + ", "+"{0:.4e}".format(end_state_probs[esname] * init_freq) + "\n")
#                 else:
#                     outfile.write(esname + ", "+"0.00" + "\n")
#     


