# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 14:21:57 2018
Updated Mon Jan 29 2018

@author: ejb

This program is intended to externally post-process results from SAPHIRE.
It reads in exported cut sets from one initiating event in the Vogtle PRA 
model, as well as an exported list of basic event probabilities, and 
quantifies end state frequencies in a way that results in less frequency 
inflation than the default SAPHIRE method.

The quantification approach here is to identify basic events that are shared 
between multiple cut sets with the same end state, quantify them without the 
shared event, and then multiply the result by the shared event probability. 
SAPHIRE, on the other hand, would include the shared event in each cut set's 
probability before combining them with minimal cut set upper bound (MCUB).

"""
import sys
import copy

 
def calc_probs_dependent(cut_sets, basic_events, cut_set_indices):
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
     
        basic_events: A dictionary of basic events, with string as key and
        containing a tuple of probability and FdT
     
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
        for event in cut_sets[cut_set_index][1]:
            if event.startswith("/"):
                ev_prb = 1.0 - basic_events[event[1:]][0]
            else:
                ev_prb = basic_events[event][0]
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
        if common_event.startswith("/"): # it's a success event
            # so set it to 1 - probability of failure
            common_prob = 1.0 - basic_events[common_event[1:]][0]
        else:
            common_prob = basic_events[common_event][0]
        
        # Make a list of which cut sets have that event, and modify
        # those cut sets' probabilities. And remove the event from them.
        common_cs = []
        other_cs = [] # list of which cut sets don't contain the common event
        for cut_set_index in cut_set_indices:
            if common_event in cut_sets[cut_set_index][1]:
                common_cs.append(cut_set_index)
                # Change its probability
                if common_prob != 0.0:
                    cs = list(cut_sets[cut_set_index])
                    cs[0] /= common_prob
                    if cs[0] > 1.0: cs[0] = 1.0
                    # And remove the common event
                    cs[1].remove(common_event)
                    cut_sets[cut_set_index] = tuple(cs)
                    # (tuples are immutable so we create a replacement one)
            else:
                other_cs.append(cut_set_index)
        
        # Find probability for the group containing the common event
        if common_prob == 0.0:
            batch_prob = 0.0
        else:
            batch_prob = calc_probs_dependent(cut_sets, basic_events, 
                                              common_cs) * common_prob
        # Find probability for the group not containing the common event
        other_prob = calc_probs_dependent(cut_sets, basic_events, other_cs)
        
        # Return the MCUB combination of those two groups' probabilities
        return 1.0 - (1.0 - batch_prob) * (1.0 - other_prob)
       
     
    else:
        # Base case: Calculate MCUB in the normal way
        probcomp = 1.0 # complement of combined prob of cut sets
        for cut_set_index in cut_set_indices:
            probcomp *= 1.0 - cut_sets[cut_set_index][0]
        return 1.0 - probcomp
    
# end of calc_probs_dependent







if __name__ == "__main__":
    """
    Note for reading this code:
        All the actual calculation happens in the for loop at the end,
        where calc_probs_dependent is called.
        Everything up to there is just setup (reading in files and creating 
        data structures for the basic events, cut sets, sequences, end states) 
        plus a little error checking.
    """
    
    if len(sys.argv) < 4:
        print("Usage: requantify_from_solve.py \"cut_set_file.txt\" " +
              " \"output_file.txt\" importances=True")
        print("If output_file is \"None\", Fussell-Vessely and RAW importances will " +
              "not be calculated.")
        exit(0)
    
    # Choose data source files:

    filename = sys.argv[1] # Cut sets exported from SAPHIRE
    
    
    # Choose whether to calculate Fussell-Vessely and RAW importances 
    if sys.argv[3] == "importances=True" or sys.argv[3] == "importances=true":
        calcfv = True
    else:
        calcfv = False
    outfilename = sys.argv[2]

    if outfilename == "None" or outfilename == "\"None\"":
        outfilename = False
        calcfv = False
 
    with open("VOGT.BEI",'r') as infile:
        lines = infile.readlines()
 
    # Make dictionary of basic events with probabilities (and calculation type)
    # This info comes from MARD output
    basic_events = {}
    for line in lines:
        pieces = line.split(",")
        if len(pieces) < 12 or line.startswith("*"): # skip irrelevant lines
            continue
        name = pieces[0].strip()
        FdT = pieces[1].strip() # I don't currently use this for anything.
        Prob = float(pieces[5].strip())
        CalcProb = float(pieces[12].strip())
        basic_events[name] = (CalcProb, FdT)
 
    cut_sets = [] # use a list here because we will reference cut sets by index
    sequences = {} # sequences and end states use dictionaries for fast access
    end_states = {}
 
#    # Make a mapping of sequences to end states:
#    with open(mapfilename, 'r') as mapfile:
#        maplines = mapfile.readlines()
#    seqmap = {}
#    linenum = 0
#    for line in maplines:
#        linenum += 1
#        if linenum < 7 or len(line) < 5:
#            continue
#        pieces = line.split(",")
#        seq = pieces[0]
#        es = pieces[1].strip()
#        try:
#            es = float(es) # if this works, it was the point estimate
#            es = pieces[2].strip() # so use the next value from the line
#            if len(es) == 0: es = pieces[3].strip()
#        except:
#            pass
#        seqmap[seq] = es
    
    with open(filename,'r') as cutsetfile:
        cutsetlines = cutsetfile.readlines()
 
    # find initiator frequency for the first cut set
    # we'll check later to make sure they all have the same initiator
    initiator = cutsetlines[2].split(",")[4]
    init_freq = basic_events[initiator][0]
    
    numbad = 0 # Number of inconsistent cut sets found
    freqbad = 0
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
        # Correct double-counted success of containment failure
        # This is introduced by PP rules in SAPHIRE. The effect is small.
        # Not needed with latest version of Vogtle model but I'm leaving it in.
        if "/1-STRC-CD-EQ6" in events and "/1-RPS-SYS-EQ6-CONT" in events:
            events.remove("/1-RPS-SYS-EQ6-CONT")
        # Look for bad cut sets
        if "1-RPS-SYS-EQ6-CONT" in events and "/1-STRC-CD-EQ6" in events:
            print("Error! Inconsistent cut set", pieces[0])
            if int(pieces[0]) < 5000: print(events)
            numbad += 1
            freqbad += freq
        
#        # Read in the name of the sequence this cut set is associated with
#        seqnum = pieces[-1].split(" ")[-1].strip()
#        seqtree = pieces[-1].split("<")[-1].split(">")[0]
#        seq = seqtree + ":" + seqnum
        end_state = pieces[-1].split(">")[-1].strip() # get end state from file
        # (Show Origin must be checked when the cut sets are exported.)
        
        # Add this cut set's index to the list for its sequence and end state
#        if seq in sequences:
#            sequences[seq].append(index)
#        else:
#            sequences[seq] = [index]
        if end_state in end_states:
            end_states[end_state].append(index)
        else:
            end_states[end_state] = [index]
        
        # Add this cut set to the list of cut sets
        cut_sets.append((prob,events,"",end_state))
 
    # Calculate cut set probabilities (imitating SAPHIRE)
    # The main reason for this is to check that none of the basic event 
    # probabilities have been modified by flag sets; if they have, SAPHIRE's 
    # cut set probability will come out different than what we calculate here.
    calc_probs = [] # cut set probabilities that I will calculate here
    for cut_set in cut_sets:
        prb = 1.0
        for event in cut_set[1]:
            if event.startswith("/"):
                prb *= (1 - basic_events[event[1:]][0])
            else:
                prb *= basic_events[event][0]
        calc_probs.append(prb)
        # Check for match with SAPHIRE.
        # There's often a small difference, so report if > 1% and > 1e-4
        if abs(prb - cut_set[0]) > 1.e-4 and (prb == 0 or 
              abs(prb - cut_set[0]) / prb > 1.e-2):
            print("Probability mismatch in cut set ", len(calc_probs)-1, 
                  "SAPHIRE:",cut_set[0],"My calc:",prb)
 
    # Calculate a probability and frequency for each end state, taking into
    # account that the cut sets are not independent (SAPHIRE doesn't do this).
    end_state_probs = {}
    end_state_fv_importances = {}
    end_state_raw_importances = {}
    tot_freq = 0
    
    fv_list = ["1-L2-OP-SCG1-1", "1-L2-BE-PZRVSTUCK-SRV", "1-L2-OP-SAG2-1",
               "1-L2-OP-SAG1", "1-L2-BE-PZRVSTUCK-PORV", "1-L2-OP-SCG1-4",
               "/1-L2-BE-CCI-DISP","1-L2-BE-INDHLF-MP","1-L2-BE-INDSGTR-HDL",
               "/1-L2-BE-IVREC", "1-L2-BE-IVSE-LP", "1-L2-BE-IVSE-NLP",
               "1-L2-BE-RCP480GPM-DEP", "1-L2-BE-BMT-CHR", "1-L2-BE-BMT-NCHR",
               "1-L2-MCR-MOV-EQ5", "1-L2-MCR-PUMP-EQ5",
               "1-L2-MCR-MOV-EQ6", "1-L2-MCR-MOV-EQ4", "1-L2-MCR-PUMP-EQ4",
               "1-L2-MCR-PUMP-EQ6"]
    
    for end_state in end_states.items():
        
        # Make a list of basic events used in this end state's cut sets.
        es_bes = []
        for cs_index in end_state[1]:
            for event in cut_sets[cs_index][1]:
                #if not event.startswith("/") and not event in es_bes:
                if event in fv_list and not event in es_bes:
                    es_bes.append(event)
#        print(end_state[0],es_bes)
        
        # This line does the meat of the calculations:
        esprob = calc_probs_dependent(copy.deepcopy(cut_sets), basic_events, 
                                      end_state[1])

        end_state_probs[end_state[0]] = esprob
        tot_freq += esprob*init_freq
        
        # Calculate Fussell-Vessely & RAW importances for all basic events used
        if calcfv:
            fv_imps = {}
            raw_imps = {}
            for be in es_bes:
                # Fussell-Vessely
                # get lists of cut set indices containing and not containing be
                cs_with_be = []
                cs_wout_be = []
                for cs_ind in end_state[1]:
                    if be in cut_sets[cs_ind][1]:
                        cs_with_be.append(cs_ind)
                    else:
                        cs_wout_be.append(cs_ind)
                p_with_be = calc_probs_dependent(copy.deepcopy(cut_sets), 
                                                 basic_events, cs_with_be)
                p_wout_be = calc_probs_dependent(copy.deepcopy(cut_sets), 
                                                 basic_events, cs_wout_be)
                fv_standard = p_with_be / esprob
                fv_revised = 1.0 - p_wout_be / esprob
                if be == "1-L2-BE-PZRVSTUCK-SRV" and end_state[0] == "1-REL-NOCF":
                    pass
                fv_imps[be] = (fv_standard, fv_revised)
                
#                # RAW
#                be_old = basic_events[be]
#                basic_events[be] = (1.0, "") # set the basic event prob to 1.0
#                cuts_mod = copy.deepcopy(cut_sets)
#                for cs_ind in end_state[1]:
#                    if be in cuts_mod[cs_ind][1]:
#                        newcs = list(cuts_mod[cs_ind])
#                        newcs[0] /= be_old[0]   # change the cut set prob too
#                        cuts_mod[cs_ind] = tuple(newcs)
#                p1 = calc_probs_dependent(cuts_mod, basic_events, end_state[1])
#                raw_imps[be] = p1 / esprob
#                basic_events[be] = be_old

            end_state_fv_importances[end_state[0]] = fv_imps
#            end_state_raw_importances[end_state[0]] = raw_imps
                
        
    # Print results
    for end_state in sorted(end_state_probs, key=end_state_probs.__getitem__, 
                            reverse=True):
        print(end_state, "Frequency:", 
              "{0:.4e}".format(end_state_probs[end_state] * init_freq))
    
    if outfilename:
        with open(outfilename,"a") as outfile:
            outfile.write("\n" + filename)
            outfile.write("\nEnd State Frequencies\n")
            es_order = ["BMT","CIF","CIF-SC","ECF","ICF-BURN","ICF-BURN-SC",
                        "ISGTR","LCF","LCF-SC","NOCF","SGTR-C","SGTR-O",
                        "SGTR-O-SC","V","V-F","V-F-SC"]
            for es in es_order:
                esname = "1-REL-" + es
                if esname in end_state_probs:
                    outfile.write(esname + ", "+"{0:.4e}".format(end_state_probs[esname] * init_freq) + "\n")
                else:
                    outfile.write(esname + ", "+"0.00" + "\n")
    
    if calcfv:
        with open(outfilename,"a") as outfile:
            # Write Fussell-Vessely importances to file
            outfile.write("\nFussell-Vessely importances   , standard"+" "*16+", alternate\n")
            for es in end_state_fv_importances.items():
                outfile.write("\n" + es[0] + "\n")
                for ev in [(x, es[1][x]) for x in sorted(es[1], key=es[1].get, 
                           reverse=True)]:
                    fvst = str(ev[1][0])+" "*(24-len(str(ev[1][0])))+", "+str(ev[1][1])
                    outfile.write(ev[0]+" "*(30 - len(ev[0]))+", "+fvst+"\n")
#            # Write RAW importances to file
#            outfile.write("\nRAW importances\n")
#            for es in end_state_raw_importances.items():
#                outfile.write("\n" + es[0] + "\n")
#                for ev in [(x, es[1][x]) for x in sorted(es[1], key=es[1].get, 
#                           reverse=True)]:
#                    outfile.write(ev[0]+" "*(30 - len(ev[0]))+", "+str(ev[1])+"\n")

    print("Initiator Frequency:", "{0:.4e}".format(init_freq))
    print("Total Frequency:", "{0:.4e}".format(tot_freq))
    print("Inflation:", "{0:.3g}".format((tot_freq / init_freq - 1)*100), 
          "% from initiator")
    # For bins 6 and 7, CDF is similar to initiator frequency, so     
    # this inflation number is reasonably accurate.



