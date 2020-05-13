# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 16:25:16 2018

@author: ejb
"""

end_states =["BMT","CIF","CIF-SC","ECF","ICF-BURN","ICF-BURN-SC","ISGTR",
             "LCF","LCF-SC","NOCF","SGTR-C","SGTR-O","SGTR-O-SC","V","V-F","V-F-SC"]

esfreq = []
esfreq.append({}) # First entry is the combined end state frequencies
                    # That way the others have index equal to their bin num
fvs = []    # Fussell-Vessely importances
fvs.append({})
for es in end_states:
    fvs[0]["1-REL-" + es] = {}

raws = []
raws.append({})
for es in end_states:
    raws[0]["1-REL-" + es] = {}

# Go through each bin's output file and read in the end state frequencies,
    # and the Fussel-Vessely importance and RAW importance for each basic
    # event used by that bin's cut sets
for binnum in range(1,8):
    filename = "bin" + str(binnum) + ".txt"
    with open(filename,'r') as file:
        lines = file.readlines()
    esfreq.append({})
    
    importance_type = None
    category = None
    fvcats = {}
    for es in end_states:
        fvcats["1-REL-" + es] = {}
    fvs.append(fvcats)
    rawcats = {}
    for es in end_states:
        rawcats["1-REL-" + es] = {}
    raws.append(rawcats)
    for line in lines:
        pieces = line.split(",")
        if len(pieces) > 1 and pieces[0].startswith("1-REL-"):
            esfreq[-1][pieces[0]] = float(pieces[1].strip())
        elif pieces[0].startswith("1-REL-"):
            category = pieces[0].strip()
        elif pieces[0].startswith("Fussell-Vessely"):
            importance_type = "FV"
        elif pieces[0].startswith("RAW importances"):
            importance_type = "RAW"
        elif len(pieces) > 1 and category is not None:
            if importance_type == "FV":
                fvs[binnum][category][pieces[0].strip()] = float(pieces[1].strip())
            if importance_type == "RAW":
                raws[binnum][category][pieces[0].strip()] = float(pieces[1].strip())
            
# Find the total release frequency in each bin
bin_totals = []
bin_totals.append(0)
for binnum in range(1,8):
    bin_totals.append(sum(item[1] for item in esfreq[binnum].items()))
bin_totals[0] = sum(bin_totals[1:]) # And the total for all bins

be_set = {}
basic_events = {}

combined_fvs = {}
combined_raws = {}

for relcat in ["1-REL-" + es for es in end_states]:
    # Find the release frequency in each release category for all bins combined
    esfreq[0][relcat] = sum([esfreq[binnum][relcat] for binnum in range(1,8)])
    
    # Make a list of basic events used by each end state (combining all bins)
    basic_events[relcat] = []
    be_set[relcat] = set()
    for binnum in range(1,8):
        for entry in fvs[binnum][relcat]:
            basic_events[relcat].append(entry)     
            be_set[relcat].add(entry)
    
    # Calculate importances for each of those basic events
    combined_fvs[relcat] = {}
    combined_raws[relcat] = {}
    for be in be_set[relcat]:
        # Calculate combined F-V importance
        denominator = esfreq[0][relcat]
        numerator = 0
        for binnum in range(1,8):
            if be in fvs[binnum][relcat]:
                numerator += esfreq[binnum][relcat] * fvs[binnum][relcat][be]
        fv_combined = numerator / denominator
        combined_fvs[relcat][be] = fv_combined
        
        #Calculate combined RAW importance
        denominator = esfreq[0][relcat]
        numerator = 0
        for binnum in range(1,8):
            if be in raws[binnum][relcat]:
                numerator += esfreq[binnum][relcat] * raws[binnum][relcat][be]
        raw_combined = numerator / denominator
        combined_raws[relcat][be] = raw_combined
    # Calculate Fussel-Vessely importances for the initiating events, too
    for binnum in range(1,8):
        if esfreq[0][relcat] > 0.0:
            combined_fvs[relcat]["1-IE-EQK-BIN-" + str(binnum)] = esfreq[binnum][relcat] / esfreq[0][relcat]
    # No RAW importances for initiating events, that would be nonsensical.
    # It would mean setting the iniator frequency to 1.0/rcy.
        
with open("out.txt","w") as outfile:
    outfile.write("\nFussel-Vessely importances, all bins combined:\n\n")
    for relcat in ["1-REL-" + es for es in end_states]:
        outfile.write("\n" + relcat + "\n")
        for ev in [(x, combined_fvs[relcat][x]) for x in sorted(combined_fvs[relcat], 
                           key=combined_fvs[relcat].get, reverse=True)]:
            outfile.write(ev[0]+" "*(30 - len(ev[0]))+", "+str(ev[1])+"\n")
    
    outfile.write("\nRAW importances, all bins combined:\n\n")
    for relcat in ["1-REL-" + es for es in end_states]:
        outfile.write("\n" + relcat + "\n")
        for ev in [(x, combined_raws[relcat][x]) for x in sorted(combined_raws[relcat], 
                           key=combined_raws[relcat].get, reverse=True)]:
            outfile.write(ev[0]+" "*(30 - len(ev[0]))+", "+str(ev[1])+"\n")
            