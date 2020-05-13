# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 12:23:39 2018

@author: ejb
"""
import pickle
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

max_samples = 4375

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


end_samples = [0, 547, 1094, 1641, 2188, 2735, 3282, 3829, 4375]
for index in range(1,len(end_samples)):
    start_sample = end_samples[index - 1] + 1
    end_sample = end_samples[index]
    file_end = file_end = str(start_sample) + "-" + str(end_sample)
    point_es_freqs = pickle.load(open("seismic_results/point_es_freqs" + file_end + ".p","rb"))
    normal = pickle.load(open("seismic_results/normalized_es_freqs" + file_end + ".p","rb"))
    ft = pickle.load(open("seismic_results/ft_probs" + file_end + ".p","rb"))
    es = pickle.load(open("seismic_results/es_freqs" + file_end + ".p","rb"))
    cd = pickle.load(open("seismic_results/cd_freqs" + file_end + ".p","rb"))
    for sample in range(start_sample, end_sample+1):
        normalized_es_freqs[sample] = normal[sample]
        ft_probs[sample] = ft[sample]
        es_freqs[sample] = es[sample]
        cd_freqs[sample] = cd[sample]




if cd_freqs[1][1] == cd_freqs[2][1]:
    print("Problem with output files! Samples are duplicated!")

# Plot the CDF distribution 
for s in range(1, max_samples+1):
    for i in range(1,8):
        cd_freqs[s][0] += cd_freqs[s][i]
cd_dist = [cd_freqs[s][0] for s in range(1, max_samples+1)]
cdmax = max(cd_dist)
cdmin = min(cd_dist)

fig = plt.figure(1)
ax = fig.add_subplot(121) # nrows, ncols, position of this figure
cdhist = plt.hist(cd_dist, bins = np.geomspace(cdmin, cdmax, num=50), density=True, log=True, cumulative=False)
plt.yscale('linear')
plt.xscale('log')

ax2 = fig.add_subplot(122)
cdhistcum = plt.hist(cd_dist, bins = np.geomspace(cdmin, cdmax, num=50), density=True, log=True, cumulative=True)
plt.yscale('linear')
plt.xscale('log')


binfig = plt.figure(2)
for binnum in range(1,8):
    ax = binfig.add_subplot(7,2,binnum*2-1)
    bin_cd = [cd_freqs[s][binnum] for s in range(1, max_samples+1)]
    bin_cdhist = plt.hist(bin_cd, bins = np.geomspace(cdmin, cdmax, num=50), density=True, log=True, cumulative=False)
    plt.yscale('linear')
    plt.xscale('log')
    ax2 = binfig.add_subplot(7,2,binnum*2)
    print("Bin", binnum, min(bin_cd), "to", max(bin_cd))
    bin_cdcum = plt.hist(bin_cd, bins = np.geomspace(cdmin, cdmax, num=50), density=True, log=True, cumulative=True)
    plt.yscale('linear')
    plt.xscale('log')
plt.tight_layout()

fig = plt.figure(3)
ftvecs = {}
i = 0
for ft in ft_probs[1].keys():
    i += 1
    ax = fig.add_subplot(1,6,i)
    ftvecs[ft] = [ft_probs[s][ft] for s in range(1, max_samples+1)]
    ftmin = min(ftvecs[ft])
    ftmax = max(ftvecs[ft])
    print(ft, ftmin, "to", ftmax)
    fthist = plt.hist(ftvecs[ft], bins = np.geomspace(ftmin, ftmax, num=50), density=True, log=True)
    plt.yscale('linear')
    plt.xscale('log')
plt.tight_layout()

#fig = plt.figure(4)
esvecs = {}
pct5 = {}
pct50 = {}
pct95 = {}
means = {}
stdevs = {}
skewness = {}
kurtoses = {}
esnum = 0
for es in es_freqs[1][0].keys():
    esnum += 1
    fig = plt.figure(frameon=False,figsize=(6.4,4.8))
    ax = fig.add_subplot(1,1,1)
#    ax = fig.add_subplot(3,4,esnum)
    esvecs[es] = [es_freqs[s][0][es] for s in range(1, max_samples+1)]
    esmin = min(esvecs[es])
    esmax = max(esvecs[es])
    print(es, esmin, "to", esmax)
    eshist = plt.hist(esvecs[es], bins = np.geomspace(esmin, esmax, num=50), density=True, log=True)
    plt.yscale('linear')
    plt.xscale('log')
    plt.title(es)
    plt.ylabel("Probability Density")
    plt.xlabel("1/rcy")
#    plt.savefig("figs/" + es + "_fig.svg", bbox_inches="tight")
    pct5[es] = np.percentile(esvecs[es], 5)
    pct50[es] = np.percentile(esvecs[es], 50)
    pct95[es] = np.percentile(esvecs[es], 95)
    means[es] = np.mean(esvecs[es])
    stdevs[es] = np.std(esvecs[es])
    skewness[es] = sp.stats.skew(esvecs[es])
    kurtoses[es] = sp.stats.kurtosis(esvecs[es])

pct5["cd"] = np.percentile(cd_dist, 5)
pct50["cd"] = np.percentile(cd_dist, 50)
pct95["cd"] = np.percentile(cd_dist, 95)
means["cd"] = np.mean(cd_dist)
stdevs["cd"] = np.std(cd_dist)
skewness["cd"] = sp.stats.skew(cd_dist)
kurtoses["cd"] = sp.stats.kurtosis(cd_dist)

allvec = [sum(x) for x in zip(*esvecs.values())]
pct5["all"] = np.percentile(allvec, 5)
pct50["all"] = np.percentile(allvec, 50)
pct95["all"] = np.percentile(allvec, 95)
means["all"] = np.mean(allvec)
stdevs["all"] = np.std(allvec)
skewness["all"] = sp.stats.skew(allvec)
kurtoses["all"] = sp.stats.kurtosis(allvec)

lerf_rcs = ["CIF","CIF-SC", "ECF", "ISGTR", "SGTR-O", "SGTR-O-SC"]

lerf = [sum(x) for x in zip(*[esvecs["1-REL-" + rc] for rc in lerf_rcs])]
pct5["lerf"] = np.percentile(lerf, 5)
pct50["lerf"] = np.percentile(lerf, 50)
pct95["lerf"] = np.percentile(lerf, 95)
means["lerf"] = np.mean(lerf)
stdevs["lerf"] = np.std(lerf)
skewness["lerf"] = sp.stats.skew(lerf)
kurtoses["lerf"] = sp.stats.kurtosis(lerf)

pct5f = {}
pct50f = {}
pct95f = {}
meansf = {}
stdevsf = {}
skewnessf = {}
kurtosesf = {}
devpct = {}
for key in pct5.keys():
    pct5f[key] = "{0:.3g}".format(pct5[key])
    pct50f[key] = "{0:.3g}".format(pct50[key])
    pct95f[key] = "{0:.3g}".format(pct95[key])
    meansf[key] = "{0:.3g}".format(means[key])
    stdevsf[key] = "{0:.3g}".format(stdevs[key])
    skewnessf[key] = "{0:.3g}".format(skewness[key])
    kurtosesf[key] = "{0:.3g}".format(kurtoses[key])
    devpct[key] = int(round(stdevs[key] / means[key] * 100,0))
    
# this part requires L2_end_states to be left over from running montecarlomulti.py
#cut_set_count = {}
#for es in es_freqs[1][0].keys():
#    count = 0
#    for i in range(1,8):
#        if es in L2_end_states[i]:
#            count += len(L2_end_states[i][es])
#    cut_set_count[es] = count

point_ests = {"all": 1.18e-05, "cd": 1.18e-05, "BMT": 1.20e-7, "CIF": 6.88e-7,
              "CIF-SC": 2.32e-8, "ECF": 1.90e-9, "ICF-BURN":2.14e-6, "ICF-BURN-SC":2.13e-8,
              "ISGTR": 1.74e-7, "LCF": 6.51e-6, "LCF-SC": 6.25e-8, "NOCF": 1.89e-6,
              "SGTR-O":4.79e-8, "SGTR-O-SC": 1.21e-7}

point_ests["lerf"] = sum([point_ests[rc] for rc in lerf_rcs])


#plt.figure()
rel_cats = ["BMT","CIF","CIF-SC","ECF","ICF-BURN","ICF-BURN-SC",
            "ISGTR","LCF","LCF-SC","NOCF","SGTR-O",
            "SGTR-O-SC"]
vectors = [esvecs["1-REL-" + es] for es in rel_cats if es in rel_cats]
confs = [(pct5["1-REL-" + es],pct95["1-REL-" + es]) for es in rel_cats]
fifths = [pct5["1-REL-" + es] for es in rel_cats]
ninetyfifths = [pct95["1-REL-" + es] for es in rel_cats]
#plt.violinplot(vectors, showmeans=True, showmedians=True, points=1000)
quartile1, medians2, quartile3 = np.percentile(vectors, [25, 50, 75], axis=1)

fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
ax1.boxplot(vectors, whis=0.0, labels=rel_cats, showmeans=False, notch=False, showfliers=False)#, sym="_")
inds = np.arange(1, len(medians2) + 1)
means2 = [means["1-REL-" + es] for es in rel_cats]
#ax1.scatter(inds, [means["1-REL-" + es] for es in rel_cats], marker="D", color="r", zorder=3)
for i in range(len(medians2)):
    ax1.add_patch(plt.Rectangle([i+1-.15,medians2[i]],0.3,means2[i] - medians2[i],color="red"))
ax1.scatter(inds, [point_ests[es] for es in rel_cats], marker="P", color="#0070FF", zorder=3)
ax1.vlines(inds, quartile3, ninetyfifths, lw=2, color="k")
ax1.vlines(inds, fifths, quartile1, lw=2, color="k")
for i in range(len(medians2)):
    ax1.hlines([fifths[i], ninetyfifths[i]], i+1-.1, i+1+.1)
plt.yscale("log")
for tick in ax1.get_xticklabels():
    tick.set_rotation(15)
ax1.yaxis.grid(which="major", color='k', linestyle='-.', linewidth=0.5)
plt.title("Seismic release category frequency uncertainty ranges, 5000-sample Monte Carlo")

