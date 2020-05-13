# Cutset_Requantification
An algorithm for getting more accurate PRA results than min-cut upper bound (for SAPHIRE output post-processing).

Motivation for external quantification
This process is an attempt to mitigate a problem that occurs in SAPHIRE release category frequency results for the Level 3 PRA project (specifically, for seismic bin 5, 6, and 7 initiators), namely, that the summed release category frequency results are much higher than the core damage frequency for the same initiator. In a perfect model, every core damage incident corresponds to one plant damage state and one release category, so the PDS frequency and release frequency should be exactly equal to the CDF. 
One significant source of this “frequency inflation” is the deletion of success terms in SAPHIRE’s default quantification—it assumes that successes occur with frequency approximately 1.0, which is often not the case following a severe seismic event. That problem can be partially resolved by using the I or W process flag to add success events to the cut set results, thereby reducing the calculated probability of cut sets that include successes. However, a reduction in cut set probability does not always sufficiently reduce the corresponding end state probability, due to the minimal cut set upper bound (MCUB) approximation SAPHIRE uses to calculate the combined probability of the cut sets.
As an example, take four cut sets in two end states (assume initiating event frequency of 1.0):
1.	A*/C -> ES1
2.	B*/C -> ES1
3.	A*C -> ES2
4.	B*C -> ES2
Each of the basic events A, B, and C has probability 0.5.
The total CCDP, combining all end states, should be 1-(1-A)(1-B) = 0.75, and the end state probabilities should be
PES1 = /C*(1-(1-A)(1-B)) = 0.375
PES2 = C*(1-(1-A)(1-B)) = 0.375
In SAPHIRE’s calculation, each cut set will be assigned probability 0.5 * 0.5 = 0.25, and when gathered each end state will have probability 1 – (1 – 0.25)(1 – 0.25) = 0.75, twice the correct value.
The BDD solver does not have this problem but has so far proven difficult to use on large models. Instead, the current process uses a script on SAPHIRE’s cut set output to quantify the conditional probability of each end state represented in those cut sets. 
The input for the script can be either solve or gather results from SAPHIRE; currently I am using the solve results. In some cases, using cut sets from SAPHIRE’s solve gives a lower end state probability than using cut sets from SAPHIRE’s gather, probably due to the presence of more shared events between cut sets. This is a point that could use more exploration.

 
Description of the external quantification process
The script’s quantification method is a variant of MCUB in which basic events shared between multiple cut sets are used as multipliers for the group of cut sets, rather than for each cut set individually. For example, to find the probability for ES2 in the example above, we first factor out the shared basic event C and calculate the combined probability of the two cut sets conditional on C:
PA or B = 1 – (1 – A)(1 – B) = 0.75
PES2 = C * PA or B = 0.5 * 0.75 = 0.375
In this simple case, that gives us the exact answer (and likewise for ES1).

Now consider a case where there are multiple shared events in the cut sets belonging to a particular end state:
1.	A*B*E -> ES1
2.	A*C*E -> ES1
3.	D*E -> ES1
In this case, there are at first two options: we can factor out the E, which is in all three cut sets, or factor out the A, which is in just two cut sets.
If we factor out A, we can calculate the probability of the first two cut sets as
P1 or 2 = A*(1 – B*E)(1 – C*E)
And the total probability, using MCUB, would be
PES1 = 1 – (1 – D*E)*(1 – P1 or 2) = 1 – (1 – D*E)(1 – A*(1 – B*E)(1 – C*E))
But, this would not be the exact answer. If instead we first factor out the E, then we can factor out the A as well for just the group of cut sets that share that event as well:
PES1 = E * (1 – (1 – D)*(1 – PAB or AC))
where PAB or AC = A*(1 – (1 – B)*(1 – C)). 
This factoring can be performed recursively any number of times, starting with the largest groups of cut sets sharing a basic event and then operating on subgroups that share another event, as long as some shared event remains within the group of cut sets after factoring. When the remaining cut sets are all independent, their probability is calculated by MCUB as it would be in SAPHIRE.

In some cases, it is not obvious which event to factor out first. For example,
1.	A*B -> ES3
2.	A*C -> ES3
3.	A*D -> ES3
4.	D*E -> ES3
5.	E*F -> ES3
6.	E*G -> ES3
Here, there are three shared events: A (cut sets 1, 2, and 3), D (cut sets 3 and 4), and E (cut sets 4, 5, and 6). It is not possible to factor out all of them simultaneously; factoring out any one of those three creates a group with no shared events. So we can either factor out the A and the E and get
PES3 = 1 – (1 – A*PB or C or D)(1 – E*PD or F or G)
Or, we can factor out the D and get
PES3 = 1 – (1 – D*PA or E)(1 – A*PB or C)(1 – E*PF or G)
The choice of which common event to account for is somewhat arbitrary, since none of them will give the exact answer. In the current implementation, the first event to be factored out is the one for which the sum of the probabilities of the cut sets that contain it is greatest; so, the second option would only be taken if
P3 + P4 > P1 + P2 + P3 	(cut sets containing D have greater probability than those containing A)
and
P3 + P4 > P4 + P5 + P6 	(cut sets containing D have greater probability than those containing E)
The ranking of shared events uses the rare event approximation. Another approach would be to just start with the event that appears in the greatest number of cut sets; this variation generally results in a slightly higher end state frequency estimate, because factoring out events that are common among low-frequency cut sets may prevent factoring out other events that are in a smaller number of higher-frequency cut sets.

Note: Some strange behavior seems to occur with compound events, since MARD doesn’t necessarily output the correct probabilities for them. This is generally not important.

