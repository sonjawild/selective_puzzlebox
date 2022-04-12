# selective_puzzlebox

This repository contains data generated and code for replicating all analyses for:

**'Manipulating actions: a selective two-option device for cognitive experiments in wild animals'**

Authors: **Sonja Wild, Gustavo Alarcón-Nieto, Michael Chimento, Lucy Aplin
**
Journal of Animal Ecology

**age.data.2021.txt**
- raw data to generate individual-level variables for birds with the following columns:
- Date: date (yyyy-mm-dd) of capture 
- Ring: individual metal EU-Ring that birds were equipped with
- Species: GRETI (=great tit); BLUTI (=blue tit); MARTI (=marsh tit); NUTHA (=nuthatch)
- Sex: M (male) or F (female), based on plumage    
- Pit: 10-digit alphanumeric code of the passive integrated transponder (PIT) tag of each individual bird
- Num.caught: Number of times each bird was captured
- Pit.year: Year within which the bird was equipped with a PIT tag (could be more than one year for the same bird, if PIT tags malfuncitoned and had to be replaced)
- hatch.year: Minimum hatch year calculated based on age at capture
- Age: age as either 'adult' or 'first.year' for each individual bird

**attributes_McDonalds_temp.txt** and **attributes_Mettnau_temp.txt**
- raw data file containing on birds present at experimental locations
- rownames: individual PIT tags
- species: GRETI (=great tit); BLUTI (=blue tit); MARTI (=marsh tit)
- solver: indicating whether the bird had solved prior to the restriction date (see manuscript); 
'no'= no solves; 'solver_left' or 'solver_right' if minimum of 10 solves prior to the restrition date; 'tutor' if knowledgeable bird from previous captive experiment
- num.solves: total number of solves up until the restriction date

**gmm.McDonalds.all.RData** and **gmm.Mettnau.all.RData**
RDA file containing raw data collected from RFID network feeders
- $gbi: group by individual matrix (1 for presence, 0 for absence)
- $metadata: details on start and endtime and location for each group
- $B: contains number of observations for each individual (same structure as gbi)

**latency_to_switch.R**
R code for replicating analyses on latencies to solve the puzzle box on the less preferred side after restriction and for reaching the solving rate from before
The large part at the start of the script that is commented out (#) was used to generate input files for survival analyses using the data described above and below.
This compiled data used for the survival analysis can be found in **latencies.ILVs.txt**

**latencies.ILVs.txt**
txt file containing the compiled individual level data needed for survival analyses:
- PIT: alphanumeric PIT code for each individual
- total_left: total number of solves (sliding door left)
- total_right: total number of solves (sliding door right)
- total_solves: number of solves in total across entire experiments
- visits_10_days_prior: number of visits (registered on receiving antenna on puzzle box) in the 10 experimental days prior to the restriction date
- solves_10_days_prior: number of solves in the 10 experimental days prior to the restriction date
- present_after: 'yes'for being present (registered on puzzle box antenna) after the restriction date
- num_visits_after: number of visits after the restriction date (registered on antenna)
- solves_after: total number of solves after the restriction date
- tutor: 'yes' if part of the independent captive experiment, 'no' if learned in the wild during this experiment
- preferred_side_prior: side preference prior to the restriction ('left' or 'right')
- pref_strength: strength of the prior side preference as proportion of solves on the preferred side
- latency.first.solve: number of visits after the restriction until the first solve occurred on the less preferred side
- latency.rate.reached: number of visits after the restriction until birds reached their solving rate (solves per visit) from prior to the restriction
- censored.first.solve: 1 if birds managed to produce a solution on their less preferred side, 0 if they failed to produce a solution
- censored.rate: 1 if birds reached their prior solving rate, 0 if not
- solve.rate.prior: solve rate as solves per visits in the 10 days prior to the restriciton
- solves.rate.after.total: solve rate as solves per visit after the restriction

**puzzle.data.RDA**
RDA file containing two slots $McDonalds and $Mettnau with the raw data collected at puzzle boxes
- PIT: individual PIT code
- event: 'arrived' if PIT tag picked up by antenna, 'left' or 'right' for solves (sliding the door); 'departed' if PIT tag no longer read by antenna; 
'displacement' if a bird displaced an existing bird on the antenna
- date: yyyy-mm-dd
- time: HH:MM:SS
- location: experimental puzzle box site ("McDonalds_back"; "McDonalds_front"; "Mettnau_back"; "Mettnau_front")
- date.time: yyyymmddHHMMSS
- week: number of experimental week since start of experiment
