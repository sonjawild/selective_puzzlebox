
# #### Manipulating actions: a selective two-action device for cognitive experiments in wild animals #### ------------------------------
# Sonja Wild, Gustavo Alarcon-Nieto, Michael Chimento, Lucy Aplin

# 1) Data preparation ------------------------------------------------

# 1.1 (Install) and load necessary libraries ------------------------------------------------
library(data.table)
library(ggplot2)
library(asnipe)
library(car)
library(MuMIn)
library(coxme)
library(egg)

# set working directory
setwd("C:/Users/swild/Desktop/Konstanz/Restricted puzzle box/git/restricted.puzzle.box")


# THE FOLLOWING CODE EXTRACTS DATA FOR ANALYSIS FROM RAW DATA
# BUT IS NOT NEEDED TO RUN THE COX PROPORTIONAL HAZARD MODELS
# USERS CAN DIRECTLY JUMP TO LINE 712 AND IMPORT THE EXTRACTED DATA DIRECTLY

# # upload species and age list
# species.list <- read.delim("age.data.2021.txt", sep=",")
# # contains information on birds' ring number, PIT tag, sex and age
# GT.list <- subset(species.list$Pit, species.list$Species=="GRETI")
# 
# # load data from puzzle boxes
# load("puzzle.data.RDA")
# 
# # reassign puzzle box data to each site (Mettnau ME and McDonalds MC)
# Mettnau.puzzle.data <- puzzle.data$Mettnau
# McDonalds.puzzle.data <- puzzle.data$McDonalds
# 
# 
# # 2) Create a data summary -----------------------------------------------------
# 
# # 2.1. Extract how many great tits have visited in total -------------------------
# 
# # create a function that adds visits to the data frame (starting with arrival of bird and ending with departure)
# add.visits <- function(puzzle.data){
#   puzzle.data <- subset(puzzle.data, puzzle.data$PIT %in% GT.list)
#   puzzle.data$visit <- NA
#   puzzle.data[1, "visit"] <- 1 # assign a 1 to the first visit
#   for(i in 2:length(puzzle.data[,1])){
#     event <- puzzle.data[i, "event"]
#     if(event %in% c("departed", "right", "left", "scrounging") &
#        puzzle.data[i, "PIT"] == puzzle.data[i-1, "PIT"] &
#        puzzle.data[i, "location"] == puzzle.data[i-1, "location"]){
#       puzzle.data[i, "visit"] <- puzzle.data[i-1, "visit"]
#     } else {
#       puzzle.data[i, "visit"] <- puzzle.data[i-1, "visit"]+1
#     }
#   }
#   return(puzzle.data)
# }
# 
# # run the function on ME
# Mettnau.puzzle.data <- add.visits(puzzle.data = Mettnau.puzzle.data)
# 
# # run the function on MC
# McDonalds.puzzle.data <- add.visits(puzzle.data = McDonalds.puzzle.data)
# 
# max(Mettnau.puzzle.data$visit)
# # 24'429 visits by great tits in Mettnau
# 
# max(McDonalds.puzzle.data$visit)
# # 26'732 visits by great tits at McDonalds
# 
# 
# # 2.2.  Number of great tits -------------------------------------------------------------------
# # how many birds were registered on the puzzle box at each site?
# num.birds.Mettnau <- length(unique(Mettnau.puzzle.data$PIT[Mettnau.puzzle.data$PIT %in% GT.list]))
# num.birds.Mettnau
# # 87
# 
# num.birds.McDonalds <- length(unique(McDonalds.puzzle.data$PIT[McDonalds.puzzle.data$PIT %in% GT.list]))
# num.birds.McDonalds
# # 74
# 
# 
# # 2.3. Number of solves  -------------------------------------
# 
# # how many solves were performed in total
# length(subset(Mettnau.puzzle.data$PIT, Mettnau.puzzle.data$event %in% c("left", "right")& Mettnau.puzzle.data$PIT %in% GT.list))
# length(subset(Mettnau.puzzle.data$PIT, Mettnau.puzzle.data$event %in% c("scrounging")))
# # 11'649 solves at Mettnau and 7577 scrounges
# 
# length(subset(McDonalds.puzzle.data$PIT, McDonalds.puzzle.data$event %in% c("left", "right") & McDonalds.puzzle.data$PIT %in% GT.list))
# length(subset(McDonalds.puzzle.data$PIT, McDonalds.puzzle.data$event %in% c("scrounging")))
# # 13'431 solves at McDonalds and 6807 scrounges
# 
# 
# # 2.4. How many solvers? --------------------------------------------------
# 
# # how many learned to solve across the entire experiment (with a minimum of 10 solves)
# 
# # Mettnau
# learners.Mettnau <- NULL
# 
# for(i in unique(Mettnau.puzzle.data$PIT[Mettnau.puzzle.data$PIT %in% GT.list])){
#   # extract the number of solves
#   num.solves.right <- length(subset(Mettnau.puzzle.data$event, Mettnau.puzzle.data$event %in% c("right") & Mettnau.puzzle.data$PIT==i))
#   num.solves.left <- length(subset(Mettnau.puzzle.data$event, Mettnau.puzzle.data$event %in% c("left") & Mettnau.puzzle.data$PIT==i))
#   num.scrounges <- length(subset(Mettnau.puzzle.data$event, Mettnau.puzzle.data$event %in% c("scrounging") & Mettnau.puzzle.data$PIT==i))
#   num.solves.total <- num.solves.left + num.solves.right
# 
#   i.data <- c(i, num.solves.right, num.solves.left, num.scrounges, num.solves.total)
#   learners.Mettnau <- rbind.data.frame(learners.Mettnau, i.data)
# }
# 
# colnames(learners.Mettnau) <- c("PIT", "solves.right", "solves.left", "scrounges", "solves.total")
# head(learners.Mettnau)
# 
# learners.Mettnau.with.10.solves <- subset(learners.Mettnau$PIT, learners.Mettnau$solves.total>=10)
# length(learners.Mettnau.with.10.solves)
# # [1] 36 solvers with a minimum of 10 solves (great tits only)
# 
# 
# # McDonalds
# learners.McDonalds <- NULL
# 
# for(i in unique(McDonalds.puzzle.data$PIT[McDonalds.puzzle.data$PIT %in% GT.list])){
#   # extract the number of solves
#   num.solves.right <- length(subset(McDonalds.puzzle.data$event, McDonalds.puzzle.data$event %in% c("right") & McDonalds.puzzle.data$PIT==i))
#   num.solves.left <- length(subset(McDonalds.puzzle.data$event, McDonalds.puzzle.data$event %in% c("left") & McDonalds.puzzle.data$PIT==i))
#   num.scrounges <- length(subset(McDonalds.puzzle.data$event, McDonalds.puzzle.data$event %in% c("scrounging") & McDonalds.puzzle.data$PIT==i))
#   num.solves.total <- num.solves.left + num.solves.right
# 
#   i.data <- c(i, num.solves.right, num.solves.left, num.scrounges, num.solves.total)
#   learners.McDonalds <- rbind.data.frame(learners.McDonalds, i.data)
# }
# 
# colnames(learners.McDonalds) <- c("PIT", "solves.right", "solves.left", "scrounges", "solves.total")
# head(learners.McDonalds)
# 
# learners.McDonalds.with.10.solves <- subset(learners.McDonalds$PIT, learners.McDonalds$solves.total>=10)
# length(learners.McDonalds.with.10.solves)
# # [1] 33 solvers with a minimum of 10 solves (great tits only)
# 
# 
# # how many were demonstrators (that learned how to solve during a previous independent experiment in captivity)
# 
# demos <- c(
#   "0700ED9C1D",
#   "0700EDCDCB",
#   "0700EE1990",
#   "0700EDAC82",
#   "0700EE41D3",
#   "0700EDB3F2",
#   "0700EE261C",
#   "0700EE2D2A",
#   "0700EE3E71",
#   "0700EDA852",
#   "0700EE3D79",
#   "0700EDB3D5"
# )
# 
# # for Mettnau
# length(which(learners.Mettnau.with.10.solves %in% demos))
# # 3 of the 36 were demonstrators (which leaves 33 learners)
# 
# length(which(learners.McDonalds.with.10.solves %in% demos))
# # and 6 of the 33 in McDonalds (which leaves 27 learners)
# 
# 
# 
# # 3) Extract latencies to start solving after restriction --------
# 
# # first load the list of the birds that had been restricted on the 22nd of February
# McD.restr <- read.delim("attributes_McDonalds_temp.txt", sep=",")
# Mettnau.restr <- read.delim("attributes_Mettnau_temp.txt", sep=",")
# # all those with the attribute 'tutor' or 'solver_left' or 'solver'right' were restricted to their less preferred side
# 
# # remove those that have not learned to solve (<10 solves)
# # and remove all others but great tits
# restr.PIT.McD <- rownames(subset(McD.restr, McD.restr$num.solves>=10 & McD.restr$species=="GRETI"))
# restr.PIT.Mettnau <- rownames(subset(Mettnau.restr, Mettnau.restr$num.solves>=10 & Mettnau.restr$species=="GRETI"))
# 
# # NOTE: two weeks prior corresponds to 10 experimental days (as boxes were removed for 2x48 hours)
# two.weeks.prior.Mettnau <-
#   subset(
#     Mettnau.puzzle.data,
#     Mettnau.puzzle.data$date > as.POSIXct("2021-02-06") &
#       as.POSIXct("2021-02-22") > Mettnau.puzzle.data$date &
#       Mettnau.puzzle.data$event %in% c("right", "left")
#   )
# # check whether the birds in the restricted file were present in the last 10 days before placing the restriction
# # and were present after the restriction
# which(restr.PIT.Mettnau %in% two.weeks.prior.Mettnau$PIT==FALSE )
# # 1 bird (number 8) was not solving in the two weeks before the restriction was put in place
# # exclude those from the analysis
# restr.PIT.Mettnau <- restr.PIT.Mettnau[-which(restr.PIT.Mettnau %in% two.weeks.prior.Mettnau$PIT==FALSE)]
# length(restr.PIT.Mettnau)
# # at Mettnau, 12 knowledgeable birds were solving in the two weeks prior to the restriction
# 
# 
# # repeat for McDonalds
# two.weeks.prior.McD <-
#   subset(
#     McDonalds.puzzle.data,
#     McDonalds.puzzle.data$date > as.POSIXct("2021-02-06") &
#       as.POSIXct("2021-02-22") > McDonalds.puzzle.data$date &
#       McDonalds.puzzle.data$event %in% c("right", "left")
#   )
# # check whether the birds in the restricted file were present in the last 10 days before placing the restriction
# which(restr.PIT.McD %in% two.weeks.prior.McD$PIT==FALSE)
# # integer(0) means that all birds were present
# length(restr.PIT.McD)
# # at McD, 18 birds were solving in the two weeks prior to restriction being put in place
# 
# # Analysis will be restricted to those 12+18 birds
# 
# # as a next step, we need to examine the block size of visits that are best to use (see SI):
# # we here extract two latencies: the latency until the first solve on the less preferred side after the restriction
# # and the latency until the bird reaches the same solving rate from prior to the restriction
# # for the latter, we need to determine how we best calculate the solving rate
# # to do so, we analyse the solving rates from before the restriction
# # we split up each bird's visits into 'blocks' of different sizes (10, 20, 30, 40, 50)
# # and calculate the solving rate within these blocks (as # solves per visits). We then extract the variance
# # in solving rate for each of these block sizes, and compare whether it significantly differs from
# # one block size to the next. (see below for a visual representation)
# # in small block sizes, we expect a lot of noise in solving rate (which increases variance)
# # in larger block sizes, we expect low variances. The optimal block size is where the variance
# # remains stable from one to the next block size
# 
# # we have written a function that calculates solving rates and variances for different block sizes automatically
# inv.block.sizes <- function(puzzle.data, block.sizes, restr.PIT){
#   var.all <- NULL
#   for(i in unique(puzzle.data$PIT)){
#     # we extract the data for individual i in the 10 days prior to the restriction
#     sub.10.days.prior.i <-
#       subset(
#         puzzle.data,
#         puzzle.data$date > as.POSIXct("2021-02-06") &
#           as.POSIXct("2021-02-22") > puzzle.data$date & puzzle.data$PIT == i
#       )
#     solve.rate.10.days.prior <- length(subset(sub.10.days.prior.i$event, sub.10.days.prior.i$event %in% c("left", "right")))/length(unique(sub.10.days.prior.i$visit))
#     if(length(sub.10.days.prior.i[,1])>0 & i %in% restr.PIT){ # if the bird was present in the 10 days before restriction
#       var.j <- NULL
# 
#       for(j in block.sizes){
#         # we try different block sizes from 10 to 50 visits and calculate the variance in solving rate for each individual i and each block size
#         # ensure that the data is ordered by date and time
# 
#         sub.10.days.prior.i <- sub.10.days.prior.i[order(sub.10.days.prior.i$date.time),]
#         visits.i <- unique(sub.10.days.prior.i$visit)
#         how.many.blocks.i <- floor(length(visits.i)/j) # extract how many times we can partition the data for i given the block size j
#         solve.rates.i <- NULL
#         for(z in 1:how.many.blocks.i){
#           visits.z <- visits.i[(z*j-j+1):(z*j)]
#           solve.rate.i.j.z <- length(subset(sub.10.days.prior.i$event, sub.10.days.prior.i$visit %in% visits.z & sub.10.days.prior.i$event %in% c("right", "left")))/length(visits.z)
#           solve.rates.i[z] <- solve.rate.i.j.z
#         }
#         # extract the proportion of solve rates that fall within 10% of the overall solve rate of individual i
#         solve.rate.deviance.i <- 100/solve.rate.10.days.prior*(solve.rates.i-solve.rate.10.days.prior)/100
#         prop.within.10.percent <- length(subset(solve.rate.deviance.i, abs(solve.rate.deviance.i)<=0.1))/length(solve.rate.deviance.i)
#         # then calculate the variance of those values
# 
#         var.i.j <- cbind.data.frame(i, round(j, 0), solve.rate.10.days.prior, var(solve.rates.i), prop.within.10.percent, z)
#         var.j <- rbind(var.j, var.i.j)
#       }
# 
#       var.all <- rbind(var.all, var.j)
#     }
# 
#   }
#     colnames(var.all) <- c("PIT", "block.size", "solve.rate.prior", "variance", "prop.within.10.percent", "num.blocks")
#     return(var.all)
# }
# 
# # run on all puzzle data together (MC and ME)
# # NOTE: this function takes quite a while to run,
# # we have added the RDA object to be loaded directly
# # block.size.analysis <- inv.block.sizes(puzzle.data = rbind(Mettnau.puzzle.data, McDonalds.puzzle.data), restr.PIT = c(restr.PIT.Mettnau, restr.PIT.McD), block.sizes = c(10, 20, 30,40, 50))
# # save(block.size.analysis, file="block.size.analysis.RDA")
# load("block.size.analysis.RDA")
# head(block.size.analysis)
# # each line represent a different bird for a different block size (e.g. line one bird 0700EE261C for a block size of 10)
# 
# # plot as boxplot
# # we can see that the variance of the solving rate gradually decreases
# # but then does not change from 40 to 50.
# # png(filename="block.sizes.png", width = 500, height = 400)
# 
# boxplot( block.size.analysis$variance ~ block.size.analysis$block.size, ylab="Variance (solving rate)", xlab="block size (# visits)")
# 
# # dev.off()
# 
# # we do a more quantitative approach and compare the different block sizes by paired t-tests
# # if they are not significantly different, we have found the optimum block size
# 
# # write a function to make t-tests to compare different block sizes
# compare.t.test <- function(block.sizes){
#   sub1 <- subset(block.size.analysis, block.size.analysis$block.size==block.sizes[1] & !is.na(block.size.analysis$variance) & block.size.analysis$num.blocks>=3)
#   sub2 <- subset(block.size.analysis, block.size.analysis$block.size==block.sizes[2] & !is.na(block.size.analysis$variance) & block.size.analysis$num.blocks>=3)
#   # subset to only those PITs who have variances in both
#   PITs.to.include <- intersect(sub1$PIT, sub2$PIT)
# 
#   t1.variance <- subset(block.size.analysis, block.size.analysis$block.size==block.sizes[1] & block.size.analysis$PIT %in% PITs.to.include)
#   t2.variance <- subset(block.size.analysis, block.size.analysis$block.size==block.sizes[2] & block.size.analysis$PIT %in% PITs.to.include)
# 
#   t.test.out <- t.test(t1.variance$variance, t2.variance$variance, paired=TRUE, alternative="two.sided")
#   if(t.test.out$p.value<0.05){
#     paste("The variances for block sizes ", block.sizes[1], " and ", block.sizes[2], " are significantly different (p=", t.test.out$p.value, ", t=", round(t.test.out$statistic, 3), ", effect size=", round(t.test.out$estimate, 3), ", N=", length(t1.variance$PIT), ").", sep="")
# 
#   } else {
#     paste("The variances for block sizes ", block.sizes[1], " and ", block.sizes[2], " do not differ significantly (p=", t.test.out$p.value, ", t=", round(t.test.out$statistic, 3), ", effect size=", round(t.test.out$estimate, 3), ", N=", length(t1.variance$PIT), ").", sep="")
#     }
#   }
# 
# # 10 vs 20
# compare.t.test(block.sizes = c(10,20))
# # [1] "The variances for block sizes 10 and 20 are significantly different (p=5.74921783514934e-06, t=5.577, effect size=0.019, N=29)."
# 
# compare.t.test(block.sizes = c(20,30))
# # [1] "The variances for block sizes 20 and 30 are significantly different (p=0.00100229146755767, t=3.744, effect size=0.009, N=25)."
# 
# compare.t.test(block.sizes = c(30,40))
# # "The variances for block sizes 30 and 40 are significantly different (p=0.00614362180715136, t=3.031, effect size=0.012, N=23)."
# 
# compare.t.test(block.sizes = c(40,50))
# # [1] "The variances for block sizes 40 and 50 do not differ significantly (p=0.0581850679852026, t=2.016, effect size=0.005, N=20)."
# 
# # the variance is not significantly different from 40 to 50 (block size), therefore, we continue with a block size of 40
# # (which maximizes precision and minimizes variance)
# 
# # We go on to extract various individual-level variables:
# # - PIT = individual PIT code
# # - total_left = total solutions 'left' over the entire duration of the experiment
# # - total_right = total solutions 'right' over the entire duration of the experiment
# # - total_solves = the sum of total_left and total_right
# # - visits_10_days_prior = number of visits in the 10 days before the restriction was put in place
# # - solves_10_days_prior = number of solves (both sides) in the 10 days before the restriction was put in place
# # - present_after = yes/no whether the bird was registered on the puzzle box after the restriction
# # - num_visits_after = the number of visits after the restriction was in place
# # - solves_after = the number of solves (both sides) after the restriction
# # - tutor = whether the bird was a tutor in a previous experiment (yes, no)
# # - preferred_side_prior = which side it preferred prior to the restriction (left, right)
# # - pref_strength = the proportion of solves on the preferred side (considering only solves before the restriction)
# # - latency_first_solve = the minutes from first landing on the box after the restriction until the first solve occurred
# # - latency.rate.reached = the minutes from the first landing on the box after the restriction until the solving rate was reached from before (-10%), always considering 50 visits at the time
# # - censored.first.solve = whether the bird managed to produce a solve on the opposite side after the restriction 1/0
# # - censored.rate = whether the bird reached the solving rate at all after the restriction (-10%) - 1/0
# # - solve.rate.prior = proportion of solves/ total number of visits in the 10 days prior to restriction
# # - solve.rate.after.total = proportion of solves / total number of visits after the restriction (can be above 1 if bird solved several times per visit)
# # - solve.rate.in.block = the solving rate in the block in which the bird either reached the solving rate from prior or if it never reached it, it reports the solve rate from the last block (not useful for analysis)
# # - solve.rate.after.rate.reached = reports the solving rate after the bird reached the solving rate from prior for the rest of the experiments
# 
# 
# # this is all automated in a function
# extract <- function(puzzle.data, restriction.date, restr.PIT, latency.rate=FALSE, block){
#   restricted <- NULL
#   for(i in restr.PIT){
#     # get the data for i after the restriction was put in place
#     sub.i.after <- subset(puzzle.data, puzzle.data$PIT ==i & puzzle.data$date > restriction.date )
#     # extract if the bird was present after the restriction (yes or no)
#     if(length(sub.i.after[,1])!=0){
#       present.after <- "yes"
#     } else {
#       present.after <- "no"
#     }
# 
#     # extract the data 10 (experimental) days prior to the restriction
#     sub.10.days.prior.i <-
#       subset(
#         puzzle.data,
#         puzzle.data$date > as.POSIXct("2021-02-06") &
#           as.POSIXct("2021-02-22") > puzzle.data$date & puzzle.data$PIT == i
#       )
#     # get the data for i across the entire experiment
#     sub.i.all <- subset(puzzle.data, puzzle.data$PIT ==i )
# 
#     # extract the number of solves over the entire time of the experiment
#     num.solves.total <-
#       length(sub.i.all[sub.i.all$event %in% c("right", "left"), ][, 1])
#     num.right.all <- length(subset(sub.i.all, sub.i.all$event %in% "right")[,1])
#     num.left.all <- length(subset(sub.i.all, sub.i.all$event %in% "left")[,1])
# 
#     # extract the number of solves in the 10 days prior to the experiment
#     num.solves.10.days.prior <-
#       length(sub.10.days.prior.i[sub.10.days.prior.i$event %in% c("right", "left"), ][, 1])
# 
# 
#     # extract the solves after the restriction
#     solves.right <- length(subset(sub.i.after$event, sub.i.after$event =="right"))
#     solves.left <- length(subset(sub.i.after$event, sub.i.after$event =="left"))
#     solves.total <- solves.right+solves.left
# 
#     if(solves.total==0){
#       censored.first.solve <- 0
#     } else {
#       censored.first.solve <- 1
#     }
# 
#     # extract which was preferred up until the restriction
#     sub.i.prior <- subset(puzzle.data, as.POSIXct("2021-02-22") > puzzle.data$date & puzzle.data$PIT ==i)
#     prior.left <- length(subset(sub.i.prior, sub.i.prior$event %in% "left")[,1])
#     prior.right <- length(subset(sub.i.prior, sub.i.prior$event %in% "right")[,1])
# 
# 
#     if(prior.right>prior.left){
#       preferred.side <- "right"
#     } else {
#       preferred.side <- "left"
#     }
# 
#     # extract the strength of the preferred side
#     prop.preferred <- max(prior.right, prior.left)/(prior.right + prior.left)
# 
#     # extract whether bird had been a tutor (meaning that it was restricted for most of the experiment)
# 
#     if(i %in% demos){
#       tutor <- "yes"
#     } else {
#       tutor <- "no"
#     }
# 
#     # extract when it first arrived on the puzzle box after the restriction
#     first.arrived <- min(sub.i.after$date.time) # extract the time it first interacted with the puzzle box after restriction
# 
#     # extract when it first solved after the restriction was put in place
#     first.solve <- min(sub.i.after$date.time [sub.i.after$event %in% c("right", "left")])
# 
#     # extract the latency between landing after the restriction and solving for the first time
# 
#     latency.to.first.solve <-
#       difftime(
#         as.POSIXct(first.solve, format = "%y%m%d%H%M%S"),
#         as.POSIXct(first.arrived, format = "%y%m%d%H%M%S"),
#         units = "mins"
#       )
#     # latency is in minutes
#     latency.to.first.solve <- as.numeric(latency.to.first.solve)
# 
#     # also extract on which visit after the first visit it made the first solve
#     if(is.na(latency.to.first.solve)){
#       latency.visits <- length(unique(sub.i.after$visit))
#     } else {
#       latency.visits <- length(unique(subset(sub.i.after$visit, sub.i.after$date.time>=first.arrived & sub.i.after$date.time<=first.solve)))
#     }
# 
# 
#     # extract the rate of solving (as solves per visits) in the 10 days before the restriction
#     num.visits.10.prior <- length(unique(sub.10.days.prior.i$visit))
# 
#     solving.rate.prior <-
#       length(subset(
#         sub.10.days.prior.i$event,
#         sub.10.days.prior.i$event %in% c("right", "left")
#       )) / num.visits.10.prior
#     # this is the proportion of solves per visit to the puzzle box
# 
#     # extract the solving rate after the restriction
#     # across all of the 10 days after
#     num.visits.after <- length(unique(sub.i.after$visit))
# 
#     solving.rate.after <-
#       length(subset(sub.i.after$event, sub.i.after$event %in% c("right", "left"))) /
#       num.visits.after
# 
#     # extract the moment at which the bird reached the solving rate from before the restriction
# 
#     if(latency.rate==TRUE & present.after=="yes" & solves.total!=0 ){ # if the birds was present after and has managed to solve after the restriction
# 
#       # calculate the rate of solving after first landing on the box
# 
#       j <- block
# 
#       count <- 1
# 
#       repeat{
#         # extract the solving rate in the first 10 visits after the first solve
#         # if it has reached the solving rate from prior to the restriction (-10%)
#         # then we consider the solving rate reached and we extract the latency
#         # as the number of visits since the restriction until the solving rate from prior was reached
#         # if the solving rate is not reached, we take the first 20 visits after the first solve etc.
# 
# 
#         # ensure that the data is ordered according to date and time
#         sub.i.after <- sub.i.after[order(sub.i.after$date.time),]
# 
#         # extract the visits after the restriction during which the bird solved
#         solving.visits <- unique(subset(sub.i.after$visit, sub.i.after$event %in% c("right", "left")))
# 
#         solving.after.first.solve.i <- subset(sub.i.after, sub.i.after$date.time>=first.solve)
#         # extract all visits after the restriction
#         all.visits.after.i <- unique(sub.i.after$visit)
# 
#         start.visit <- solving.visits[count]
# 
#         # extract the first 40 visits (block size) from the start visit
#         which.visit.start <- which(all.visits.after.i==start.visit)
#         which.visits.j <- all.visits.after.i[which.visit.start:(which.visit.start+block-1)]
# 
# 
#         # extract the first j visits after the first solve and calculate the solving rate (per visits)
#         solving.rate.j <-
#           length(subset(
#             solving.after.first.solve.i$event,
#             solving.after.first.solve.i$visit %in% which.visits.j &
#               solving.after.first.solve.i$event %in% c("right", "left")
#           )) / length(which.visits.j)
# 
#         # if the solving rate has been reached, then we extract the latency of the number of visits that were needed to reach this rate
# 
#         if(solving.rate.j >= solving.rate.prior-(0.1*solving.rate.prior)){
# 
#           # extract the date.time of the first solve from the batch in which the solving rate has been reached
#         date.time.solve.rate.reached <- min(subset(solving.after.first.solve.i$date.time, solving.after.first.solve.i$event %in% c("right", "left")
#                & solving.after.first.solve.i$visit %in% which.visits.j))
# 
#         latency.visits.rate.reached <- length(unique(subset(sub.i.after$visit, sub.i.after$date.time>=first.arrived & sub.i.after$date.time<=date.time.solve.rate.reached)))
# 
#         # also extract the solve rate after the solve rate from prior has been reached
#         solve.rate.after.switch <- length(
#           subset(
#             sub.i.after$event,
#             sub.i.after$date.time >= date.time.solve.rate.reached &
#               sub.i.after$event %in% c("right", "left")
#           )
#         ) / length(unique(
#           subset(
#             sub.i.after$visit,
#             sub.i.after$date.time >= date.time.solve.rate.reached
#           )
#         ))
# 
#         censored <- 1
#         break
#         } else if(length(which.visits.j)<block | is.na(which.visits.j[block])){
#           latency.visits.rate.reached <- num.visits.after
#           censored <- 0
#           solve.rate.after.switch <- 0
# 
#           break
#         } else {
#           j <- j+block
#           count <- count+1
#         }
# 
#       }
# 
#     } else {
#       latency.visits.rate.reached <- num.visits.after
#       solve.rate.after.switch <- 0
#       censored <- 0
#       solving.rate.j <- 0
#     }
# 
# 
# 
#     data.i <- c(
#       i,
#       num.left.all,
#       num.right.all,
#       num.left.all + num.right.all,
#       num.visits.10.prior,
#       num.solves.10.days.prior,
#       present.after,
#       num.visits.after,
#       solves.total,
#       tutor,
#       preferred.side,
#       round(prop.preferred, 2),
#       latency.visits,
#       latency.visits.rate.reached,
#       censored.first.solve,
#       censored,
#       round(solving.rate.prior,2),
#       round(solving.rate.after,2),
#       round(solving.rate.j, 2),
#       round(solve.rate.after.switch, 2)
#     )
#     restricted <- rbind.data.frame(restricted, data.i)
# 
#     colnames(restricted) <-
#       c(
#         "PIT",
#         "total_left",
#         "total_right",
#         "total_solves",
#         "visits_10_days_prior",
#         "solves_10_days_prior",
#         "present_after",
#         "num_visits_after",
#         "solves_after",
#         "tutor",
#         "preferred_side_prior",
#         "pref_strength",
#         "latency.first.solve",
#         "latency.rate.reached",
#         "censored.first.solve",
#         "censored.rate",
#         "solve.rate.prior",
#         "solve.rate.after.total",
#         "solve.rate.in.block",
#         "solve.rate.after.rate.reached"
#       )
# 
# 
#     # remove those that weren't present after the restriction
#     restricted <- subset(restricted, restricted$present_after!="no")
# 
#   }
#   return(restricted)
# }
# 
# # specify the restriction date
# restriction.date <- as.POSIXct("210222", format="%y%m%d")
# 
# # run on Mettnau data (warnings can be safely ignored)
# Mettnau.restr.data <-
#   extract(
#     puzzle.data = Mettnau.puzzle.data,
#     restriction.date = restriction.date,
#     restr.PIT = restr.PIT.Mettnau,
#     latency.rate = TRUE,
#     block = 40
#   )
# 
# # run on McD data (warnings can be safely ignored)
# McDonalds.restr.data <-
#   extract(
#     puzzle.data = McDonalds.puzzle.data,
#     restriction.date = restriction.date,
#     restr.PIT = restr.PIT.McD,
#     latency.rate = TRUE,
#     block = 40
#   )
# 
# Mettnau.restr.data
# 
# # 4) Extract network position ---------------------------------------------
# 
# # 4.1. Load network data --------------------------------------------------
# 
# load("gmm.Mettnau.all.Rdata")
# load("gmm.McDonalds.all.Rdata")
# 
# 
# # 4.2. Restrict the network data to the last 4 weeks before the restriction --------
# # dates are between 27.01. - 19.2.21
# 
# # Mettnau
# # create a network (SRI) for individuals with at least 10 sightings
# Mettnau.network <- get_network(association_data = gmm.Mettnau$gbi[which(gmm.Mettnau$metadata$Start>"210126000000"),which(colSums(gmm.Mettnau$gbi)>=10)], data_format = "GBI", association_index = "SRI")
# # restrict the network to great tits
# Mettnau.network <- Mettnau.network[rownames(Mettnau.network) %in% GT.list, colnames(Mettnau.network) %in% GT.list]
# # ensure that the birds with restrictions are part of the network
# restr.PIT.Mettnau %in% rownames(Mettnau.network)
# 
# # McDonalds
# # create a network (SRI) for individuals with at least 10 sightings
# McD.network <- get_network(association_data = gmm.McD$gbi[which(gmm.McD$metadata$Start>"210126000000"),which(colSums(gmm.McD$gbi)>=10)], data_format = "GBI", association_index = "SRI")
# # restrict the network to great tits
# McD.network <- McD.network[rownames(McD.network) %in% GT.list, colnames(McD.network) %in% GT.list]
# # ensure that the birds with restrictions are part of the network
# restr.PIT.McD %in% rownames(McD.network)
# 
# # 4.3. Extract the summed strength to knowledgeable individuals which were solving in the 10 days prior --------
# Mettnau.summed.strength <- NULL
# 
# for(i in Mettnau.restr.data$PIT){
#   summed.strengt.w.informed.ind <- sum(Mettnau.network[rownames(Mettnau.network)==i, colnames(Mettnau.network) %in% restr.PIT.Mettnau])
#   Mettnau.summed.strength[which(Mettnau.restr.data$PIT==i)] <- summed.strengt.w.informed.ind
# }
# 
# 
# McD.summed.strength <- NULL
# 
# for(i in McDonalds.restr.data$PIT){
#   summed.strengt.w.informed.ind <- sum(McD.network[rownames(McD.network)==i, colnames(McD.network) %in% restr.PIT.McD])
#   McD.summed.strength[which(McDonalds.restr.data$PIT==i)] <- summed.strengt.w.informed.ind
# }
# 
# # add this to the data frame with the individual level variables
# 
# Mettnau.restr.data$summed.strength <- Mettnau.summed.strength
# Mettnau.restr.data$site <- rep("Mettnau", length(Mettnau.restr.data[,1]))
# McDonalds.restr.data$summed.strength <- McD.summed.strength
# McDonalds.restr.data$site <- rep("McDonalds", length(McDonalds.restr.data[,1]))
# 
# # combine data from both sites into one data frame
# restr.data.combined <- rbind(Mettnau.restr.data, McDonalds.restr.data)
# 
# head(restr.data.combined)
# # same data frame as before, but with two additional columns (summed strength and site)
# 
# # 4.4. Add age and species to the data frame ------------------------------
# 
# age.vec <- NULL
# species.vec <- NULL
# 
# for(i in restr.data.combined$PIT){
#   age <- unique(subset(species.list$Age.in.2021, species.list$Pit==i))
#   species <- unique(subset(species.list$Species, species.list$Pit==i))
# 
#   age.vec[which(restr.data.combined$PIT==i)] <- age
#   species.vec[which(restr.data.combined$PIT==i)] <- species
# 
# }
# 
# restr.data.combined$age <- age.vec
# restr.data.combined$species <- species.vec
# head(restr.data.combined)
# length(restr.data.combined[,1])

# write.table(restr.data.combined, file="latencies.ILVs.txt")

# this contains the extracted data in a ready-to-go data frame
restr.data.combined <- read.delim("latencies.ILVs.txt", sep=" ")

# 5) Analysis -------------------------------------------------------------
# we restrict analyses to birds that have solved at least 5 times during the 10 days prior to restriction
# as calculating solving rates for low frequency solvers is very noisy
restr.data.combined.sub <- subset(restr.data.combined, as.numeric(restr.data.combined$solves_10_days_prior)>=5)

# 5.1. Summary on switching behaviour  -----------------------------------

# How many managed to switch to the other side
paste("Out of " , length(restr.data.combined.sub$PIT), " birds, ", length(restr.data.combined.sub$PIT)-length(restr.data.combined.sub$PIT[restr.data.combined.sub$solves_after==0]), " managed to switch to the other side after the restriction was put in place.", sep="")

# How long was the latency to the first solve on average?

mean(na.omit(as.numeric(restr.data.combined.sub$latency.first.solve)))
min(na.omit(as.numeric(restr.data.combined.sub$latency.first.solve)))
max(na.omit(as.numeric(restr.data.combined.sub$latency.first.solve)))

# on average, they produced their first solve on their 12th visit (range 1-43)

# How many reached their prior solving rate [-10%]

paste("Out of ", length(restr.data.combined.sub$PIT), " birds, ", length(restr.data.combined.sub[restr.data.combined.sub$censored.rate==1,"PIT"]), " reached their solving rate from prior.", sep="")

# How long did it take them to reach the solving rate
mean(na.omit(as.numeric(restr.data.combined.sub$latency.rate.reached)))
min(na.omit(as.numeric(restr.data.combined.sub$latency.rate.reached)))
max(na.omit(as.numeric(restr.data.combined.sub$latency.rate.reached)))


# for those that managed to switch:
# Was the solve rate before and after significantly different?


solve.rate <- cbind.data.frame(subset(restr.data.combined.sub$PIT, restr.data.combined.sub$censored.rate!=0), 
                               as.numeric(subset(restr.data.combined.sub$solve.rate.prior, restr.data.combined.sub$censored.rate!=0)),
       as.numeric(subset(restr.data.combined.sub$solve.rate.after.rate.reached, restr.data.combined.sub$censored.rate!=0)))
colnames(solve.rate) <- c("PIT", "before", "after.switch")
solve.rate

# check for normality
shapiro.test(solve.rate$before)
shapiro.test(solve.rate$after.switch)

qqplot(solve.rate$before, solve.rate$after.switch)

# data does not appear to be normally distributed
# we use a non-parametric test:

wilcox.test(as.numeric(solve.rate$before), as.numeric(solve.rate$after.switch), paired=TRUE)

# Wilcoxon signed rank test with continuity correction
# 
# data:  solve.rate$before and solve.rate$after.switch
# V = 66.5, p-value = 0.4204
# alternative hypothesis: true location shift is not equal to 0

# Create a figure that shows the solving rate before and after

solve.rate.to.plot <- cbind.data.frame(c(solve.rate$before, solve.rate$after.switch),
c(rep("prior", length(solve.rate$before)), rep("after", length(solve.rate$after.switch))),
rep(solve.rate$PIT, 2))

colnames(solve.rate.to.plot) <- c("solve.rate", "event", "PIT")

solve.rate.to.plot$event <- factor(solve.rate.to.plot$event,
                       levels = c('prior','after'),ordered = TRUE)


ggplot(solve.rate.to.plot, aes(x=event, y=as.numeric(solve.rate), fill=event))+
  geom_boxplot(fill=c('#74f8e5', "#a4d4af"))+
  geom_point()+
  theme_classic()+
  geom_line(aes(group=PIT))+
  ylab("solving rate")+
  xlab("restriction")+
  theme(legend.position = "none")
  
  

# 5.2 Check for multicollinearity between covariates--------
# using variance inflation factors (VIFs)

# we include all birds - whether they managed to switch or not
restr.data.combined.for.cox <- subset(restr.data.combined.sub, restr.data.combined.sub$solves_after>=0)
str(restr.data.combined.for.cox)

# assign correct classes to variables
restr.data.combined.for.cox$pref_strength <- as.numeric(restr.data.combined.for.cox$pref_strength)
restr.data.combined.for.cox$latency.first.solve <- as.numeric(restr.data.combined.for.cox$latency.first.solve)
restr.data.combined.for.cox$latency.rate.reached <- as.numeric(restr.data.combined.for.cox$latency.rate.reached)
restr.data.combined.for.cox$solve.rate.prior <- as.numeric(restr.data.combined.for.cox$solve.rate.prior)
restr.data.combined.for.cox$summed.strength <- as.numeric(restr.data.combined.for.cox$summed.strength)
restr.data.combined.for.cox$solves_10_days_prior <- as.numeric(restr.data.combined.for.cox$solves_10_days_prior)
length(restr.data.combined.for.cox$PIT)


# check for multi-collinearity among variables
fit.multi.collinearty <-
  lm(
    latency.first.solve ~ solve.rate.prior + summed.strength + age + tutor +
      pref_strength + solves_10_days_prior,
    data = restr.data.combined.for.cox
  )
summary(fit.multi.collinearty)
vif(fit.multi.collinearty)

# solve.rate.prior      summed.strength        age                tutor        pref_strength        solves_10_days_prior 
# 4.098490             1.373441             1.242394             2.246181             1.482034             3.289813 

# we get an inflated VIF for the solve rate prior (4.1) - hence we drop it and recalculate the VIFs
fit.multi.collinearty <-
  lm(
    latency.first.solve ~  summed.strength + age + tutor +
      pref_strength + solves_10_days_prior,
    data = restr.data.combined.for.cox
  )
summary(fit.multi.collinearty)
vif(fit.multi.collinearty)
# summed.strength                  age                tutor        pref_strength solves_10_days_prior 
# 1.244140             1.159382             1.523147             1.464692             1.359561 

# this reduces VIFs to below 2.5


# 5.3. Latencies - cox proportional hazard models -------------------------

# 5.3.1. Latency until first solve -------------------------------------------

# create a 'survival' object with the time it took to produce the first solution and a vector stating whether birds managed to produce a solution or not
survival.latency.first <-
  Surv(
    time = restr.data.combined.for.cox$latency.first.solve,
    event = as.numeric(restr.data.combined.for.cox$censored.first.solve)
  )
survival.latency.first
# this shows the number of visits it took until the bird produced its first solution on the less preferred side 
# the + next to the number means that it never managed to produce a solution until the experiment ended (censored date)
# the number in this case indicates the total number of visits after the restriction until the end of the experiment


# we test for the assumptions for a cox prop hazard model (ignoring the random effect (site) for now)

test.ph.latency.first <-   coxph(
  survival.latency.first ~   summed.strength + age + tutor +
    pref_strength + solves_10_days_prior,
  data = restr.data.combined.for.cox
)

out.test.ph.latency.first <-   cox.zph(
test.ph.latency.first
)

out.test.ph.latency.first
# 
# chisq df    p
# summed.strength      0.1316  1 0.72
# age                  0.1361  1 0.71
# tutor                1.5333  1 0.22
# pref_strength        0.0231  1 0.88
# solves_10_days_prior 0.2825  1 0.60
# GLOBAL               4.4000  5 0.49

# The proportional hazard assumption is supported by a non-significant relationship between residuals and time, and refuted by a significant relationship
# here, we have no significant values, indicating the the hazard assumption is supported and we can use the cox models for our data


# we now run the cox proportional hazard model including the random effect (for site)
# we first run the gloabl model (including all explanatory variables)
# we include an interaction term for preferred strength and the number of solves in the 10 days prior to the restriction
m.first.solve <-
  coxme(
    survival.latency.first ~   summed.strength + age + tutor +
      pref_strength + solves_10_days_prior + (1|site),
    data = restr.data.combined.for.cox
  )


print(m.first.solve)
# p -values show whether Beta (the slope) is significantly different from 0 (with a value <0.05) - not the case for the global model
# coef are the regression coefficients. A positive value means that the hazard is higher - meaning that solves are likely to occur sooner
# a negative value means that the hazard is lower, meaning that the solves occur later
# exp(coef) are the hazard ratios - they give the effect of the covariates. 
# For example, a hazard rate of 0.72 (for tutor) indicates that tutors were slower at switching (negative value) than non-tutors by a factor of 0.72

# we can use the 'dredge' function to run all candidate models
options(na.action = "na.fail")
dredge.m.first.solve <- dredge(m.first.solve)
dredge.m.first.solve # contains all possible models sorted according to Akaike weights
# we can see that the top model includes the preference strength and the solve rate prior
# but it is closely followed by alternative models
options(na.action = "na.omit")

# we can visualize the cumulative Akaike weights for covariates
plot(dredge.m.first.solve)

# to calculate support for covariates, we can get the summed Akaike weights
sw(subset(model.sel(dredge.m.first.solve)))

#                     solves_10_days_prior pref_strength summed.strength tutor age 
# Sum of weights:      0.91                 0.42          0.34            0.34  0.21
# N containing models:   16                   16            16              16    16

# we consider weights with more support than against (>0.5) as important

# extract model average estimate based on models with deltaAIc of <4 (subset - as 4 is default)
mle.first.solve <- model.avg(dredge.m.first.solve, subset = delta <= 4)
mle.first.solve
# Coefficients: 
#          pref_strength solves_10_days_prior   tutoryes summed.strength agefirst.year
# full       -1.145801          0.003052886 -0.2125277      -0.2642769  -0.006703916
# subset     -2.863728          0.003052886 -0.6651375      -0.8534901  -0.039290190

# we can extract the hazard ratios by exponentiation from the 'subset'
# which only considers models within a deltaAICc of 4
exp(mle.first.solve$coefficients[2,])

# get 95% CI intervals, again considering models within a deltaAICc of 4
confset.95p <- get.models(dredge.m.first.solve, subset = delta < 4)
avgmod.95p <- model.avg(confset.95p)
summary(avgmod.95p)
c.first <- confint(avgmod.95p)
# we get lower and upper bounds of the hazard ratios by exponentiation
exp(c.first[,1])
exp(c.first[,2])


# 5.3.2 Latency till rate from prior reached ------------------------------

### we repeat this for the latency to reach their prior rate

# create a 'survival' object with the time it took to produce the first solution and a vector stating whether birds managed to produce a solution or not
survival.latency.rate <-
  Surv(
    time = restr.data.combined.for.cox$latency.rate.reached,
    event = as.numeric(restr.data.combined.for.cox$censored.rate)
  )
survival.latency.rate
# this shows the number of visits it took until the bird reached the solving rate from prior
# the + next to the number means that it never managed to reach the prior rate (censored date)
# the number in this case indicates the total number of visits after the restriction until the end of the experiment


# we test for the assumptions for a cox prop hazard model (ignoring the random effect (site) for now)

test.ph.latency.rate <-   coxph(
  survival.latency.rate ~   summed.strength + age + tutor +
    pref_strength + solves_10_days_prior ,
  data = restr.data.combined.for.cox
)

out.test.ph.latency.rate <-   cox.zph(
  test.ph.latency.rate
)

out.test.ph.latency.rate

#                         chisq df     p
# summed.strength      2.73e-05  1 0.996
# age                  1.25e-01  1 0.724
# tutor                1.83e+00  1 0.176
# pref_strength        3.27e+00  1 0.071
# solves_10_days_prior 2.89e+00  1 0.089
# GLOBAL               1.02e+01  5 0.071

# again, no p-values are significant, meaning that the assumptions are met

# we now run the cox proportional hazard model including the random effect (for site)
# we rate run the global model (including all explanatory variables)
# we include an interaction term for preferred strength and the number of solves in the 10 days prior to the restriction
m.rate.solve <-
  coxme(
    survival.latency.rate ~  summed.strength  + solves_10_days_prior + pref_strength  + tutor + age + (1|site),
    data = restr.data.combined.for.cox
  )


print(m.rate.solve)
# p -values show whether Beta (the slope) is significantly different from 0 (with a value <0.05) - not the case for the global model
# coef are the regression coefficients. A positive value means that the hazard is higher - meaning that solves are likely to occur sooner
# a negative value means that the hazard is lower, meaning that the solves occur later
# exp(coef) are the hazard ratios - they give the effect of the covariates. 

# we can use the 'dredge' function to run all candidate models
options(na.action = "na.fail")
dredge.m.rate.solve <- dredge(m.rate.solve)
dredge.m.rate.solve # contains all possible models sorted according to Akaike weights
# we can see that the top model includes the preference strength and the number of solves 10 days prior
options(na.action = "na.omit")

# we can visualize the cumulative Akaike weights for covariates
plot(dredge.m.rate.solve)

# to calculate support for covariates, we can get the summed Akaike weights
sw(subset(model.sel(dredge.m.rate.solve)))

# pref_strength summed.strength tutor age  solves_10_days_prior
# Sum of weights:      0.67          0.30            0.28  0.27 0.24                
# N containing models:   16            16              16    16   16                      

# extract model average estimate based on models with deltaAIc of <4 (subset - as 4 is default)
mle.rate.solve <- model.avg(dredge.m.rate.solve, subset = delta <= 4)
mle.rate.solve
# Coefficients: 
#            pref_strength agefirst.year summed.strength    tutoryes solves_10_days_prior
# full       -3.388433    -0.1024626      -0.2857158 -0.09634687         0.0001062727
# subset     -4.540927    -0.5052755      -1.0823228 -0.62101714         0.0008444579
exp(mle.rate.solve$coefficients[2,])


# get 95% CI
confset.95p <- get.models(dredge.m.rate.solve, subset = delta < 4)
avgmod.95p <- model.avg(confset.95p)
summary(avgmod.95p)
mle.rate.c <- confint(avgmod.95p)
exp(mle.rate.c[,1])
exp(mle.rate.c[,2])

# 5.4. Plotting hazard ratios ---------------------------------------------

# plot hazard ratios for each covariate (model averaged estimate) with the 95% CI.

cov.first <- as.vector(colnames(mle.first.solve$coefficients))
mle.first <- as.vector(mle.first.solve$coefficients[2,])
hazard.first <- exp(mle.first)
confint.lower.first <- as.vector(exp(confint(mle.first.solve)[,1]))
confint.upper.first <- as.vector(exp(confint(mle.first.solve)[,2]))

# combine coefficients and confidence intervals into data frame
# make reference levels for the binary variables
comb.first <- cbind.data.frame(cov.first, mle.first, hazard.first, confint.lower.first, confint.upper.first)
comb.first <- rbind.data.frame(comb.first, c("tutor_ref", NA, 1+(1-as.numeric(comb.first$hazard.first[which(comb.first$cov.first=="tutoryes")])), NA, NA))
comb.first <- rbind.data.frame(comb.first, c("age_ref", NA, 1+(1-as.numeric(comb.first$hazard.first[which(comb.first$cov.first=="agefirst.year")])), NA, NA))

class(comb.first)
comb.first$hazard.first <- as.numeric(comb.first$hazard.first)
comb.first$confint.lower.first <- as.numeric(comb.first$confint.lower.first)
comb.first$confint.upper.first <- as.numeric(comb.first$confint.upper.first)
comb.first$names.first <-
  c(
    "Strength of side preference",
    "Num. solves prior",
    "Tutor (yes)",
    "Summed strength with informed inds.",
    "Age (first year)",
    "Tutor (no - reference)",
    "Age (adults - reference)"
  )

# reorder:
comb.first <- comb.first[order(comb.first$cov.first, decreasing=TRUE),]
# make sure it keeps the order
comb.first$cov.first <- factor(comb.first$cov.first, levels = comb.first$cov.first)

str(comb.first)
first.p <- ggplot(data=comb.first, aes(x=as.numeric(hazard.first), y=cov.first))+
  geom_point(shape=15, size=4)+
  scale_x_log10(limits = c(0.0001,10), breaks=c(0.001, 0.01, 0.1, 1, 10), labels=c(0.001, 0.01, 0.1, 1, 10))+
  geom_vline(xintercept = 1, linetype="dashed")+
  geom_errorbarh(aes(xmax=confint.upper.first, xmin=confint.lower.first), height=0.3, na.rm=TRUE)+
  scale_y_discrete(labels=comb.first$names.first)+
  xlab("Hazard ratio (log scale)")+
  ylab(NULL)+
  ggtitle("A) Latency until first solve")+
  theme(text = element_text(size=14))+
  theme(plot.title = element_text(hjust = 0))

## now plot the same for the latency till the rate was reached

cov.rate <- as.vector(colnames(mle.rate.solve$coefficients))
mle.rate <- as.vector(mle.rate.solve$coefficients[2,])
hazard.rate <- exp(mle.rate)
confint.lower.rate <- as.vector(exp(confint(mle.rate.solve)[,1]))
confint.upper.rate <- as.vector(exp(confint(mle.rate.solve)[,2]))


comb.rate <- cbind.data.frame(cov.rate, mle.rate, hazard.rate, confint.lower.rate, confint.upper.rate)
comb.rate <- rbind.data.frame(comb.rate, c("tutor_ref", NA, 1+(1-as.numeric(comb.rate$hazard.rate[which(comb.rate$cov.rate=="tutoryes")])), NA, NA))
comb.rate <- rbind.data.frame(comb.rate, c("age_ref", NA, 1+(1-as.numeric(comb.rate$hazard.rate[which(comb.rate$cov.rate=="agefirst.year")])), NA, NA))

class(comb.rate)
comb.rate$hazard.rate <- as.numeric(comb.rate$hazard.rate)
comb.rate$confint.lower.rate <- as.numeric(comb.rate$confint.lower.rate)
comb.rate$confint.upper.rate <- as.numeric(comb.rate$confint.upper.rate)
comb.rate$names.rate <-
  c(
    "Strength of side preference",
    "Age (first year)",
    "Summed strength with informed inds.",
    "Tutor (yes)",
      "Num. solves prior",
    "Tutor (no - reference)",
    "Age (adults - reference)"
  )



# reorder:
comb.rate <- comb.rate[order(comb.rate$cov.rate, decreasing=TRUE),]
# make sure it keeps the order
comb.rate$cov.rate <- factor(comb.rate$cov.rate, levels = comb.rate$cov.rate)

str(comb.rate)
rate.p <- ggplot(data=comb.rate, aes(x=as.numeric(hazard.rate), y=cov.rate))+
  geom_point(shape=15, size=4)+
  scale_x_log10(limits = c(0.0001,10), breaks=c(0.001, 0.01, 0.1, 1, 10), labels=c(0.001, 0.01, 0.1, 1, 10))+
  geom_vline(xintercept = 1, linetype="dashed")+
  geom_errorbarh(aes(xmax=confint.upper.rate, xmin=confint.lower.rate), height=0.3, na.rm=TRUE)+
  scale_y_discrete(labels=NULL)+
  xlab("Hazard ratio (log scale)")+
  ylab(NULL)+
  ggtitle("B) Latency until prior rate reached")+
  theme(text = element_text(size=14))+
  theme(plot.title = element_text(hjust = 0, size=14))

require(gridExtra)

png(filename="hazard.ratios.png", width = 800, height = 400)

ggarrange(first.p, rate.p, ncol=2)


dev.off()

