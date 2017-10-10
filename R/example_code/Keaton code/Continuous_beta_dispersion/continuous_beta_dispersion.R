# continuous_beta_dispersion.R

####################################### DESCRIPTION #######################################
# The point of this script is to take a continuous variable, in this case, Factor1.min0 and
# bin it into a specified interval (0.05 below), and then calculate the mean community
# distance for each person to all other people within that bin. As examples I have provided
# a distance matrix (shuar.dist) and sample metadata (sol.factors.data)
###########################################################################################

# Set path to directory below and un-comment line
# setwd(".../Continuous_beta_dispersion/")

# Required libraries, sources and objects
library(doBy)
library(ggplot2)
theme_set(theme_bw())
source("generate_distance_data.R")
load("shuar_dist_obj", verbose=TRUE)
load("sol_factors_data_obj", verbose=TRUE)

interval <- 0.05 # Set the interval you want to bin the continuous variable by

# The funciton below takes a distance matrix and converts from the matrix format to a data
# frame with the following structure:

# smpl1 | smpl2 | distance
# -------------------------
#     A |     B |    0.543
#     A |     C |    0.023
#     B |     A |    0.543
#     B |     C |    0.172
#     C |     A |    0.023
#     C |     A |    0.172

# NOTICE there are redundant lines A <=> B and B <=> A. This is necessary to do the proper
# computation for each sample.
# The script will also take metadata from a your sample metadata file and add those
# comparisions (for categorical data) or differences in measure (for continuous data) to
# the data.frame, thus:

# smpl1 | smpl2 | distance | Com_rename.comp | Factor1.min0
# ----------------------------------------------------------
#     A |     B |    0.543 |        CC1<>CC2 |        0.012
#     A |     C |    0.023 |        CC1<>CC3 |        0.287
#     B |     A |    0.543 |        CC1<>CC2 |        0.012
#     B |     C |    0.172 |        CC2<>CC3 |        0.325
#     C |     A |    0.023 |        CC1<>CC3 |        0.287
#     C |     B |    0.172 |        CC2<>CC3 |        0.325

shuar.dist.data <- generate_distance_data(dist.obj=shuar.dist,
                                          smpl.data=sol.factors.data,
                                          measures.of.interest=c("Com_rename",
                                                                 "time.to.sucua",
                                                                 "Factor1.min0",
                                                                 "Factor2.min0",
                                                                 "Factor3.min0"),
                                          verbose=TRUE)

# If this script takes a long time to run, it may be easier to save the object so you
# don't have to keep re-running it.
save(shuar.dist.data, file="shuar_dist_data_obj")

# Calculate the bin size for the measure of interest
fact1Diff.interval <- interval * max(shuar.dist.data$Factor1.min0.diff)

# Subset only those comparisons that fall within the interval
fact1.diffBin.interval <- subset(shuar.dist.data, Factor1.min0.diff <= fact1Diff.interval)

# The steps below make sure all the comparisons for each sample will be made
fact1.diffBin.interval1 <- fact1.diffBin.interval[, -2]
names(fact1.diffBin.interval1)[1] <- "smpl"
fact1.diffBin.interval2 <- fact1.diffBin.interval[, -1]
names(fact1.diffBin.interval2)[1] <- "smpl"
fact1.diffBin <- rbind(fact1.diffBin.interval1, fact1.diffBin.interval2)

# This step finds the mean (and sd and number) of the community distance for each sample
# to all the other within the interval
fact1.diffBin.sum <- summaryBy(distance ~ smpl,
                               data=fact1.diffBin,
                               FUN=c(mean, sd, length))

# Merge the distance means with the original sample metadata
fact1.diffBin.data <- merge(sol.factors.data, fact1.diffBin.sum, by.x="row.names", by.y=1)

# Linear model testing
fact1.betadisp.lm <- lm(distance.mean ~ Factor1.min0, data=fact1.diffBin.data)
(fact1.betadisp.lm.sum <- summary(fact1.betadisp.lm))

# Plot results
tt.cols <- c("black", "gray50", "gray65", "gray80") # colors for plotting

fact1.diffBin.betadisp <- ggplot(fact1.diffBin.data, aes(x=Factor1.min0, y=distance.mean)) +
    geom_point(aes(fill=time.to.sucua), shape=21, color='white', size=1.2) +
    stat_smooth(method="lm", size=0.5, color='black', alpha=0.5) +
    scale_fill_manual(values=tt.cols,
                      guide=FALSE) +
    labs(x="House Modernity", y="β-dispersion") +
    theme(plot.background=element_blank(),
          plot.margin=unit(c(4,4,4,4), "pt"),
          panel.grid=element_blank())
fact1.diffBin.betadisp
