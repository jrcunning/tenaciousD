# LOAD DATA AND LIBRARIES -----
library(lme4)
library(scales)
library(lsmeans)
library(reshape2)
library(plyr)
library(lattice)
library(mgcv)
library(MASS)

# Define function to use for plotting confidence intervals
addpoly <- function(x, y1, y2, col=alpha("lightgrey", 0.8), ...){
  ii <- order(x)
  y1 <- y1[ii]
  y2 <- y2[ii]
  x <- x[ii]
  polygon(c(x, rev(x)), c(y1, rev(y2)), col=col, border=NA, ...)
}

# Import data
data <- read.csv("tenaciousD_data.csv")

# Adjust and transform data -----
# Create factor version of time
data$timef <- factor(data$time)
# Calculate proportion D in each sample
data$propD <- data$D.SH / (data$C.SH + data$D.SH)
# Categorize C- and D-dominated corals based on community composition at time zero
dom <- with(data[data$time==0, ], na.omit(data.frame(sample=sample, 
                                                     dom=ifelse(is.na(propD), NA, ifelse(propD > 0.5, "D", "C")))))
# Merge dominant symbiont classification with rest of data
data <- merge(data, dom, by="sample", all.x=T)
table(data[data$time==0, "dom"])  # 27 C-dominant and 88 D-dominant corals = 115 corals with qPCR data at t0
data <- data[with(data, order(sample, time)), ]
# Keep only those corals with data at time zero
data <- data[!is.na(data$dom),]
# Replace zeros with detection limits (just below minimum detected value)
table(data$C.SH==0) # 23% of samples had no detectable clade C
table(data$D.SH==0) # 11% of samples had no detectable clade D
min(data[data$C.SH!=0,"C.SH"], na.rm=T) # detection limit for C is ~1e-6
min(data[data$D.SH!=0,"D.SH"], na.rm=T) # detection limit for D is ~1e-4
data[data$C.SH==0 & !is.na(data$C.SH), "C.SH"] <- 1e-6
data[data$D.SH==0 & !is.na(data$D.SH), "D.SH"] <- 1e-4
data$tot.SH <- data$C.SH + data$D.SH
data[data$tot.SH<0.000101 & !is.na(data$tot.SH), "tot.SH"] <- 0.000101
# Calculate relative change in Fv/Fm
for (sample in data$sample) {
  sdata <- data[which(data$sample==sample), ]
  for (time in sdata$time) {
    sdata[which(sdata$time==time), "rfvfm"] <- sdata[which(sdata$time==time), "fvfm"] / sdata[which(sdata$time==0), "fvfm"]
  }
  data[which(data$sample==sample), "rfvfm"] <- sdata[, "rfvfm"]
}


# Prepare raw data for plotting
datsumm <- rbind(expand.grid(time=seq(0,63,7), ramp=factor("cool"), dom=factor(c("C","D"))),
                 expand.grid(time=seq(0,42,7), ramp=factor("heat"), dom=factor(c("C","D"))))
datsumm$rfvfm.mean <- aggregate(data$rfvfm, by=list(interaction(data$time, data$dom, data$ramp)), FUN=mean, na.rm=T)$x
datsumm$rfvfm.sd <- aggregate(data$rfvfm, by=list(interaction(data$time, data$dom, data$ramp)), FUN=sd, na.rm=T)$x
datsumm$logtot.SH.mean <- aggregate(log10(data$tot.SH), by=list(interaction(data$time, data$dom, data$ramp)), FUN=mean, na.rm=T)$x
datsumm$logtot.SH.sd <- aggregate(log10(data$tot.SH), by=list(interaction(data$time, data$dom, data$ramp)), FUN=sd, na.rm=T)$x
datlist <- split(datsumm, f=datsumm$ramp)
datlist <- lapply(datlist, function(x) rev(split(x, f=x$dom)))