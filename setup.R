# LOAD DATA AND LIBRARIES -----
library(lme4)
library(scales)
library(lsmeans)
library(reshape2)
library(plyr)
library(lattice)
library(mgcv)
library(MASS)
library(pander)

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
data$prevBleach <- ifelse(data$history %in% c("24A'", "B'", "A'"), "NB", "B")

# Adjust and transform data -----
# Create factor version of time
data$timef <- factor(data$time)
#Based on the relative abundance of clades C and D at time zero, each coral was categorized as 
#either initially clade C-dominated or clade D-dominated (groups hereafter referred to as C corals 
#and D corals, respectively). For samples missing qPCR data at time zero, data from the previous 
#time point (Silverstein et al. 2015) were substituted in the case of non-bleached corals. For 
#bleached corals, no data were substituted but cores were categorized as C- or D-dominated based 
#on data from the subsequent time point. For samples in which no clade C or D was detected, the 
#S/H ratio for that clade was set to a value just below the detection threshold, defined as the 
#minimum S/H ratio detected for that clade across the entire dataset (1e-6 for clade C, 1e-4 for 
#clade D).
# Calculate proportion D in each sample
data$propD <- data$D.SH / (data$C.SH + data$D.SH)
# Categorize C- and D-dominated corals based on community composition at time zero
dom <- with(data[data$time==0, ], na.omit(data.frame(sample=sample, 
                                                     dom=ifelse(is.na(propD), NA, ifelse(propD > 0.5, "D", "C")))))

# Count number of cores with mixed communities
syms <- aggregate(data.frame(C=data$C.SH, D=data$D.SH), by=list(core=data$sample), FUN=mean, na.rm=T)
symstab <- addmargins(table(syms$C!=0, syms$D!=0)) # True if ever contained C (rows) or D (columns)
dimnames(symstab) <- list(C_detected=c("no", "yes", "Sum"), D_detected=c("no", "yes", "Sum"))
ftable(symstab)  # 129/158 had both detected at least once (81.6%)

# Merge dominant symbiont classification with rest of data
data <- merge(data, dom, by="sample", all.x=T)
table(data[data$time==0, "dom"])  # 27 C-dominant and 88 D-dominant corals = 115 corals with qPCR data at t0
data <- data[with(data, order(sample, time)), ]
# Assign dominant symbiont to corals missing data for t0
# 2_24 is missing data but was never heated, so assign C dominant
data[data$sample=="2_24", "dom"] <- "C"
# All other corals without data for t0 were subsequently D dominant at other times, so assign D dominance
data$dom[is.na(data$dom)] <- "D"
# Replace zeros with detection limits (just below minimum detected value)
table(data$C.SH==0) # 22% of samples had no detectable clade C
table(data$D.SH==0) # 16% of samples had no detectable clade D
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
datsumm$fvfm.mean <- aggregate(data$fvfm, by=list(interaction(data$time, data$dom, data$ramp)), FUN=mean, na.rm=T)$x
datsumm$fvfm.sd <- aggregate(data$fvfm, by=list(interaction(data$time, data$dom, data$ramp)), FUN=sd, na.rm=T)$x
datsumm$rfvfm.mean <- aggregate(data$rfvfm, by=list(interaction(data$time, data$dom, data$ramp)), FUN=mean, na.rm=T)$x
datsumm$rfvfm.sd <- aggregate(data$rfvfm, by=list(interaction(data$time, data$dom, data$ramp)), FUN=sd, na.rm=T)$x
datsumm$logtot.SH.mean <- aggregate(log10(data$tot.SH), by=list(interaction(data$time, data$dom, data$ramp)), FUN=mean, na.rm=T)$x
datsumm$logtot.SH.sd <- aggregate(log10(data$tot.SH), by=list(interaction(data$time, data$dom, data$ramp)), FUN=sd, na.rm=T)$x
datlist <- split(datsumm, f=datsumm$ramp)
datlist <- lapply(datlist, function(x) rev(split(x, f=x$dom)))