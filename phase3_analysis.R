# LOAD DATA AND LIBRARIES -----
library(effects)
library(lme4)
library(spida)
library(scales)
library(lsmeans)
library(reshape2)
library(plyr)
library(lattice)

addpoly <- function(x, y1, y2, col=alpha("lightgrey", 0.8), ...){
  ii <- order(x)
  y1 <- y1[ii]
  y2 <- y2[ii]
  x <- x[ii]
  polygon(c(x, rev(x)), c(y1, rev(y2)), col=col, border=NA, ...)
}

data <- read.csv("phase3.csv")

# Calculate proportion D in each sample -----
data$propD <- data$D.SH / (data$C.SH + data$D.SH)
hist(data$propD)

# Remove zeros and ones and recalc proportions
data[is.na(data$C.SH), "C.SH"] <- 1e-6
data[data$C.SH==0, "C.SH"] <- 1e-6
data[is.na(data$D.SH), "D.SH"] <- 1e-6
data[data$D.SH==0, "D.SH"] <- 1e-6
data$propD2 <- data$D.SH / (data$C.SH + data$D.SH)
data[is.na(data$tot.SH), "propD2"] <- NA
data[1:20,]
na.omit(data)$propD2

# Plot changes in proportion D over time
xyplot(propD ~ time | ramp, groups= ~ sample, data = na.omit(data), type="o", lty=1)
xyplot(asin(sqrt(propD)) ~ time | ramp, groups= ~ sample, data = na.omit(data), type="o", lty=1)
#xyplot(propD ~ time | ramp + mother, groups= ~ sample, data = na.omit(data), type="o", lty=1)

# Categorize C- and D-dominated corals based on community composition at time zero
dom <- with(data[data$time==0, ], na.omit(data.frame(sample=sample, dom=ifelse(propD > 0.5, "D", "C"))))
# Some corals are missing qpcr data for time zero; for these, assign based on mean propD of all samples from each unassigned coral
dom2 <- with(data[!data$sample %in% dom$sample, ], na.omit(aggregate(propD, by=list(sample=sample), FUN=mean, na.rm=T)))
dom2$dom <- ifelse(dom2$x > 0.5, "D", "C")
dom <- merge(dom, dom2, all=T)[,-3]
dom

# Categorize C- and D-dominant corals based on >90% one or the other at t0, omitting mixtures
nlevels(data$sample)
dom <- with(data[data$time==0, ], 
            data.frame(sample=sample, propD=propD, 
                       dom=ifelse(propD > 0.9, "D", ifelse(propD < 0.1, "C", ifelse(is.na(propD), NA, "mix")))))
dom <- with(data[data$time==0, ], 
            data.frame(sample=sample, propD=propD, 
                       dom=ifelse(propD > 0.9, "D", ifelse(propD < 0.1, "C", NA))))

nlevels(droplevels(dom[!is.na(dom$dom), "sample"]))
# 55 coral samples are not assigned as C or D corals.... OK?

# Merge dominant symbiont classification with rest of data
data <- merge(data, dom, by="sample", all.x=T)
table(data[data$time==0, "dom"])  # 21 C-dominant and 83 D-dominant corals = 104 corals with qPCR data
data <- data[with(data, order(sample, time)), ]
nlevels(droplevels(data$sample))

# Calculate relative change in Fv/Fm and difference in Fv/Fm (rfvfm, dfvfm)
for (sample in data$sample) {
  sdata <- data[which(data$sample==sample), ]
  for (time in sdata$time) {
    sdata[which(sdata$time==time), "rfvfm"] <- sdata[which(sdata$time==time), "fvfm"] / sdata[which(sdata$time==0), "fvfm"]
    sdata[which(sdata$time==time), "dfvfm"] <- sdata[which(sdata$time==time), "fvfm"] - sdata[which(sdata$time==0), "fvfm"]
    }
  data[which(data$sample==sample), "rfvfm"] <- sdata[, "rfvfm"]
  data[which(data$sample==sample), "dfvfm"] <- sdata[, "dfvfm"]
}
head(data, 50)

# Fv/Fm ANALYSIS ---------
# Visualize all data
xyplot(rfvfm ~ time | ramp + dom, groups= ~ sample, data = data, type="o", lty=1)

# Fit mixed model
mod.all <- lmerTest::lmer(rfvfm ~ poly(time, 2) * ramp * dom + (1|mother/sample), data=data)
anovatab <- lmerTest::anova(mod.all)
anovatab

# Test C vs. D differences at each date using model fit
rg <- ref.grid(mod.all, at=list(time=c(0,7,14,21,28,35,42,49,56,63)))
lsm <- lsmeans(rg, specs=pairwise ~ dom | ramp * time)  # CORRECT FOR MULTIPLE TESTING?
pvals <- summary(rbind(contrast(lsm, "tukey", adjust="mvt")))
fvfmsigs <- pvals[which(pvals$p.value < 0.05), ]
fvfmsigs
write.csv(round(anovatab, digits=3), file="output/Table1.csv")
# Pseudo-r2 value-- squared correlation between fitted and observed values
summary(lm(model.response(model.frame(mod.all)) ~ fitted(mod.all)))$r.squared
# Generate predictions and confidence intervals by parametric bootstrapping
pred.all <- expand.grid(time=seq(0,63,1), ramp=factor(c("cool", "heat")), dom=factor(c("C", "D")))
bootfit <- bootMer(mod.all, FUN=function(x) predict(x, pred.all, re.form=NA), nsim=100)
# Extract 95% confidence interval on predicted values
pred.all$fit <- predict(mod.all, pred.all, re.form=NA)
pred.all$lci <- apply(bootfit$t, 2, quantile, 0.025)
pred.all$uci <- apply(bootfit$t, 2, quantile, 0.975)
# Prepare data for plotting
datsumm <- data.frame(
  colsplit(as.character(levels(droplevels(interaction(data$dom, data$ramp, data$time)))),
           pattern="\\.", names=c("dom", "ramp", "time")),
  mean=aggregate(data$rfvfm, by=list(interaction(data$dom, data$ramp, data$time)), FUN=mean, na.rm=T)$x,
  sd=aggregate(data$rfvfm, by=list(interaction(data$dom, data$ramp, data$time)), FUN=sd, na.rm=T)$x,
  se=aggregate(data$rfvfm, by=list(interaction(data$dom, data$ramp, data$time)), 
               FUN=function(x) sd(x, na.rm=T)/sqrt(length(na.omit(x))))$x,
  conf95=aggregate(data$rfvfm, by=list(interaction(data$dom, data$ramp, data$time)), 
                   FUN=function(x) sd(x, na.rm=T)/sqrt(length(na.omit(x))) * qt(0.975, length(na.omit(x))-1))$x
  )
datlist <- split(datsumm, f=datsumm$ramp)
datlist <- lapply(datlist, function(x) rev(split(x, f=x$dom)))
predlist <- split(pred.all, f=pred.all$ramp)
predlist <- lapply(predlist, function(x) rev(split(x, f=x$dom))) 
predlist$heat$C <- predlist$heat$C[predlist$heat$C$time <= 42, ]
predlist$heat$D <- predlist$heat$D[predlist$heat$D$time <= 42, ]

# S/H ANALYSIS ------------
# analyse SH ratios with time as categorical, not continuous, because it's not really a good time series.
# i.e., not enough temporal resolution to be confident in quadratic effect of time....
# get SH data frame and plot raw data
shdf <- droplevels(subset(data, !is.na(tot.SH)))
shdf[shdf$tot.SH==0,"tot.SH"] <- NA
range(shdf$tot.SH, na.rm=T)
shdf[is.na(shdf$tot.SH), "tot.SH"] <- 2e-6
xyplot(log10(tot.SH) ~ time | ramp + dom, groups= ~ sample, data = shdf, type="o", lty=1)
# create list split by dom&ramp, fit mixed model for each group
shl <- split(shdf, f=interaction(shdf$ramp, shdf$dom))
mods <- llply(shl, function(df) lmerTest::lmer(log10(tot.SH) ~ factor(time) + (1|mother/sample), data=df))

# run anova and post-hoc comparisons over time for each group
anovas <- llply(mods, anova)
lsm <- llply(mods, function(mod) lsmeans(mod, specs="time", contr="dunnett"))
lsm
# plot lsm -- are these exactly same as geomeans of raw data?
# get contrasts for comparisons
contr <- ldply(sapply(lsm, "[[", 2), summary)
contr$bt <- 10^(contr$estimate)
contr$loss <- -(1-contr$bt)
shsigs <- contr[which(contr$p.value < 0.05), ]
shsigs
# get lsmeans from mixed model (slightly diff than raw data bc of random effects)
lsmeans <- ldply(sapply(lsm, "[[", 1), summary)
lsmeans
# get means from raw data
geomeans <- ldply(shl, function(df) aggregate(list(mean=log10(df$tot.SH)), by=list(time=df$time), FUN=mean))
geosds <- ldply(shl, function(df) aggregate(list(sd=log10(df$tot.SH)), by=list(time=df$time), FUN=sd))
geostats <- merge(geomeans, geosds, by=c(".id", "time"))
geostats <- split(geostats, f=geostats$.id)
geostats

# Multipanel figure -------------
  layout(mat=matrix(c(1,1,2,2,3,3,4,4), ncol=2))
  par(mgp=c(1.5,0.25,0), tck=-0.06)
  # Cooling Fv/Fm
  par(mar=c(1,3,4,1))
  plot(NA, xlim=c(0,63), ylim=c(0, 1), bty="n", tck=-0.03,
       xaxt="n", ann=F)
  title("A. Cooling", adj=0)
  mtext(side=2, text="Relative Fv/Fm", cex=0.66, line=1.5)
  with(predlist[["cool"]], {
    lapply(predlist[["cool"]], function(dom) {
      addpoly(dom$time, dom$lci, dom$uci, col=alpha(list("C"="blue", "D"="red")[[dom$dom[1]]], 0.2), xpd=NA)
      lines(dom$time, dom$fit, lty=1)
    })
  })
  # Plot raw data +/- standard deviation
  lapply(datlist[["cool"]], function(dom) {
    arrows(dom$time, dom$mean + dom$sd, dom$time, dom$mean - dom$sd, code=3, angle=90, length=0.05, xpd=NA,
           col=list("C"="blue", "D"="red")[[dom$dom[1]]])
    points(dom$mean ~ dom$time, pch=21, bg=list("C"="blue", "D"="red")[[dom$dom[1]]], ylim=c(0, 1), cex=1)
  })
  with(fvfmsigs[which(fvfmsigs$ramp=="cool"), ],
       points(time, c(1,1,1), pch="*", cex=1.5))
  # Cooling S/H
  par(mar=c(5,3,0,1))
  plot(NA, xlim=c(0,63), ylim=c(-6, 0.5), bty="n", tck=-0.03, ylab="", xlab="", xaxt="n")
  axis(side=1, at=seq(0,63,7), labels=seq(0,63,7), line=2, tck=-0.03)
  axis(side=1, at=seq(0,63,7), labels=NA, lwd=0, lwd.ticks=par("lwd"), line=2, tck=0.03)
  mtext(side=1, text="days (below) and temperature (°C, above)", cex=0.66, line=3.5)
  axis(side=1, at=seq(0,62,7)+3.5, labels=paste0(seq(23,15,-1), "°"), tick=F, line=0.5)
  mtext(side=2, text="log10 S/H ratio", cex=0.66, line=1.5)
  with(geostats$cool.C, {
    lines(time, mean, type="o", col="black", lty=2, pch=21, bg="blue")
    arrows(time, mean-sd, time, mean+sd, code=3, angle=90, length=0.05, col="blue")
    text(time, mean, labels=c("","","","-94.70%"), pos=2)
    text(time-1, mean-0.25, labels=c("","","","*"), pos=4, xpd=T, cex=1.5)})  # values in shsigs
  with(geostats$cool.D, {
    lines(time, mean, type="o", col="black", lty=2, pch=21, bg="red")
    arrows(time, mean-sd, time, mean+sd, code=3, angle=90, length=0.05, col="red")})
  # with(geostats$cool.mix, {
  #   lines(time, mean, type="o", col="black", lty=2, pch=21, bg="red")
  #   arrows(time, mean-sd, time, mean+sd, code=3, angle=90, length=0.05, col="red")})
  # Heating Fv/Fm
  par(mar=c(1,3,4,1))
  plot(NA, xlim=c(0,63), ylim=c(0, 1), bty="n", tck=-0.03, ylab="fvfm", xlab="days",
       xaxt="n", ann=F)
  title("B. Heating", adj=0)
  mtext(side=2, text="Relative Fv/Fm", cex=0.66, line=1.5)
  with(predlist[["heat"]], {
    lapply(predlist[["heat"]], function(dom) {
      addpoly(dom$time, dom$lci, dom$uci, col=alpha(list("C"="blue", "D"="red")[[dom$dom[1]]], 0.2), xpd=NA)
      lines(dom$time, dom$fit, lty=1)
    })
  })
  with(fvfmsigs[which(fvfmsigs$ramp=="heat"), ],
       points(time, c(1.1,1.1,1,0.9), pch="*", cex=1.5, xpd=T))
  # Plot raw data +/- standard deviation
  lapply(datlist[["heat"]], function(dom) {
    arrows(dom$time, dom$mean + dom$sd, dom$time, dom$mean - dom$sd, code=3, angle=90, length=0.05, xpd=NA,
           col=list("C"="blue", "D"="red")[[dom$dom[1]]])
    points(dom$mean ~ dom$time, pch=21, bg=list("C"="blue", "D"="red")[[dom$dom[1]]], ylim=c(0, 1), cex=1)
  })
  # Heating S/H
  par(mar=c(5,3,0,1))
  plot(NA, xlim=c(0,63), ylim=c(-6, 0.5), bty="n", tck=-0.03, ylab="", xlab="", xaxt="n")
  axis(side=1, at=seq(0,63,7), labels=seq(0,63,7), tck=-0.03, line=2)
  axis(side=1, at=seq(0,42,7), lwd=0, lwd.ticks=par("lwd"), labels=NA, tck=0.03, line=2)
  mtext(side=1, text="days (below) and temperature (°C, above)", cex=0.66, line=3.5)
  mtext(side=2, text="log10 S/H ratio", cex=0.66, line=1.5)
  axis(side=1, at=seq(0,41,7)+3.5, labels=paste0(seq(30,35,1), "°"), tick=F, line=0.5)
  with(geostats$heat.C, {
    lines(time, mean, type="o", col="black", lty=2, pch=21, bg="blue")
    arrows(time, mean-sd, time, mean+sd, code=3, angle=90, length=0.05, col="blue")
    text(time, mean, labels=c("","-99.83%", "-99.99%"), pos=4)
    text(time+1, mean-0.25, labels=c("","*", "*"), pos=2, cex=1.5)})  # values for %loss in shsigs
  with(geostats$heat.D, {
    lines(time, mean, type="o", col="black", lty=2, pch=21, bg="red")
    arrows(time, mean-sd, time, mean+sd, code=3, angle=90, length=0.05, col="red")
    text(time, mean, labels=c("","", "-93.82%"), pos=4)
    text(time+1, mean-0.25, labels=c("","", "*"), pos=2, cex=1.5)})
  # with(geostats$heat.mix, {
  #   lines(time, mean, type="o", col="black", lty=2, pch=21, bg="red")
  #   arrows(time, mean-sd, time, mean+sd, code=3, angle=90, length=0.05, col="red")})





# # plot contrasts from lsmeans + SE
# contrl <- split(contr, f=contr$.id)
# plot(NA, xlim=c(0,63), ylim=c(-10, 2), bty="n", tck=-0.03, ylab="log tot.SH", xlab="days")
# lapply(contrl[c(1,2)], function(x) lines(c(0, as.numeric(substr(x$contrast, 1, 2))),
#                                          c(0, x$estimate), type="o"))
# lapply(contrl[c(1,2)], function(x) arrows(as.numeric(substr(x$contrast, 1, 2)), x$estimate-x$SE, 
#                                           as.numeric(substr(x$contrast, 1, 2)), x$estimate+x$SE,
#                                           code=3, angle=90, length=0.05))
# plot(NA, xlim=c(0,63), ylim=c(-10, 2), bty="n", tck=-0.03, ylab="log tot.SH", xlab="days")
# lapply(contrl[c(3,4)], function(x) lines(c(0, as.numeric(substr(x$contrast, 1, 2))),
#                                          c(0, x$estimate), type="o"))
# lapply(contrl[c(3,4)], function(x) arrows(as.numeric(substr(x$contrast, 1, 2)), x$estimate-x$SE, 
#                                           as.numeric(substr(x$contrast, 1, 2)), x$estimate+x$SE,
#                                           code=3, angle=90, length=0.05))


# C vs. D scatterplots
shdf <- droplevels(subset(data, !is.na(tot.SH)))
shdf[which(shdf[, "C.SH"]==0), "C.SH"] <- 1e-6
shdf[which(shdf[, "D.SH"]==0), "D.SH"] <- 1e-6
head(shdf)
par(mfrow=c(1,1))

plot(log10(C.SH) ~ log10(D.SH), data=subset(shdf, ramp=="heat" & time=="0" & dom=="D"), ylim=c(-5,1), xlim=c(-5,1))
library(mixtools)
with(subset(shdf, ramp=="heat" & time=="0" & dom=="D"), {
  plot(log10(C.SH) ~ log10(D.SH), ylim=c(-6,0), xlim=c(-6,0))
  dataEllipse(log10(D.SH), log10(C.SH), add=T, levels=0.9)
})
with(subset(shdf, ramp=="heat" & time=="28" & dom=="D"), {
  points(log10(C.SH) ~ log10(D.SH), ylim=c(-6,0), xlim=c(-6,0))
  dataEllipse(log10(D.SH), log10(C.SH), add=T, levels=0.9)
})
with(subset(shdf, ramp=="heat" & time=="42" & dom=="D"), {
  points(log10(C.SH) ~ log10(D.SH), ylim=c(-6,0), xlim=c(-6,0))
  dataEllipse(log10(D.SH), log10(C.SH), add=T, levels=0.9)
})
with(subset(shdf, ramp=="heat" & time=="0" & dom=="C"), {
  points(log10(C.SH) ~ log10(D.SH), ylim=c(-6,0), xlim=c(-6,0))
  dataEllipse(log10(D.SH), log10(C.SH), add=T, levels=0.9, col="blue")
})
with(subset(shdf, ramp=="heat" & time=="28" & dom=="C"), {
  points(log10(C.SH) ~ log10(D.SH), ylim=c(-6,0), xlim=c(-6,0))
  dataEllipse(log10(D.SH), log10(C.SH), add=T, levels=0.9, col="blue")
})
with(subset(shdf, ramp=="heat" & time=="42" & dom=="C"), {
  points(log10(C.SH) ~ log10(D.SH), ylim=c(-6,0), xlim=c(-6,0))
  dataEllipse(log10(D.SH), log10(C.SH), add=T, levels=0.9, col="blue")
})

with(subset(shdf, ramp=="cool" & time=="0" & dom=="D"), {
  plot(log10(C.SH) ~ log10(D.SH), ylim=c(-6,0), xlim=c(-6,0))
  dataEllipse(log10(D.SH), log10(C.SH), add=T, levels=0.9)
})
with(subset(shdf, ramp=="cool" & time=="28" & dom=="D"), {
  points(log10(C.SH) ~ log10(D.SH), ylim=c(-6,0), xlim=c(-6,0))
  dataEllipse(log10(D.SH), log10(C.SH), add=T, levels=0.9)
})
with(subset(shdf, ramp=="cool" & time=="42" & dom=="D"), {
  points(log10(C.SH) ~ log10(D.SH), ylim=c(-6,0), xlim=c(-6,0))
  dataEllipse(log10(D.SH), log10(C.SH), add=T, levels=0.9)
})
with(subset(shdf, ramp=="cool" & time=="63" & dom=="D"), {
  points(log10(C.SH) ~ log10(D.SH), ylim=c(-6,0), xlim=c(-6,0))
  dataEllipse(log10(D.SH), log10(C.SH), add=T, levels=0.9)
})
with(subset(shdf, ramp=="cool" & time=="0" & dom=="C"), {
  points(log10(C.SH) ~ log10(D.SH), ylim=c(-6,0), xlim=c(-6,0))
  dataEllipse(log10(D.SH), log10(C.SH), add=T, levels=0.9, col="blue")
})
with(subset(shdf, ramp=="cool" & time=="28" & dom=="C"), {
  points(log10(C.SH) ~ log10(D.SH), ylim=c(-6,0), xlim=c(-6,0))
  dataEllipse(log10(D.SH), log10(C.SH), add=T, levels=0.9, col="blue")
})
with(subset(shdf, ramp=="cool" & time=="42" & dom=="C"), {
  points(log10(C.SH) ~ log10(D.SH), ylim=c(-6,0), xlim=c(-6,0))
  dataEllipse(log10(D.SH), log10(C.SH), add=T, levels=0.9, col="blue")
})
with(subset(shdf, ramp=="cool" & time=="63" & dom=="C"), {
  points(log10(C.SH) ~ log10(D.SH), ylim=c(-6,0), xlim=c(-6,0))
  dataEllipse(log10(D.SH), log10(C.SH), add=T, levels=0.9, col="blue")
})

plot(log10(C.SH) ~ log10(D.SH), data=subset(shdf, ramp=="heat" & time=="28" & dom=="D"), ylim=c(-5,1), xlim=c(-5,1))
plot(log10(C.SH) ~ log10(D.SH), data=subset(shdf, ramp=="heat" & time=="42" & dom=="D"), ylim=c(-5,1), xlim=c(-5,1))

plot(log10(C.SH) ~ log10(D.SH), data=subset(shdf, ramp=="heat" & time=="0" & dom=="C"), ylim=c(-5,1), xlim=c(-5,1))
plot(log10(C.SH) ~ log10(D.SH), data=subset(shdf, ramp=="heat" & time=="28" & dom=="C"), ylim=c(-5,1), xlim=c(-5,1))
plot(log10(C.SH) ~ log10(D.SH), data=subset(shdf, ramp=="heat" & time=="42" & dom=="C"), ylim=c(-5,1), xlim=c(-5,1))

plot(log10(C.SH) ~ log10(D.SH), data=subset(shdf, ramp=="heat" & time=="0" & dom=="mix"), ylim=c(-5,1), xlim=c(-5,1))
plot(log10(C.SH) ~ log10(D.SH), data=subset(shdf, ramp=="heat" & time=="28" & dom=="mix"), ylim=c(-5,1), xlim=c(-5,1))
plot(log10(C.SH) ~ log10(D.SH), data=subset(shdf, ramp=="heat" & time=="42" & dom=="mix"), ylim=c(-5,1), xlim=c(-5,1))

plot(log10(C.SH) ~ log10(D.SH), data=subset(shdf, ramp=="cool" & time=="0" & dom=="mix"), ylim=c(-5,1), xlim=c(-5,1))
plot(log10(C.SH) ~ log10(D.SH), data=subset(shdf, ramp=="cool" & time=="28" & dom=="mix"), ylim=c(-5,1), xlim=c(-5,1))
plot(log10(C.SH) ~ log10(D.SH), data=subset(shdf, ramp=="cool" & time=="63" & dom=="mix"), ylim=c(-5,1), xlim=c(-5,1))

plot(log10(C.SH) ~ log10(D.SH), data=subset(shdf, ramp=="cool" & time=="42" & dom=="mix"), ylim=c(-5,1), xlim=c(-5,1))

means <- aggregate(shdf[,c("C.SH","D.SH")], by=list(ramp=shdf$ramp, dom=shdf$dom, time=shdf$time), FUN=function(x) mean(log10(x)))
colnames(means) <- c("ramp", "dom", "time", "Cmean", "Dmean")
ses <- aggregate(shdf[,c("C.SH","D.SH")], by=list(ramp=shdf$ramp, dom=shdf$dom, time=shdf$time), FUN=function(x) sd(log10(x))/sqrt(length(x)))
colnames(ses) <- c("ramp", "dom", "time", "Cse", "Dse")

dat <- merge(means, ses, all=T)

par(mfrow=c(1,2), mar=c(3,3,1,1), mgp=c(1.5,0.4,0), tcl=-0.3)

    df <- dat[dat$ramp=="cool", ]
    plot(NA, xlim=c(-6,0), ylim=c(-6,0),xlab="Clade D abundance (log10 S/H)",ylab="Clade C abundance (log10 S/H)")
    abline(a=0,b=1,lty=2)
    arrows(df$Dmean, df$Cmean+df$Cse, df$Dmean, df$Cmean-df$Cse, length=0.025, angle=90, code=3, col="black")
    arrows(df$Dmean+df$Dse, df$Cmean, df$Dmean-df$Dse, df$Cmean, length=0.025, angle=90, code=3, col="black")
    for (j in 2:nrow(df)) {
      if (df$time[j] > df$time[j-1]) {
        arrows(df$Dmean[j-1], df$Cmean[j-1], df$Dmean[j], df$Cmean[j], length=0.1, code=2, lwd=2,
                     col=c("blue","red","green")[df$dom[j]])
      }
    }
    text(-6, 0, labels=expression(bold("A. Cooling")), adj=0)
    legend("bottomright", legend=c("C-dominant", "D-dominant", "C+D"), lty=1, col=c("blue","red","green"), inset=0.05)
    df <- dat[dat$ramp=="heat", ]
    plot(NA, xlim=c(-6,0), ylim=c(-6,0), xlab="Clade D abundance (log10 S/H)",ylab="Clade C abundance (log10 S/H)")
    abline(a=0,b=1,lty=2)
    arrows(df$Dmean, df$Cmean+df$Cse, df$Dmean, df$Cmean-df$Cse, length=0.025, angle=90, code=3, col="black")
    arrows(df$Dmean+df$Dse, df$Cmean, df$Dmean-df$Dse, df$Cmean, length=0.025, angle=90, code=3, col="black")
    for (j in 2:nrow(df)) {
      if (df$time[j] > df$time[j-1]) {
        arrows(df$Dmean[j-1], df$Cmean[j-1], df$Dmean[j], df$Cmean[j], length=0.1, code=2, lwd=2,
               col=c("blue","red","green")[df$dom[j]])
      }
    }
    text(-6, 0, labels=expression(bold("B. Heating")), adj=0)

# PLOT LSMEANS
Cmod <- lmer(log10(C.SH) ~ factor(time) * ramp * dom + (1|mother/sample), data=shdf)
Clsm <- lsmeans(Cmod, specs=c("time", "ramp", "dom"))
Clsmeans <- summary(Clsm)
Dmod <- lmer(log10(D.SH) ~ factor(time) * ramp * dom + (1|mother/sample), data=shdf)
Dlsm <- lsmeans(Dmod, specs=c("time", "ramp", "dom"))
Dlsmeans <- summary(Dlsm)

plot(Clsmeans$lsmean ~ Dlsmeans$lsmean)
df1 <- data.frame(time=Clsmeans$time, ramp=Clsmeans$ramp, dom=Clsmeans$dom, Cmean=Clsmeans$lsmean, Dmean=Dlsmeans$lsmean, Cse=Clsmeans$SE, Dse=Dlsmeans$SE)
df <- subset(df1, ramp=="cool")
plot(NA, xlim=c(-6,0), ylim=c(-6,0),xlab="Clade D abundance (log10 S/H)",ylab="Clade C abundance (log10 S/H)")
abline(a=0,b=1,lty=2)
arrows(df$Dmean, df$Cmean+df$Cse, df$Dmean, df$Cmean-df$Cse, length=0.025, angle=90, code=3, col="black")
arrows(df$Dmean+df$Dse, df$Cmean, df$Dmean-df$Dse, df$Cmean, length=0.025, angle=90, code=3, col="black")
for (j in 2:nrow(df)) {
  if (df$time[j] > df$time[j-1]) {
    arrows(df$Dmean[j-1], df$Cmean[j-1], df$Dmean[j], df$Cmean[j], length=0.1, code=2, lwd=2,
           col=c("blue","red","green")[df$dom[j]])
  }
}
text(-6, 0, labels=expression(bold("A. Cooling")), adj=0)
legend("bottomright", legend=c("C-dominant", "D-dominant", "C+D"), lty=1, col=c("blue","red","green"), inset=0.05)
df <- subset(df1, ramp=="heat")

plot(NA, xlim=c(-6,0), ylim=c(-6,0), xlab="Clade D abundance (log10 S/H)",ylab="Clade C abundance (log10 S/H)")
abline(a=0,b=1,lty=2)
arrows(df$Dmean, df$Cmean+df$Cse, df$Dmean, df$Cmean-df$Cse, length=0.025, angle=90, code=3, col="black")
arrows(df$Dmean+df$Dse, df$Cmean, df$Dmean-df$Dse, df$Cmean, length=0.025, angle=90, code=3, col="black")
for (j in 2:nrow(df)) {
  if (df$time[j] > df$time[j-1]) {
    arrows(df$Dmean[j-1], df$Cmean[j-1], df$Dmean[j], df$Cmean[j], length=0.1, code=2, lwd=2,
           col=c("blue","red","green")[df$dom[j]])
  }
}
text(-6, 0, labels=expression(bold("B. Heating")), adj=0)
contrast(lsm, method="pairwise", adjust="none")

