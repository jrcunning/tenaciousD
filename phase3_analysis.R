# LOAD DATA AND LIBRARIES -----
library(lme4)
library(scales)
library(lsmeans)
library(reshape2)
library(plyr)
library(lattice)

# Define function to use for plotting confidence intervals
addpoly <- function(x, y1, y2, col=alpha("lightgrey", 0.8), ...){
  ii <- order(x)
  y1 <- y1[ii]
  y2 <- y2[ii]
  x <- x[ii]
  polygon(c(x, rev(x)), c(y1, rev(y2)), col=col, border=NA, ...)
}

# Import data
data <- read.csv("phase3_2.csv")

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

# Test for effect of clade in corals with same history -----
# Do any corals in either treatment have same history but different dominant clade?
table(data$history, data$dom, data$ramp)  # Yes: history B' in heating and E' in cooling
# Test if corals with same history showed different responses based on dominant clade

# History E' (DCMU-24-ctrl-24) in cooling treatment
Edf <- subset(data, history=="E'" & ramp=="cool")
table(unique(Edf[,c("sample", "dom")])$dom)  # 7 corals with Cdom, 4 corals with Ddom
mod <- lmerTest::lmer(rfvfm ~ timef * dom + (1|mother/sample), data=Edf, na.action=na.omit)
lmerTest::anova(mod)  # clade does not impact Fv/Fm
plot(Effect(c("timef", "dom"), mod))
mod <- lmerTest::lmer(log(tot.SH) ~ timef * dom + (1|mother/sample), data=Edf, na.action=na.omit)
lmerTest::anova(mod)  # clade marginally impacts totSH
plot(Effect(c("timef", "dom"), mod))

# History B' (ctrl-29-ctrl-29) in heating treatment
Bdf <- subset(data, history=="B'" & ramp=="heat")
table(unique(Bdf[,c("sample", "dom")])$dom)  # 6 corals with Cdom, 4 corals with Ddom
mod <- lmerTest::lmer(rfvfm ~ timef * dom + (1|mother/sample), data=Bdf, na.action=na.omit)
lmerTest::anova(mod) # clade does not impact Fv/Fm
plot(Effect(c("timef", "dom"), mod))
mod <- lmerTest::lmer(log(tot.SH) ~ timef * dom + (1|mother/sample), data=Bdf, na.action=na.omit)
lmerTest::anova(mod)
plot(Effect(c("timef", "dom"), mod))  # clade strongly impacts totSH

# Test for effect of history in corals with same clade -----
# Do any corals have the same clade but different history?
table(data$dom, data$history, data$ramp) #Yes, cooling Cdom, A' (ctrl-24-ctrl-24) vs. E' (DCMU-24-ctrl-24)
df <- droplevels(subset(data, ramp=="cool" & dom=="C"))
table(unique(df[,c("sample","history")])$history) # 3 corals A', 7 corals E'
mod <- lmerTest::lmer(rfvfm ~ timef * history + (1|mother/sample), data=df, na.action=na.omit)
lmerTest::anova(mod)  # no impact of history on fv/fm
plot(Effect(c("timef", "history"), mod))
mod <- lmerTest::lmer(log(tot.SH) ~ timef * history + (1|mother/sample), data=df, na.action=na.omit)
lmerTest::anova(mod)  # no impact of history on totSH
plot(Effect(c("timef", "history"), mod))

# Analyze Fv/Fm -----
finalheat <- subset(data, ramp=="heat" & time==42)
table(finalheat$dom, finalheat$fvfm)
# Cooling treatment
gm.cool <- gamm(rfvfm ~ dom + s(time, by=dom, k=4), random=list(mother=~1, sample=~1), 
                data=subset(data, ramp=="cool"), na.action=na.omit)
pdat.cool <- expand.grid(time=seq(0,63,1), ramp=factor("cool"), dom=factor(c("C", "D")))
gm.cool.pred <- data.frame(cbind(pdat.cool, fit=predict(gm.cool$gam, pdat.cool)))
# Simulate predicted values from posterior distribution of beta 1000 times
set.seed(789)
Rbeta <- mvrnorm(n = 1000, coef(gm.cool$gam), vcov(gm.cool$gam))
Xp <- predict(gm.cool$gam, newdata = pdat.cool, type = "lpmatrix")
sim <- Xp %*% t(Rbeta)
# Extract 95% confidence intervals of simulated values for plotting
gm.cool.pred$lci <- apply(sim, 1, quantile, 0.025)
gm.cool.pred$uci <- apply(sim, 1, quantile, 0.975)
# Calculate 95% confidence interval on difference between C and D corals at select dates (=sampling dates)
sim1 <- cbind(pdat.cool, sim)
sim1 <- sim1[sim1$time %in% c(0,7,14,21,28,35,42,49,56,63), ]
sim1split <- split(sim1, f=sim1$time, drop=T)
CDdiff.95CI <- lapply(sim1split, function(x) quantile(apply(x[,4:ncol(x)], 2, diff), c(0.025, 0.975)))
# Test if 95% CI on difference between C and D corals contains zero
cool.test0 <- lapply(CDdiff.95CI, function(x) x[1] < 0 & 0 < x[2])
cool.test0[cool.test0==F]  # Shows when C and D are significantly different with p<0.05

# Heating treatment
gm.heat <- gamm(rfvfm ~ dom + s(time, by=dom, k=4), random=list(mother=~1, sample=~1), 
                data=subset(data, ramp=="heat"), na.action=na.omit)
pdat.heat <- expand.grid(time=seq(0,42,1), ramp=factor("heat"), dom=factor(c("C", "D")))
gm.heat.pred <- data.frame(cbind(pdat.heat, fit=predict(gm.heat$gam, pdat.heat)))
# Simulate predicted values from posterior distribution of beta 1000 times
set.seed(789)
Rbeta <- mvrnorm(n = 1000, coef(gm.heat$gam), vcov(gm.heat$gam))
Xp <- predict(gm.heat$gam, newdata = pdat.heat, type = "lpmatrix")
sim <- Xp %*% t(Rbeta)
# Extract 95% confidence intervals of simulated values for plotting
gm.heat.pred$lci <- apply(sim, 1, quantile, 0.025)
gm.heat.pred$uci <- apply(sim, 1, quantile, 0.975)
# Calculate 95% confidence interval on difference between C and D corals at select dates (=sampling dates)
sim1 <- cbind(pdat.heat, sim)
sim1 <- sim1[sim1$time %in% c(0,7,14,21,28,35,42,49,56,63), ]
sim1split <- split(sim1, f=sim1$time, drop=T)
CDdiff.95CI <- lapply(sim1split, function(x) quantile(apply(x[,4:ncol(x)], 2, diff), c(0.025, 0.975)))
# Test if 95% CI on difference between C and D corals contains zero
heat.test0 <- lapply(CDdiff.95CI, function(x) x[1] < 0 & 0 < x[2])
heat.test0[heat.test0==F]  # Shows when C and D are significantly different with p<0.05

# Analyze total S/H ratio ------------
initsh <- subset(data, time==0 & dom=="D")
table(initsh$propD==0) # 13 all C
table(initsh$propD==1) # 12 all D
# 115 - 25 = 90 mixed --> 78.26% mixed
hist(log10(initsh[initsh$dom=="C", "propD"]))
hist(log10(1-initsh[initsh$dom=="D", "propD"]))
boxplot(log10(1-initsh[initsh$dom=="D", "propD"]))

finalsh <- subset(data, ramp=="heat" & time==42)
table(finalsh$dom, finalsh$tot.SH)
# Fit model for cooling treatment and test for effect of past thermal history
coolSHmod <- lmerTest::lmer(log10(tot.SH) ~ timef * dom + (1|mother/sample), 
                            data=subset(data, ramp=="cool"), na.action=na.omit)
# Fit model for heating treatment and test for effect of past thermal history
heatSHmod <- lmerTest::lmer(log10(tot.SH) ~ timef * dom + (1|mother/sample), 
                            data=subset(data, ramp=="heat"), na.action=na.omit)
# Compare S/H ratio within each group to time zero
cool.dunnett <- rbind(contrast(lsmeans::lsmeans(coolSHmod, specs=c("timef", "dom")), "dunnett", by="dom"))
heat.dunnett <- rbind(contrast(lsmeans::lsmeans(heatSHmod, specs=c("timef", "dom")), "dunnett", by="dom"))
# Calculate percent change relative to time zero and extract significant contrasts
cool.dunnett <- within(data.frame(summary(cool.dunnett)), {
  loss <- paste0(sprintf("%.1f", round(-(1-(10^estimate)) * 100, 1)), "%")
  day <- substr(contrast, 1, 2)
})
cool.diff <- cool.dunnett[cool.dunnett$p.value < 0.05, ]
heat.dunnett <- within(data.frame(summary(heat.dunnett)), {
  loss <- paste0(sprintf("%.1f", round(-(1-(10^estimate)) * 100, 1)), "%")
  day <- substr(contrast, 1, 2)
})
heat.diff <- heat.dunnett[heat.dunnett$p.value < 0.05, ]

# Figure 1: Fv/Fm and total S/H under (A.) cooling and (B.) heating -------------
# Prepare raw data for plotting
datsumm <- rbind(expand.grid(time=seq(0,63,7), ramp=factor("cool"), dom=factor(c("C","D"))),
                 expand.grid(time=seq(0,42,7), ramp=factor("heat"), dom=factor(c("C","D"))))
datsumm$rfvfm.mean <- aggregate(data$rfvfm, by=list(interaction(data$time, data$dom, data$ramp)), FUN=mean, na.rm=T)$x
datsumm$rfvfm.sd <- aggregate(data$rfvfm, by=list(interaction(data$time, data$dom, data$ramp)), FUN=sd, na.rm=T)$x
datsumm$logtot.SH.mean <- aggregate(log10(data$tot.SH), by=list(interaction(data$time, data$dom, data$ramp)), FUN=mean, na.rm=T)$x
datsumm$logtot.SH.sd <- aggregate(log10(data$tot.SH), by=list(interaction(data$time, data$dom, data$ramp)), FUN=sd, na.rm=T)$x
datlist <- split(datsumm, f=datsumm$ramp)
datlist <- lapply(datlist, function(x) rev(split(x, f=x$dom)))
# Create figure
#pdf(file = "output/Figure1.pdf", width=6.85, height=3.425)
layout(mat=matrix(c(1,1,2,2,3,3,4,4), ncol=2))
par(mgp=c(1.5,0.25,0), tck=-0.06, xpd=NA)
# Cooling Fv/Fm
par(mar=c(1,3,4,1))
plot(NA, xlim=c(0,63), ylim=c(0, 1), bty="n", tck=-0.03,
     xaxt="n", ann=F)
title("A. Cooling", adj=0, cex.main=1.5, line=2)
mtext(side=2, text="Relative Fv/Fm", cex=0.75, line=1.5)
with(subset(gm.cool.pred, dom=="C"), {
  addpoly(time, uci, lci, col=alpha("blue", 0.4))
  lines(time, fit)
})
with(subset(gm.cool.pred, dom=="D"), {
  addpoly(time, uci, lci, col=alpha("red", 0.4))
  lines(time, fit)
})
points(c(21,28,35,42), c(1,0.97,0.94,0.91), pch="*", cex=2, xpd=T)
# Plot raw data +/- standard deviation
lapply(datlist[["cool"]], function(dom) {
  arrows(dom$time, dom$rfvfm.mean + dom$rfvfm.sd, dom$time, dom$rfvfm.mean - dom$rfvfm.sd, code=3, angle=90, length=0.05, xpd=NA,
         col=list("C"="blue", "D"="red")[[dom$dom[1]]])
  points(dom$rfvfm.mean ~ dom$time, pch=21, bg=list("C"="blue", "D"="red")[[dom$dom[1]]], ylim=c(0, 1), cex=1)
})
legend("bottom", legend=c("Corals initially dominated by clade C", "Corals initially dominated by clade D"),
       pch=21, pt.bg=c("blue","red"), bty="n")
# Cooling S/H
par(mar=c(5,3,0,1))
plot(NA, xlim=c(0,63), ylim=c(-4, 0.5), bty="n", tck=-0.03, ylab="", xlab="", xaxt="n")
axis(side=1, at=seq(0,63,7), labels=seq(0,63,7), line=2, tck=-0.03)
axis(side=1, at=seq(0,63,7), labels=NA, lwd=0, lwd.ticks=par("lwd"), line=2, tck=0.03)
mtext(side=1, text="days (below) and temperature (°C, above)", cex=0.75, line=3.5)
axis(side=1, at=seq(0,62,7)+3.5, labels=paste0(seq(23,15,-1), "°"), tick=F, line=0.5)
mtext(side=2, text="log10 S/H ratio", cex=0.75, line=1.5)
lapply(datlist[["cool"]], function(dom) {
  dom <- na.omit(dom)
  arrows(dom$time, dom$logtot.SH.mean + dom$logtot.SH.sd, dom$time, dom$logtot.SH.mean - dom$logtot.SH.sd, code=3, angle=90, length=0.05, xpd=NA,
         col=list("C"="blue", "D"="red")[[dom$dom[1]]])
  lines(dom$time, dom$logtot.SH.mean, type="l", col="black", lty=2, pch=21)
  points(dom$logtot.SH.mean ~ dom$time, pch=21, bg=list("C"="blue", "D"="red")[[dom$dom[1]]], ylim=c(0, 1), cex=1)
})
for (i in 1:nrow(cool.diff)) {
  with(cool.diff, {
    dom <- dom[i]
    points(as.numeric(day[i])+3, with(datlist[["cool"]][[dom]], logtot.SH.mean[time==day[i]]), pch="*", cex=2)
    text(day[i], with(datlist[["cool"]][[dom]], logtot.SH.mean[time==day[i]]), labels=loss[i], pos=c(2,2,2)[i])
  })
}

     
# Heating Fv/Fm
par(mar=c(1,3,4,1))
plot(NA, xlim=c(0,63), ylim=c(0, 1), bty="n", tck=-0.03, ylab="fvfm", xlab="days",
     xaxt="n", ann=F, xpd=NA)
title("B. Heating", adj=0, cex.main=1.5, xpd=NA, line=2)
mtext(side=2, text="Relative Fv/Fm", cex=0.75, line=1.5)
with(subset(gm.heat.pred, dom=="C" & time <=42), {
  addpoly(time, uci, lci, col=alpha("blue", 0.4))
  lines(time, fit)
})
with(subset(gm.heat.pred, dom=="D" & time <=42), {
  addpoly(time, uci, lci, col=alpha("red", 0.4))
  lines(time, fit)
})
points(c(7,14,21), c(1.1,1.08,0.99), pch="*", cex=2, xpd=T)
# Plot raw data +/- standard deviation
lapply(datlist[["heat"]], function(dom) {
  arrows(dom$time, dom$rfvfm.mean + dom$rfvfm.sd, dom$time, dom$rfvfm.mean - dom$rfvfm.sd, code=3, angle=90, length=0.05, xpd=NA,
         col=list("C"="blue", "D"="red")[[dom$dom[1]]])
  points(dom$rfvfm.mean ~ dom$time, pch=21, bg=list("C"="blue", "D"="red")[[dom$dom[1]]], ylim=c(0, 1), cex=1)
})
# Heating S/H
par(mar=c(5,3,0,1))
plot(NA, xlim=c(0,63), ylim=c(-4, 0.5), bty="n", tck=-0.03, ylab="", xlab="", xaxt="n")
axis(side=1, at=seq(0,63,7), labels=seq(0,63,7), tck=-0.03, line=2)
axis(side=1, at=seq(0,42,7), lwd=0, lwd.ticks=par("lwd"), labels=NA, tck=0.03, line=2)
mtext(side=1, text="days (below) and temperature (°C, above)", cex=0.75, line=3.5)
mtext(side=2, text="log10 S/H ratio", cex=0.75, line=1.5)
axis(side=1, at=seq(0,41,7)+3.5, labels=paste0(seq(30,35,1), "°"), tick=F, line=0.5)
lapply(datlist[["heat"]], function(dom) {
  dom <- na.omit(dom)
  arrows(dom$time, dom$logtot.SH.mean + dom$logtot.SH.sd, dom$time, dom$logtot.SH.mean - dom$logtot.SH.sd, code=3, angle=90, length=0.05, xpd=NA,
         col=list("C"="blue", "D"="red")[[dom$dom[1]]])
  lines(dom$time, dom$logtot.SH.mean, type="l", col="black", lty=2, pch=21)
  points(dom$logtot.SH.mean ~ dom$time, pch=21, bg=list("C"="blue", "D"="red")[[dom$dom[1]]], ylim=c(0, 1), cex=1)
})
for (i in 1:nrow(heat.diff)) {
  with(heat.diff, {
    dom <- dom[i]
    points(as.numeric(day[i])+3, with(datlist[["heat"]][[dom]], logtot.SH.mean[time==day[i]]), pch="*", cex=2)
    text(day[i], with(datlist[["heat"]][[dom]], logtot.SH.mean[time==day[i]]), labels=loss[i], pos=c(2,2,2)[i])
  })
}
#dev.off()

# Modeling C and D responses together -----
df <- data[!is.na(data$tot.SH), ]
df2 <- melt(df, id.vars=c(1:6,9:13), value.name="SH")
colnames(df2)[colnames(df2)=="variable"] <- "clade"
# Fit model for cooling treatment
coolmod <- lmerTest::lmer(log10(SH) ~ timef * dom * clade + (timef+dom+clade|mother/sample), 
                          data=subset(df2, ramp=="cool"))
coolmod.lsm <- lsmeans(coolmod, specs=c("clade", "timef", "dom"))
coolmod.contr <- data.frame(summary(contrast(coolmod.lsm, "consec", by=c("dom", "clade"))))
coolmod.sigs <- coolmod.contr[coolmod.contr$p.value < 0.05, ]
coolmod.sigs$loss <- -(1-10^coolmod.sigs$estimate)
coolmod.sigs
coolmod.lsm <- data.frame(summary(coolmod.lsm))

# Fit model for heating treatment
heatmod <- lmerTest::lmer(log10(SH) ~ timef * dom * clade + (timef+dom+clade|mother/sample), 
                          data=subset(df2, ramp=="heat"))
heatmod.lsm <- lsmeans(heatmod, specs=c("clade", "timef", "dom"))
heatmod.contr <- data.frame(summary(contrast(heatmod.lsm, "consec", by=c("dom", "clade"))))
heatmod.sigs <- heatmod.contr[heatmod.contr$p.value < 0.05, ]
heatmod.sigs$loss <- -(1-10^heatmod.sigs$estimate)
heatmod.sigs
heatmod.lsm <- data.frame(summary(heatmod.lsm))

# Figure 2: Mixed community dynamics under (A.) cooling and (B.) heating -----
# Plot mean abundances of C and D over time by treatment and dominant clade at start
#pdf(file="output/Figure2.pdf", width=6.85, height=3.425)
par(mfrow=c(1,2), mar=c(3,3,2,1), mgp=c(1.5,0.4,0), tcl=-0.3, xpd=F)
df <- dcast(coolmod.lsm, timef + dom ~ clade, value.var="lsmean")
df$C.se <- dcast(coolmod.lsm, timef + dom ~ clade, value.var="SE")$C.SH
df$D.se <- dcast(coolmod.lsm, timef + dom ~ clade, value.var="SE")$D.SH
df <- with(df, df[order(dom, timef), ])
plot(NA, xlim=c(-6,0), ylim=c(-6,0), xlab="Clade D (log10 S/H)", ylab="Clade C (log10 S/H)", cex.axis=0.75, cex.lab=0.75)
abline(a=0,b=1,lty=2)
arrows(df$D.SH, df$C.SH+df$C.se, df$D.SH, df$C.SH-df$C.se, length=0.025, angle=90, code=3, col="black")
arrows(df$D.SH+df$D.se, df$C.SH, df$D.SH-df$D.se, df$C.SH, length=0.025, angle=90, code=3, col="black")
for (j in 2:nrow(df)) {
  if (as.numeric(df$time)[j] > as.numeric(df$time)[j-1]) {
    arrows(df$D.SH[j-1], df$C.SH[j-1], df$D.SH[j], df$C.SH[j], length=0.07, code=2, lwd=2,
           col=c("blue","red","green")[df$dom[j]])
  }
}
with(df[4, ], text(D.SH, C.SH, labels=expression("" %down%C), adj=c(1,-4), cex=0.7))
with(df[8, ], text(D.SH, C.SH, labels=expression("" %down%C), adj=c(1.5,-1), cex=0.7))

xx <- seq(-5,-4,len=4)
arrows(x0=xx[-4], y0=-0.2, x1=xx[-1], code=2, length=0.07, lwd=1.5, col="blue")
arrows(x0=xx[-4], y0=-0.4, x1=xx[-1], code=2, length=0.07, lwd=1.5, col="red")
text(xx, -0.6, labels=c(24,20,18,15), cex=0.7, adj=c(0.5, 0.35))
text(-5.5, c(-0.2, -0.4, -0.6), labels=c("C-dom.", "D-dom.", "°C"), cex=0.75, adj=c(0.5, 0.35))
rect(-6,-0.75,-3.8,0)

text(par("usr")[1], 0.5, labels=expression(bold("A. Cooling")), adj=0, xpd=NA)


df <- dcast(heatmod.lsm, timef + dom ~ clade, value.var="lsmean")
df$C.se <- dcast(heatmod.lsm, timef + dom ~ clade, value.var="SE")$C.SH
df$D.se <- dcast(heatmod.lsm, timef + dom ~ clade, value.var="SE")$D.SH
df <- with(df, df[order(dom, timef), ])
plot(NA, xlim=c(-6,0), ylim=c(-6,0), xlab="Clade D (log10 S/H)",ylab="Clade C (log10 S/H)", cex.axis=0.75, cex.lab=0.75)
abline(a=0,b=1,lty=2)
arrows(df$D.SH, df$C.SH+df$C.se, df$D.SH, df$C.SH-df$C.se, length=0.025, angle=90, code=3, col="black")
arrows(df$D.SH+df$D.se, df$C.SH, df$D.SH-df$D.se, df$C.SH, length=0.025, angle=90, code=3, col="black")
for (j in 2:nrow(df)) {
  if (as.numeric(df$time)[j] > as.numeric(df$time)[j-1]) {
    arrows(df$D.SH[j-1], df$C.SH[j-1], df$D.SH[j], df$C.SH[j], length=0.07, code=2, lwd=2,
           col=c("blue","red","green")[df$dom[j]])
  }
}
text(par("usr")[1], 0.5, labels=expression(bold("B. Heating")), adj=0, xpd=NA)
with(df[2, ], text(D.SH, C.SH, labels=expression("" %down%C), adj=c(1.2,-12), cex=0.7))
with(df[3, ], text(D.SH, C.SH, labels=expression("" %down%C), adj=c(0.55,-4), cex=0.7))
with(df[5, ], text(D.SH, C.SH, labels=expression("" %down%C), adj=c(0.9,-4.5), cex=0.7))
with(df[6, ], text(D.SH, C.SH, labels=expression("" %down%D), adj=c(-1,-1.4), cex=0.7))

xx <- seq(-5,-4,len=3)
arrows(x0=xx[-3], y0=-0.2, x1=xx[-1], code=2, length=0.07, lwd=1.5, col="blue")
arrows(x0=xx[-3], y0=-0.4, x1=xx[-1], code=2, length=0.07, lwd=1.5, col="red")
text(xx, -0.6, labels=c(29,33,35), cex=0.7, adj=c(0.5, 0.35))
text(-5.5, c(-0.2, -0.4, -0.6), labels=c("C-dom.", "D-dom.", "°C"), cex=0.75, adj=c(0.5, 0.35))
rect(-6,-0.75,-3.8,0)
#dev.off()
