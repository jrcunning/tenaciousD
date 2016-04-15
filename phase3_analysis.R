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

# Adjust and transform data -----
# Create factor version of time
data$timef <- factor(data$time)
# Calculate proportion D in each sample
data$propD <- data$D.SH / (data$C.SH + data$D.SH)
# Categorize C- and D-dominated corals based on community composition at time zero
dom <- with(data[data$time==0, ], na.omit(data.frame(sample=sample, 
                                                     dom=ifelse(is.na(propD), NA, ifelse(propD > 0.5, "D", "C")))))
dom <- with(data[data$time==0, ], na.omit(data.frame(sample=sample,
                                                     dom=ifelse(is.na(propD), NA, ifelse(propD > 0.9, "D", ifelse(propD < 0.1, "C", NA))))))
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
data[data$tot.SH<0.000101 & !is.na(data$tot.SH), "tot.SH"] <- 0.000101
# Calculate relative change in Fv/Fm
for (sample in data$sample) {
  sdata <- data[which(data$sample==sample), ]
  for (time in sdata$time) {
    sdata[which(sdata$time==time), "rfvfm"] <- sdata[which(sdata$time==time), "fvfm"] / sdata[which(sdata$time==0), "fvfm"]
    }
  data[which(data$sample==sample), "rfvfm"] <- sdata[, "rfvfm"]
}

# Test effect of past thermal history on initial Fv/Fm and S/H ratio -----
fullmod <- lmerTest::lmer(fvfm ~ history * ramp * dom + (1|mother), data=subset(data, time==0))
fullmod.lsm <- lsmeans::lsmeans(fullmod, specs=c("history","dom", "ramp"))
tests <- data.frame(summary(pairs(fullmod.lsm, by=c("dom", "ramp"))))
tests[!is.na(tests$p.value),]  # No significant differences due to thermal history within groups

fullmod <- lmerTest::lmer(log(tot.SH) ~ history * ramp * dom + (1|mother), data=subset(data, time==0))
fullmod.lsm <- lsmeans::lsmeans(fullmod, specs=c("history","dom", "ramp"))
tests <- data.frame(summary(pairs(fullmod.lsm, by=c("dom", "ramp"))))
tests[!is.na(tests$p.value),]  # No significant differences due to thermal history within groups

# TEST EFFECT OF HISTORY ON FVFM AND SH RATIO
mod <- lmerTest::lmer(rfvfm ~ timef + history + (1|mother/sample), 
                      data=subset(data, ramp=="heat" & dom=="C")); lmerTest::anova(mod) # no effect
lmerTest::step(mod, reduce.random=F)
mod <- lmerTest::lmer(rfvfm ~ timef * history + (1|mother/sample), 
                      data=subset(data, ramp=="heat" & dom=="D")); lmerTest::anova(mod) # no effect
lmerTest::step(mod, reduce.random=F)
mod <- lmerTest::lmer(rfvfm ~ timef * history + (1|mother/sample), 
                      data=subset(data, ramp=="cool" & dom=="C")); lmerTest::anova(mod) # no effect
lmerTest::step(mod, reduce.random=F)
mod <- lmerTest::lmer(rfvfm ~ timef * history + (1|mother/sample), 
                      data=subset(data, ramp=="cool" & dom=="D")); lmerTest::anova(mod) # YES effect
mod <- lmerTest::step(mod, reduce.random=F)$model
modlsm <- lsmeans::lsmeans(mod, specs=c("history"), contr="pairwise")
# group A* has Fv/Fm values ~0.06 higher than C' and C* under cooling treatment
with(subset(data, ramp=="cool" & dom=="D"), aggregate(propD, by=list(history), FUN=mean, na.rm=T))
mod <- lmerTest::lmer(log(tot.SH) ~ timef * history + (1|mother/sample), 
                      data=subset(data, ramp=="heat" & dom=="C")); lmerTest::anova(mod) # no effect
lmerTest::step(mod, reduce.random=F)
mod <- lmerTest::lmer(log(tot.SH) ~ timef * history + (1|mother/sample), 
                      data=subset(data, ramp=="heat" & dom=="D")); lmerTest::anova(mod) # no effect
lmerTest::step(mod, reduce.random=F)
mod <- lmerTest::lmer(log(tot.SH) ~ timef * history + (1|mother/sample), 
                      data=subset(data, ramp=="cool" & dom=="C")); lmerTest::anova(mod) # no effect
lmerTest::step(mod, reduce.random=F)
mod <- lmerTest::lmer(log(tot.SH) ~ timef * history + (1|mother/sample), 
                      data=subset(data, ramp=="cool" & dom=="D")); lmerTest::anova(mod) # no effect
lmerTest::step(mod, reduce.random=F)

table(data$dom, data$history, data$ramp)

Edf <- subset(data, history=="E'")
mod <- lmerTest::lmer(log(tot.SH) ~ timef * dom + (1|mother/sample), data=Edf)
lmerTest::anova(mod)
lmerTest::step(mod)
plot(Effect(c("timef","dom"), mod))

Bdf <- subset(data, history=="B'")
mod <- lmerTest::lmer(log(tot.SH) ~ timef * dom + (1|mother/sample), data=Bdf)
lmerTest::anova(mod)
plot(Effect(c("timef","dom"), mod))

# Analyze Fv/Fm -----
# Visualize all data
# xyplot(rfvfm ~ time | ramp + dom, groups= ~ sample, data = data, type="o", lty=1)



# Fit Generalized Additive Mixed Model (GAMM) to Fv/Fm data
gm <- gamm(rfvfm ~ ramp + dom + s(time, by=interaction(ramp, dom), k=5), random=list(mother=~1, sample=~1), data=data)
pdat <- rbind(expand.grid(time=seq(0,63,1), ramp=factor("cool"), dom=factor(c("C", "D"))),
              expand.grid(time=seq(0,42,1), ramp=factor("heat"), dom=factor(c("C","D"))))
# Get predicted values
gmpreds <- data.frame(cbind(pdat, predict(gm$gam, pdat, re.form=NA, se.fit=T)))
# Simulate predicted values from posterior distribution of beta 1000 times
set.seed(789)
Rbeta <- mvrnorm(n = 1000, coef(gm$gam), vcov(gm$gam))
Xp <- predict(gm$gam, newdata = pdat, type = "lpmatrix")
sim <- Xp %*% t(Rbeta)
# Extract 95% confidence intervals of simulated values for plotting
gmpreds$lci <- apply(sim, 1, quantile, 0.025)
gmpreds$uci <- apply(sim, 1, quantile, 0.975)

# Calculate 95% confidence interval on difference between C and D corals at select dates (=sampling dates)
sim1 <- cbind(pdat, sim)
sim1 <- sim1[sim1$time %in% c(0,7,14,21,28,35,42,49,56,63), ]
sim1split <- split(sim1, f=interaction(sim1$time, sim1$ramp, drop=T))
CDdiff.95CI <- lapply(sim1split, function(x) quantile(apply(x[,4:ncol(x)], 2, diff), c(0.025, 0.975)))
# Test if 95% CI on difference between C and D corals contains zero
test0 <- lapply(CDdiff.95CI, function(x) x[1] < 0 & 0 < x[2])
test0[test0==F]  # Shows when C and D are significantly different with p<0.05

# Analyze total S/H ratio ------------
# Fit model for cooling treatment and test for effect of past thermal history
coolSHmod <- lmerTest::lmer(log10(tot.SH) ~ timef * dom * history + (1|mother/sample), 
                            data=subset(data, ramp=="cool"))
lmerTest::anova(coolSHmod)
coolSHlsm <- lsmeans::lsmeans(coolSHmod, specs=c("timef", "history", "dom"), by=c("dom", "timef"))
hist.pw <- data.frame(summary(pairs(coolSHlsm))); na.omit(hist.pw)
hist.pw[which(hist.pw$p.value < 0.01), ]  # No contrasts are significant at alpha=0.01
coolSHmod <- update(coolSHmod, formula=~timef * dom + (1|mother/sample)) # remove history from model
# Fit model for heating treatment and test for effect of past thermal history
heatSHmod <- lmerTest::lmer(log10(tot.SH) ~ timef * dom * history + (1|mother/sample), 
                            data=subset(data, ramp=="heat"))
lmerTest::anova(heatSHmod)
heatSHlsm <- lsmeans::lsmeans(heatSHmod, specs=c("timef", "history", "dom"), by=c("dom", "timef"))
hist.pw <- data.frame(summary(pairs(heatSHlsm))); na.omit(hist.pw)
hist.pw[which(hist.pw$p.value < 0.01), ]  # No contrasts are significant at alpha=0.01
heatSHmod <- update(heatSHmod, formula=~timef * dom + (1|mother/sample)) # remove history from model
# Compare S/H ratio within each group to time zero
cool.dunnett <- rbind(contrast(lsmeans::lsmeans(coolSHmod, specs=c("timef", "dom")), "dunnett", by="dom"))
heat.dunnett <- rbind(contrast(lsmeans::lsmeans(heatSHmod, specs=c("timef", "dom")), "dunnett", by="dom"))
# Calculate percent change relative to time zero and extract significant contrasts
cool.dunnett <- within(data.frame(summary(cool.dunnett)), {
  loss <- paste0(round(-(1-(10^estimate)) * 100, 1), "%")
  day <- substr(contrast, 1, 2)
})
cool.diff <- cool.dunnett[cool.dunnett$p.value < 0.01, ]
heat.dunnett <- within(data.frame(summary(heat.dunnett)), {
  loss <- paste0(round(-(1-(10^estimate)) * 100, 1), "%")
  day <- substr(contrast, 1, 2)
})
heat.diff <- heat.dunnett[heat.dunnett$p.value < 0.01, ]

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
title("A. Cooling", adj=0)
mtext(side=2, text="Relative Fv/Fm", cex=0.75, line=1.5)
with(subset(gmpreds, ramp=="cool" & dom=="C"), {
  addpoly(time, uci, lci, col=alpha("blue", 0.4))
  lines(time, fit)
})
with(subset(gmpreds, ramp=="cool" & dom=="D"), {
  addpoly(time, uci, lci, col=alpha("red", 0.4))
  lines(time, fit)
})
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
title("B. Heating", adj=0)
mtext(side=2, text="Relative Fv/Fm", cex=0.75, line=1.5)
with(subset(gmpreds, ramp=="heat" & dom=="C" & time <=42), {
  addpoly(time, uci, lci, col=alpha("blue", 0.4))
  lines(time, fit)
})
with(subset(gmpreds, ramp=="heat" & dom=="D" & time <=42), {
  addpoly(time, uci, lci, col=alpha("red", 0.4))
  lines(time, fit)
})
points(c(7,14,21), c(1.11,1.08,0.99), pch="*", cex=2, xpd=T)
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

# Analyze C and D abundances over time by treatment and dominant clade at start -----
# Clade C
df <- data[!is.na(data$tot.SH), ]
Cmod <- lmer(log10(C.SH) ~ timef * ramp * dom + (1|mother/sample), data=df)
Clsm <- lsmeans::lsmeans(Cmod, specs=c("timef", "ramp", "dom"))
Clsmeans <- summary(Clsm)
# Clade D
Dmod <- lmer(log10(D.SH) ~ timef * ramp * dom + (1|mother/sample), data=df)
Dlsm <- lsmeans::lsmeans(Dmod, specs=c("timef", "ramp", "dom"))
Dlsmeans <- summary(Dlsm)
# Combine into single data frame
df1 <- data.frame(time=Clsmeans$timef, ramp=Clsmeans$ramp, dom=Clsmeans$dom, Cmean=Clsmeans$lsmean, Dmean=Dlsmeans$lsmean, Cse=Clsmeans$SE, Dse=Dlsmeans$SE)

# Figure 2: Community dynamics under (A.) cooling and (B.) heating -----
# Plot mean abundances of C and D over time by treatment and dominant clade at start
#pdf(file="output/Figure2.pdf", width=6.85, height=3.425)
par(mfrow=c(1,2), mar=c(3,3,2,1), mgp=c(1.5,0.4,0), tcl=-0.3, xpd=F)
df <- subset(df1, ramp=="cool")
plot(NA, xlim=c(-6,0), ylim=c(-6,0), xlab="Clade D (log10 S/H)", ylab="Clade C (log10 S/H)", cex.axis=0.75, cex.lab=0.75)
abline(a=0,b=1,lty=2)
arrows(df$Dmean, df$Cmean+df$Cse, df$Dmean, df$Cmean-df$Cse, length=0.025, angle=90, code=3, col="black")
arrows(df$Dmean+df$Dse, df$Cmean, df$Dmean-df$Dse, df$Cmean, length=0.025, angle=90, code=3, col="black")
for (j in 2:nrow(df)) {
  if (as.numeric(df$time)[j] > as.numeric(df$time)[j-1]) {
    arrows(df$Dmean[j-1], df$Cmean[j-1], df$Dmean[j], df$Cmean[j], length=0.1, code=2, lwd=2,
           col=c("blue","red","green")[df$dom[j]])
  }
}
text(par("usr")[1], 0.5, labels=expression(bold("A. Cooling")), adj=0, xpd=NA)
legend("bottomright", legend=c("C-dominated", "D-dominated"), 
       lty=1, lwd=2, col=c("blue","red","green"), inset=0.05, cex=0.75, bty="n")
df <- subset(df1, ramp=="heat")

plot(NA, xlim=c(-6,0), ylim=c(-6,0), xlab="Clade D (log10 S/H)",ylab="Clade C (log10 S/H)", cex.axis=0.75, cex.lab=0.75)
abline(a=0,b=1,lty=2)
arrows(df$Dmean, df$Cmean+df$Cse, df$Dmean, df$Cmean-df$Cse, length=0.025, angle=90, code=3, col="black")
arrows(df$Dmean+df$Dse, df$Cmean, df$Dmean-df$Dse, df$Cmean, length=0.025, angle=90, code=3, col="black")
for (j in 2:nrow(df)) {
  if (as.numeric(df$time)[j] > as.numeric(df$time)[j-1]) {
    arrows(df$Dmean[j-1], df$Cmean[j-1], df$Dmean[j], df$Cmean[j], length=0.1, code=2, lwd=2,
           col=c("blue","red","green")[df$dom[j]])
  }
}
text(par("usr")[1], 0.5, labels=expression(bold("B. Heating")), adj=0, xpd=NA)
#dev.off()

# -----


mod <- gamm(rfvfm ~ s(time, by=history, k=3) + history, data=subset(data, ramp=="heat" & dom=="D"))
plot(mod)
summary(mod$gam)

mod <- lmerTest::lmer(rfvfm ~ poly(time, 2) * history + (1|mother/sample), data=subset(data, ramp=="heat" & dom=="D"))
plot(Effect(c("time", "history"), mod), x.var="time", multiline=T)
lmerTest::anova(mod, test="F")


