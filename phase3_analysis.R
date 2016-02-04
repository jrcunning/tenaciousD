library(effects)
library(lme4)
library(spida)
library(scales)
library(lsmeans)

addpoly <- function(x, y1, y2, col=alpha("lightgrey", 0.8), ...){
  ii <- order(x)
  y1 <- y1[ii]
  y2 <- y2[ii]
  x <- x[ii]
  polygon(c(x, rev(x)), c(y1, rev(y2)), col=col, border=NA, ...)
}


data <- read.csv("phase3.csv")

# Calculate proportion D in each sample
data$propD <- data$D.SH / (data$C.SH + data$D.SH)
hist(data$propD)

# Plot changes in proportion D over time
library(lattice)
xyplot(propD ~ time | ramp, groups= ~ sample, data = na.omit(data), type="o", lty=1)
xyplot(asin(sqrt(propD)) ~ time | ramp, groups= ~ sample, data = na.omit(data), type="o", lty=1)
#xyplot(propD ~ time | ramp + mother, groups= ~ sample, data = na.omit(data), type="o", lty=1)

# # Categorize C- and D-dominated corals based on community composition at time zero
# dom <- with(data[data$time==0, ], na.omit(data.frame(sample=sample, dom=ifelse(propD > 0.5, "D", "C"))))
# # Some corals are missing qpcr data for time zero; for these, assign based on mean propD of all samples from each unassigned coral
# dom2 <- with(data[!data$sample %in% dom$sample, ], na.omit(aggregate(propD, by=list(sample=sample), FUN=mean, na.rm=T)))
# dom2$dom <- ifelse(dom2$x > 0.5, "D", "C")
# dom <- merge(dom, dom2, all=T)[,-3]

# Categorize C- and D-dominant corals based on >90% one or the other at t0, omitting mixtures
nlevels(data$sample)
dom <- with(data[data$time==0, ], 
            data.frame(sample=sample, propD=propD, 
                       dom=ifelse(propD > 0.9, "D", ifelse(propD < 0.1, "C", NA))))

nlevels(droplevels(dom[!is.na(dom$dom), "sample"]))
# 55 coral samples are not assigned as C or D corals.... OK?
# Alternative possibility: assign as C or D based on average propD < or > 0.5 over time?

# Merge dominant symbiont classification with rest of data
data <- merge(data, dom[,-2], by="sample", all.x=T)
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

# Fv/Fm ANALYSIS
#fdat <- data[which(!is.na(data$fvfm) & !is.na(data$dom)), ]
# Visualize all data
xyplot(fvfm ~ time | ramp + dom, groups= ~ sample, data = data, type="o", lty=1)
xyplot(rfvfm ~ time | ramp + dom, groups= ~ sample, data = data, type="o", lty=1)
xyplot(dfvfm ~ time | ramp + dom, groups= ~ sample, data = data, type="o", lty=1)

# Fit mixed model
mod.all <- lmerTest::lmer(fvfm ~ poly(time, 2) * ramp * dom + (poly(time, 2)|mother/sample), data=data)
mod.all <- lmerTest::lmer(dfvfm ~ poly(time, 2) * ramp * dom + (1|mother/sample), data=data)
mod.all <- lmerTest::lmer(rfvfm ~ poly(time, 2) * ramp * dom + (1|mother/sample), data=data)
anovatab <- lmerTest::anova(mod.all)
anovatab
#write.csv(round(anovatab, digits=3), file="output/Table1.csv")
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
# Plot figure
par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=c(1.2, 0.4, 0))
for (ramp in c("cool", "heat")) {
  with(datlist[[ramp]], {
    # Create plot frame for each reef
    plot(NA, xlim=c(0,63), ylim=c(0, 1.2), bty="n", tck=-0.03, ylab="fvfm", xlab="days")
    if (ramp=="cool") text(x=seq(0,63,7), y=1.1, labels = c(24, 23, 22, 21, 20, 19, 18, 17, 16, 15), xpd=T, cex=0.6)
    if (ramp=="heat") text(x=seq(0,63,7), y=1.1, labels = c(29, 30, 31, 32, 33, 34, 35, "", "", ""), xpd=T, cex=0.6)
    title(paste("", ramp), line=-0.9, adj=0, outer=F)
    legend("bottomleft", col=c("blue", "red"), legend=c("C", "D"), pch=21, lty=1, cex=0.75, bty="n")
    # Plot model fit line and shaded CI for C and D corals
    with(predlist[[ramp]], {
      lapply(predlist[[ramp]], function(dom) {
        addpoly(dom$time, dom$lci, dom$uci, col=alpha(list("C"="blue", "D"="red")[[dom$dom[1]]], 0.2), xpd=NA)
        lines(dom$time, dom$fit, lty=1)
      })
    })
    # Plot raw data +/- standard deviation
    lapply(datlist[[ramp]], function(dom) {
      arrows(dom$time, dom$mean + dom$sd, dom$time, dom$mean - dom$sd, code=3, angle=90, length=0.05, xpd=NA,
             col=list("C"="blue", "D"="red")[[dom$dom[1]]])
      points(dom$mean ~ dom$time, pch=21, bg=list("C"="blue", "D"="red")[[dom$dom[1]]], ylim=c(0, 1), cex=0.75)
    })
  })
}

# S/H ------------
# Calculate relative change in S/H
for (sample in data$sample) {
  sdata <- data[which(data$sample==sample), ]
  for (time in sdata$time) {
    sdata[which(sdata$time==time), "rSH"] <- sdata[which(sdata$time==time), "tot.SH"] / sdata[which(sdata$time==0), "tot.SH"]
    sdata[which(sdata$time==time), "plost"] <- 1 - sdata[which(sdata$time==time), "rSH"]
  }
  data[which(data$sample==sample), "rSH"] <- sdata[, "rSH"]
  data[which(data$sample==sample), "plost"] <- sdata[, "plost"]
}
head(data, 50)
boxplot(data$plost)
hist(log(data$rSH))

# TOTAL S/H ANALYSIS
# Visualize all data
xyplot(tot.SH ~ time | ramp + dom, groups= ~ sample, data = na.omit(data), type="o", lty=1)
xyplot(log10(tot.SH) ~ time | ramp + dom, groups= ~ sample, data = na.omit(data), type="o", lty=1)
xyplot(log(rSH) ~ time | ramp + dom, groups= ~ sample, data = na.omit(data), type="o", lty=1)
xyplot(rSH ~ time | ramp + dom, groups= ~ sample, data = na.omit(data), type="o", lty=1, ylim=c(0,1))

# # Plot geometric means
# shdat <- subset(data, !is.na(tot.SH))
# shdat
# gmeans <- aggregate(log(shdat$tot.SH), by=list(interaction(shdat$ramp, shdat$dom, shdat$time)), FUN=mean, na.rm=T)
# gmeans$bt <- exp(gmeans$x)
# gmeans
# gmeans <- split(gmeans, f=substr(gmeans$Group.1, 1, 6))
# plot(bt ~ seq(1,4), data=gmeans$cool.C, type="o")


# Fit mixed model
#sp <- function(x) gsp(x, knots=0, degree=2)
shmod.all.full <- lmerTest::lmer(log(tot.SH) ~ poly(time, 2) * ramp * dom + (1|mother/sample), data=data)
# model will not fit because there are samples with zero as total, cannot log transform
# change zeros to just below minimum detected symbionts
min(data$tot.SH[which(data$tot.SH!=0)], na.rm=T)  # min. is 1.08e-5, so change zeros to 1e-5
data[which(data$tot.SH==0), "tot.SH"] <- 1e-5
# refit model
xyplot(log10(tot.SH) ~ time | ramp + dom, groups= ~ sample, data = na.omit(data[,-10]), type="o", lty=1)
shmod.all.full <- lmerTest::lmer(log(tot.SH) ~ time * ramp * dom + (1|mother/sample), data=data)


# Test significance of fixed effects by backwards selection
shmodselect <- lmerTest::step(shmod.all.full, lsmeans.calc=F, difflsmeans.calc=F, alpha.fixed=0.05)
shmodselect$anova.table
# Rebuild model omitting non-significant fixed effects
shmod.all <- update(shmod.all.full, formula(shmodselect$model))
# Identify outliers with standardized residuals > 2.5
#shout <- abs(residuals(shmod.all)) > sd(residuals(shmod.all)) * 2.5
#data[shout, ]  # outlying data points
# Refit model without outliers
#shmod.all <- lmerTest::lmer(log(tot.SH) ~ poly(time, 2) * ramp * dom + (1|mother/sample), data=data[!shout, ])
# Print and save ANOVA table for model
shanovatab <- lmerTest::anova(shmod.all)
shanovatab
#write.csv(round(shanovatab, digits=3), file="output/Table2.csv")
# Pseudo-r2 value-- squared correlation between fitted and observed values
summary(lm(model.response(model.frame(shmod.all)) ~ fitted(shmod.all)))$r.squared
# Generate predictions and confidence intervals by parametric bootstrapping
shpred.all <- expand.grid(time=seq(0,63,1), ramp=factor(c("cool", "heat")), dom=factor(c("C", "D")))
shbootfit <- bootMer(shmod.all, FUN=function(x) predict(x, shpred.all, re.form=NA), nsim=10)
# Extract 95% confidence interval on predicted values
shpred.all$fit <- predict(shmod.all, shpred.all, re.form=NA)
shpred.all$lci <- apply(shbootfit$t, 2, quantile, 0.025)
shpred.all$uci <- apply(shbootfit$t, 2, quantile, 0.975)
# Prepare data for plotting
shdatsumm <- data.frame(
  colsplit(as.character(levels(droplevels(interaction(data$dom, data$ramp, data$time)))),
           pattern="\\.", names=c("dom", "ramp", "time")),
  mean=aggregate(log(data$tot.SH), by=list(interaction(data$dom, data$ramp, data$time)), FUN=mean, na.rm=T)$x,
  sd=aggregate(log(data$tot.SH), by=list(interaction(data$dom, data$ramp, data$time)), FUN=sd, na.rm=T)$x,
  se=aggregate(log(data$tot.SH), by=list(interaction(data$dom, data$ramp, data$time)), 
               FUN=function(x) sd(x, na.rm=T)/sqrt(length(na.omit(x))))$x,
  conf95=aggregate(log(data$tot.SH), by=list(interaction(data$dom, data$ramp, data$time)), 
                   FUN=function(x) sd(x, na.rm=T)/sqrt(length(na.omit(x))) * qt(0.975, length(na.omit(x))-1))$x
)
shdatsumm <- na.omit(shdatsumm)
shdatlist <- split(shdatsumm, f=shdatsumm$ramp)
shdatlist <- lapply(shdatlist, function(x) rev(split(x, f=x$dom)))
shpredlist <- split(shpred.all, f=shpred.all$ramp)
shpredlist <- lapply(shpredlist, function(x) rev(split(x, f=x$dom))) 
shpredlist$heat$C <- shpredlist$heat$C[shpredlist$heat$C$time <= 42, ]
shpredlist$heat$D <- shpredlist$heat$D[shpredlist$heat$D$time <= 42, ]
# Plot figure
par(mfrow=c(2,1), mar=c(2,3,1,1), mgp=c(1.75, 0.4, 0))
for (ramp in c("cool", "heat")) {
  with(shdatlist[[ramp]], {
    # Create plot frame for each reef
    plot(NA, xlim=c(0,63), ylim=c(-12, 3), bty="n", tck=-0.03, ylab="log tot.SH", xlab="days")
    if (ramp=="cool") text(x=seq(0,63,7), y=0.6, labels = c(24, 23, 22, 21, 20, 19, 18, 17, 16, 15), xpd=T, cex=0.6)
    if (ramp=="heat") text(x=seq(0,63,7), y=0.6, labels = c(29, 30, 31, 32, 33, 34, 35, "", "", ""), xpd=T, cex=0.6)
    title(paste("", ramp), line=-0.9, adj=0, outer=F)
    legend("bottomleft", col=c("blue", "red"), legend=c("C", "D"), pch=21, lty=1, cex=0.75, bty="n")
    # Plot model fit line and shaded CI for C and D corals
    with(shpredlist[[ramp]], {
      lapply(shpredlist[[ramp]], function(dom) {
        addpoly(dom$time, dom$lci, dom$uci, col=alpha(list("C"="blue", "D"="red")[[dom$dom[1]]], 0.2), xpd=NA)
        lines(dom$time, dom$fit, lty=1)
      })
    })
    # Plot raw data +/- standard deviation
    lapply(shdatlist[[ramp]], function(dom) {
      arrows(dom$time, dom$mean + dom$sd, dom$time, dom$mean - dom$sd, code=3, angle=90, length=0.05, xpd=NA,
             col=list("C"="blue", "D"="red")[[dom$dom[1]]])
      points(dom$mean ~ dom$time, pch=21, bg=list("C"="blue", "D"="red")[[dom$dom[1]]], ylim=c(0, 1), cex=0.75)
    })
  })
}

# Convert model predictions to relative values and percent loss
for (ramp in c("cool", "heat")) {
  for (dom in c("C", "D")) {
    for (time in shpredlist[[ramp]][[dom]]$time) {
      #df <- shpredlist[[ramp]][[dom]]
      shpredlist[[ramp]][[dom]][which(shpredlist[[ramp]][[dom]]$time==time), "rfit"] <- exp(shpredlist[[ramp]][[dom]][which(shpredlist[[ramp]][[dom]]$time==time), "fit"]) / exp(shpredlist[[ramp]][[dom]][which(shpredlist[[ramp]][[dom]]$time=="0"), "fit"])
    }
  }
}

# Plot figure with percent remaining
par(mfrow=c(2,1), mar=c(2,3,1,1), mgp=c(1.75, 0.4, 0))
for (ramp in c("cool", "heat")) {
  with(shdatlist[[ramp]], {
    # Create plot frame for each reef
    plot(NA, xlim=c(0,63), ylim=c(0, 2.5), bty="n", tck=-0.03, ylab="log tot.SH", xlab="days")
    if (ramp=="cool") text(x=seq(0,63,7), y=0.6, labels = c(24, 23, 22, 21, 20, 19, 18, 17, 16, 15), xpd=T, cex=0.6)
    if (ramp=="heat") text(x=seq(0,63,7), y=0.6, labels = c(29, 30, 31, 32, 33, 34, 35, "", "", ""), xpd=T, cex=0.6)
    title(paste("", ramp), line=-0.9, adj=0, outer=F)
    legend("bottomleft", col=c("blue", "red"), legend=c("C", "D"), pch=21, lty=1, cex=0.75, bty="n")
    # Plot model fit line and shaded CI for C and D corals
    with(shpredlist[[ramp]], {
      lapply(shpredlist[[ramp]], function(dom) {
        #addpoly(dom$time, dom$lci, dom$uci, col=alpha(list("C"="blue", "D"="red")[[dom$dom[1]]], 0.2), xpd=NA)
        lines(dom$time, dom$rfit, lty=1)
      })
    })
#     # Plot raw data +/- standard deviation
#     lapply(shdatlist[[ramp]], function(dom) {
#       arrows(dom$time, dom$mean + dom$sd, dom$time, dom$mean - dom$sd, code=3, angle=90, length=0.05, xpd=NA,
#              col=list("C"="blue", "D"="red")[[dom$dom[1]]])
#       points(dom$mean ~ dom$time, pch=21, bg=list("C"="blue", "D"="red")[[dom$dom[1]]], ylim=c(0, 1), cex=0.75)
#     })
  })
}

# analyse SH ratios with time as categorical, not continuous, because it's not really a good time series.
# i.e., not enough temporal resolution to be confident in quadratic effect of time....
# get SH data frame and plot raw data
shdf <- droplevels(subset(data, !is.na(tot.SH)))
xyplot(log(tot.SH) ~ time | ramp + dom, groups= ~ sample, data = shdf, type="o", lty=1)
# create list split by dom&ramp, fit mixed model for each group
shl <- split(shdf, f=interaction(shdf$ramp, shdf$dom))
mods <- llply(shl, function(df) lmerTest::lmer(log(tot.SH) ~ factor(time) + (1|mother/sample), data=df))

# run anova and post-hoc comparisons over time for each group
anovas <- llply(mods, anova)
lsm <- llply(mods, function(mod) lsmeans(mod, specs="time", contr="dunnett"))
# plot lsm -- are these exactly same as geomeans of raw data?
# get contrasts for comparisons
contr <- ldply(sapply(lsm, "[[", 2), summary)
contr$bt <- exp(contr$estimate)
contr$loss <- -(1-contr$bt)
contr
# get lsmeans from mixed model (slightly diff than raw data bc of random effects)
lsmeans <- ldply(sapply(lsm, "[[", 1), summary)
lsmeans
# get means from raw data
geomeans <- ldply(shl, function(df) aggregate(list(mean=log(df$tot.SH)), by=list(time=df$time), FUN=mean))
geosds <- ldply(shl, function(df) aggregate(list(sd=log(df$tot.SH)), by=list(time=df$time), FUN=sd))
geostats <- merge(geomeans, geosds, by=c(".id", "time"))
geostats <- split(geostats, f=geostats$.id)
geostats
# plot raw data + sd
par(mfrow=c(1,2))
plot(NA, xlim=c(0,63), ylim=c(-12, 1), bty="n", tck=-0.03, ylab="log tot.SH", xlab="days")
with(geostats$cool.C, {
  lines(time, mean, type="o", col="blue", lty=2, pch=21, bg="blue")
  arrows(time, mean-sd, time, mean+sd, code=3, angle=90, length=0.05, col="blue")})
with(geostats$cool.D, {
  lines(time, mean, type="o", col="red", lty=2, pch=21, bg="red")
  arrows(time, mean-sd, time, mean+sd, code=3, angle=90, length=0.05, col="red")})
plot(NA, xlim=c(0,63), ylim=c(-12, 1), bty="n", tck=-0.03, ylab="log tot.SH", xlab="days")
with(geostats$heat.C, {
  lines(time, mean, type="o", col="blue", lty=2, pch=21, bg="blue")
  arrows(time, mean-sd, time, mean+sd, code=3, angle=90, length=0.05, col="blue")})
with(geostats$heat.D, {
  lines(time, mean, type="o", col="red", lty=2, pch=21, bg="red")
  arrows(time, mean-sd, time, mean+sd, code=3, angle=90, length=0.05, col="red")})
# plot contrasts from lsmeans + SE
contrl <- split(contr, f=contr$.id)
plot(NA, xlim=c(0,63), ylim=c(-10, 2), bty="n", tck=-0.03, ylab="log tot.SH", xlab="days")
lapply(contrl[c(1,2)], function(x) lines(c(0, as.numeric(substr(x$contrast, 1, 2))),
                                         c(0, x$estimate), type="o"))
lapply(contrl[c(1,2)], function(x) arrows(as.numeric(substr(x$contrast, 1, 2)), x$estimate-x$SE, 
                                          as.numeric(substr(x$contrast, 1, 2)), x$estimate+x$SE,
                                          code=3, angle=90, length=0.05))
plot(NA, xlim=c(0,63), ylim=c(-10, 2), bty="n", tck=-0.03, ylab="log tot.SH", xlab="days")
lapply(contrl[c(3,4)], function(x) lines(c(0, as.numeric(substr(x$contrast, 1, 2))),
                                         c(0, x$estimate), type="o"))
lapply(contrl[c(3,4)], function(x) arrows(as.numeric(substr(x$contrast, 1, 2)), x$estimate-x$SE, 
                                          as.numeric(substr(x$contrast, 1, 2)), x$estimate+x$SE,
                                          code=3, angle=90, length=0.05))
contrl





# calculate lsmeans relative to group t0 (% remaining)
for (dom in c("C", "D")) {
  for (time in c(0,28,42,63)) {
    lsm[which(lsm$dom==dom & lsm$time==time), "prem"] <- exp(lsm[which(lsm$dom==dom & lsm$time==time), "lsmean"]) / exp(lsm[which(lsm$dom==dom & lsm$time==0), "lsmean"])
  }
}
lsm
# backtransform CL relative to t0
for (dom in c("C", "D")) {
  for (time in c(0,28,42,63)) {
    lsm[which(lsm$dom==dom & lsm$time==time), "lower"] <- exp(lsm[which(lsm$dom==dom & lsm$time==time), "lower.CL"]) / exp(lsm[which(lsm$dom==dom & lsm$time==0), "lsmean"])
    lsm[which(lsm$dom==dom & lsm$time==time), "upper"] <- exp(lsm[which(lsm$dom==dom & lsm$time==time), "upper.CL"]) / exp(lsm[which(lsm$dom==dom & lsm$time==0), "lsmean"])
  }
}
lsm
bars <- barplot(lsm$prem)
arrows(bars, lsm$lower, bars, lsm$upper)
arrows()
lsm
barplot(lsm$lsmean)

?ref.grid
na.omit(heat)

?lsmeans





















# FvFm vs. propD over time
graphics.off()
plot(fvfm ~ propD, data)
modc <- lmer(fvfm ~ propD * time + (time|mother/sample), data[which(data$ramp=="cool"), ])
plot(Effect(c("propD", "time"), mod), multiline=T, x.var="propD")
modh <- lmer(fvfm ~ propD * time + (time|mother/sample), data[which(data$ramp=="heat"), ])
plot(Effect(c("propD", "time"), mod), multiline=T, x.var="propD")
pdf(file = "fvfmpropD.pdf")
par(mfrow=c(2,1))
layout(mat=matrix(c(1,1,2,2), 2, 2, byrow = TRUE))
plot(Effect(c("propD", "time"), modc), multiline=T, x.var="propD", ylim=c(0,0.6))
plot(Effect(c("propD", "time"), modh), multiline=T, x.var="propD", ylim=c(0,0.6))
dev.off()

# identify some corals that switched/shuffled during exposure
hist(aggregate(data$propD, by=list(sample=data$sample), FUN=sd, na.rm=T)$x)
propDsd <- aggregate(data$propD, by=list(sample=data$sample), FUN=sd, na.rm=T)
shuff <- propDsd[which(propDsd$x > 0.1), "sample"]
shuff
shuff <- data[data$sample %in% shuff, ]
shuff
plot(fvfm ~ propD, shuff)

propDmod <- lmer(asin(sqrt(propD)) ~ ramp * dom * time + (time|mother/sample), data=data)
ind <- ranef(propDmod)
hist(ind$`sample:mother`$time)
ind$`sample:mother`[which(ind$`sample:mother`$time > )]

# divide heated C-corals into ones that switched to D and ones that didn't
Cheat <- na.omit(subset(data, data$ramp=="heat" & data$dom=="C"))
length(unique(Cheat$sample))  # 17 heated samples started C-dominant
xyplot(propD ~ time | ramp, groups= ~ sample, data = na.omit(Cheat), type="o", lty=1)
# i think the switches from D dom back to C dom between last two time points are artifacts from 
# putting in placeholder data when no symbionts at all were detected with qPCR

b <- Cheat[which(Cheat$time==28), c("sample", "propD")]
b$CtoD <- ifelse(b$propD > 0.5, T, F)

div <- merge(Cheat, b[,c(1,3)], by="sample")
div <- div[with(div, order(sample, time)), ]
div
# plot fvfm for CtoD true vs false corals
xyplot(fvfm ~ time | CtoD, groups= ~ sample, data = div, type="o", lty=1)

# ANIMATION OF FVFM AND PROPD OVER TIME
for time()

plot(fvfm ~ propD, data=subset(data, ramp=="cool" & time==0), ylim=c(0,0.6))
plot(fvfm ~ propD, data=subset(data, ramp=="cool" & time==28), ylim=c(0,0.6))
segments(x0=subset(data, ramp=="cool" & time==0)[,"propD"],
         y0=subset(data, ramp=="cool" & time==0)[,"fvfm"],
         x1=subset(data, ramp=="cool" & time==28)[,"propD"],
         y1=subset(data, ramp=="cool" & time==28)[,"fvfm"])
plot(fvfm ~ propD, data=subset(data, ramp=="cool" & time==28), ylim=c(0,0.6))
plot(fvfm ~ propD, data=subset(data, ramp=="cool" & time==42), ylim=c(0,0.6))
segments(x0=subset(data, ramp=="cool" & time==28)[,"propD"],
         y0=subset(data, ramp=="cool" & time==28)[,"fvfm"],
         x1=subset(data, ramp=="cool" & time==42)[,"propD"],
         y1=subset(data, ramp=="cool" & time==42)[,"fvfm"])
plot(fvfm ~ propD, data=subset(data, ramp=="cool" & time==42), ylim=c(0,0.6))
plot(fvfm ~ propD, data=subset(data, ramp=="cool" & time==63), ylim=c(0,0.6))
segments(x0=subset(data, ramp=="cool" & time==42)[,"propD"],
         y0=subset(data, ramp=="cool" & time==42)[,"fvfm"],
         x1=subset(data, ramp=="cool" & time==63)[,"propD"],
         y1=subset(data, ramp=="cool" & time==63)[,"fvfm"])
plot(fvfm ~ propD, data=subset(data, ramp=="cool" & time==63), ylim=c(0,0.6))

plot(fvfm ~ propD, data=subset(data, ramp=="heat" & time==0), ylim=c(0,0.6))
plot(fvfm ~ propD, data=subset(data, ramp=="heat" & time==28), ylim=c(0,0.6))
segments(x0=subset(data, ramp=="heat" & time==0)[,"propD"],
         y0=subset(data, ramp=="heat" & time==0)[,"fvfm"],
         x1=subset(data, ramp=="heat" & time==28)[,"propD"],
         y1=subset(data, ramp=="heat" & time==28)[,"fvfm"])
plot(fvfm ~ propD, data=subset(data, ramp=="heat" & time==28), ylim=c(0,0.6))
plot(fvfm ~ propD, data=subset(data, ramp=="heat" & time==42), ylim=c(0,0.6))
segments(x0=subset(data, ramp=="heat" & time==28)[,"propD"],
         y0=subset(data, ramp=="heat" & time==28)[,"fvfm"],
         x1=subset(data, ramp=="heat" & time==42)[,"propD"],
         y1=subset(data, ramp=="heat" & time==42)[,"fvfm"])
plot(fvfm ~ propD, data=subset(data, ramp=="heat" & time==42), ylim=c(0,0.6))


