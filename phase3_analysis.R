library(effects)
library(lmer)
library(spida)
library(scales)

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
xyplot(propD ~ time | ramp + mother, groups= ~ sample, data = na.omit(data), type="o", lty=1)

# Categorize C- and D-dominated corals based on community composition at time zero
dom <- with(data[data$time==0, ], na.omit(data.frame(sample=sample, dom=ifelse(propD > 0.5, "D", "C"))))
# Some corals are missing qpcr data for time zero; for these, assign based on mean propD of all samples from each unassigned coral
dom2 <- with(data[!data$sample %in% dom$sample, ], na.omit(aggregate(propD, by=list(sample=sample), FUN=mean, na.rm=T)))
dom2$dom <- ifelse(dom2$x > 0.5, "D", "C")
dom <- merge(dom, dom2, all=T)[,-3]
# Merge dominant symbiont with data
data <- merge(data, dom, by="sample", all.x=T)
table(data[data$time==0, "dom"])  # 37 C-dominant and 113 D-dominant corals = 150 corals with qPCR data
data <- data[with(data, order(sample, time)), ]


# Fv/Fm ANALYSIS
# Visualize all data
xyplot(fvfm ~ time | ramp + dom, groups= ~ sample, data = data, type="o", lty=1)
# Fit mixed model
#sp <- function(x) gsp(x, knots=0, degree=2)
mod.all.full <- lmerTest::lmer(fvfm ~ poly(time, 2) * ramp * dom + (poly(time, 2)|mother/sample), data=data)
# Test significance of fixed effects by backwards selection
modselect <- lmerTest::step(mod.all.full, lsmeans.calc=F, difflsmeans.calc=F, alpha.fixed=0.05)
modselect$anova.table
# Rebuild model omitting non-significant fixed effects
mod.all <- update(mod.all.full, formula(modselect$model))
# Identify outliers with standardized residuals > 2.5
out <- abs(residuals(mod.all)) > sd(residuals(mod.all)) * 2.5
data[out, ]  # outlying data points
# Refit model without outliers
mod.all <- lmerTest::lmer(fvfm ~ poly(time, 2) * ramp * dom + (poly(time, 2)|mother/sample), data=data[!out, ])
# Print and save ANOVA table for model
anovatab <- lmerTest::anova(mod.all)
write.csv(round(anovatab, digits=3), file="output/Table1.csv")
# Pseudo-r2 value-- squared correlation between fitted and observed values
summary(lm(model.response(model.frame(mod.all)) ~ fitted(mod.all)))$r.squared
# Generate predictions and confidence intervals by parametric bootstrapping
pred.all <- expand.grid(time=seq(0,63,1), ramp=factor(c("cool", "heat")), dom=factor(c("C", "D")))
bootfit <- bootMer(mod.all, FUN=function(x) predict(x, pred.all, re.form=NA), nsim=10)
# Extract 90% confidence interval on predicted values
pred.all$fit <- predict(mod.all, pred.all, re.form=NA)
pred.all$lci <- apply(bootfit$t, 2, quantile, 0.05)
pred.all$uci <- apply(bootfit$t, 2, quantile, 0.95)
# Prepare data for plotting
datsumm <- data.frame(
  expand.grid(dom=levels(dat$dom), ramp=levels(dat$ramp), time=seq(0,63,7)),
  mean=aggregate(data$fvfm, by=list(interaction(data$dom, data$ramp, data$time)), FUN=mean, na.rm=T)$x,
  sd=aggregate(data$fvfm, by=list(interaction(data$dom, data$ramp, data$time)), FUN=sd, na.rm=T)$x,
  se=aggregate(data$fvfm, by=list(interaction(data$dom, data$ramp, data$time)), 
               FUN=function(x) sd(x, na.rm=T)/sqrt(length(na.omit(x))))$x,
  conf95=aggregate(data$fvfm, by=list(interaction(data$dom, data$ramp, data$time)), 
                   FUN=function(x) sd(x, na.rm=T)/sqrt(length(na.omit(x))) * qt(0.975, length(na.omit(x))-1))$x
  )
datsumm <- na.omit(datsumm)
datlist <- split(datsumm, f=datsumm$ramp)
datlist <- lapply(datlist, function(x) rev(split(x, f=x$dom)))
predlist <- split(pred.all, f=pred.all$ramp)
predlist <- lapply(predlist, function(x) rev(split(x, f=x$dom))) 
predlist$heat$C <- predlist$heat$C[predlist$heat$C$time <= 42, ]
predlist$heat$D <- predlist$heat$D[predlist$heat$D$time <= 42, ]
# Plot figure
par(mfrow=c(2,1), mar=c(2,3,1,1), mgp=c(1.75, 0.4, 0))
for (ramp in c("cool", "heat")) {
  with(datlist[[ramp]], {
    # Create plot frame for each reef
    plot(NA, xlim=c(0,63), ylim=c(0, 0.7), bty="n", tck=-0.03, ylab="fvfm")
    title(paste("", ramp), line=-0.9, adj=0, outer=F)
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


