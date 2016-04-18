source('setup.R')

# Test for effect of clade in corals with same history -----
# Do any corals in either treatment have same history but different dominant clade?
table(data$history, data$dom, data$ramp)  # Yes: history B' in heating and E' in cooling
# History E' (DCMU-24-ctrl-24) in cooling treatment
Edf <- subset(data, history=="E'" & ramp=="cool")
table(unique(Edf[,c("sample", "dom")])$dom)  # 7 C-dom, 4 D-dom
mod <- lmerTest::lmer(rfvfm ~ timef * dom + (1|mother/sample), data=Edf, na.action=na.omit)
lmerTest::anova(mod)  # clade does not impact Fv/Fm
mod <- lmerTest::lmer(log(tot.SH) ~ timef * dom + (1|mother/sample), data=Edf, na.action=na.omit)
lmerTest::anova(mod)  # clade marginally impacts totSH
# History B' (ctrl-29-ctrl-29) in heating treatment
Bdf <- subset(data, history=="B'" & ramp=="heat")
table(unique(Bdf[,c("sample", "dom")])$dom)  # 6 C-dom, 4 D-dom
mod <- lmerTest::lmer(rfvfm ~ timef * dom + (1|mother/sample), data=Bdf, na.action=na.omit)
lmerTest::anova(mod) # clade does not impact Fv/Fm
mod <- lmerTest::lmer(log(tot.SH) ~ timef * dom + (1|mother/sample), data=Bdf, na.action=na.omit)
lmerTest::anova(mod) # clade strongly impacts totSH
plot(Effect(c("timef", "dom"), mod))  

# Test for effect of history in corals with same clade -----
# Do any corals have the same clade but different history?
table(data$dom, data$history, data$ramp) #Yes, cooling Cdom, A' (ctrl-24-ctrl-24) vs. E' (DCMU-24-ctrl-24)
df <- droplevels(subset(data, ramp=="cool" & dom=="C"))
table(unique(df[,c("sample","history")])$history) # 3 corals A', 7 corals E'
mod <- lmerTest::lmer(rfvfm ~ timef * history + (1|mother/sample), data=df, na.action=na.omit)
lmerTest::anova(mod)  # no impact of history on fv/fm
mod <- lmerTest::lmer(log(tot.SH) ~ timef * history + (1|mother/sample), data=df, na.action=na.omit)
lmerTest::anova(mod)  # no impact of history on totSH
