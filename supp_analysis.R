source('setup.R')

# To evaluate the potential influence of previous bleaching history on the response of Fv/Fm and 
# total S/H over time, we tested for a significant symbiont*time interaction in corals with the same 
# history but different dominant symbionts (highly significant), and for a history*time interaction 
# in corals with the same dominant symbiont but different histories (not significant; see Supplementary 
# Information). While these tests were only possible for a small subset of corals (since most D corals 
# had previously bleached and most C corals had not), they support dominant symbiont rather than 
# bleaching history as the primary driver of coral response.

# Test for clade*time effect in corals that had previously bleached in cooling treatment
table(data$prevBleach, data$dom, data$ramp)
df <- subset(data, prevBleach=="B" & ramp=="cool")
table(unique(df[,c("sample", "dom")])$dom)  # 7 C-dom, 35 D-dom
mod <- lmerTest::lmer(rfvfm ~ timef * dom + (1|mother/sample), data=df, na.action=na.omit)
lmerTest::anova(mod)  # clade*time significant
mod <- lmerTest::lmer(log(tot.SH) ~ timef * dom + (1|mother/sample), data=df, na.action=na.omit)
lmerTest::anova(mod)  # clade*time significant  # CLADE IS IMPORTANT

# Test for clade*time effect in corals that had not previously bleached in heating treatment
table(data$prevBleach, data$dom, data$ramp)
df <- subset(data, prevBleach=="NB" & ramp=="heat")
table(unique(df[,c("sample", "dom")])$dom)  # 19 C-dom, 4 D-dom
mod <- lmerTest::lmer(rfvfm ~ timef * dom + (1|mother/sample), data=df, na.action=na.omit)
lmerTest::anova(mod)  # clade*time not significant
mod <- lmerTest::lmer(log(tot.SH) ~ timef * dom + (1|mother/sample), data=df, na.action=na.omit)
lmerTest::anova(mod)  # clade*time significant  # CLADE IS IMPORTANT



# Test for prevBleach*time effect in corals with clade C in cooling treatment
table(data$prevBleach, data$dom, data$ramp)
df <- subset(data, dom=="C" & ramp=="cool")
table(unique(df[,c("sample", "prevBleach")])$prevBleach)  # 7 B, 12 NB
mod <- lmerTest::lmer(rfvfm ~ timef * prevBleach + (1|mother/sample), data=df, na.action=na.omit)
lmerTest::anova(mod)  # prevBleach*time not significant
mod <- lmerTest::lmer(log(tot.SH) ~ timef * prevBleach + (1|mother/sample), data=df, na.action=na.omit)
lmerTest::anova(mod)  # prevBleach*time not significant  # PREVBLEACH IS NOT IMPORTANT

# Test for prevBleach*time effect in corals with clade D in heating treatment
table(data$prevBleach, data$dom, data$ramp)
df <- subset(data, dom=="D" & ramp=="heat")
table(unique(df[,c("sample", "prevBleach")])$prevBleach)  # 55 B, 4 NB
mod <- lmerTest::lmer(rfvfm ~ timef * prevBleach + (1|mother/sample), data=df, na.action=na.omit)
lmerTest::anova(mod)  # prevBleach*time is significant
mod <- lmerTest::lmer(log(tot.SH) ~ timef * prevBleach + (1|mother/sample), data=df, na.action=na.omit)
lmerTest::anova(mod)  # prevBleach*time not significant  # PREVBLEACH IS NOT IMPORTANT




