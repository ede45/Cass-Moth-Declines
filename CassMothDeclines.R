# 4 July 2022
# by Elizabeth de Jongh
# R script for "New Zealand endemic moth fauna declines 72% in abundance over 60 years despite stable land use"

##### Research aims #####

#### Measure changes in total moth catch from 1961-2021
#### Measure changes in moth catch for macro-moths and micro-moths from 1961-2021
#### Test whether more common or uncommon moth species had larger declines in abundance (abundance class analysis)
#### Measure changes moth community composition across sample seasons and decades (ordination)
#### Measure changes in average temperature near Cass
#### Measure changes in vegetation cover from 1990-2022

##### Import data #####

## data for total moth, macro-moth, micro-moth, and species catches
moth.dat <- read.csv("moth.dat.csv", header = T, sep = ",", fileEncoding = "UTF-8-BOM")

moth.dat$decade <- as.factor(moth.dat$decade)
moth.dat$season <- as.factor(moth.dat$season)
moth.dat$month <- factor(moth.dat$month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Oct", "Nov", "Dec"))
moth.dat$site2 <- factor(moth.dat$site, levels = c("S", "R"))

moth.dat$obs.effect <- as.factor(1:nrow(moth.dat)) # create object to account for overdispersion

# data for graphing 
# (file containing final model back transformed glmm coefficients and SE for total moths, macro-moths, and micro-moths across decades)
glmm.coef.decade <- read.csv("glmm.decade.csv", header = T, sep = ",", fileEncoding = "UTF-8-BOM")
glmm.coef.decade$decade <- as.factor(glmm.coef.decade$decade)

## data for abundance class analysis
abun.class.dat <- read.csv("abun.class.dat.csv", header = T, sep = ",", fileEncoding = "UTF-8-BOM")

abun.class.dat$decade2 <- as.factor(abun.class.dat$decade)
abun.class.dat$abun.class.log <- as.factor(abun.class.dat$abun.class.log)

# data for graphing 
# (file containing final model back transformed glmm coefficients and SE for abundance classes across decades)
glmm.coef.abun.class <- read.csv("glmm.abun.class.csv", header = T, sep = ",", fileEncoding = "UTF-8-BOM")

glmm.coef.abun.class$abun.class.log <- as.factor(glmm.coef.abun.class$abun.class.log)

## data for ordination

# data file where Nov 2022 is lumped together with 2020 season and 1989 season is removed
ord.dat <- read.csv("Edj_moth_ord2.csv", header = TRUE, sep = ",", fileEncoding = "UTF-8-BOM")

# data file with species of interest for plot (spp with at least 100 individuals caught in total)
interest.dat <- read.csv("interest.spp.csv", header = TRUE, sep = ",", fileEncoding = "UTF-8-BOM")

# data file with species not of interest for plot (all other spp)
not.interest.dat <- read.csv("not.interest.spp.csv", header = TRUE, sep = ",", fileEncoding = "UTF-8-BOM")

## estimated Cass annual temps 1964 - 2021
temp.dat <- read.csv("cass.annual.temp.dat.csv", header = T, sep = ",", fileEncoding = "UTF-8-BOM")

# data for graphing (estimated Cass annual temps 1964 - 2021)
temp.dat2 <- read.csv("cass.annual.temp.dat2.csv", header = T, sep = ",", fileEncoding = "UTF-8-BOM")


### theme for ggplot

library(ggplot2)
library(ggprism)
library(patchwork)

theme <- theme(panel.background = element_blank(), 
               axis.line = element_line(colour = "black", size = 1), legend.key = element_blank(), 
               axis.text = element_text(colour = "black", size = 15), 
               axis.title = element_text(colour = "black", size = 17), 
               axis.ticks = element_line(size = 1), 
               legend.text = element_text(size = 13), 
               legend.title = element_text(size = 13))




##### Measure changes in total moth catch from 1961-2021 #####

#### model selection

## poisson GLMM with season within site as random term
# fixed effects = site as factor, light type, month as factor, decade as factor
# random effects = season within site, observation effect
# offset = total trapping time out of 3 hr to account for diff trapping times

library(glmmTMB)

moth.glmm1 <- glmmTMB(total.moth ~ site2 + light + month + decade + (1 | season.site) + (1 | obs.effect) + 
                        offset(log(ttime.out.of.3hr)), family = "poisson", data = moth.dat)

glmmTMB:::Anova.glmmTMB(moth.glmm1)

summary(moth.glmm1)

# post hoc comparison test (to make pairwise comparisons between decades)

library(emmeans)

emmeans(moth.glmm1, list(pairwise ~ decade))

## poisson GLMM with nested season within decade as random effect
# fixed effects = site as factor, light type, month as factor, decade as factor
# random effects = season nested within decade, observation effect
# offset = total trapping time out of 3 hr to account for diff trapping times

moth.glmm2 <- glmmTMB(total.moth ~ site2 + light + month + decade + (1 | decade/season) + (1 | obs.effect) 
                      + offset(log(ttime.out.of.3hr)), family = "poisson", data = moth.dat)

glmmTMB:::Anova.glmmTMB(moth.glmm2)

summary(moth.glmm2)

# post hoc comparison test 

emmeans(moth.glmm2, list(pairwise ~ decade))

### same as moth.glmm1 but with negative binomial distribution (FINAL MODEL)

moth.glmm3 <- glmmTMB(total.moth ~ site2 + light + month + decade + (1 | season.site) + (1 | obs.effect) + 
                        offset(log(ttime.out.of.3hr)), family = "nbinom2", data = moth.dat)

glmmTMB:::Anova.glmmTMB(moth.glmm3)

summary(moth.glmm3)

# post hoc comparison test 

emmeans(moth.glmm3, list(pairwise ~ decade))

### compare AIC of glmm1, 2, and 3

AIC(moth.glmm1, moth.glmm2, moth.glmm3)

### compare BIC of glmm1, 2, and 3

BIC(moth.glmm1, moth.glmm2, moth.glmm3)

# graph fitted means and se (log10)
# with lines connecting points

glmm3.coef.decade <- subset(glmm.coef.decade, moth.type == "Total moths") # subset data

ggplot(data = glmm3.coef.decade, aes(x = decade.q, y = log10.est)) + 
  geom_point(size = 5) + 
  geom_errorbar(aes(ymin = log10.lower, ymax = log10.upper), width = 1.5, size = 1.1) + 
  geom_line(size = 1.25) + 
  theme + 
  scale_x_continuous(breaks = c(1962, 1988, 2021), labels = c("1961-62", "1987-89", "2020-2021")) + 
  labs(x = "Decade", y = "Log10(mean moth catch per 3 hours)")




##### Measure changes in moth catch for macro-moths and micro-moths from 1961-2021 #####

#### macro-moths

## negbin GLMM with season within site as random term (FINAL MODEL)
# same as final model for total moths

macro.glmm <- glmmTMB(macromoth ~ site2 + light + month + decade + (1 | season.site) + (1 | obs.effect) + 
                        offset(log(ttime.out.of.3hr)), family = "nbinom2", data = moth.dat)

glmmTMB:::Anova.glmmTMB(macro.glmm)

summary(macro.glmm)

# post hoc comparison test 

emmeans(macro.glmm, list(pairwise ~ decade))

#### micro-moths

## negbin GLMM with season within site as random term (FINAL MODEL)
# same analysis as total moths and macro-moths

micro.glmm <- glmmTMB(micromoth ~ site2 + light + month + decade + (1 | season.site) + (1 | obs.effect) + 
                        offset(log(ttime.out.of.3hr)), family = "nbinom2", data = moth.dat)

glmmTMB:::Anova.glmmTMB(micro.glmm)

summary(micro.glmm)

# post hoc comparison test 

emmeans(micro.glmm, list(pairwise ~ decade))

#### graph fitted means and se (log10) for macro- and micro-moths
# with line connecting points

macromicro.coef.decade2 <- subset(glmm.coef.decade, 
                                  moth.type != "Total moths")

ggplot(data = macromicro.coef.decade2, aes(x = decade.q, y = log10.est, shape = moth.type)) + 
  geom_point(size = 7) + 
  geom_errorbar(aes(ymin = log10.lower, ymax = log10.upper), width = 1.5, size = 1.1) + 
  geom_line(aes(x = decade.q, y = log10.est), size = 1.25) + 
  scale_x_continuous(breaks = c(1962, 1988, 2021), labels = c("1961-62", "1987-89", "2020-2021")) + 
  labs(x = "Decade", y = "Log10(mean moth catch per 3 hours)") + 
  theme(legend.title = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 1), legend.key = element_blank(), 
        axis.text = element_text(colour = "black", size = 15), 
        axis.title = element_text(colour = "black", size = 17), 
        axis.ticks = element_line(size = 1), 
        legend.text = element_text(size = 13), 
        legend.position = c(0.8, 0.9))




##### Test whether more common or uncommon moth species had larger declines in abundance (abundance class analysis) #####
### replicates = species

### negbin GLMM with season.site as random effect to account for nested data structure (FINAL MODEL)
# fixed effects = decade*abundance class interaction
# random effects = season within site, observation effect
# offset = total trapping time out of 3 hr to account for diff trapping times

abun.class.dat$obs.effect <- as.factor(1:nrow(abun.class.dat)) # create object to account for overdispersion

abun.class.glmm1 <- glmmTMB(count ~ decade2*abun.class.log + (1 | season.site) + (1 | obs.effect) + 
                              offset(log(num.nights.per.season)), data = abun.class.dat, family = "nbinom2")

glmmTMB:::Anova.glmmTMB(abun.class.glmm1)

summary(abun.class.glmm1)

## post hoc comparison test 

emmeans(abun.class.glmm1, list(pairwise ~ decade2 | abun.class.log))

# graph fitted means and se (log 10)

ggplot(data = glmm.coef.abun.class, aes(x = decade.q, y = log10.est, shape = abun.class.log)) + 
  geom_point(size = 6) + 
  scale_shape_manual(values=c(17, 18, 16, 15)) + 
  geom_line(size = 1) + 
  geom_errorbar(aes(ymin = log10.lower, ymax = log10.upper), width = 1.5, size = 1.1) + 
  scale_x_continuous(breaks = c(1962, 1988, 2021), labels = c("1961-62", "1987-89", "2020-2021")) + 
  labs(x = "Decade", y = "Log10(mean moth catch per species per night)", colour = "Abundance class") + 
  theme(legend.title = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 1), legend.key = element_blank(), 
        axis.text = element_text(colour = "black", size = 15), 
        axis.title = element_text(colour = "black", size = 17), 
        axis.ticks = element_line(size = 1), 
        legend.text = element_text(size = 13), 
        legend.position = c(0.8, 0.9))




##### Measure changes moth community composition across sample seasons and decades (ordination) #####

#### run ordination

library(vegan)

ord <- metaMDS(ord.dat[,4:219], k = 2)
# metaMDS = does a NMDS ordination
# [,4:219] = only include columns 4-219
# k = number of axes

#### plot ordination

site.scores <- as.data.frame(scores(ord, "sites")) # create a dataframe of the site scores to plot

site.scores$site <- ord.dat$site # create a column of site names, from the row names of data.scores

write.csv(site.scores, file = "site.scores.csv") # export site scores csv to add season and site2 column

# added season, season number, full site name, and decade to site.scores.csv and saved as site.scores2.csv

# re-import site scores for plotting
site.scores <- read.csv("site.scores2.csv", header = T, sep = ",", fileEncoding = "UTF-8-BOM")

site.scores$season <- as.factor(site.scores$season)
site.scores$decade <- as.factor(site.scores$decade)

# create a dataframe of the species scores to plot
species.scores <- as.data.frame(scores(ord, "species"))

# create a column of species, from the rownames of species.scores
species.scores$species <- rownames(species.scores)

## make plot with species of interest as text (common spp) and spp of noninterest as crosses 
## common moth species = species with at least 100 individuals caught in total (36 spp)

library(lemon)

# subset site scores to group sample seasons by site and decade
site.scoresR <- subset(site.scores, site2 %in% c("Ribbonwood"))
site.scoresR12 <- subset(site.scoresR, season.num <= 2)
site.scoresR23 <- subset(site.scoresR, season.num %in% c(2, 3))
site.scoresR34 <- subset(site.scoresR, season.num %in% c(3, 4))
site.scoresR45 <- subset(site.scoresR, season.num %in% c(4, 5))

site.scoresS <- subset(site.scores, site2 %in% c("Sugarloaf"))
site.scoresS12 <- subset(site.scoresS, season.num <= 2)
site.scoresS23 <- subset(site.scoresS, season.num %in% c(2, 3))
site.scoresS34 <- subset(site.scoresS, season.num %in% c(3, 4))
site.scoresS45 <- subset(site.scoresS, season.num %in% c(4, 5))

ggplot() + 
  geom_point(data = not.interest.dat, aes(x = NMDS1, y = NMDS2), size = 4, shape = 3, colour = "gray30") + 
  geom_text(data = interest.dat, aes(x = NMDS1, y = NMDS2, label = species), size = 4, colour = "black") + 
  geom_pointline(data = site.scoresR12, aes(x = NMDS1, y = NMDS2),  colour = "#F8766D", 
                 arrow = arrow(length = unit(0.25,"cm"), ends = "last", type = "closed"), distance = 8.5) + 
  geom_pointline(data = site.scoresR23, aes(x = NMDS1, y = NMDS2),  colour = "#F8766D", 
                 arrow = arrow(length = unit(0.25,"cm"), ends = "last", type = "closed"), distance = 8.5) + 
  geom_pointline(data = site.scoresR34, aes(x = NMDS1, y = NMDS2),  colour = "#F8766D", 
                 arrow = arrow(length = unit(0.25,"cm"), ends = "first", type = "closed"), distance = 8.5) + 
  geom_pointline(data = site.scoresR45, aes(x = NMDS1, y = NMDS2),  colour = "#F8766D", 
                 arrow = arrow(length = unit(0.25,"cm"), ends = "first", type = "closed"), distance = 8.5) + 
  geom_pointline(data = site.scoresS12, aes(x = NMDS1, y = NMDS2),  colour = "#00BFC4", 
                 arrow = arrow(length = unit(0.25,"cm"), ends = "first", type = "closed"), distance = 8.5) + 
  geom_pointline(data = site.scoresS23, aes(x = NMDS1, y = NMDS2),  colour = "#00BFC4", 
                 arrow = arrow(length = unit(0.25,"cm"), ends = "last", type = "closed"), distance = 8.5) + 
  geom_pointline(data = site.scoresS34, aes(x = NMDS1, y = NMDS2),  colour = "#00BFC4", 
                 arrow = arrow(length = unit(0.25,"cm"), ends = "first", type = "closed"), distance = 8.5) + 
  geom_pointline(data = site.scoresS45, aes(x = NMDS1, y = NMDS2),  colour = "#00BFC4", 
                 arrow = arrow(length = unit(0.25,"cm"), ends = "first", type = "closed"), distance = 8.5) + 
  geom_point(data = site.scores, aes(x = NMDS1, y = NMDS2, shape = site2, colour = site2), size = 8) + 
  geom_text(data = site.scores, aes(x = NMDS1, y = NMDS2, label = season.num)) + 
  coord_equal() + 
  labs(colour = "Site", shape = "Site") + 
  theme(panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", size = 0.75, fill = NA), 
        axis.line = element_line(colour = "black", size = 0.75), 
        legend.key = element_blank(), 
        axis.text = element_text(colour = "black", size = 13), 
        axis.title = element_text(colour = "black", size = 15), 
        axis.ticks = element_line(colour = "black", size = 1), 
        legend.position = c(0.9, 0.9), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 13))




##### Measure changes in average temperature near Cass #####

# linear regression for average annual temps from 1964 to 2013
annual.temp.lm <- lm(Annual ~ Year, data = temp.dat)

anova(annual.temp.lm, test = "F")

summary(annual.temp.lm)

# linear regression for average summer temps (October to April) from 1964 to 2013
summer.temp.lm <- lm(Oct.Apr ~ Year, data = temp.dat)

anova(summer.temp.lm, test = "F")

summary(summer.temp.lm)

# linear regression for average winter temps (May to September) from 1964 to 2013
winter.temp.lm <- lm(May.Sep ~ Year, data = temp.dat)

anova(winter.temp.lm, test = "F")

summary(winter.temp.lm)

# plot avg summer and winter temps from 1964 to 2021

summer.dat <- subset(temp.dat2, Months %in% c("Oct-Apr"))
winter.dat <- subset(temp.dat2, Months %in% c("May-Sep"))

ggplot(data = temp.dat2) + 
  geom_line(aes(x = Year, y = Temp, colour = Months), size = 1.5) + 
  scale_y_continuous(name = "Cass mean air temperature (°C)", limits = c(0, 14), 
                     breaks = c(0, 2, 4, 6, 8, 10, 12, 14)) + 
  geom_smooth(data = summer.dat, aes(x = Year, y = Temp), method = "lm", formula = "y ~ x", 
              colour = "#D55E00", size = 1.5) + 
  geom_smooth(data = winter.dat, aes(x = Year, y = Temp), method = "lm", formula = "y ~ x", 
              colour = "#0072B2", size = 1.5) + 
  scale_colour_manual(values=c("Oct-Apr" = "#D55E00", "May-Sep" = "#0072B2")) + 
  scale_x_continuous(name = "Year", breaks = c(1960, 1970, 1980, 1990, 2000, 2010, 2020)) + 
  theme(legend.title = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 1), 
        legend.key = element_blank(), 
        axis.text = element_text(colour = "black", size = 15), 
        axis.title = element_text(colour = "black", size = 17), 
        axis.ticks = element_line(size = 1), 
        legend.position = "top", 
        legend.text = element_text(size = 13))















