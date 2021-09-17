## ---- message = FALSE------------------------------------------------------------------------------------------------------------------
# Required packages
library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)

# Load the simulated data set
load('data_lung.Rda')

dim(df)

# A snipet of the data frame
head(df)

# Histogram of the survival times
hist(df$time, breaks = 30, probability = TRUE, xlab = "Time", main = "Simulated data")
box()

# Summaries

# Deprivation level
table(df$dep)

# Age at diagnosis
summary(df$age)

# Proportion of observed times to event
mean(df$status)

# Kaplan Meier estimators for each deprivation level
km_dep_fit <- survfit(Surv(time, status) ~ dep, data=df)
autoplot(km_dep_fit, conf.int = FALSE, xlab = "Time", ylab = "Survival")


## ---- message=FALSE--------------------------------------------------------------------------------------------------------------------
#############################################################
#
# GJRM model fitting 
# Data: Women diagnosed with lung cancer, England
# The Simulacrum data set
#
#############################################################


# Load packages:
library(GJRM)
library(knitr)

# Select background mortality rate only for uncensored observations (this is our convention)
hrate.select = df$hrate[df$status == 1]


#################
# Model fitting
#################

# Using only model equation:

eq_GJRM = list(time ~ s(log(time), bs = "mpi") + s(agec, bs='cr') + ti(log(time), agec) +
               dep.2 + dep.3 + dep.4 + dep.5)


# Fitting three different hazard structures and retaining the best model based on smallest AIC:

# Generalised proportional hazards model
out_GJRM_1 = gamlss(eq_GJRM, data = df, surv = TRUE, margin = 'PH',
                    cens = status, type.cens = "R", hrate = hrate.select)

# Generalised proportional odds model
out_GJRM_2 = gamlss(eq_GJRM, data = df, surv = TRUE, margin = 'PO',
                    cens = status, type.cens = "R", hrate = hrate.select)

# Generalised probit model
out_GJRM_3 = gamlss(eq_GJRM, data = df, surv = TRUE, margin = 'probit',
                    cens = status, type.cens = "R", hrate = hrate.select)

best <- which.min(c(AIC(out_GJRM_1), AIC(out_GJRM_2), AIC(out_GJRM_3)))

if(best == 1) out_GJRM <- out_GJRM_1
if(best == 2) out_GJRM <- out_GJRM_2
if(best == 3) out_GJRM <- out_GJRM_3

# Best model
best

# AIC
c(AIC(out_GJRM_1), AIC(out_GJRM_2), AIC(out_GJRM_3))




## ---- message=FALSE, results = "hide"--------------------------------------------------------------------------------------------------
# Estimate net survival (and confidence intervals) for the whole cohort 
# at times t = 1, 3 , 5 (years):

net.surv <- hazsurv.plot(out_GJRM, type = 'surv', t.range = c(1,5), newdata = df, ls = 3,
                         intervals = TRUE, n.sim = 1000, plot.out = FALSE)



## --------------------------------------------------------------------------------------------------------------------------------------
net.surv$s
net.surv$CIs



## ---- message=FALSE, results = "hide"--------------------------------------------------------------------------------------------------
# Estimate net survival (and confidence intervals) for deprivation level 1 (least deprived patients) 
# at times t = 1, 3 , 5 (years):

# Define sub-population index for deprivation level 1
pop.indices.dep1 <- (df$dep.1 == 1)

# Create sub-population data using previous index
data.dep1 <- df[pop.indices.dep1, ]

# Estimate net survival
net.surv.dep1 <- hazsurv.plot(out_GJRM, type = 'surv', newdata = data.dep1, t.range = c(1,5),
                              ls = 3, intervals = TRUE, n.sim = 1000, plot.out = FALSE)


## --------------------------------------------------------------------------------------------------------------------------------------
net.surv.dep1$s
net.surv.dep1$CIs


## ---- message=FALSE, results = "hide"--------------------------------------------------------------------------------------------------
# Estimate net survival (and confidence intervals) for deprivation level 5 (most deprived patients) 
# at times t = 1, 3 , 5 (years):

# Define sub-population index for deprivation level 5
pop.indices.dep5 <- (df$dep.5 == 1)

# Create sub-population data using previous index
data.dep5 <- df[pop.indices.dep5, ]

# Estimate net survival
net.surv.dep5 <- hazsurv.plot(out_GJRM, type = 'surv', newdata = data.dep5, t.range = c(1,5),
                              ls = 3, intervals = TRUE, n.sim = 1000, plot.out = FALSE)


## --------------------------------------------------------------------------------------------------------------------------------------
net.surv.dep5$s
net.surv.dep5$CIs


## --------------------------------------------------------------------------------------------------------------------------------------
# Export results to a table (use library(xtable) to export to LaTeX)
MNSG <- cbind(c(1,3,5), net.surv$s, net.surv$CIs[,1], net.surv$CIs[,2], 
              net.surv.dep1$s, net.surv.dep1$CIs[,1], net.surv.dep1$CIs[,2],
              net.surv.dep5$s, net.surv.dep5$CIs[,1], net.surv.dep5$CIs[,2])

colnames(MNSG) <- c('F-up (yrs)','NS pop','95% LL', '95% UL',
                    'NS dep 1','95% LL', '95% UL',
                    'NS dep 5','95% LL', '95% UL')
kable(MNSG, digits = 3)


## ---- message=FALSE, results = "hide"--------------------------------------------------------------------------------------------------
# Plot of net survival by deprivation (least deprived versus most deprived):


net.surv.dep1.plot <- hazsurv.plot(out_GJRM, type = 'surv', newdata = data.dep1, t.range = c(0,5),
                              ls = 100, intervals = TRUE, n.sim = 1000, plot.out = TRUE)

net.surv.dep5.plot <- hazsurv.plot(out_GJRM, type = 'surv', newdata = data.dep5, t.range = c(0,5),
                                   ls = 100, intervals = TRUE, n.sim = 1000, plot.out = TRUE)



## --------------------------------------------------------------------------------------------------------------------------------------
# Comparison 
plot(seq(0,5, length.out = 100), net.surv.dep1.plot$s, las=1, type = 'l', lwd = 1.5, ylim = c(0,1),
     ylab = 'Net survival', xlab = 'Follow-up time (years)', cex.axis=1.25, cex.lab=1.25)

polygon(c(seq(0,5, length.out = 100),rev(seq(0,5, length.out = 100))), 
        c(net.surv.dep1.plot$CIs[,2],rev(net.surv.dep1.plot$CIs[,1])),col="grey", border=NA)

polygon(c(seq(0,5, length.out = 100),rev(seq(0,5, length.out = 100))),
        c(net.surv.dep5.plot$CIs[,2],rev(net.surv.dep5.plot$CIs[,1])),col="grey", border=NA)

lines(seq(0,5, length.out = 100), net.surv.dep1.plot$s, lwd = 1.5, lty = 1)
lines(seq(0,5, length.out = 100), net.surv.dep5.plot$s, lwd = 1.5, lty = 2)

legend('topright', legend = c('Least deprived', 'Most deprived'), lwd = 1.2, lty = c(1,5))

