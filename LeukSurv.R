## ---- message = FALSE------------------------------------------------------------------------------------------------------------------
#############################################################
#
# GJRM spatial model fitting 
# Data: Leukemia data set
#############################################################

rm(list=ls())

# Load packages:
library(GJRM)
library(xtable)
library(knitr)
library(survival)
library(spBayesSurv)
library(BayesX)
library(R2BayesX)
library(sp)
library(fields)
library(surveillance)

library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)

#------------------------------------------------------------------------
# Life tables
#------------------------------------------------------------------------
LT <- as.data.frame(read.table("England_LT_2010_2015_dep_gor.txt",header = T))
head(LT)

LT$sex <- LT$sex-1


#----------------------------------------------------------------------------------------------
# Leukemia data
#---------------------------------------------------------------------------------------------
data(LeukSurv)
?LeukSurv
head(LeukSurv)
dim(LeukSurv)

# Sample size
n <- nrow(LeukSurv)

# Age
age <- as.numeric(LeukSurv$age)

# Sex
sex <- as.numeric(LeukSurv$sex)

# white blood cell count at diagnosis
wbc <- as.numeric(LeukSurv$wbc)

# Townsend score 
tpi <- as.numeric(LeukSurv$tpi)

# Survival time
time <- LeukSurv$time/365.24

# Age at the end of follow up
age.out <- round(age + time)

# Creating deprivation quintiles based on the tpi
quintiles <- quantile(LeukSurv$tpi,c(0.2,0.4,0.6,0.8))

dep0 <- vector()

for(i in 1:n){
  if(LeukSurv$tpi[i] <= quintiles[1]) dep0[i] <- 1
  if(quintiles[1] < LeukSurv$tpi[i] & LeukSurv$tpi[i] <= quintiles[2]) dep0[i] <- 2
  if(quintiles[2] < LeukSurv$tpi[i] & LeukSurv$tpi[i] <= quintiles[3]) dep0[i] <- 3
  if(quintiles[3] < LeukSurv$tpi[i] & LeukSurv$tpi[i] <= quintiles[4]) dep0[i] <- 4
  if(quintiles[4] < LeukSurv$tpi[i] ) dep0[i] <- 5
}

# Population hazard rates
# gor = 2 for all patients as this seems to correspond to North West England
# https://www.ons.gov.uk/methodology/geography/ukgeographies/administrativegeography/england
MAT.ID <- cbind(1:n, rep(2010,n), age.out,sex, dep0, rep(2,n))
MAT.ID <- as.data.frame(MAT.ID)
colnames(MAT.ID) <- c("index", "X_year", "age", "sex", "dep","gor")

hrates <- merge(x = MAT.ID, y =  LT, all.x = TRUE, by = c("X_year", "age", "sex", "dep","gor"), sort = FALSE)
hrates <- hrates[order(hrates$index),]
pop.rates <- hrates$rate


# Create data frame for model fitting
df <- data.frame(time = time, age = age, agec = scale(age), 
                 hrate = pop.rates, wbc = wbc, tpi = tpi, sex = sex,
                 status = as.logical(LeukSurv$cens), district= as.factor(LeukSurv$district))

# Select background mortality rate only for uncensored observations (this is our convention)
hrate.select = df$hrate[df$status == 1]

# Scale age at diagnosis
df$agec <- as.numeric(scale(df$age))

# Scale wbc at diagnosis
df$wbcc <- as.numeric(scale(df$wbc))

# Scale tpi at diagnosis
df$tpic <- as.numeric(scale(df$tpi))

# Set-up the polygons -----------------------------------------------------------------------------------------------------

# Map object
nwengland <- read.bnd(system.file("otherdata/nwengland.bnd",
                                  package = "spBayesSurv"))

nwsp = bnd2sp(nwengland)

nwspd = SpatialPolygonsDataFrame(Sr=nwsp, data=data.frame(NAME_0=paste0("X",1:length(nwsp))))

class(nwspd)

poldistrict = polys.setup(nwspd)

# Save the polygons
final.poldistrict = poldistrict$polys

# The final polygons will be:
xt <- list(polys = final.poldistrict)



## ---- message=FALSE--------------------------------------------------------------------------------------------------------------------
#################
# Model fitting
#################

eq_GJRM = list(time ~ s(log(time), bs = "mpi") + s(agec, bs='cr') + ti(log(time), agec) + sex +
                 wbcc + tpic +                
                 s(district, bs = 'mrf', xt = xt) ) 

               
# Fitting three different hazard structures and retaining the best model based on smallest AIC:
start_time <- Sys.time()
# Generalised proportional hazards model
out_GJRM_1 = gamlss(eq_GJRM, data = df, surv = TRUE, margin = 'PH',
                    cens = status, type.cens = "R", hrate = hrate.select)
end_time <- Sys.time()
run_time <- end_time - start_time
run_time

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

conv.check(out_GJRM)

best
# 2
c(AIC(out_GJRM_1), AIC(out_GJRM_2), AIC(out_GJRM_3))



## ---- message=FALSE, results = "hide"--------------------------------------------------------------------------------------------------
# Estimate net survival (and confidence intervals) for the whole cohort 
# at times t = 1, 3, 5 years

net.surv <- matrix(0, ncol = 3, nrow = 24)

times = c(1,5)

for(i in 1:24){
  df.temp <- df[which(df$district == i),]
  
  net.surv[i,] <- hazsurv.plot(out_GJRM, type = 'surv', t.range = times, ls = 3, newdata = df.temp,
                               intervals = FALSE, n.sim = 1000, plot.out = FALSE)$s
  }


## --------------------------------------------------------------------------------------------------------------------------------------

# One year net survival by district
nsd1 <- cbind(1:24, net.surv[,1])
colnames(nsd1) <- cbind("District", "1-year Net Survival")
kable(nsd1, digits = 3)

#One-year net survival, LeukSurv
polys.map(final.poldistrict, net.surv[,1], zlim=c(0,0.51), las=1, lab = "",  
          cex.main = 1.2, rev.col=FALSE, main = "", scheme = "heat")


#Three-year net survival, LeukSurv
polys.map(final.poldistrict, net.surv[,2], zlim=c(0,0.51), las=1, lab = "",  
          cex.main = 1.2, rev.col=FALSE, main = "", scheme = "heat")


#Five-year net survival, LeukSurv
polys.map(final.poldistrict, net.surv[,3], zlim=c(0,0.51), las=1, lab = "",  
          cex.main = 1.2, rev.col=FALSE, main = "", scheme = "heat")

