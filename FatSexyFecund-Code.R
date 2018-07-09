
# Worked example 1: Fat, sexy and fecund

# Before running any code, import the data file FatSexyFecund-Data.csv from the Supporting Information
# (uncomment next line and change path as necessary)
# PTraits <- read_csv("~/Desktop/FatSexyFecund-Data.csv")

library(lavaan)
library(readr)

# Define variables body size (Z), mating success (M), reproductive success (R) and fecundity-per-mate (F)

Z <- PTraits$Z
M <- PTraits$M
R <- PTraits$R
F <- R/M

# Calculate relative values for M, R and F by dividing each value by the mean

Mrel <- M/mean(M)
Rrel <- R/mean(R)
Frel <- F/mean(F)

# Calculate the path coefficients as given in Table 2. Variable names of the form betaAB.C are partial
# regression coefficients of A on B, given additional covariate C (i.e. based on the multiple regression
# A~B+C+D). Simple regression coefficients of A on B are written betaAB or betaAB.0. The Bateman gradient
# is betass. Variables names of the form beta...raw are based on raw values of M, R and F, whereas
# variables of the form beta...rel are based on relative values of these variables.

betaMZraw <- summary(lm(M~Z))$coefficients["Z",1]
betaFZraw <- summary(lm(F~Z))$coefficients["Z",1]
betaRZraw <- summary(lm(R~Z+M))$coefficients["Z",1]
betassraw <- summary(lm(R~M))$coefficients["M",1]
betaRMraw <- summary(lm(R~Z+M))$coefficients["M",1]
betaRM.Fraw <- summary(lm(R~M+F))$coefficients["M",1]
betaRF.0raw <- summary(lm(R~F))$coefficients["F",1]
betaRF.Mraw <- summary(lm(R~M+F))$coefficients["F",1]

betaMZrel <- summary(lm(Mrel~Z))$coefficients["Z",1]
betaFZrel <- summary(lm(Frel~Z))$coefficients["Z",1]
betaRZrel <- summary(lm(Rrel~Z+Mrel))$coefficients["Z",1]
betassrel <- summary(lm(Rrel~Mrel))$coefficients["Mrel",1]
betaRMrel <- summary(lm(Rrel~Z+Mrel))$coefficients["Mrel",1]
betaRM.Frel <- summary(lm(Rrel~Mrel+Frel))$coefficients["Mrel",1]
betaRF.0rel <- summary(lm(Rrel~Frel))$coefficients["Frel",1]
betaRF.Mrel <- summary(lm(Rrel~Mrel+Frel))$coefficients["Frel",1]

# Reproduce the estimates of sexual selection, remaining selection and total selection in Table 2, based on
# the models of Arnold (1994), Conner (1996) and Jones (2009), as well as our model. We also reproduce
# our model in the lavaan package (see below).

# Arnold (1994)

data.frame("Arnold.raw"=c(var(Z)*betaRM.Fraw*betaMZraw,var(Z)*betaRF.Mraw*betaFZraw,cov(Z,R)),
           "Arnold.relative"=c(var(Z)*betaRM.Frel*betaMZrel,var(Z)*betaRF.Mrel*betaFZrel,cov(Z,Rrel)),
           row.names=c("sexual.selection","remaining.selection","total.selection"))

# Conner (1996)

data.frame("Conner.raw"=c(var(Z)*betassraw*betaMZraw,var(Z)*betaRF.0raw*betaFZraw,cov(Z,R)),
           "Conner.relative"=c(var(Z)*betassrel*betaMZrel,var(Z)*betaRF.0rel*betaFZrel,cov(Z,Rrel)),
           row.names=c("sexual.selection","remaining.selection","total.selection"))

# Jones (2009)

data.frame("Jones.raw"=c(var(Z)*betassraw*betaMZraw,cov(Z,R)-var(Z)*betassraw*betaMZraw,cov(Z,R)),
           "Jones.relative"=c(var(Z)*betassrel*betaMZrel,cov(Z,Rrel)-var(Z)*betassrel*betaMZrel,cov(Z,Rrel)),
           row.names=c("sexual.selection","remaining.selection","total.selection"))

# Our model

data.frame("our.model.raw"=c(var(Z)*betaRMraw*betaMZraw,var(Z)*betaRZraw,cov(Z,R)),
           "our.model.relative"=c(var(Z)*betaRMrel*betaMZrel,var(Z)*betaRZrel,cov(Z,Rrel)),
           row.names=c("sexual.selection","remaining.selection","total.selection"))

# Define our path model in the syntax of the lavaan package (shown here for raw data). Note that * is the
# syntax for naming path coefficients (i.e. it does not represent an iteraction term). This model is based
# on a single trait (body size), but if multiple traits are included, then their covariances must also be
# modelled (see worked example 'Cheating vs Caring' in Supporting Information).

fsf.model <- 'Z~~varZ*Z   # assign name varZ to the variance in Z
M ~ bMZ*Z                 # model for mating success
R ~ bRZ*Z + bRM*M         # model for reproductive success
xZ := varZ*bRM*bMZ        # sexual selection
rZ := varZ*bRZ            # remaining selection
sZ := xZ + rZ'            # total selection

# Fit the path model using the lavaan package. For consistency with equations (5) and (6) in the main text,
# we use estimator="GLS", which is equivalent to ordinary least  squares in this case. By default, lavaan
# uses maximum likelihood estimation. Many alternative estimation techniques are supported, including
# robust variants (see tutorial at http://lavaan.ugent.be/tutorial/).

fsf.fit <- sem(fsf.model,data=PTraits, estimator="GLS")

# Summarise the fitted path model

summary(fsf.fit)


