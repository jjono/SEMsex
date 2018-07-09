
# Worked example 2: Size selection in male pipefish

# Many thanks to Kenyon Mobley and colleagues for providing the dataset for this worked example
# (see main text for references)

# Before running any code, import the data file Pipefish-Data.csv from the Supporting Information
# (uncomment next line and change path as necessary)
# PipefishData <- read_csv("~/Desktop/Pipefish-Data.csv")

library(lavaan)
library(readr)

# Change code below to analyse data from populations at any of the six sites.
# Possible site values are (see main text): "VA"  "NC"  "TB"  "SJ"  "TX1" "TX3"

PipefishSiteData <- PipefishData[which(PipefishData$Site=="NC"),]

# Variance-standardise trait values and relativise (i.e. mean-standardise) mating and reproductive success
# within the chosen population

PipefishSiteData$Z <- (PipefishSiteData$Z-mean(PipefishSiteData$Z))/sd(PipefishSiteData$Z)
PipefishSiteData$M <- PipefishSiteData$M/mean(PipefishSiteData$M)
PipefishSiteData$R <- PipefishSiteData$R/mean(PipefishSiteData$R)

# Define our path model in the syntax of the lavaan package. Note that * is the syntax for naming path
# coefficients (i.e. it does not represent an iteraction term). This model is based on a single trait
# (body length), but if multiple traits are included, then their covariances must also be modelled
# (see worked example 'Cheating vs Caring' in Supporting Information).

pipefish.model <- 'Z~~varZ*Z    # assign name varZ to the variance in Z
M ~ bMZ*Z                       # model for mating success
R ~ bRM*M + bRZ*Z               # model for reproductive success
xZ := varZ*bMZ*bRM              # sexual selection
rZ := varZ*bRZ                  # remaining selection
sZ := varZ*bMZ*bRM + varZ*bRZ'  # total selection

# Fit the path model using the lavaan package. For consistency with equations (5) and (6) in the main text,
# we use estimator="GLS", which is equivalent to ordinary least  squares in this case. By default, lavaan
# uses maximum likelihood estimation. Many alternative estimation techniques are supported, including
# robust variants (see tutorial at http://lavaan.ugent.be/tutorial/).

pipefish.fit <- sem(pipefish.model,data=PipefishSiteData, estimator="GLS")

# Summarise the fitted path model

summary(pipefish.fit)

# The code below reproduces the estimates of sexual selection, remaining selection and total selection
# in Table 3 for our model and that of Jones (2009).

Z <- PipefishSiteData$Z
M <- PipefishSiteData$M
R <- PipefishSiteData$R

# Calculate the mating gradient on Z (betaMZ), the simple Bateman gradient (betass), the partial Bateman
# gradient controlling for Z (betaRM.Z), and the selection gradient on Z controlling for M (betaRZ.M)

betaMZ <- summary(lm(M~Z))$coefficients["Z",1]
betass <- summary(lm(R~M))$coefficients["M",1]
betaRM.Z <- summary(lm(R~M+Z))$coefficients["M",1]
betaRZ.M <- summary(lm(R~M+Z))$coefficients["Z",1]

# Selection estimates

data.frame("Jones.2009"=c(var(Z)*betaMZ*betass,cov(Z,R)-var(Z)*betaMZ*betass,cov(Z,R)),
           "our.model"=c(var(Z)*betaMZ*betaRM.Z,var(Z)*betaRZ.M,cov(Z,R)),
           row.names=c("sexual.selection","remaining.selection","total.selection"))


