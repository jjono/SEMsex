
# Worked exampe 3: Cheating vs caring (from Supporting Information)

# This code is in two sections. In Section 1, we calculate phenotypic and quantitative genetic path
# coefficients using the package lavaan. In Section 2, we use the package ASReml to estimate the additive
# genetic variances and covariances that are needed for the quantitative genetic analysis in Section 1.
# The sections are presented in this order because ASReml is proprietary and some users may wish to use
# alternative means of estimating additive genetic variances and covariances.

# Before running any code, import the data files
# CheatingVsCaring-RelTraitData.csv
# CheatingVsCaring-Pedigree.csv
# from the Supporting Information (uncomment next two lines and change path as necessary)
# traits <- read_csv("~/Desktop/CheatingVsCaring-RelTraitData.csv")
# ped<-read_csv("~/Desktop/CheatingVsCaring-Pedigree.csv")

# Section 1: Path models

library(lavaan)
library(readr)

# Define our path model in the syntax of the lavaan package. The path model is the same for the phenotypic
# and quantitative genetic analyses. Note that * is the lavaan syntax for naming path coefficients
# (i.e. it does not represent an iteraction term).

cvc.model <- 'O~~varO*O               # assign name varO to the variance in O
C~~varC*C                             # assign name varC to the variance in C
O~~covCO*C                            # assign name covCO to the covariance between C and O
M ~ bMO*O + bMC*C                     # model for mating success
R ~ bRO*O + bRC*C + bRM*M             # model for reproductive success

xO := varO*bRM*bMO + covCO*bRM*bMC    # sexual selection on O
rO := varO*bRO + covCO*bRC            # remaining selection on O
sO := xO + rO                         # total selection on O

xC := varC*bRM*bMC + covCO*bRM*bMO    # sexual selection on C
rC := varC*bRC + covCO*bRO            # remaining selection on C
sC := xC + rC'                        # total selection on C

# Fit the phenotypic model using the lavaan package

pheno.fit <- sem(cvc.model,
                 data=traits)

# Summarise the phenotypic model

summary(pheno.fit)

# The lower triangle of the additive genetic variance-covariance matrix for the variables (C,O,M,R), which
# is estimated in Section 2.

lower <-'
0.5043 
-0.2134   0.5302
-0.0862   0.1861   0.0945
0.0829  0.0065  0.0055  0.0192'

# The full additive genetic variance-covariance matrix

cov.matrix <- 
  getCov(lower, names = c("C","O","M","R"))

# Fit the quantitative genetic model using the lavaan package

quant.gen.fit <- sem(cvc.model,
            sample.cov = cov.matrix,
            sample.nobs = 10000,
            se="none")

# Summarise the quantitative genetic model

summary(quant.gen.fit)



# Section 2: Estimation of additive genetic variance-covariance matrix using the package ASReml

library(asreml)
library('pedantics')

traits$IndividualID<-as.factor(traits$IndividualID)
head(traits)
dim(traits)

for (i in 1:3) {ped[,i]<-as.factor(ped[,i])}
str(ped)
head(ped)

ainv<-asreml.Ainverse(ped)$ginv
dim(ainv)

splom(traits[,2:5])

ped2<-ped
colnames(ped2)<-c('id','dam','sire')

# Multivariate animal model
m2<-asreml(cbind(C,O,M,R) ~ trait,
                 random=~ped(IndividualID):us(trait),ginverse=list(IndividualID=ainv),
                 rcov=~units:us(trait),
                 data=traits, maxiter=40000,
                 na.method.X="exclude", na.method.Y="include", singular.ok=TRUE,trace=T)
summary(m2,nice=T)$nice

# Define function to create matrix with SEs from varcomp list
# object<-summary(model)$varcom[wanted.rows]  # list of VCs
# VC.matrix<-as.matrix(data.frame(summary(model,nice=T)$nice)$ped)  # matrix of VCs
# ntrait<-number of traits

make.matrix.withSE<-function (object,VC.matrix,ntraits,dp) {
  k<-length(object[,1])
  list.VC<-object[,2]
  list.SE<-object[,3]
  list.z<-object[,4]
  SE.matrix<-matrix(rep(NA,ntraits^2),ntraits)
  full.matrix<-SE.matrix
  # now turn the column of SEs into a matrix as well, and concatenate with estimates
  for (i in 1:(ntraits^2)) {
    cell.i<-match(VC.matrix[i],list.VC)
    SE.matrix[i]<-list.SE[cell.i]
    z.i<-list.z[cell.i]
#    star<-if(abs(z.i)>2) "*"
        full.matrix[i]<-paste(round(VC.matrix[i],dp)," (",round(SE.matrix[i],dp),")",sep="")
  }
 as.data.frame(full.matrix)
}

# Multivariate phenotypic model
m3<-asreml(cbind(C,O,M,R) ~ trait,
           rcov=~units:us(trait),
           data=traits, maxiter=40000,
           na.method.X="exclude", na.method.Y="include", singular.ok=TRUE,trace=T)
summary(m3)$varcomp

m3.P.list<-summary(m3)$varcomp[2:11,]
VC.P.matrix<-as.matrix(data.frame(summary(m3,nice=T)$nice))
make.matrix.withSE(m3.list,VC.P.matrix,4,4)

# Now the animal model
summary(m2)$varcomp
summary(m2,nice=T)$nice

m2.G.list<-summary(m2)$varcomp[1:10,]
VC.G.matrix<-as.matrix(data.frame(summary(m2,nice=T)$nice$ped))
make.matrix.withSE(m2.G.list,VC.G.matrix,4,4)

m2.R.list<-summary(m2)$varcomp[12:21,]
VC.R.matrix<-as.matrix(data.frame(summary(m2,nice=T)$nice$R))
make.matrix.withSE(m2.R.list,VC.R.matrix,4,4)







