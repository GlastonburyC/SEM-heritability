
## FORMAT FILE FOR HERITABILITY ESTIMATION

allDup <- function (value) 
{ 
    duplicated(value) | duplicated(value, fromLast = TRUE) 
} 

cell1=NULL
cell2=NULL
bmi1=NULL
bmi2=NULL
zyg=NULL

for(i in seq(from=1, to=dim(final)[1], by=2)){
  cell1=append(cell1,final$dermis.craig[i])
  
  zyg=append(zyg,final$Zyg[i])
  
}

for(i in seq(from=2, to=dim(final)[1], by=2)){
  cell2=append(cell2,final$dermis.craig[i])
  
}


cell1<-as.data.frame(cell1)
cell2<-as.data.frame(cell2)
res<-cbind(cell1,cell2,zyg)

require(OpenMx)
# Load Library
# -----------------------------------------------------------------------------

# Load Data
data(twinData)

# Select Variables for Analysis
selVars   <- c('cell1','cell2')
aceVars   <- c("A1","C1","E1","A2","C2","E2")

# Select Data for Analysis
mzData    <- subset(res, zyg=='MZ', selVars)
dzData    <- subset(res, zyg=='DZ', selVars)

# Generate Descriptive Statistics
colMeans(mzData,na.rm=TRUE)
colMeans(dzData,na.rm=TRUE)
cov(mzData,use="complete")
cov(dzData,use="complete")
# Prepare Data
# -----------------------------------------------------------------------------

# Path objects for Multiple Groups
manifestVars=selVars
latentVars=aceVars
# variances of latent variables
latVariances <- mxPath( from=aceVars, arrows=2, 
                        free=FALSE, values=1 )
# means of latent variables
latMeans     <- mxPath( from="one", to=aceVars, arrows=1, 
                        free=FALSE, values=0 )
# means of observed variables
obsMeans     <- mxPath( from="one", to=selVars, arrows=1, 
                        free=TRUE, values=20, labels="mean" )
# path coefficients for twin 1
pathAceT1    <- mxPath( from=c("A1","C1","E1"), to="cell1", arrows=1, 
                        free=TRUE, values=.5,  label=c("a","c","e") )
# path coefficients for twin 2
pathAceT2    <- mxPath( from=c("A2","C2","E2"), to="cell2", arrows=1, 
                        free=TRUE, values=.5,  label=c("a","c","e") )
# covariance between C1 & C2
covC1C2      <- mxPath( from="C1", to="C2", arrows=2, 
                        free=FALSE, values=1 )
# covariance between A1 & A2 in MZ twins
covA1A2_MZ   <- mxPath( from="A1", to="A2", arrows=2, 
                        free=FALSE, values=1 )
# covariance between A1 & A2 in DZ twins
covA1A2_DZ   <- mxPath( from="A1", to="A2", arrows=2, 
                        free=FALSE, values=.5 )

# Data objects for Multiple Groups
dataMZ       <- mxData( observed=mzData, type="raw" )
dataDZ       <- mxData( observed=dzData, type="raw" )

# Combine Groups
paths        <- list( latVariances, latMeans, obsMeans, 
                      pathAceT1, pathAceT2, covC1C2 )
modelMZ      <- mxModel(model="MZ", type="RAM", manifestVars=selVars, 
                        latentVars=aceVars, paths, covA1A2_MZ, dataMZ )
modelDZ      <- mxModel(model="DZ", type="RAM", manifestVars=selVars, 
                        latentVars=aceVars, paths, covA1A2_DZ, dataDZ )
minus2ll     <- mxAlgebra( expression=MZ.fitfunction + DZ.fitfunction, 
                           name="minus2loglikelihood" )
obj          <- mxFitFunctionAlgebra( "minus2loglikelihood" )
twinACEModel     <- mxModel(model="ACE", modelMZ, modelDZ, minus2ll, obj )

# Run Model
twinACEFit   <- mxRun(twinACEModel)
twinACESum   <- summary(twinACEFit)
twinACESum
# Fit ACE Model with RawData and Path-Style Input
# -----------------------------------------------------------------------------

# Generate & Print Output
# mean
M  <- mxEval(mean, twinACEFit)
# additive genetic variance, a^2
A  <- mxEval(a*a, twinACEFit)
# shared environmental variance, c^2
C  <- mxEval(c*c, twinACEFit)
# unique environmental variance, e^2
E  <- mxEval(e*e, twinACEFit)
# total variance
V  <- (A+C+E)
# standardized A
a2 <- A/V
# standardized C
c2 <- C/V
# standardized E
e2 <- E/V
# table of estimates
estACE <- rbind(cbind(A,C,E),cbind(a2,c2,e2))
# likelihood of ACE model
LL_ACE <- mxEval(fitfunction, twinACEFit)
# Get Model Output
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------

# Change Path Objects
# path coefficients for twin 1
pathAceT1    <- mxPath( from=c("A1","C1","E1"), to="bmi1", arrows=1, 
                        free=c(T,F,T), values=c(.6,0,.6),  label=c("a","c","e") )
# path coefficients for twin 2
pathAceT2    <- mxPath( from=c("A2","C2","E2"), to="bmi2", arrows=1, 
                        free=c(T,F,T), values=c(.6,0,.6),  label=c("a","c","e") )

# Combine Groups
paths        <- list( latVariances, latMeans, obsMeans, 
                      pathAceT1, pathAceT2, covC1C2 )
modelMZ      <- mxModel(model="MZ", type="RAM", manifestVars=selVars, 
                        latentVars=aceVars, paths, covA1A2_MZ, dataMZ )
modelDZ      <- mxModel(model="DZ", type="RAM", manifestVars=selVars, 
                        latentVars=aceVars, paths, covA1A2_DZ, dataDZ )
twinAEModel      <- mxModel(model="AE", modelMZ, modelDZ, minus2ll, obj )

# Run Model
twinAEFit    <- mxRun(twinAEModel)
twinAESum    <- summary(twinAEFit)
# Fit AE model
# -----------------------------------------------------------------------------

# Generate & Print Output
M  <- mxEval(mean, twinAEFit)
A  <- mxEval(a*a, twinAEFit)
C  <- mxEval(c*c, twinAEFit)
E  <- mxEval(e*e, twinAEFit)
V  <- (A+C+E)
a2 <- A/V
c2 <- C/V
e2 <- E/V
estAE <- rbind(cbind(A, C, E),cbind(a2, c2, e2))
LL_AE <- mxEval(fitfunction, twinAEFit)
LRT_ACE_AE <- LL_AE - LL_ACE
estACE
estAE
LRT_ACE_AE
# Get Model Output
# -----------------------------------------------------------------------------

estACE
estAE
LRT_ACE_AE
#Print relevant output
# -----------------------------------------------------------------------------
