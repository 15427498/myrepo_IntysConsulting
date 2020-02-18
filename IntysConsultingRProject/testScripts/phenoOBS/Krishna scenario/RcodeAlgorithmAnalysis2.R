
library(BCS.JAMES3) 
library(MASS)
library(BCS.Base)
library(DiGGer)
library(xlsx)
sessionInfo()

#Note: We use the same scenario as in "RcodeAlexCheck"
datadir <-"/gpfs/gssgpfs1/biogrid/workspace/gehzz/projects/Pheno_OBS/USE_2016"
setwd(datadir)

#-----------------------------------------------------------------------------
# Read the .csv data
#-----------------------------------------------------------------------------

dat <- read.csv("phenotypes.csv", header=TRUE)
dat <- clean(dat)
geno <- read.csv("genotypes.csv", header=TRUE)
mhead(geno)
rownames(geno) <- geno$X
geno <- geno[, -which(colnames(geno) == "X")]
geno <- as.matrix(geno)

#NOTE:
#In the previous analyses, we run some experiments to test the three search techniques in a fixed scenario. 
#In this script, we will look more closely at the impact of the SUBSET SIZE on the problem complexity and, analogously, on the algorithmic performance.

###############################################################################################################

#-----------------------------------------------------------------------------
# (1) Generate and run the phenoOBSDataManagement data
#-----------------------------------------------------------------------------

parentsPhenotypeDataFrame <- data.frame(parentNames = dat$Material.Name,
                                        YIELD = dat$YIELD_EST,
                                        LINT = dat$LINT_EST)
parentsGeneticSimilarityMatrix <- geno
noOfTraits <- ncol(parentsPhenotypeDataFrame)-1
constraintParameters <- matrix(nrow = noOfTraits, ncol = 2)
constraintParameters[1,1] <- 0
traitObjectiveWeights <- matrix(nrow = noOfTraits, ncol = 2)
traitObjectiveWeights[1,] <- c(0.7, 1)  # YIELD
traitObjectiveWeights[2,] <- c(0.3, 1)  # LINT
javaInputData <- phenoObsDataManagement(parentsPhenotypeDataFrame,
                                        parentsGeneticSimilarityMatrix,
                                        constraintParameters,
                                        traitObjectiveWeights)

###############################################################################################################

#-----------------------------------------------------------------------------
# (2) Set preliminary fixed arguments
#-----------------------------------------------------------------------------

#I.e. for consistency we use the exact same algorithm termination conditions across the analysis scripts (yet differing domain and objective spaces in this analysis).
offspringPhenotype <- javaInputData[[1]]
parentsGeneticSimilarityArray <- javaInputData[[2]]
parentsPhenotypeDataFrame <- javaInputData[[3]]
geneticSimilarityWeight <- 0.4
penaltyReplicatedIndividuals <- 10
inclusionSet <- c(647)
maximumSearchTime <- 20000
timeWithoutImprovement <- 10000
maximumNumberMoves <- 80000000

###############################################################################################################

#--------------------------------------------------------------------------------------------
# (3) Conduct analysis with different subset sizes for the RANDOM DESCENT algorithm
#--------------------------------------------------------------------------------------------

#Experiment 1

solutionSize <- 20
numberSearches <- 1
searchIds <- c("Search 1")
numberRuns <- 10
searchMethod <- "RD"

phenoObsAnalysisOutput1 <- phenoObsAlgorithmPerformance(parentsPhenotypeDataFrame,
                                                        offspringPhenotype,
                                                        parentsGeneticSimilarityArray,
                                                        geneticSimilarityWeight,
                                                        penaltyReplicatedIndividuals,
                                                        inclusionSet,
                                                        solutionSize,
                                                        maximumSearchTime,
                                                        timeWithoutImprovement,
                                                        maximumNumberMoves,
                                                        numberSearches,
                                                        numberRuns,
                                                        searchIds,
                                                        searchMethod,
                                                        analysisParametersMS,
                                                        analysisParametersPT)

bestFoundValues1 <- phenoObsAnalysisOutput1[[1]]
print(bestFoundValues1)
print(mean(bestFoundValues1))
boxplot(t(bestFoundValues1),
        data=bestFoundValues1, 
        main="RD best found solutions for n = 20",
        ylab="Objective value",
        ylim=c(0.27, 0.282)) 
searchId <- 1
updatedBestFoundValues1 <- phenoObsAnalysisOutput1[[3]]
updatedBestFoundTimes1 <- phenoObsAnalysisOutput1[[4]]
plot(updatedBestFoundTimes1[searchId,1,], 
     updatedBestFoundValues1[searchId,1,], 
     main="RD convergence for n = 20",
     xlab="Time", 
     ylab="Objective value",
     xlim=c(0, 500),
     ylim=c(0.15, 0.282),
     type='l', 
     col='black', 
     lwd=0.5)
for(l in 2:numberRuns){
lines(updatedBestFoundTimes1[searchId,l,], 
      updatedBestFoundValues1[searchId,l,],
      type='l', 
      col='black', 
      lwd=0.5)}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Experiment 2

solutionSize <- 40
numberSearches <- 1
searchIds <- c("Search 1")
numberRuns <- 10
searchMethod <- "RD"

phenoObsAnalysisOutput12 <- phenoObsAlgorithmPerformance(parentsPhenotypeDataFrame,
                                                        offspringPhenotype,
                                                        parentsGeneticSimilarityArray,
                                                        geneticSimilarityWeight,
                                                        penaltyReplicatedIndividuals,
                                                        inclusionSet,
                                                        solutionSize,
                                                        maximumSearchTime,
                                                        timeWithoutImprovement,
                                                        maximumNumberMoves,
                                                        numberSearches,
                                                        numberRuns,
                                                        searchIds,
                                                        searchMethod,
                                                        analysisParametersMS,
                                                        analysisParametersPT)

bestFoundValues12 <- phenoObsAnalysisOutput12[[1]]
print(bestFoundValues12)
print(mean(bestFoundValues12))
boxplot(t(bestFoundValues12),
        data=bestFoundValues12, 
        main="RD best found solutions for n = 40",
        ylab="Objective value",
        ylim=c(0.25, 0.26)) 
searchId <- 1
updatedBestFoundValues12 <- phenoObsAnalysisOutput12[[3]]
updatedBestFoundTimes12 <- phenoObsAnalysisOutput12[[4]]
plot(updatedBestFoundTimes12[searchId,1,], 
     updatedBestFoundValues12[searchId,1,], 
     main="RD convergence for n = 40",
     xlab="Time", 
     ylab="Objective value",
     xlim=c(0, 500),
     ylim=c(0.15, 0.282),
     type='l', 
     col='black', 
     lwd=0.5)
for(l in 2:numberRuns){
  lines(updatedBestFoundTimes12[searchId,l,], 
        updatedBestFoundValues12[searchId,l,],
        type='l', 
        col='black', 
        lwd=0.5)}

#--------------------------------------------------------------------------------------------
# (4) Conduct analysis with different subsets for the METROPOLIS SEARCH algorithm
#--------------------------------------------------------------------------------------------

#Experiment 1

solutionSize <- 20
numberSearches <- 6
numberRuns <- 8
searchIds <- c("Search 1","Search 2","Search 3","Search 4","Search 5","Search 6")
searchMethod <- "MS"
analysisParametersMS <- c(0.001,0.0005,0.000075,0.000025,0.000005,0.0000005) 

phenoObsAnalysisOutput2 <- phenoObsAlgorithmPerformance(parentsPhenotypeDataFrame,
                                                        offspringPhenotype,
                                                        parentsGeneticSimilarityArray,
                                                        geneticSimilarityWeight,
                                                        penaltyReplicatedIndividuals,
                                                        inclusionSet,
                                                        solutionSize,
                                                        maximumSearchTime,
                                                        timeWithoutImprovement,
                                                        maximumNumberMoves,
                                                        numberSearches,
                                                        numberRuns,
                                                        searchIds,
                                                        searchMethod,
                                                        analysisParametersMS,
                                                        analysisParametersPT)

bestFoundValues2 <- phenoObsAnalysisOutput2[[1]]
print(bestFoundValues2)
boxplot(t(bestFoundValues2),
        data=bestFoundValues2, 
        main="MS best found solutions for n = 20",
        names=analysisParametersMS,
        xlab="Temperature", 
        ylab="Objective value",
        ylim=c(0.27, 0.282)) 
searchId <- 3
updatedBestFoundValues2 <- phenoObsAnalysisOutput2[[3]]
updatedBestFoundTimes2 <- phenoObsAnalysisOutput2[[4]]
plot(updatedBestFoundTimes2[searchId,1,], 
     updatedBestFoundValues2[searchId,1,], 
     main="MS convergence for n = 20 and T = 0.000075",
     xlab="Time", 
     ylab="Objective value",
     xlim=c(0, 10000),
     ylim=c(0.15, 0.282),
     type='l', 
     col='black', 
     lwd=0.5)
for(l in 2:numberRuns){
  lines(updatedBestFoundTimes2[searchId,l,], 
        updatedBestFoundValues2[searchId,l,],
        type='l', 
        col='black', 
        lwd=0.5)}
