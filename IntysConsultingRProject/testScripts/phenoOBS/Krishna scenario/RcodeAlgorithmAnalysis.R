
library(MASS)
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

#I.e. for consistency we use the exact same domain space, objective space and algorithm termination conditions across these experiments.
offspringPhenotype <- javaInputData[[1]]
parentsGeneticSimilarityArray <- javaInputData[[2]]
parentsPhenotypeDataFrame <- javaInputData[[3]]
geneticSimilarityWeight <- 0.4
penaltyReplicatedIndividuals <- 10
inbreedingCoefficientMetric <- FALSE
numberOffspring <- nrow(parentsGeneticSimilarityMatrix )
Gmatrix <- parentsGeneticSimilarityMatrix 
inbreedingCoefficientBoundsEstimatorSampleSize <- 300
inclusionSet <- 1
exclusionSet <- NULL
solutionSize <- 10
maximumSearchTime <- 20000
timeWithoutImprovement <- 1000
maximumNumberMoves <- 80000000

###############################################################################################################

#--------------------------------------------------------------------------------------------
# (3) Conduct analysis for the METROPOLIS SEARCH algorithm
#--------------------------------------------------------------------------------------------

#Experiment 1

numberSearches <- 7
numberRuns <- 5
searchIds <- c("Search 1","Search 2","Search 3","Search 4","Search 5","Search 6","Search 7")
searchMethod <- "MS"
analysisParametersMS <- c(20,5,1,0.05,0.001,0.0005,0.000001) 

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
boxplot(t(bestFoundValues1),
        data=bestFoundValues1, 
        main="MS best found solutions",
        names=analysisParametersMS,
        xlab="Temperature", 
        ylab="Objective value") 

#--------------------------------------------------------------------------------------------------------------

#Experiment 2

numberSearches <- 10
numberRuns <- 5
searchIds <- c("Search 1","Search 2","Search 3","Search 4","Search 5","Search 6","Search 7","Search 8","Search 9","Search 10")
searchMethod <- "MS"
analysisParametersMS <- c(0.01,0.005,0.0025,0.001,0.00075,0.0005,0.00025,0.0001,0.00005,0.000025) 

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
        main="MS best found solutions",
        names=analysisParametersMS,
        xlab="Temperature", 
        ylab="Objective value") 

#--------------------------------------------------------------------------------------------------------------

#Experiment 3

numberSearches <- 1
numberRuns <- 1
searchIds <- c("Search 1")
searchMethod <- "RD"
analysisParametersMS <- c(100) 

phenoObsAnalysisOutput3 <- phenoObsAlgorithmPerformance(parentsPhenotypeDataFrame,
                                                        offspringPhenotype,
                                                        parentsGeneticSimilarityArray,
                                                        geneticSimilarityWeight,
                                                        penaltyReplicatedIndividuals,
                                                        inbreedingCoefficientMetric,
                                                        Gmatrix,
                                                        inbreedingCoefficientBoundsEstimatorSampleSize,
                                                        inclusionSet,
                                                        exclusionSet,
                                                        solutionSize,
                                                        maximumSearchTime,
                                                        timeWithoutImprovement,
                                                        maximumNumberMoves,
                                                        numberSearches,
                                                        numberRuns,
                                                        searchIds = as.character(c(1:numberSearches)),
                                                        searchMethod,
                                                        analysisParametersMS,
                                                        analysisParametersPT)
  
bestFoundObjectiveValues3 <- phenoObsAnalysisOutput3[[1]]
bestFoundIds3 <- phenoObsAnalysisOutput3[[2]]
bestFoundUpdatedObjectiveValues3 <- phenoObsAnalysisOutput3[[3]]
bestFoundUpdatedTimesObjectiveValues3 <- phenoObsAnalysisOutput3[[4]]
print(bestFoundObjectiveValues3)
print(bestFoundIds3)
print(bestFoundUpdatedObjectiveValues3)
print(bestFoundUpdatedTimesObjectiveValues3)


plot(bestFoundUpdatedTimesObjectiveValues3,bestFoundUpdatedObjectiveValues3)

print(mean(bestFoundValues3))
boxplot(t(bestFoundValues3),
        data=bestFoundValues3, 
        main="MS best found solutions with temperature T --> 0",
        names=analysisParametersMS,
        ylab="Objective value") 

###############################################################################################################

#--------------------------------------------------------------------------------------------
# (4) Conduct analysis for the RANDOM DESCENT algorithm
#--------------------------------------------------------------------------------------------

numberSearches <- 1
searchIds <- c("Search 1")
numberRuns <- 1
searchMethod <- "RD"

phenoObsAnalysisOutput4 <- phenoObsAlgorithmPerformance(parentsPhenotypeDataFrame,
                                                        offspringPhenotype,
                                                        parentsGeneticSimilarityArray,
                                                        geneticSimilarityWeight,
                                                        penaltyReplicatedIndividuals,
                                                        inbreedingCoefficientMetric,
                                                        Gmatrix,
                                                        inbreedingCoefficientBoundsEstimatorSampleSize,
                                                        inclusionSet,
                                                        exclusionSet,
                                                        solutionSize,
                                                        maximumSearchTime,
                                                        timeWithoutImprovement,
                                                        maximumNumberMoves,
                                                        numberSearches,
                                                        numberRuns,
                                                        searchIds = as.character(c(1:numberSearches)),
                                                        searchMethod,
                                                        analysisParametersMS,
                                                        analysisParametersPT)

bestFoundObjectiveValues4 <- phenoObsAnalysisOutput4[[1]]
bestFoundIds4 <- phenoObsAnalysisOutput4[[2]]
bestFoundUpdatedObjectiveValues4 <- phenoObsAnalysisOutput4[[3]]
bestFoundUpdatedTimesObjectiveValues4 <- phenoObsAnalysisOutput4[[4]]
print(bestFoundObjectiveValues4)
print(bestFoundIds4)
print(bestFoundUpdatedObjectiveValues4)
print(bestFoundUpdatedTimesObjectiveValues4)

plot(bestFoundUpdatedTimesObjectiveValues4,bestFoundUpdatedObjectiveValues4)





boxplot(t(bestFoundObjectiveValues4),
        data=bestFoundObjectiveValues4, 
        main="RD best found solutions",
        ylab="Objective value") 

###############################################################################################################

#--------------------------------------------------------------------------------------------
# (5) Conduct analysis for the parallel tempering algorithm
#--------------------------------------------------------------------------------------------

#Experiment 1

numberSearches <- 6
numberRuns <- 2
searchIds <- c("Search 1","Search 2","Search 3","Search 4","Search 5","Search 6")
searchMethod <- "PT"
analysisParametersPT <- matrix(nrow=numberSearches,ncol=3)
analysisParametersPT[,1] <- c(3,3,6,6,9,9) 
analysisParametersPT[,2] <- c(0.0001,0.0000000001,0.0001,0.0000000001,0.0001,0.0000000001) 
analysisParametersPT[,3] <- c(1,0.001,1,0.001,1,0.001) 

phenoObsAnalysisOutput5 <- phenoObsAlgorithmPerformance(parentsPhenotypeDataFrame,
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

bestFoundValues5 <- phenoObsAnalysisOutput5[[1]]
print(bestFoundValues5)
boxplot(t(bestFoundValues5),
        data=bestFoundValues5, 
        main="PT best found solutions",
        #names=analysisParametersMS,
        xlab="Search id", 
        ylab="Objective value") 
updatedBestFoundTimes5 <- phenoObsAnalysisOutput5[[4]]
print(updatedBestFoundTimes5[2,,])

#--------------------------------------------------------------------------------------------------------------

#Experiment 2

numberSearches <- 10
numberRuns <- 5
searchIds <- c("Search 1","Search 2","Search 3","Search 4","Search 5","Search 6","Search 7","Search 8","Search 9","Search 10")
searchMethod <- "PT"
analysisParametersPT <- matrix(nrow=numberSearches,ncol=3)
analysisParametersPT[,1] <- c(1,2,3,4,5,6,7,8,9,10) 
analysisParametersPT[,2] <- c(0.00000001,0.00000001,0.00000001,0.00000001,0.00000001,0.00000001,0.00000001,0.00000001,0.00000001,0.00000001) 
analysisParametersPT[,3] <- c(0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005) 

phenoObsAnalysisOutput6 <- phenoObsAlgorithmPerformance(parentsPhenotypeDataFrame,
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

bestFoundValues6 <- phenoObsAnalysisOutput6[[1]]
print(bestFoundValues6)
boxplot(t(bestFoundValues6),
        data=bestFoundValues6, 
        main="PT best found solutions",
        #names=analysisParametersMS,
        xlab="Search id", 
        ylab="Objective value") 
searchId <- 6
updatedBestFoundValues6 <- phenoObsAnalysisOutput6[[3]]
print(updatedBestFoundValues6[searchId,1,])
updatedBestFoundTimes6 <- phenoObsAnalysisOutput6[[4]]
print(updatedBestFoundTimes6[searchId,1,])
plot(updatedBestFoundTimes6[searchId,1,], 
     updatedBestFoundValues6[searchId,1,], 
     main="Search 10, run 1: convergence graph",
     xlab="Time", 
     ylab="Objective value",
     xlim=c(0, 500),
     ylim=c(min(updatedBestFoundValues6[searchId,1,])-0.02, 0.32),
     type='l', 
     col='black', 
     lwd=1)
globalOptimum <- c(rep(0.3133231,1000000))
lines(globalOptimum, 
      type='l', 
      col='blue', 
      lty=2,
      lwd=1)
searchMean <- c(rep(bestFoundValues6[searchId,],1000000))
lines(searchMean, 
      type='l', 
      col='dark green', 
      lty=2,
      lwd=1)

#--------------------------------------------------------------------------------------------------------------

#Experiment 3 (check parallel tempering with MS configuration)

numberSearches <- 1
numberRuns <- 15
searchIds <- c("Search 1")
searchMethod <- "PT"
analysisParametersPT <- matrix(nrow=numberSearches,ncol=3)
analysisParametersPT[,1] <- c(1) 
analysisParametersPT[,2] <- c(0.000000000000001) 
analysisParametersPT[,3] <- c(0.0005) 

phenoObsAnalysisOutput7 <- phenoObsAlgorithmPerformance(parentsPhenotypeDataFrame,
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

bestFoundValues7 <- phenoObsAnalysisOutput7[[1]]
print(bestFoundValues7)
boxplot(t(bestFoundValues7),
        data=bestFoundValues7, 
        main="PT best found solutions",
        #names=analysisParametersMS,
        xlab="Search id", 
        ylab="Objective value") 

