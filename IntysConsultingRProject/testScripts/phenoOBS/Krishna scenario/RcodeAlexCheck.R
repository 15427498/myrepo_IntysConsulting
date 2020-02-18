
library(BCS.JAMES3) 
library(MASS)
library(BCS.Base)
library(DiGGer)
library(xlsx)
sessionInfo()

setwd("C:/Users/Gebruiker/Documents/Intys/RStudio/IntysConsultingRProject/testScripts/phenoOBS/Input")

#-----------------------------------------------------------------------------
# Read the .csv data
#-----------------------------------------------------------------------------

dat <- read.csv("phenotypes.csv", header=TRUE)
#dat <- clean(dat)
geno <- read.csv("genotypes.csv", header=TRUE)
#mhead(geno)
rownames(geno) <- geno$X
geno <- geno[, -which(colnames(geno) == "X")]
geno <- as.matrix(geno)

###############################################################################################################

#-----------------------------------------------------------------------------
#(1)(a) Generate the phenoOBSDataManagement data
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

#-----------------------------------------------------------------------------
#(1)(b) Test the phenoOBSDataManagement function 
#-----------------------------------------------------------------------------

javaInputData <- phenoObsDataManagement(parentsPhenotypeDataFrame,
                                        parentsGeneticSimilarityMatrix,
                                        constraintParameters,
                                        traitObjectiveWeights)
print(javaInputData)

###############################################################################################################

#-----------------------------------------------------------------------------
#(2)(a) Generate the phenoOBSOptimisation data
#-----------------------------------------------------------------------------

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
searchMethod <- "PT"
searchParametersMS <- 5
searchParametersPT <- c(6,0.001,9)

#-----------------------------------------------------------------------------
#(2)(b) Test the phenoOBSOptimisation function 
#-----------------------------------------------------------------------------

phenoObsOptimsationResults <- phenoObsOptimisation(parentsPhenotypeDataFrame,
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
                                                   searchMethod,
                                                   searchParametersMS,
                                                   searchParametersPT)
print(phenoObsOptimsationResults)

###############################################################################################################

#-----------------------------------------------------------------------------
#(3)(a) Generate the phenoOBSOutput data
#-----------------------------------------------------------------------------

parentsPhenotypeDataFrame <- javaInputData[[3]]
optimalSolutionIdsR <- phenoObsOptimsationResults[[1]]
bestFoundSolutionGeneticDiversity <- phenoObsOptimsationResults[[2]]
overallOffspringInformation <- javaInputData[[4]]
solutionSize <- phenoObsOptimsationResults[[3]]
parentsGeneticSimilarityMatrix <- javaInputData[[5]]

#-----------------------------------------------------------------------------
#(3)(b) Test the phenoOBSOutput function 
#-----------------------------------------------------------------------------

phenoObsSummaryOutput <- phenoObsOutput(parentsPhenotypeDataFrame,
                                        optimalSolutionIdsR,
                                        bestFoundSolutionGeneticDiversity,
                                        overallOffspringInformation,
                                        solutionSize,
                                        parentsGeneticSimilarityMatrix)
print(phenoObsSummaryOutput)














