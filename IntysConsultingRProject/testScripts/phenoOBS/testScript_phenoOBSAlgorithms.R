library(IntysConsulting) 
library(MASS)
library(james.analysis)
library(rJava)

.jinit()
.jaddClassPath("C:\\Users\\Gebruiker\\Documents\\Intys\\RStudio\\IntysConsultingRProject\\pkg\\inst\\java\\BCSJames3.jar")

#-----------------------------------------------------------------------------
#(1)(a) Generate the phenoOBSDataManagement data
#-----------------------------------------------------------------------------

parentsPhenotypeDataFrame <- data.frame(
  parentNames = c("A","B","C","D"),
  Yield = c(1200,1300,1250,1150),
  Lint = c(50,40,45,55),
  BollType = c(4.5,4,3,4.5))
parentsGeneticSimilarityMatrix <- matrix(nrow = 4, ncol = 4)
parentsGeneticSimilarityMatrix[1,] = c(1,0.60,0.40,0.30)
parentsGeneticSimilarityMatrix[2,] = c(0.60,1,0.50,0.35)
parentsGeneticSimilarityMatrix[3,] = c(0.40,0.50,1,0.55)
parentsGeneticSimilarityMatrix[4,] = c(0.30,0.35,0.55,1)
constraintParameters <- matrix(nrow = 3, ncol = 2)
constraintParameters[1,1] <- 1190
traitObjectiveWeights <- matrix(nrow = 3, ncol = 2)
traitObjectiveWeights[1,] = c(0.7,1)
traitObjectiveWeights[2,] = c(0.3,1)
traitObjectiveWeights[3,] = c(0,1)

#-----------------------------------------------------------------------------
#(1)(b) Test the phenoOBSDataManagement function at once
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
penaltyReplicatedIndividuals <- 1.08
inbreedingCoefficientMetric <- FALSE
Gmatrix <- matrix(nrow = 4, ncol = 4)
Gmatrix[1,] = c(1,0.60,0.40,0.30)
Gmatrix[2,] = c(0.60,1,0.50,0.35)
Gmatrix[3,] = c(0.40,0.50,1,0.55)
Gmatrix[4,] = c(0.30,0.35,0.55,1)
inbreedingCoefficientBoundsEstimatorSampleSize <- 300
inclusionSet <- 1
exclusionSet <- NULL
solutionSize <- 2
maximumSearchTime <- 500
timeWithoutImprovement <- 200
maximumNumberMoves <- 400
searchMethod <- "MS"
searchParametersMS <- 5
searchParametersPT <- c(6,0.001,9)

#-----------------------------------------------------------------------------
#(2)(b) Test the phenoOBSOptimisation function at once
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
#(3)(b) Test the phenoOBSOutput function at once
#-----------------------------------------------------------------------------

phenoObsSummaryOutput <- phenoObsOutput(parentsPhenotypeDataFrame,
                                        optimalSolutionIdsR,
                                        bestFoundSolutionGeneticDiversity,
                                        overallOffspringInformation,
                                        solutionSize,
                                        parentsGeneticSimilarityMatrix)
print(phenoObsSummaryOutput)


#-----------------------------------------------------------------------------
#(4)(a) Generate the phenoOBSAlgorithmPerformance data (test for PT)
#-----------------------------------------------------------------------------

offspringPhenotype <- javaInputData[[1]]
parentsGeneticSimilarityArray <- javaInputData[[2]]
parentsPhenotypeDataFrame <- javaInputData[[3]]
geneticSimilarityWeight <- 0.4
penaltyReplicatedIndividuals <- 1.08
inclusionSet <- c(4)
solutionSize <- 2
maximumSearchTime <- 5000
timeWithoutImprovement <- 2000
maximumNumberMoves <- 40000
numberSearches <- 4
numberRuns <- 5
searchIds <- c("Search 1","Search 2","Search 3","Search 4")
searchMethod <- "RD"
analysisParametersMS <- c(200,2,0.02,0.0002)
analysisParametersPT <- matrix(nrow=numberSearches, ncol=3)
analysisParametersPT[1,] <- c(5,0.1,6)
analysisParametersPT[2,] <- c(5,0.001,2)

#-----------------------------------------------------------------------------
#(4)(b) Test the phenoOBSAlgorithmPerformance function step by step
#-----------------------------------------------------------------------------

#Some pre-calculations
numberParents <- nrow(parentsPhenotypeDataFrame)
offspringPopulationSize <- (numberParents*(numberParents-1))/2

#Perform several parametric checks
if(geneticSimilarityWeight < 0 || geneticSimilarityWeight > 1)
  stop("The genetic similarity sub-objective Weight should be selected between 0 and 1")
if(!(penaltyReplicatedIndividuals >= 1))
  stop("The penalty on replicated individuals has a minimum value of 1")
if(!(is.null(inclusionSet))){
  if(length(inclusionSet) > solutionSize)
    stop("The size of the inclusion set is larger than the specified solution size")
  if(any(!(inclusionSet %in% offspringPhenotype[,1])))
    stop("At least one id elements from the inclusion set is infeasible")
}
if(solutionSize < 2 || solutionSize >= nrow(offspringPhenotype))
  stop("Infeasible solution size")
if(!(maximumSearchTime > 0))
  stop("Infeasible or unspecified maximum search time value")
if(!(timeWithoutImprovement > 0) && !is.null(timeWithoutImprovement))
  stop("Infeasible time without improvement value")
if(!(maximumNumberMoves > 0) && !is.null(maximumNumberMoves))
  stop("Infeasible maximum number of moves selected")
if(numberSearches != length(searchIds))
  stop("The length of searchIds does not match the numberSearches argument")
if(searchMethod != "MS" && searchMethod  != "RD" && searchMethod != "PT")
  stop("Search methodology not recognised")
if(searchMethod == "MS"){
  analysisParametersMSJ <- .jarray(analysisParametersMS, dispatch=TRUE)
  if(length(analysisParametersMS) != numberSearches)
    stop("Incorrect number of entries in analysisParametersMS")
  if(any(analysisParametersMS < 0))
    stop("The MS initial temperature settings must be larger than zero")}
if(searchMethod == "PT"){
  analysisParametersPTJ <- .jarray(as.matrix(analysisParametersPT, nrow=numberSearches, ncol=3), dispatch=TRUE)
  if(nrow(analysisParametersPT) != numberSearches || ncol(analysisParametersPT) != 3)
    stop("Incorrect number of entries in searchParametersPT")
  if(any(!(as.double(analysisParametersPT[,1]) - as.integer(analysisParametersPT[,1]) == 0)) || analysisParametersPT[,1] < 1)
    stop("The PT number of parallel processors must be a positive integer")
  if(any(analysisParametersPT[,2] >= analysisParametersPT[,3]))
    stop("The PT initial temperature setting must be larger than the final temperature setting")
  if(any(analysisParametersPT[,2] <= 0))
    stop("The PT initial temperature setting must be larger than zero")}

#Convert the (feasible) inclusion set ids to the correct Java ids
inclusionSetJava <- c()
for(k in 1:length(inclusionSet)){
  inclusionSetJava[k] <- which(offspringPhenotype[,1] == inclusionSet[k])
}

#Dummy value for redundant argument (see assumptions)
problemIds <- c("Problem 1")
numberProblems <- 1

#Translate the required R objects into Java Objects.
breedingPopulationSizeJ <- as.integer(numberParents, dispatch=TRUE)
offspringPhenotypeJ <- .jarray(as.matrix(offspringPhenotype, nrow=nrow(offspringPhenotype), ncol=4), dispatch=TRUE)
parentsGeneticSimilarityArrayJ <- .jarray(parentsGeneticSimilarityArray, dispatch=TRUE)
geneticSimilarityWeightJ <- as.double(geneticSimilarityWeight, dispatch=TRUE)
penaltyReplicatedIndividualsJ <- as.double(penaltyReplicatedIndividuals, dispatch=TRUE)
if(is.null(inclusionSetJava)){
  inclusionSetJavaJ <- .jnull()
} else {
  inclusionSetJavaJ <- .jarray(inclusionSetJava, dispatch=TRUE)
}
solutionSizeJ <- as.integer(solutionSize, dispatch=TRUE)
maximumSearchTimeJ <- as.integer(maximumSearchTime, dispatch=TRUE)
timeWithoutImprovementJ <- as.integer(timeWithoutImprovement, dispatch=TRUE)
maximumNumberMovesJ <- as.integer(maximumNumberMoves, dispatch=TRUE)
numberProblemsJ <- as.integer(numberProblems, dispatch=TRUE)
numberSearchesJ <- as.integer(numberSearches, dispatch=TRUE)
numberRunsJ <- as.integer(numberRuns, dispatch=TRUE)
problemIdsJ <- .jarray(problemIds, dispatch=TRUE)
searchIdsJ <- .jarray(searchIds, dispatch=TRUE)

#Create the java input data object 
inputData <- .jrcall("phenoObs/ObsData", method="createObsData", 
                     breedingPopulationSizeJ,
                     offspringPhenotypeJ, 
                     parentsGeneticSimilarityArrayJ,
                     geneticSimilarityWeightJ,
                     penaltyReplicatedIndividualsJ,
                     inbreedingCoefficientMetricJ,
                     GmatrixJ,
                     inbreedingCoefficientBoundsEstimatorSampleSizeJ,
                     solutionSizeJ,
                     inclusionSetJavaJ)

#Create the java problem parameters object
problemParameters <- .jrcall("jamesgeneral/SingleSubsetProblemParameters", method="createSingleSubsetProblemParameters",
                             inputData, solutionSizeJ)

#Create the java search termination criteria object
searchTerminationCriteria <- .jrcall("jamesgeneral/SearchParameters", method="createSearchParameters",
                                     maximumSearchTimeJ,
                                     timeWithoutImprovementJ,
                                     maximumNumberMovesJ)

#Create the analysis object
analysisParameters <- .jrcall("jamesgeneral/AnalysisParameters", method="createAnalysisParameters",
                              numberProblemsJ,
                              numberSearchesJ,
                              numberRunsJ,
                              problemIdsJ,
                              searchIdsJ)

#Configure the analysis objects as per the specified search method
if(searchMethod == "MS"){
  initialTemperatureJ <- as.double(10)  #This dummy temperature is used only in launching the Java process; it is redundant in the analysis
  searchConfiguration <- .jrcall("jamesgeneral/MSSearchParameters", method="createMSSearchParameters",
                                 searchTerminationCriteria,
                                 initialTemperatureJ)
  analysisConfiguration <- .jrcall("jamesgeneral/MSAnalysisParameters", method="createMSAnalysisParameters",
                                   analysisParameters,
                                   analysisParametersMSJ)
  analysis <- .jrcall("phenoObs/ObsAnalysis", method="getAnalysis",
                      problemParameters,
                      searchConfiguration,
                      analysisConfiguration)

} else if (searchMethod == "RD"){
  searchConfiguration <- .jrcall("jamesgeneral/RDSearchParameters", method="createRDSearchParameters",
                                 searchTerminationCriteria)
  analysisConfiguration <- .jrcall("jamesgeneral/RDAnalysisParameters", method="createRDAnalysisParameters",
                                   analysisParameters)
  analysis <- .jrcall("phenoObs/ObsAnalysis", method="getAnalysis",
                      problemParameters,
                      searchConfiguration,
                      analysisConfiguration)

} else{
  numberProcessorsJ <- as.integer(5) #Dummy arguments, same as for MS
  finalTemperatureJ <- as.double(8)
  initialTemperatureJ <- as.double(1)
  searchConfiguration <- .jrcall("jamesgeneral/PTSearchParameters", method="createPTSearchParameters",
                                 searchTerminationCriteria,
                                 numberProcessorsJ,
                                 finalTemperatureJ,
                                 initialTemperatureJ)
  analysisConfiguration <- .jrcall("jamesgeneral/PTAnalysisParameters", method="createPTAnalysisParameters",
                                   analysisParameters,
                                   analysisParametersPTJ)
  analysis <- .jrcall("phenoObs/ObsAnalysis", method="getAnalysis",
                      problemParameters,
                      searchConfiguration,
                      analysisConfiguration)}

#get the analysis results object and extract the relevant information
results <- .jrcall("phenoObs/ObsAnalysis", method="getAnalysisResults", analysis)
numberSearches <- .jrcall(results, method="getNumSearches", problemIds[1])
numberRuns <- .jrcall(results, method="getNumRuns", problemIds[1], searchIds[1])

#extract and store the runs information
# - TEST (with search1 -> run1)
run1 <- .jrcall(results, method="getRun", problemIds[1], searchIds[1], as.integer(1))
run1BestFoundValue <- .jrcall("phenoObs/ObsAnalysis", method="getBestFoundValue", run1)
print(run1BestFoundValue)
run1BestFoundSolutionJ <- .jrcall(run1, method="getBestSolution")
run1BestFoundSolutionIdsJ <- .jrcall(run1BestFoundSolutionJ, "getSelectedIDs")
run1BestFoundSolutionIdsR <- .jrcall(inputData, "getObsIdsR", run1BestFoundSolutionIdsJ)
print(run1BestFoundSolutionIdsR)
run1UpdatedValues <- .jrcall("phenoObs/ObsAnalysis", method="getRunValues", run1)
print(run1UpdatedValues)
run1UpdatedTimes <- .jrcall("phenoObs/ObsAnalysis", method="getRunTimes", run1)
print(run1UpdatedTimes)

# - GENERIC
searchRunBestFoundValues <- matrix(nrow=numberSearches, ncol=numberRuns)
searchRunBestFoundSolutions <- array(rep(NA, numberSearches*numberRuns*solutionSize), dim=c(numberSearches, numberRuns, solutionSize))
searchRunUpdatedValuesLengths <- matrix(nrow=numberSearches, ncol=numberRuns)
for(s in 1:numberSearches){
  for(r in 1:numberRuns){
    run <- .jrcall(results, method="getRun", problemIds[1], searchIds[s], as.integer(r-1))
    searchRunBestFoundValues[s,r] <- .jrcall("phenoObs/ObsAnalysis", method="getBestFoundValue", run)
    searchRunBestFoundSolutionJ <- .jrcall(run, method="getBestSolution")
    searchRunBestFoundSolutionIdsJ <- .jrcall(searchRunBestFoundSolutionJ, "getSelectedIDs")
    searchRunBestFoundSolutions[s,r,] <- .jrcall(inputData, "getObsIdsR", searchRunBestFoundSolutionIdsJ)
    searchRunUpdatedValuesLengths[s,r] <- length(.jrcall("phenoObs/ObsAnalysis", method="getRunValues", run))
  }
}
maximumSearchRunUpdatedValuesLengths <- max(searchRunUpdatedValuesLengths)
searchRunUpdatedValues <- array(rep(NA, numberSearches*numberRuns*maximumSearchRunUpdatedValuesLengths), dim=c(numberSearches, numberRuns, maximumSearchRunUpdatedValuesLengths))
searchRunUpdatedTimes <- array(rep(NA, numberSearches*numberRuns*maximumSearchRunUpdatedValuesLengths), dim=c(numberSearches, numberRuns, maximumSearchRunUpdatedValuesLengths))
for(s in 1:numberSearches){
  for(r in 1:numberRuns){
    run <- .jrcall(results, method="getRun", problemIds[1], searchIds[s], as.integer(r-1))
    vectorLength <- searchRunUpdatedValuesLengths[s,r]
    if(vectorLength < maximumSearchRunUpdatedValuesLengths){
      searchRunUpdatedValues[s,r,1:vectorLength] <- .jrcall("phenoObs/ObsAnalysis", method="getRunValues", run)
      searchRunUpdatedValues[s,r,(vectorLength+1):maximumSearchRunUpdatedValuesLengths] <- searchRunUpdatedValues[s,r,vectorLength]
      searchRunUpdatedTimes[s,r,1:vectorLength] <- .jrcall("phenoObs/ObsAnalysis", method="getRunTimes", run)
      searchRunUpdatedTimes[s,r,(vectorLength+1):maximumSearchRunUpdatedValuesLengths] <- searchRunUpdatedTimes[s,r,vectorLength]

    } else{
      searchRunUpdatedValues[s,r,] <- .jrcall("phenoObs/ObsAnalysis", method="getRunValues", run)
      searchRunUpdatedTimes[s,r,] <- .jrcall("phenoObs/ObsAnalysis", method="getRunTimes", run)
    }
  }
}

#-----------------------------------------------------------------------------
#(4)(c) Test the phenoOBSAlgorithmPerformance function at once
#-----------------------------------------------------------------------------

phenoObsAnalysisOutput <- phenoObsAlgorithmPerformance(parentsPhenotypeDataFrame,
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
print(phenoObsAnalysisOutput)





