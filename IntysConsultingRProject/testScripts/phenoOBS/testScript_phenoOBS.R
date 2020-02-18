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
#(1)(b) Test the phenoOBSDataManagement function step by step
#-----------------------------------------------------------------------------

#Some pre-calculations
numberParents <- nrow(parentsPhenotypeDataFrame)
offspringPopulationSize <- (numberParents*(numberParents-1))/2
numberTraits <- ncol(parentsPhenotypeDataFrame)-1

#Perform several data configuration and parametric checks
if(any(is.na(parentsPhenotypeDataFrame)))
  stop("Missing data in parentsPhenotypeDataFrame")
if(any(is.na(parentsGeneticSimilarityMatrix)))
  stop("Missing data in parentsGeneticSimilarityMatrix")
if(any(parentsGeneticSimilarityMatrix < 0))
  stop("All data entries in parentsGeneticSimilarityMatrix must be greater than or equal to 0")
if(ncol(parentsGeneticSimilarityMatrix) != ncol(parentsGeneticSimilarityMatrix))
  stop("The relationship matrix must be square")
for(r in 1:numberTraits){
  if(traitObjectiveWeights[r,1] < 0 && !is.na(traitObjectiveWeights[r,1]))
    stop("Sub-objective weights may not take on negative values")
  if(traitObjectiveWeights[r,2] != 1 && traitObjectiveWeights[r,2] != 0 && !is.na(traitObjectiveWeights[,2]))
    stop("Sub-objectives categorisation may only take on binary values 0 or 1")}

#Fill in missing input data values 
for(r in 1:numberTraits){
  if(is.na(constraintParameters[r,1])){
    constraintParameters[r,1] <- -1000000000}
  if(is.na(constraintParameters[r,2])){
    constraintParameters[r,2] <- 1000000000}
  if(is.na(traitObjectiveWeights[r,1])){ 
    traitObjectiveWeights[r,1] <- 0  
    traitObjectiveWeights[r,2] <- 0}}

#Construct the overall offspring data frame
overallOffspringInformation <-  matrix(nrow = offspringPopulationSize, ncol = 3 + numberTraits)
index <- 1
for(i in 1:(numberParents-1)){
  for(j in (i+1):numberParents){
    overallOffspringInformation[index,1] <- index
    overallOffspringInformation[index,2] <- i
    overallOffspringInformation[index,3] <- j
    for(t in 1:numberTraits){
      overallOffspringInformation[index,3+t] <- (parentsPhenotypeDataFrame[i,t+1] + parentsPhenotypeDataFrame[j,t+1])/2
    }
    index <- index + 1
  }
}

#Construct the feasible offspring data frame and normalise the trait sub-objective scores
maxTraitValues <- c()
minTraitValues <- c()
for(t in 1:numberTraits){
  maxTraitValues[t] <- max(overallOffspringInformation[,3+t])
  minTraitValues[t] <- min(overallOffspringInformation[,3+t])
}
feasibleOffspringInformation <- matrix(nrow = offspringPopulationSize, ncol = 3 + numberTraits)
count <- 1
for(i in 1:offspringPopulationSize){
  infeasible <- 0
  trait <- 1
  while(infeasible == 0 && trait <= numberTraits){
    if(overallOffspringInformation[i,3+trait] < constraintParameters[trait,1] || overallOffspringInformation[i,3+trait] > constraintParameters[trait,2]){
      infeasible <- 1
    } 
    trait <- trait + 1
  }
  if(infeasible == 0){
    feasibleOffspringInformation[count,1] <- overallOffspringInformation[count,1]
    feasibleOffspringInformation[count,2] <- overallOffspringInformation[count,2]
    feasibleOffspringInformation[count,3] <- overallOffspringInformation[count,3]
    for(t in 1:numberTraits){
      feasibleOffspringInformation[count,t+3] <- (overallOffspringInformation[count,t+3] - minTraitValues[t]) / (maxTraitValues[t] - minTraitValues[t])
    } 
  } 
  count <- count + 1
}
feasibleOffspringInformation <- feasibleOffspringInformation[rowSums(is.na(feasibleOffspringInformation))!=ncol(feasibleOffspringInformation), ]

#Construct the feasible offspring data frame with offspring Ids, corresponding parents Ids and summed trait sub-objective evaluations
numberFeasibleIds <- nrow(feasibleOffspringInformation)
offspringInformationSimplified <- matrix(nrow = numberFeasibleIds, ncol = 4)
traitObjectiveWeightsAdjusted <- traitObjectiveWeights[1:numberTraits,1]*traitObjectiveWeights[1:numberTraits,2]
for(i in 1:numberFeasibleIds){
  offspringInformationSimplified[i,1] <- feasibleOffspringInformation[i,1]
  offspringInformationSimplified[i,2] <- feasibleOffspringInformation[i,2]
  offspringInformationSimplified[i,3] <- feasibleOffspringInformation[i,3]
  traitVectorValues <- feasibleOffspringInformation[i,4:(numberTraits+3)]
  offspringInformationSimplified[i,4] <- sum(traitObjectiveWeightsAdjusted * traitVectorValues)
}

#Convert the upper triangular part of the relationship matrix into an array and, again, normalise the corresponding values
maxGeneticRelationship <- max(parentsGeneticSimilarityMatrix[parentsGeneticSimilarityMatrix!=max(parentsGeneticSimilarityMatrix)])
minGeneticRelationship <- min(parentsGeneticSimilarityMatrix)
parentsGeneticSimilarityArray <- c()
count <- 1
for(i in 1:(numberParents-1)){
  for(j in (i+1):numberParents){
    parentsGeneticSimilarityArray[count] <- (parentsGeneticSimilarityMatrix[i,j] - minGeneticRelationship) / (maxGeneticRelationship - minGeneticRelationship)
    count <- count + 1
  }
}

#check data configuration display
print(parentsPhenotypeDataFrame)
print(parentsGeneticSimilarityMatrix)
print(constraintParameters)
print(traitObjectiveWeights)
print(overallOffspringInformation)
print(maxTraitValues)
print(minTraitValues)
print(feasibleOffspringInformation)
print(offspringInformationSimplified)
print(parentsGeneticSimilarityArray)

#-----------------------------------------------------------------------------
#(1)(c) Test the phenoOBSDataManagement function at once
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
inclusionSet <- NULL
exclusionSet <- NULL
solutionSize <- 2
maximumSearchTime <- 500
timeWithoutImprovement <- 200
maximumNumberMoves <- 400
searchMethod <- "MS"
searchParametersMS <- 5
searchParametersPT <- c(6,0.001,9)

#CHANGE name of GMatrix to pedigree-based matrix values mapped on the range [0,2]?

#-----------------------------------------------------------------------------
#(2)(b) Test the phenoOBSOptimisation function step by step
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
if(searchMethod != "MS" && searchMethod  != "RD" && searchMethod != "PT")
  stop("Search methodology not recognised")
if(length(searchParametersMS) != 1 || length(searchParametersPT) != 3)
  stop("Incorrect number of arguments provided for the selected search methodology")
if(searchParametersMS < 0)
  stop("The MS initial temperature setting must be larger than zero")
if(!(as.double(searchParametersPT[1]) - as.integer(searchParametersPT[1]) == 0) || searchParametersPT[1] < 1)
  stop("The PT number of parallel processors must be a positive integer")
if(searchParametersPT[2] >= searchParametersPT[3])
  stop("The PT initial temperature setting must be larger than the final temperature setting")
if(searchParametersPT[2] <= 0)
  stop("The PT initial temperature setting must be larger than zero")

#Dismiss the offspring ids listed in the exclusion set 
for(i in 1:nrow(offspringPhenotype)){
  if(is.element(i,exclusionSet)){
    offspringPhenotype <- offspringPhenotype[-which(offspringPhenotype[,1] == i),]
  }
}

#Convert the R (feasible) inclusion set ids to Java ids
inclusionSetJava <- c()
for(k in 1:length(inclusionSet)){
  inclusionSetJava[k] <- which(offspringPhenotype[,1] == inclusionSet[k]) - 1
}

#Translate the required R objects into Java Objects.
breedingPopulationSizeJ <- as.integer(numberParents, dispatch=TRUE)
offspringPhenotypeJ <- .jarray(as.matrix(offspringPhenotype, nrow=nrow(offspringPhenotype), ncol=4), dispatch=TRUE)
parentsGeneticSimilarityArrayJ <- .jarray(parentsGeneticSimilarityArray, dispatch=TRUE)  
geneticSimilarityWeightJ <- as.double(geneticSimilarityWeight, dispatch=TRUE)
penaltyReplicatedIndividualsJ <- as.double(penaltyReplicatedIndividuals, dispatch=TRUE)
if(isTRUE(inbreedingCoefficientMetric)){
  inbreedingCoefficientMetricJ <- new(J("java.lang.String"), "TRUE")
} else {
  inbreedingCoefficientMetricJ <- new(J("java.lang.String"), "FALSE")
}
GmatrixJ <- .jarray(as.matrix(Gmatrix, nrow=nrow(Gmatrix), ncol=ncol(Gmatrix)), dispatch=TRUE)
if(is.null(inclusionSetJava)){
  inclusionSetJavaJ <- .jnull()
} else {
  inclusionSetJavaJ <- .jarray(as.integer(inclusionSetJava), dispatch=TRUE)
}
inbreedingCoefficientBoundsEstimatorSampleSizeJ <- as.integer(inbreedingCoefficientBoundsEstimatorSampleSize, dispatch=TRUE)
solutionSizeJ <- as.integer(solutionSize, dispatch=TRUE)
maximumSearchTimeJ <- as.integer(maximumSearchTime, dispatch=TRUE)
timeWithoutImprovementJ <- as.integer(timeWithoutImprovement, dispatch=TRUE)
maximumNumberMovesJ <- as.integer(maximumNumberMoves, dispatch=TRUE)

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

#Configure the search object
if(searchMethod == "MS"){
  print("Configuring MS search")
  initialTemperatureJ <- as.double(searchParametersMS, dispatch=TRUE)
  searchConfiguration <- .jrcall("jamesgeneral/MSSearchParameters", method="createMSSearchParameters",
                                 searchTerminationCriteria,
                                 initialTemperatureJ)
} else if (searchMethod == "RD"){
  print("Configuring RD search")
  searchConfiguration <- .jrcall("jamesgeneral/RDSearchParameters",method="createRDSearchParameters",
                                 searchTerminationCriteria)
} else{
  print("Configuring PT search")
  numberProcessorsJ <- as.integer(searchParametersPT[1], dispatch=TRUE)
  finalTemperatureJ <- as.double(searchParametersPT[2], dispatch=TRUE)
  initialTemperatureJ <- as.double(searchParametersPT[3], dispatch=TRUE)
  searchConfiguration <- .jrcall("jamesgeneral/PTSearchParameters", method="createPTSearchParameters",
                                 searchTerminationCriteria,
                                 numberProcessorsJ,
                                 finalTemperatureJ,
                                 initialTemperatureJ)}

#Create the search object
search <- .jrcall("phenoObs/ObsFactory", method="createSearch",
                  problemParameters, searchConfiguration)

#Launch the search 
.jrcall(search, method="start")

#Extract the best found solution and algorithm run details
runTime <- .jrcall(search, method="getRuntime") 
convergenceTime <- runTime - .jrcall(search, method="getTimeWithoutImprovement")
solution <- .jrcall(search, method="getBestSolution")
bestFoundSolutionIdsJ <- .jrcall(solution, "getSelectedIDs")
bestFoundSolutionIdsR <- .jrcall(inputData, "getObsIdsR", bestFoundSolutionIdsJ)
bestFoundSolutionGeneticDiversity <- .jrcall(inputData,  method = "getObsGeneticSimilarity", bestFoundSolutionIdsJ)
solutionEvaluation <- .jrcall(search, method="getBestSolutionEvaluation")

#-----------------------------------------------------------------------------
#(2)(c) Test the phenoOBSOptimisation function at once
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
#(3)(b) Test the phenoOBSOutput function step by step
#-----------------------------------------------------------------------------

#Construct a data frame containing information on the sub-objective score attainments of the best found solution
parentNames <- parentsPhenotypeDataFrame[,1]
numberTraits <- ncol(overallOffspringInformation) - 3
traitNames <- colnames(parentsPhenotypeDataFrame)
traitNames <- traitNames[2:(numberTraits+1)]
bestFoundSolutionObjectiveScores <- matrix(nrow = 1, ncol = numberTraits+1)
for(t in 1:numberTraits){
  score <- 0
  for(i in 1:solutionSize){
    score <- score + overallOffspringInformation[optimalSolutionIdsR[i],3+t]
  }
  bestFoundSolutionObjectiveScores[t] <- score/solutionSize
} 
bestFoundSolutionObjectiveScores[numberTraits+1] <- bestFoundSolutionGeneticDiversity
DFbestFoundSolutionObjectiveScores = data.frame(bestFoundSolutionObjectiveScores)   
colnames(DFbestFoundSolutionObjectiveScores) <- c(traitNames, "Gen sim")

#Finally, construct a data frame summarising information on the individual ids in the best found solution
bestFoundSolutionIdsInformation <- matrix(nrow = solutionSize, ncol = 6+numberTraits)
for(i in 1:solutionSize){
  bestFoundSolutionIdsInformation[i,1:3] <- overallOffspringInformation[optimalSolutionIdsR[i],1:3]
  bestFoundSolutionIdsInformation[i,4] <- toString(parentNames[overallOffspringInformation[optimalSolutionIdsR[i],2]])
  bestFoundSolutionIdsInformation[i,5] <- toString(parentNames[overallOffspringInformation[optimalSolutionIdsR[i],3]])
  bestFoundSolutionIdsInformation[i,6:(5+numberTraits)] <- as.numeric(overallOffspringInformation[optimalSolutionIdsR[i],4:(3+numberTraits)])
  bestFoundSolutionIdsInformation[i,6+numberTraits] <- as.numeric(parentsGeneticSimilarityMatrix[overallOffspringInformation[optimalSolutionIdsR[i],2],overallOffspringInformation[optimalSolutionIdsR[i],3]])
}
colnames(bestFoundSolutionIdsInformation) <- c("offspring_id", "parent1_id", "parent2_id", "parent1_name", "parent2_name", traitNames, "Gen sim")
DFbestFoundSolutionIdsInformation = data.frame(bestFoundSolutionIdsInformation)  
colnames(DFbestFoundSolutionIdsInformation) <- c("offspring_id", "parent1_id", "parent2_id", "parent1_name", "parent2_name", traitNames, "Gen sim")
print(DFbestFoundSolutionObjectiveScores)
print(DFbestFoundSolutionIdsInformation)

#-----------------------------------------------------------------------------
#(3)(c) Test the phenoOBSOutput function at once
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
maximumSearchTime <- 500
timeWithoutImprovement <- 200
maximumNumberMoves <- 400
numberSearches <- 2
numberRuns <- 2
searchIds <- c("Search 1","Search 2")
searchMethod <- "PT"
analysisParametersMS <- c(8,0.5)
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

#Translate the remaining required R objects into Java Objects.
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





