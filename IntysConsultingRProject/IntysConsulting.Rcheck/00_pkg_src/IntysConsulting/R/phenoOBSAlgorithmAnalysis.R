#' Algorithm performances analysis of pheno OBS.
#'
#' Run an analysis to investigate the parameter tuning sensitivity and convergence properties of the studied algorithms.
#' 
#' @param parentsPhenotypeDataFrame (required) A data frame containing the parentNames in the first column, the trait names as headings in the remianing columns, and the parent phenotype information in the remaining matrix entries.
#' See example below for more clarity.
#' This data frame may also be obtained from the phenoOBSDataManagement output.
#' @param offspringPhenotype (required) A matrix containing:
#' * The feasible offspring Ids in the first column,
#' * The corresponding offspring parents in the second and third columns, and
#' * The weighted, normalised trait objective evaluation of the offspring.
#' This Matrix is obtained from the phenoOBSDataManagement output.
#' @param parentsGeneticSimilarityArray (required) An array containing the normalised upper triangular entries of the parents genetic similarity matrix.
#' This Matrix is obtained from the phenoOBSDataManagement output.
#' @param geneticSimilarityWeight (required) A value in the range [0,1] for the importance of solution similarity relative to the traits.
#' For example, if specified to be 0.3, then the algorithm will automatically allocate 0.7 toward finding solutions that are strong across the entirety of the traits sub-objective.
#' Note, however, that the individual trait importances are in phenoOBSDataManagement.
#' @param penaltyReplicatedIndividuals (default, double) A soft-constraint value penalising the presence of replicated individuals in a solution.
#' Note that this value is not normalised alongside the set of upper triangular matrix entries, and so it must have at least a value of one.
#' @param inclusionSet (default) An array containing the integer ids of offsping crosses that the user desires to observe in the returned solution.
#' The default value is a null vector, so that no id inclusion constraint is imposed on the selection of crosses.
#' @param solutionSize (required, integer) The number of crosses to be extracted.
#' @param maximumSearchTime (required, integer) The total duration of each search run.
#' This parameter is configured in milliseconds.
#' @param timeWithoutImprovement (default, integer) An upper bound on the total duration of the optimisation process.
#' This parameter is configured in milliseconds.
#' @param maximumNumberMoves (default, integer) An upper bound on the total number of search space solutions visited and assessed during the course of the optimisation process.
#' This parameter is configured in milliseconds.
#' @param numberSearches (required, integer) The total number of search method configurations probed in the analysis.
#' @param numberRuns (required, integer) The total number of algorithm runs performed for each search.
#' @param searchIds (default, String array) A [1*S] vector of strings used as identifyers for the various search method configurations.
#' The default argument consists of the string array ("1", "2", ... , "S").
#' @param searchMethod (required, string) Summon either one of these three proposed algorithms to analyse: 
#' 1. Random descent (enter "RD"); 
#' 2. Metropolis search (enter "MS") or 
#' 3. Parallel  tempering (enter "PT").
#' @param analysisParametersMS (default, N.A.) A [1*S] vector specifying the initial temperatures of the Metropolis search algorithm.
#' @param analysisParametersPT (default, N.A.) A [S*3] data vector specifying the search parameters of the parallel tempering search algorithm.
#' In the first, second and third columns (respectively):
#' 1. processors (integer) : The number of processors running in parallel.
#' 2. minTemp (double) : Minimum (i.e. final) temperature.
#' 3. maxTemp (double) : Maximum (i.e. initial) temperature.
#' 
#' @details Symbols notation:
#' * S : The total number of searches in the analysed scenario.
#' * R : The total number of runs in each search.
#' * n : The size of the solution.
#' Assumptions: (any of these may be altered at a later stage if necessary, as per user request) 
#' 1. This function operates for one problem at a time .
#' 2. The number of runs is the same across all searches (this may be altered at a later stage, as per user request).
#' 3. Each run is assigned the same termination criteria.
#' 
#' @return A list containing the following four objects (respectively): 
#' 1. A 2D array containing the best found objective values for each run (columns dimension) and for each search (rows dimension),
#' 2. A 3D array containing the best found solution ids (depth dimension) for each run (columns dimension) and for each search (rows dimension),
#' 3. A 3D array containing the updated new best found solution objective values throughout the algorithm (depth dimension) for each run (columns dimension) and for each search (rows dimension), and
#' 4. A 3D array containing the updated new best found solution objective times throughout the algorithm (depth dimension) for each run (columns dimension) and for each search (rows dimension).
#' Note: the last value across the depth dimension of this array in any given search and run therefore corresponds to the convergence time.
#' 
#' @examples 
#' parentsPhenotypeDataFrame <- data.frame(
#' parentNames = c("A","B","C","D"),
#' Yield = c(1200,1300,1250,1150),
#' Lint = c(50,40,45,55),
#' BollType = c(4.5,4,3,4.5))
#' offspringPhenotype <- matrix(nrow = 5, ncol = 4)
#' offspringPhenotype[1,] <- c(1,1,2,0.6)
#' offspringPhenotype[2,] <- c(2,1,3,0.5)
#' offspringPhenotype[3,] <- c(4,2,3,0.7)
#' offspringPhenotype[4,] <- c(5,2,4,0.5)
#' offspringPhenotype[5,] <- c(6,3,4,0.4)
#' parentsGeneticSimilarityArray <- c(1,0.33333333,0,0.66666667,0.166666667,0.8333333333)
#' geneticSimilarityWeight <- 0.4
#' penaltyReplicatedIndividuals <- 1.08
#' inclusionSet <- c(4)
#' solutionSize <- 2
#' maximumSearchTime <- 500
#' timeWithoutImprovement <- 200
#' maximumNumberMoves <- 400
#' numberSearches <- 2
#' numberRuns <- 2
#' searchIds <- c("Search 1","Search 2")
#' searchMethod <- "MS"
#' analysisParametersMS <- c(8,0.5) 
#' analysisParametersPT <- matrix(nrow=numberSearches, ncol=3)
#' analysisParametersPT[1,] <- c(5,0.1,6)
#' analysisParametersPT[2,] <- c(5,0.001,2)
#' 
#' @author A.Colmant and G.De Meyer
#' @export
#' 
phenoObsAlgorithmPerformance <- function(parentsPhenotypeDataFrame,
                                 offspringPhenotype,
                                 parentsGeneticSimilarityArray,
                                 geneticSimilarityWeight,
                                 penaltyReplicatedIndividuals,
                                 inclusionSet = NULL,
                                 solutionSize,
                                 maximumSearchTime = 60000,
                                 timeWithoutImprovement = 10000,
                                 maximumNumberMoves = 90000000,
                                 numberSearches,
                                 numberRuns,
                                 searchIds = as.character(c(1:numberSearches)),
                                 searchMethod,
                                 analysisParametersMS = rep(NA, numberSearches),
                                 analysisParametersPT = matrix(NA, nrow=numberSearches, ncol=3)){
  
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
    if(any(analysisParametersPT[,3] <= analysisParametersPT[,2]))
      stop("The PT initial temperature setting must be larger than the final temperature setting")
    if(any(analysisParametersPT[,3] <= 0))
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
  
  analysisOutput <- list(searchRunBestFoundObjectiveValues = searchRunBestFoundValues, 
                         searchRunBestFoundSolutionsIds = searchRunBestFoundSolutions, 
                         searchRunUpdatedObjectiveValues = searchRunUpdatedValues,
                         searchRunUpdatedTimes = searchRunUpdatedTimes)
  
  return(analysisOutput)
  
}
