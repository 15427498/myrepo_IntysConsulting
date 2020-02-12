#' Optimisation of pheno OBS.
#'
#' Run a JAMES subset selection algorithm to attemptively find an optimal solution to a pheno OBS problem.
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
#' @param inbreedingCoefficientMetric (default, boolean) Set to TRUE if the inbreeding coefficient ought to be used to measure the diversity in the population, and FALSE if the more basic genetic similarity metric ought to be used.
#' The default value is TRUE.
#' @param Gmatrix (default) A symmetrical matrix containing values on the range [0,2] depicting the relationship between the parents.
#' This matrix is for use only in the inbreeding coefficient metric and thus only needs to be provided if the argument inbreedingCoefficientMetric is TRUE.
#' @param inbreedingCoefficientBoundsEstimatorSampleSize (default, integer) A value specifying the number of samples employed in estimating the lower and upper bound values of the inbreeding coefficient metric.
#' Thse bounds are used toward the normalisation process of the diversity sub-objective function.
#' The higher the complexity of the problem, the higher this value ought to be.
#' The default value is set to a conservative 10000.
#' @param inclusionSet (default, integer array) A vector containing the integer ids of offsping crosses that the user desires to include a priori in the returned solution.
#' The default value is a null vector, so that no id inclusion constraint is imposed on the selection of crosses.
#' @param exclusionSet (default, integer array) A vector containing the integer ids of offsping crosses that the user desires to exclude a priori from the returned solution.
#' The default value is a null vector, so that no id exclusion constraint is imposed on the selection of crosses.
#' @param solutionSize (required, integer) The number of crosses to be extracted.
#' @param maximumSearchTime (required, integer) An upper bound on the total duration of the optimisation process.
#' This parameter is configured in milliseconds.
#' @param timeWithoutImprovement (default, integer) An upper bound on the total duration of the optimisation process.
#' This parameter is configured in milliseconds.
#' @param maximumNumberMoves (default, integer) An upper bound on the total number of search space solutions visited and assessed during the course of the optimisation process.
#' @param searchMethod (default) Summon either one of these three proposed algorithms to run: 
#' 1. Random descent (enter "RD"); 
#' 2. Metropolis search (enter "MS") or 
#' 3. Parallel  tempering (enter "PT").
#' @param searchParametersMS (default, double) A [1*1] data vector specifying the initial temperature of the Metropolis search algorithm.
#' @param searchParametersPT (default) A [1*3] data vector specifying the search parameters of the parallel tempering search algorithm.
#' In the first, second and third slot (respectively):
#' 1. processors (integer) : The number of processors running in parallel.
#' 2. minTemp (double) : Minimum (i.e. final) temperature.
#' 3. maxTemp (double) : Maximum (i.e. initial) temperature.
#'
#' @return A list containing the following five objects (respectively): 
#' 1. The best found ids in R format,
#' 2. The genetic diversity score of the best found solution,
#' 3. The size of the extracted solution (for use in the phenoOBSOutput function),
#' 4. The objective function value of the best found solution, and 
#' 5. The convergence time (in milliseconds) of the algorithm (i.e. how long it took to reach the vicinity of the best found solution).
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
#' inbreedingCoefficientMetric <- TRUE
#' Gmatrix <- matrix(nrow = 4, ncol = 4)
#' Gmatrix[1,] = c(1,0.60,1.40,0.30)
#' Gmatrix[2,] = c(0.60,1,0.50,0.35)
#' Gmatrix[3,] = c(1.40,0.50,1,1.55)
#' Gmatrix[4,] = c(0.30,0.35,1.55,1)
#' inbreedingCoefficientBoundsEstimatorSampleSize <- 300
#' inclusionSet <- c(4)
#' exclusionSet <- c(6)
#' solutionSize <- 2
#' maximumSearchTime <- 500
#' timeWithoutImprovement <- 200
#' maximumNumberMoves <- 400
#' searchMethod <- c("MS")
#' searchParametersMS <- c(5)
#' searchParametersPT <- c(6,0.001,9)
#' 
#' @author A.Colmant, K.Baert and G.De Meyer
#' @export
#' 
phenoObsOptimisation <- function(parentsPhenotypeDataFrame,
                                 offspringPhenotype,
                                 parentsGeneticSimilarityArray,
                                 geneticSimilarityWeight,
                                 penaltyReplicatedIndividuals,
                                 inbreedingCoefficientMetric = TRUE,
                                 Gmatrix = matrix(nrow=1,ncol=1),
                                 inbreedingCoefficientBoundsEstimatorSampleSize = 10000,
                                 inclusionSet = NULL,
                                 exclusionSet = NULL,
                                 solutionSize,
                                 maximumSearchTime,
                                 timeWithoutImprovement = 10000,
                                 maximumNumberMoves = 90000000,
                                 searchMethod = "MS",
                                 searchParametersMS = c(0.00005),
                                 searchParametersPT = c(6,0.0000001,0.001)){
  
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
  if(any(Gmatrix < 0) || any(Gmatrix > 2))
    stop("The entries in Gmatrix should all be between 0 and 2")
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
  
  optimisationOutput <- list(bestFoundOffspringIds = bestFoundSolutionIdsR, 
                             bestFoundSolutionGeneticDiversity = bestFoundSolutionGeneticDiversity, 
                             numberOfCrosses = solutionSize,
                             bestFoundSolutionObjectiveValue = solutionEvaluation,
                             algorithmConvergenceTime = convergenceTime)
  
  return(optimisationOutput)
  
}
