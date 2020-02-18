#' Output information of pheno OBS.
#'
#' Present the output results obtained from a pheno OBS search.
#'  
#' @param parentsPhenotypeDataFrame (required) A data frame containing the parentNames in the first column, the trait names as headings in the remianing columns, and the parent phenotype information in the remaining matrix entries.
#' See example below for more clarity.
#' This data frame may also be obtained from the phenoOBSDataManagement output.
#' @param optimalSolutionIdsR (required, integer) A [1*n] vector containing the best found Ids, expressed in R format,
#  and obtained from the phenoOBSOptimisation output.
#' @param bestFoundSolutionGeneticDiversity (required, double) The genetic diversity of the best found solution.
#' @param overallOffspringInformation (required) A matrix containing information on all offspring combinations.
#' @param solutionSize (required, integer) The number of crosses to be extracted.
#' @param parentsGeneticSimilarityMatrix (required, double) An [N*N] data frame, also known as relationship (triangular) matrix, containing the degree of genetic similarity accross all pairs of parents (i.e. accross all potential crosses).
#' The entries ought to be mapped on the range [0,1].
#' @details symbols notation:
#' * n : The size of the solution
#' 
#' @return A list containing the following two objects (respectively): 
#' 1. A [1*(T+1)] doubles vector containing the solution sub-objective scores (respectively): 
#'     * Trait scores, in the same order as the parentPhenotype data in vector entries 1 to T, and
#'     * Solution genetic similarity score in vector entry T+1.
#' 2. A [n*(5+T+1)] string matrix containing:
#'     * Best found offspring IDs in the first column,
#'     * Corresponding parent IDs in the second and third columns,
#'     * Corresponding parent IDs in the fourth and fifth columns,
#'     * Corresponding trait scores in the next T columns, and
#'     * Corresponding relationship matrix entry in the last column.
#' 
#' @examples 
#' parentsPhenotypeDataFrame <- data.frame(
#' parentNames = c("A","B","C","D"),
#' Yield = c(1200,1300,1250,1150),
#' Lint = c(50,40,45,55),
#' BollType = c(4.5,4,3,4.5))
#' optimalSolutionIdsR <- c(1,4)
#' bestFoundSolutionGeneticDiversity <- 0.80167
#' overallOffspringInformation <- matrix(nrow = 6, ncol = 6)
#' overallOffspringInformation [1,] <- c(1,1,2,1250,45,4.25)
#' overallOffspringInformation [2,] <- c(2,1,3,1225,47.5,3.75)
#' overallOffspringInformation [3,] <- c(3,1,4,1175,52.5,4.5)
#' overallOffspringInformation [3,] <- c(4,2,3,1275,42.5,3.5)
#' overallOffspringInformation [4,] <- c(5,2,4,1225,47.5,4.25)
#' overallOffspringInformation [5,] <- c(6,3,4,1200,50,3.75)
#' solutionSize <- 2
#' parentsGeneticSimilarityMatrix <- matrix(nrow = length(parentsPhenotypeDataFrame$parentNames), ncol = length(parentsPhenotypeDataFrame$parentNames))
#' parentsGeneticSimilarityMatrix[1,] = c(1,0.60,0.40,0.30)
#' parentsGeneticSimilarityMatrix[2,] = c(0.60,1,0.50,0.35)
#' parentsGeneticSimilarityMatrix[3,] = c(0.40,0.50,1,0.55)
#' parentsGeneticSimilarityMatrix[4,] = c(0.30,0.35,0.55,1)
#' 
#' @author A.Colmant, K.Baert and G.De Meyer
#' @export
#' 
phenoObsOutput <- function(parentsPhenotypeDataFrame,
                           optimalSolutionIdsR,
                           bestFoundSolutionGeneticDiversity,
                           overallOffspringInformation,
                           solutionSize,
                           parentsGeneticSimilarityMatrix){
  
  #Construct a data object containing information on the sub-objective score attainments of the best found solution
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

  phenoObsOutput <- list(bestFoundSolutionSubObjectiveScores = DFbestFoundSolutionObjectiveScores,
                         bestFoundSolutionCrossesInformation = DFbestFoundSolutionIdsInformation)
  
  return(phenoObsOutput)
  
}




