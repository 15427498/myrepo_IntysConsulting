#' #' Objective function of pheno OBS.
#' #'
#' #' Evaluate the appreciation of a solution to a given pheno OBS problem using the proposed model objective function.
#' #' 
#' #' @param solutionOffspringIdsR (required, integer) A [1*n] vector containing a set of offspring ids expressed in R format.
#' #' @param numberParents (required, integer) The number of candidate parents used in this selection problem.
#' #' @param offspringPhenotype (required) A matrix containing:
#' #' * The feasible offspring Ids in the first column,
#' #' * The corresponding offspring parents in the second and third columns, and
#' #' * The weighted, normalised trait objective evaluation of the offspring.
#' #' @param parentsGeneticSimilarityArray (required) An array containing the normalised upper triangular entries of the parents genetic similarity matrix.
#' #' @param geneticSimilarityWeight (required) A value in the range [0,1] for the importance of solution diversity relative to the traits.
#' #' For example, if specified to be 0.3, then the algorithm will automatically allocate 0.7 toward finding solutions that are strong across the entirety of the traits sub-objective.
#' #' Note, however, that the individual trait importances are in phenoOBSDataManagement.
#' #' @param penaltyReplicatedIndividuals (default, double) A soft-constraint value penalising the presence of replicated individuals in a solution.
#' #' Note that this value is not normalised alongside the set of upper triangular matrix entries, and so it must have at least a value of one.
#' #' @details symbols notation:
#' #' * n : The size of the solution.
#' #' 
#' #' @return 
#' #' * 1: The model objective function evaluation of the given solution.
#' #' * 2: The model trait sub-objective evaluation of the given solution (non weighted).
#' #' * 3: The model genetic similarity sub-objective evaluation of the given solution (non weighted).
#' #' 
#' #' @examples 
#' #' solutionOffspringIdsR <- c(1,4)
#' #' numberParents <- 4
#' #' offspringPhenotype <- matrix(nrow = 5, ncol = 4)
#' #' offspringPhenotype[1,] <- c(1,1,2,0.6)
#' #' offspringPhenotype[2,] <- c(2,1,3,0.5)
#' #' offspringPhenotype[3,] <- c(4,2,3,0.7)
#' #' offspringPhenotype[4,] <- c(5,2,4,0.5)
#' #' offspringPhenotype[5,] <- c(6,3,4,0.4)
#' #' parentsGeneticSimilarityArray <- c(1,0.33333333,0,0.66666667,0.166666667,0.8333333333)
#' #' geneticSimilarityWeight <- 0.4
#' #' penaltyReplicatedIndividuals <- 1.08
#' #' 
#' #' @author A.Colmant and F.Pita
#' #' @export
#' #' 
#' phenoOBSObjectiveFunction <- function(solutionOffspringIdsR,
#'                                  numberParents,
#'                                  offspringPhenotype,
#'                                  parentsGeneticSimilarityArray,
#'                                  geneticSimilarityWeight,
#'                                  penaltyReplicatedIndividuals){
#'   
#'   #Some pre-calculations
#'   solutionSize <- length(solutionOffspringIdsR)
#'   traitsWeight <- 1 - geneticSimilarityWeight 
#'   
#'   #Convert the array of offspring ids to one of parent ids
#'   solutionParentIdsR <- c(2*solutionSize)
#'   for(i in 1:solutionSize){
#'     offspringId <- solutionOffspringIdsR[i]
#'     feasibleOffspringId <- which(offspringPhenotype[,1] == offspringId)
#'     solutionParentIdsR[2*i-1] = offspringPhenotype[feasibleOffspringId,2]
#'     solutionParentIdsR[2*i] = offspringPhenotype[feasibleOffspringId,3]
#'   }
#'   
#'   #Calculate the solution traits (collective) sub-objective value 
#'   traitsObjectiveValue <- 0
#'   for(i in 1:solutionSize){
#'     offspringId <- solutionOffspringIdsR[i]
#'     feasibleOffspringId <- which(offspringPhenotype[,1] == offspringId)
#'     traitsObjectiveValue <- traitsObjectiveValue + offspringPhenotype[feasibleOffspringId,4]
#'   }
#'   
#'   #Calculate the solution genetic similarity sub-objective value 
#'   # -Part1: Assess the genetic similarity accross the parent individuals (i.e. short-term genetic similarity benefits)
#'   geneticSimilarityObjectiveValue1 <- 0
#'   for(i in 1:solutionSize){
#'     offspringId <- solutionOffspringIdsR[i]
#'     geneticSimilarityObjectiveValue1 <- geneticSimilarityObjectiveValue1 + parentsGeneticSimilarityArray[offspringId]
#'   }
#'   
#'   # -Part2: Assess the genetic similarity accross the offspring individuals (i.e. long-term genetic similarity benefits)
#'   geneticSimilarityObjectiveValue2 <- 0
#'   firstPair <- 1
#'   penalty <- penaltyReplicatedIndividuals
#'   while(firstPair < solutionSize){
#'     secondPair <- firstPair + 1
#'     while(secondPair < solutionSize + 1){
#'       parentIds <- c(solutionParentIdsR[2*firstPair-1],solutionParentIdsR[2*secondPair-1])
#'       parent1IdAdj <- min(parentIds) - 1
#'       parent2IdAdj <- max(parentIds) - 1
#'       if(parent1IdAdj != parent2IdAdj){
#'         offspringId <- parent1IdAdj * (numberParents-1) + parent2IdAdj - ((parent1IdAdj * (parent1IdAdj + 1)) / 2)
#'         geneticSimilarityObjectiveValue2 <- geneticSimilarityObjectiveValue2 + parentsGeneticSimilarityArray[offspringId]
#'       } else{
#'         geneticSimilarityObjectiveValue2 <- geneticSimilarityObjectiveValue2 + penalty
#'       }
#'       parentIds <- c(solutionParentIdsR[2*firstPair-1],solutionParentIdsR[2*secondPair])
#'       parent1IdAdj <- min(parentIds) - 1
#'       parent2IdAdj <- max(parentIds) - 1
#'       if(parent1IdAdj != parent2IdAdj){
#'         offspringId <- parent1IdAdj * (numberParents-1) + parent2IdAdj - ((parent1IdAdj * (parent1IdAdj + 1)) / 2)
#'         geneticSimilarityObjectiveValue2 <- geneticSimilarityObjectiveValue2 + parentsGeneticSimilarityArray[offspringId]
#'       } else{
#'         geneticSimilarityObjectiveValue2 <- geneticSimilarityObjectiveValue2 + penalty
#'       }
#'       parentIds <- c(solutionParentIdsR[2*firstPair],solutionParentIdsR[2*secondPair-1])
#'       parent1IdAdj <- min(parentIds) - 1
#'       parent2IdAdj <- max(parentIds) - 1
#'       if(parent1IdAdj != parent2IdAdj){
#'         offspringId <- parent1IdAdj * (numberParents-1) + parent2IdAdj - ((parent1IdAdj * (parent1IdAdj + 1)) / 2)
#'         geneticSimilarityObjectiveValue2 <- geneticSimilarityObjectiveValue2 + parentsGeneticSimilarityArray[offspringId]
#'       } else{
#'         geneticSimilarityObjectiveValue2 <- geneticSimilarityObjectiveValue2 + penalty
#'       }
#'       parentIds <- c(solutionParentIdsR[2*firstPair],solutionParentIdsR[2*secondPair])
#'       parent1IdAdj <- min(parentIds) - 1
#'       parent2IdAdj <- max(parentIds) - 1
#'       if(parent1IdAdj != parent2IdAdj){
#'         offspringId <- parent1IdAdj * (numberParents-1) + parent2IdAdj - ((parent1IdAdj * (parent1IdAdj + 1)) / 2)
#'         geneticSimilarityObjectiveValue2 <- geneticSimilarityObjectiveValue2 + parentsGeneticSimilarityArray[offspringId]
#'       } else{
#'         geneticSimilarityObjectiveValue2 <- geneticSimilarityObjectiveValue2 + penalty
#'       }
#'       secondPair <- secondPair + 1
#'     }
#'     firstPair <- firstPair + 1
#'   }
#'   geneticSimilarityObjectiveValue2 <- geneticSimilarityObjectiveValue2 / (2*(solutionSize - 1))
#'   
#'   #Combine the two sub-objective scores then normalise and weigh them 
#'   traitsObjectiveValue <- traitsObjectiveValue/solutionSize
#'   geneticSimilarityObjectiveValue <- (geneticSimilarityObjectiveValue1 + geneticSimilarityObjectiveValue2) / (2*solutionSize)
#'   objectiveValue <- (traitsWeight * traitsObjectiveValue) - (geneticSimilarityWeight * geneticSimilarityObjectiveValue) 
#'        
#'   output <- list(objectiveValue = objectiveValue, traitsSubObjectiveValue = traitsObjectiveValue, geneticSimilaritySubObjectiveValue = geneticSimilarityObjectiveValue)
#'   return(output)
#'                               
#' }             


#TEST
# library(BCS.JAMES3)
# library(MASS)
# library(BCS.Base)
# library(DiGGer)
# sessionInfo()
# solutionOffspringIdsR <- c(1,5)
# numberParents <- 4
# offspringPhenotype <- matrix(nrow = 5, ncol = 4)
# offspringPhenotype[1,] <- c(1,1,2,0.6)
# offspringPhenotype[2,] <- c(2,1,3,0.5)
# offspringPhenotype[3,] <- c(4,2,3,0.7)
# offspringPhenotype[4,] <- c(5,2,4,0.5)
# offspringPhenotype[5,] <- c(6,3,4,0.4)
# parentsGeneticSimilarityArray <- c(1,0.33333333,0,0.66666667,0.166666667,0.8333333333)
# geneticSimilarityWeight <- 0.4
# penaltyReplicatedIndividuals <- 1.08
# phenoOBSObjectiveFunction(solutionOffspringIdsR,
#                                       numberParents,
#                                       offspringPhenotype,
#                                       parentsGeneticSimilarityArray,
#                                       geneticSimilarityWeight,
#                                       penaltyReplicatedIndividuals)
                                   