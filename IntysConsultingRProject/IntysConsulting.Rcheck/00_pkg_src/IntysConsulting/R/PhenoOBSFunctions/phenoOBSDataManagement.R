#' Data management of pheno OBS.
#'
#' Manage, convert and simplify the pheno OBS user input data to be presented as Java input.
#' 
#' @param parentsPhenotypeDataFrame (required) A data frame containing the parentNames in the first column, the trait names as headings in the remaining columns, and the parent phenotype information in the remaining matrix entries.
#' See example below for more clarity.
#' @param parentsGeneticSimilarityMatrix (required, double) An [N*N] data frame, also known as relationship (triangular) matrix, containing the degree of genetic similarity accross all pairs of parents (i.e. accross all potential crosses).
#' The entries ought to be mapped on the range [0,1].
#' @param constraintParameters (required, double matrix) A [T*2] data frame delimiting the feasible domain space with 
#' the use of threshold operational constraints. In the first and second column (respectively):
#' * lbvalue : Lower numerical threshold.
#' * ubvalue : Upper numerical threshold.
#' These numerical thresholds impose tolerance bounds on the expected trait values that the resulting offspring corsses are allowed to undertake.
#' Leave blank fields for unconstrained traits.
#' Imposing both upper and lower thresholds on a specific trait is permitted.
#' Check the constriction level of your resulting search space with the complementary function "phenoOBSCheckSearchSpace".
#' The default setting is NULL.
#' @param traitObjectiveWeights (required, double matrix) A [T*2] data frame containing the trait sub-objective weights.
#' In the first and second column (respectively):
#' * weight (double) : Sub-objective weight of the criterion.
#' * is.maximisation? (boolean) : Enter 1 if the sub-objective is to be maximised, and 0 otherwise.
#' Leave blank fields (or zeros) for traits not considered in the configuration of the objective space.
#' Although the weights need not necessarily add up to one, they need to be proportional to one another as a direct reflection of their relative importance.
#' @details Symbols notation:
#' * T : The total number of input traits,
#' * N : The total number of parents,
#' * n : The size of the solution.
#' 
#' @return A list containing the following five objects (respectively): 
#' 1. A matrix containing:
#'     * The feasible offspring Ids in the first column,
#'     * The corresponding offspring parent ids in the second and third columns, and
#'     * The weighted, normalised trait objective evaluation of the offspring.
#' 2. An array containing the normalised upper triangular entries of the parents genetic similarity matrix.
#' 3. The parents phenotype input data frame (same as described earlier) for use in the other pheno OBS functions. 
#' 4. A matrix containing the overall offspring information (for use in the output function).
#' 5. The parents genetic similarity matrix described above (for use in the output function).
#'
#' @examples 
#' parentsPhenotypeDataFrame <- data.frame(
#' parentNames = c("A","B","C","D"),
#' Yield = c(1200,1300,1250,1150),
#' Lint = c(50,40,45,55),
#' BollType = c(4.5,4,3,4.5))
#' parentsGeneticSimilarityMatrix <- matrix(nrow = 4, ncol = 4)
#' parentsGeneticSimilarityMatrix[1,] = c(1,0.60,0.40,0.30)
#' parentsGeneticSimilarityMatrix[2,] = c(0.60,1,0.50,0.35)
#' parentsGeneticSimilarityMatrix[3,] = c(0.40,0.50,1,0.55)
#' parentsGeneticSimilarityMatrix[4,] = c(0.30,0.35,0.55,1)
#' constraintParameters <- matrix(nrow = 3, ncol = 2)
#' constraintParameters[1,1] <- 1190
#' traitObjectiveWeights <- matrix(nrow = 3, ncol = 2)
#' traitObjectiveWeights[1,] = c(0.7,1)
#' traitObjectiveWeights[2,] = c(0.3,1)
#' traitObjectiveWeights[3,] = c(0,1)
#' 
#' @author A.Colmant, K.Baert and G.De Meyer
#' @export
#' 
phenoObsDataManagement <- function(parentsPhenotypeDataFrame,
                                   parentsGeneticSimilarityMatrix,
                                   constraintParameters = NULL,
                                   traitObjectiveWeights){
  
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
  
  #Fill in missing values 
  if(all(is.na(constraintParameters))){
    constraintParameters <- matrix(nrow=numberTraits,ncol=2)
  }
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
  
  dataManagementOutput <- list(condensedOffspringInformation = offspringInformationSimplified, 
                               parentsGeneticSimilarityArray = parentsGeneticSimilarityArray, 
                               parentsPhenotypeDataFrame = parentsPhenotypeDataFrame, 
                               completeOffspringInformation = overallOffspringInformation, 
                               parentsGeneticSimilarityMatrix = parentsGeneticSimilarityMatrix)
  
  return(dataManagementOutput)
  
}


