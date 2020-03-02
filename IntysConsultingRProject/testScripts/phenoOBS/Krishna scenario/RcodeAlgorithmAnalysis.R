
library(MASS)
library(xlsx)
library(ggplot2)
library(gifski)
library(av)
library(gganimate)
library(chron)
library(transformr)
sessionInfo()

#Note: We use the same scenario as in "RcodeAlexCheck"
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
solutionSize <- 20
maximumSearchTime <- 20000
timeWithoutImprovement <- 3000
maximumNumberMoves <- 80000000
analysisParametersMS = rep(NA, 1)
analysisParametersPT = matrix(NA, nrow=1, ncol=3)

###############################################################################################################

#--------------------------------------------------------------------------------------------
# (3) Conduct analysis for the RANDOM DESCENT algorithm
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

mergedObjectRD <- cbind(bestFoundUpdatedTimesObjectiveValues4[,1,],bestFoundUpdatedObjectiveValues4[,1,])
mergedObjectRD[1,1] <- 0
colnames(mergedObjectRD) <- c("Time.milliseconds", "Objective.value")
plot(mergedObjectRD)
write.xlsx(mergedObjectRD,"C:/Users/Gebruiker/Documents/Intys/PowerBI/RStudioExportedData/RDAlgorithmProgression.xlsx", row.names = FALSE)
mergedObjectRD <- data.frame("Time.milliseconds" = mergedObjectRD[,1], "Objective.value" = mergedObjectRD[,2])

graphStatic <- ggplot(mergedObjectRD,aes(x=Time.milliseconds, y=Objective.value)) +
  geom_point() +
  ylab("Objective value") + xlab("Time (milliseconds)") +
  #scale_x_continuous(breaks=seq(0, 6, 1)) +
  #scale_y_continuous(breaks=seq(0, 45, 10)) +
  theme(legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.ticks.length = unit(5, "pt"),
        axis.title.y = element_text(size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 11))
graphStatic

lineGraphDynamic <- graphStatic + geom_point(aes(group = seq_along(Time.milliseconds))) + 
  transition_states(Time.milliseconds,
                    #transition_length = 2,
                    #state_length = 2,
                    wrap = TRUE) + 
  transition_reveal(Time.milliseconds)
animate(lineGraphDynamic, fps=20)





###############################################################################################################

#--------------------------------------------------------------------------------------------
# (4) Conduct analysis for the METROPOLIS SEARCH algorithm
#--------------------------------------------------------------------------------------------

#Experiment 1

numberSearches <- 5
numberRuns <- 1
searchIds <- c("Search 1","Search 2","Search 3","Search 4","Search 5")
searchMethod <- "MS"
analysisParametersMS <- c(20,5,0.1,0.0001,0.000001) 

phenoObsAnalysisOutput1 <- phenoObsAlgorithmPerformance(parentsPhenotypeDataFrame,
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

bestFoundObjectiveValues1 <- phenoObsAnalysisOutput1[[1]]
bestFoundIds1 <- phenoObsAnalysisOutput1[[2]]
bestFoundUpdatedObjectiveValues1 <- phenoObsAnalysisOutput1[[3]]
bestFoundUpdatedTimesObjectiveValues1 <- phenoObsAnalysisOutput1[[4]]
print(bestFoundObjectiveValues1)
print(bestFoundIds1)
print(bestFoundUpdatedObjectiveValues1)
print(bestFoundUpdatedTimesObjectiveValues1)

numberObservations <- length(bestFoundUpdatedObjectiveValues1)/numberSearches
matrix <- cbind(bestFoundUpdatedTimesObjectiveValues1[s,1,],bestFoundUpdatedObjectiveValues1[s,1,])
mergedObjectsMS <- data.frame("Temperature.initial" = analysisParametersMS[1], "Time.milliseconds" = matrix[,1], "Objective.value" = matrix[,2])
for(s in 2:numberSearches){
  matrix <- cbind(bestFoundUpdatedTimesObjectiveValues1[s,1,],bestFoundUpdatedObjectiveValues1[s,1,])
  mergedObjectsMS <- rbind(mergedObjectsMS,data.frame("Temperature.initial" = analysisParametersMS[s], "Time.milliseconds" = matrix[,1], "Objective.value" = matrix[,2]))
}

graphStatic <- ggplot(mergedObjectsMS,aes(x=Time.milliseconds, 
                                          y=Objective.value, 
                                          #group=Temperature.initial, 
                                          colour=as.factor(Temperature.initial)),
                                          #size = pop,  frame = year)
                                          ) +
  geom_line() +
  #geom_point(mergedObjectsMS[[2]], aes(x=Time.milliseconds, y=Objective.value), color = "darkred") +
  ylab("Objective value") + xlab("Time (milliseconds)") +
  labs(fill='NEW LEGEND TITLE') +
  #scale_x_continuous(breaks=seq(0, 6, 1)) +
  #scale_y_continuous(breaks=seq(0, 45, 10)) +
  scale_colour_manual(values=c("red","green","orange","blue","purple")) +
  theme(legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.ticks.length = unit(5, "pt"),
        axis.title.y = element_text(size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 11))
graphStatic

lineGraphDynamic <- graphStatic + geom_point(aes(group = seq_along(Time.milliseconds))) + 
  transition_states(Time.milliseconds,
                    #transition_length = 2,
                    #state_length = 2,
                    wrap = TRUE) + 
  transition_reveal(Time.milliseconds)
animate(lineGraphDynamic, fps=20,
        height = 800, width =800)
size(grap)





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

