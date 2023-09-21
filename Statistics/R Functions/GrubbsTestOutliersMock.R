## Load the R libraries ####
library(biotools)
library(tidyverse)
library(rstatix)
library(emmeans)
library(dplyr)
library(outliers)

##Load Dataset
#Read The CSV file with the Variables 
samples <- read.csv(commandArgs(TRUE)[1]) 
samples <- as.data.frame(samples)

#Retrieve Variables 
GroupBy <- as.character(commandArgs(TRUE)[2])
Variable <- as.character(commandArgs(TRUE)[3])

#Using rename() change column name 
samples <- samples %>% 
  rename("Groups" = GroupBy)

#Find Number of Tones 
tones <- c(startsWith(colnames(samples), Variable))
num_tones <- sum(tones, na.rm = TRUE) + 2

samples <- samples %>%
  gather(key = "Tones", value = "Score", as.character(colnames(samples)[3: num_tones])) %>%
  convert_as_factor(Tones)


## Check For Outliars 

#Function to recursively find Outliars with Grubbs test, until p.val < 0.05
grubbs.flag <- function(x) {
  outliers <- NULL
  test <- x
  grubbs.result <- grubbs.test(test)
  pv <- grubbs.result$p.value
  pvalues <- c(pv)
  # throw an error if there are too few values for the Grubb's test
  if (length(test) < 3 ) stop("Grubb's test requires > 2 input values")
  while(pv < 0.05) {
    outliers <- c(outliers,as.numeric(strsplit(grubbs.result$alternative," ")[[1]][3]))
    test <- x[!x %in% outliers]
    # stop if all but two values are flagged as outliers
    if (length(test) < 3 ) {
      warning("All but two values flagged as outliers")
      break
    }
    grubbs.result <- grubbs.test(test)
    pv <- grubbs.result$p.value
    pvalues <- c(pvalues, pv)
  }
  out <- list()
  out$OutliarsTable <- data.frame(X=x,Outlier=(x %in% outliers))
  pvalues <- head(pvalues, - 1) 
  out$Out_pVal <- pvalues
  return(out)
}
#Create Outliars table TRUE/FALSE
ResultsOut <- grubbs.flag(samples[,4])
outliars_table <- ResultsOut$OutliarsTable
pvalues <- ResultsOut$Out_pVal
#Retrieve Only Outliars
outliars <- c()
index <- c()
for (x in row_number(outliars_table)){
  if (outliars_table$Outlier[x] == TRUE) {
    outliars <- append(outliars, outliars_table$X[x])
    index <- append(index, x)
  }
}

# # Remove Outliars From Sample Dataset 
# OutliarsRow <- list()
# if (is.null(outliars)== FALSE) {
#   to_remove <- c()
#   for (x in 1:length(outliars)){
#     to_remove <- append(to_remove, as.integer(row.names(samples[samples$Score == outliars[x], ])))
#     OutliarsRow <- append(OutliarsRow, list(samples[samples$Score == outliars[x], ]))
#   }
# }
# OutliarsRow <- reduce(OutliarsRow, full_join)
# animals <- c(as.integer(OutliarsRow$Animal))
# for (x in 1:length(animals)) {
#   samples <- samples[samples$Animal != animals[x], ] 
# }

# Create Table to Put in The CSV 
if (length(outliars) != 0) {
  outliars <- data.frame(Index = index, Outliers = outliars, P_Values = pvalues)
} else {
  outliars$Result <- "There are 0 outliars"
}

write.csv(outliars, "OutliarsTable.csv", row.names = FALSE)
