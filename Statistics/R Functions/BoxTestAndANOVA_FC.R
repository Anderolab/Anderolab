## Load the R libraries ####
library(biotools)
library(rstatix)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(emmeans)
library(pwr)
library(CGPfunctions)
library(ggthemes)
library(corrplot)

#For Trials 
#samples <- read.csv("FreezData.csv")
#GroupBy <- as.character("Sexes")
#Variable <- as.character("Tone")
#Id <- as.character("Animal")

#Get dataset
samples <- read.csv(commandArgs(TRUE)[1])

#Retrieve Variables 
GroupBy <- as.character(commandArgs(TRUE)[2])
Variable <- as.character(commandArgs(TRUE)[3])
Id <- as.character(commandArgs(TRUE)[4])

#Retrieve New Folder Path 
Folder <- as.character(commandArgs(TRUE)[5])

## Loading Data ####
#To Remove NA: 
samples <- na.omit(samples) 

#Using rename() change column name 
samples <- samples %>% 
  rename("Groups" = GroupBy, "Animals" = Id)
#Find Number of Tones 
tones <- c(startsWith(colnames(samples), Variable))
num_tones <- sum(tones, na.rm = TRUE) + 2

#Directory to save file
dir.create(Folder) #make new folder for Statistical Results 

## Box's M-test ####

# Perform Box's M-test
BoxPlotRes <- boxM(samples[,3:num_tones], samples[, "Groups"])
BoxPlotResTable <- data_frame(BoxPlotRes$statistic,BoxPlotRes$parameter,BoxPlotRes$p.value)

## Fix Dataset for ANOVA ####

samples <- samples %>%
  gather(key = "Tones", value = "Score", as.character(colnames(samples[,3:num_tones]))) %>%
  convert_as_factor(Animals, Groups)

## Check Mean and SD and Plot Data ####

#Check mean and sd 
SummaryStatMeanSD <- get_summary_stats(group_by(samples, Groups, Tones), 
                                       Score, type = "mean_sd")
SummaryStatMeanSD <- as.data.frame(SummaryStatMeanSD)

#Box plot of data 
ggboxplot(
  samples, x = "Tones", y = "Score",
  color = "Groups", palette = "jco"
)
ggsave("BoxPlot.png", plot = last_plot(), path = Folder, width = 400, height = 400, units = "mm")

## Check ANOVA Assumptions ####

# Check Normality 
model  <- lm(Score ~ Groups*Tones, data = samples) #Build the linear model
ggqqplot(residuals(model)) #Create a QQ plot of residuals
ggsave("QQPlot.png", plot = last_plot(), path = Folder, width = 400, height = 400, units = "mm")

# Compute Shapiro-Wilk test of normality
ShapiroTest <- as.data.frame(shapiro_test(residuals(model))) #If p value NOT significant is GOOD
#Check normality assumption by groups. 
#Computing Shapiro-Wilk test for each combinations of factor levels:
ShapiroTestGroups <- as.data.frame(shapiro_test(
  group_by(samples, Groups, Tones), Score)
)
#Plot 
ggqqplot(samples, "Score", ggtheme = theme_bw()) +
  facet_grid(Groups ~ Tones)
ggsave("QQPlotGroups.png", plot = last_plot(), path = Folder, width = 400, height = 400, units = "mm")

#Check for Homogneity of variance assumption with Levene Test
LeveneTest <- as.data.frame(levene_test(samples, Score ~ Groups*Tones)) 
#If not significant GOOD

## Perform Two Way Repeated Measure ANOVA ####
model.aov <- aov(Score ~ 
                   Groups * Tones + 
                   Error(Animals/(Groups * Tones)), data = samples)
summary_aov <- summary(model.aov)
aov_1 <- as.data.frame(lapply(summary_aov$`Error: Animals`, tidy))
aov_2 <- as.data.frame(lapply(summary_aov$`Error: Animals:Tones`, tidy))

## Plot ANOVA ####

## Don't Know if it's the correct one :(
PlotAnova <- (Plot2WayANOVA(formula = Score ~ Tones*Groups, dataframe = samples,
                            confidence = .95, mean.label = TRUE,
                            posthoc.method = "bonf", 
                            show.dots = TRUE,
                            interact.line.size = 1,
                            ggtheme = ggplot2::theme_light(),
                            ci.line.size = 1)
)

ggsave("ANOVAPlot.png", plot = last_plot(), path = Folder, width = 200, height = 200, units = "mm")

## Posthoc Analisis ####
# Bonferroni post hoc analysis 

if (typeof(PlotAnova$PosthocTable) == "list") {
  if (length(PlotAnova$PosthocTable) == 1) {
    PostHoc_1 <- PlotAnova$PosthocTable$Groups
    PostHoc_2 <- "No significant effects"
  } else {
    PostHoc_1 <- PlotAnova$PosthocTable$Tones
    PostHoc_2 <- PlotAnova$PosthocTable$`Tones:Groups`
  }
} else {
  PostHoc_1 <- "No significant effects"
  PostHoc_2 <- "No significant effects"
}

#Tukey post hoc analysis 

Tukey <- TukeyHSD(aov(Score ~ Groups*Tones,
             data = samples))
TukeyGroups <- Tukey$Groups
TukeyTones <- Tukey$Tones
TukeyGroups_Tones <- Tukey$`Groups:Tones`

# Make Final CSV ####

# Start a sink file with a CSV extension

sink_path <- paste(Folder,"/SummaryStatistics.csv", sep = "")

sink(sink_path)
cat('BoxPlotResTable')
write.csv(BoxPlotResTable)
cat('\n')
cat('\n')
cat('\n')
cat('SummaryStatMeanSD')
write.csv(SummaryStatMeanSD)
cat('\n')
cat('\n')
cat('\n')
cat('ShapiroTest')
write.csv(ShapiroTest)
cat('\n')
cat('\n')
cat('\n')
cat('ShapiroTestGroups')
write.csv(ShapiroTestGroups)
cat('\n')
cat('\n')
cat('\n')
cat('LeveneTest')
write.csv(LeveneTest)
cat('\n')
cat('\n')
cat('\n')
cat('Error_Animals')
write.csv(aov_1)
cat('\n')
cat('\n')
cat('\n')
cat('Error_Animals_Tones')
write.csv(aov_2)
cat('\n')
cat('\n')
cat('\n')
cat('PosthocTones')
write.csv(PostHoc_1)
cat('\n')
cat('\n')
cat('\n')
cat('PosthocTonesGroups')
write.csv(PostHoc_2)
cat('\n')
cat('\n')
cat('\n')
cat('TukeyGroups')
write.csv(TukeyGroups)
cat('\n')
cat('\n')
cat('\n')
cat('TukeyTones')
write.csv(TukeyTones)
cat('\n')
cat('\n')
cat('\n')
cat('TukeyGroups_Tones')
write.csv(TukeyGroups_Tones)
cat('\n')
cat('\n')

# Close the sink
sink()