
#############################################################
#
#  Code to compute calculations presented in the hons project
#  Revision 11/23
#  Emma Morrison: 2409764@dundee.ac.uk
#  Davide Bulgarelli: d.bulgarelli@dundee.ac.uk 
#############################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################
#required packages
library("ggplot2")
library("PMCMRplus")

#############################################################
#set working directory-Davide CPU
setwd("/cluster/db/R_shared/Emma//")
#set working directory-Emma CPU
#TO BE COMPLETED
#############################################################

##################################################################
#Import the data
#################################################################

#import the dataset
EM_info <-read.delim("EM_results_1123.txt", row.names = 1, header = T)

#inspect the files
EM_info

#check the type of columns: independent variables, e.g., genotypes, should be factor while dependent variables, e.g., biomass, should be numeric 
#Independent variables
class(EM_info$Genotype)
class(EM_info$Operator)
#convert to factor
EM_info$Genotype <- as.factor(EM_info$Genotype)
EM_info$Operator <- as.factor(EM_info$Operator)
class(EM_info$Genotype)
class(EM_info$Operator)
#dependent variables
class(EM_info$Biomass)
class(EM_info$Zadok)
class(EM_info$Qubit)
class(EM_info$Purified.DNA)
                           
###################################################################################################
#Let's re-order factors and colors
###################################################################################################

#Genotype
levels(EM_info$Genotype)
EM_info$Genotype <- factor(EM_info$Genotype, levels=c("Bulk", "SL17", "Barke", "F3.31.65", "F3.41.78.3", "SL52"))

#color coding
DB_cols <- c("#000000", "#56B4E9", "#0072B2", "#009E73", "#F179A7", "#F0E442")

###################################################################################################
#Let's look at plant's characteristics
###################################################################################################

#subset for plant samples using a negative selection '!='
EM_info_plant <-EM_info[which(EM_info$Genotype !='Bulk'), ]
EM_info_plant

###################################################################################################
#Inspect independent variables and data distribution 
###################################################################################################

hist(EM_info_plant$Zadok)
hist(EM_info_plant$Biomass)

#a bit of work on the figures...

hist(EM_info_plant$Zadok, main = paste("Recombinant lines Zadok's distributin"), xlab = paste("Stages"), ylab = paste ("number of samples"), breaks = 1, ylim = c(0,36))
hist(EM_info_plant$Biomass, main = paste("Recombinant lines Biomass' distribution"), xlab = paste("Aboveground biomass"), ylab = paste ("number of samples"), ylim = c(0,15))

###################################################################################################
#Biomass data visualization using boxplot 
###################################################################################################
#useful website http://www.sthda.com/english/wiki/ggplot2-box-plot-quick-start-guide-r-software-and-data-visualization 

#Genotype
p <- ggplot(EM_info_plant, aes(x=Genotype, y=Biomass)) + 
  geom_boxplot()
p
#something more elaborate...
p <- ggplot(EM_info_plant, aes(x=Genotype, y=Biomass, color=Genotype)) + 
  geom_boxplot()
p
#adding individual replicates
p + geom_jitter(shape=16, position=position_jitter(0.2))

dev.off()
#four colors
DB_cols_plants <- c("#56B4E9", "#0072B2", "#009E73", "#F179A7", "#F0E442")
DB_cols_plants_2 <- c("#56B4E9", "#0072B2", "#009E73", "#F179A7", "#F0E442")
p <- ggplot(EM_info_plant, aes(x=Genotype, y=Biomass), color=Genotype) + 
  geom_boxplot() + scale_color_manual(values = DB_cols_plants)
p
#adding individual replicates
p + geom_jitter(shape=16, position=position_jitter(0.2))
p

#formal assessment of data distribution to identify the appropriate test 
#Shapiro test: https://en.wikipedia.org/wiki/Shapiro%E2%80%93Wilk_test
shapiro.test(EM_info_plant$Biomass)
# p value below .05 denotes data distribution deviating from a normal one, use non- parametric tests

#statistical test
#general formula before the sign goes the name fo the column of dependent variable, after the name of the column of independent you want to test
#Microhabitat effect
kruskal.test(Biomass ~ Genotype, data = EM_info_plant)
#the output tells us that at least one group of samples is different from another one, but which one?

#To answer this question we require what is called post-hoc test
kwAllPairsDunnTest(EM_info_plant$Biomass, EM_info_plant$Genotype, p.adjust.method="BH")

#Zadok (no graphs)

#formal assessment of data distribution to identify the appropriate test 
#Shapiro test: https://en.wikipedia.org/wiki/Shapiro%E2%80%93Wilk_test
shapiro.test(EM_info_plant$Zadok)

#Microhabitat effect
kruskal.test(Zadok ~ Genotype, data = EM_info_plant)
#no statistical differences

###################################################################################################
#DNA extraction (note that we include bulk soil here)
###################################################################################################

hist(EM_info$Qubit)
hist(EM_info$Purified.DNA)

#a bit of work on the figures...
hist(EM_info_plant$Qubit, main = paste("Extracted DNA distribution"), xlab = paste("Concentration"), ylab = paste ("number of samples"), breaks = 5, ylim = c(0,15))
hist(EM_info_plant$Purified.DNA, main = paste("Purified DNA distribution"), xlab = paste("Concentration"), ylab = paste ("number of samples"), breaks = 10, xlim = c(0,200), ylim = c(0,10))

#here we could ask the question whether input DNA may predict the efficiency of amplification (i.e., is the extracted DNA comparable with the amplified one?)
ggplot(EM_info, aes(x=Qubit, y=Purified.DNA)) + 
  geom_point(color='#2980B9', size = 4) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color='#2C3E50')
#it seems clear that for two samples we do have sufficient DNA but either the amplification or the purification did not work well..

#There is a statistical parameter called Spearman rank's correlation that will inform us as to whether the two dataset are signficantly linked
corr <- cor.test(x=EM_info$Qubit, y=EM_info$Purified.DNA, method = 'spearman')
#it is a rank correlation and similar to the KW ties may be problematic
corr
#The p value is above 0.05, meaning that the two dataset aren't significantly linked. This means that either (or both) efficiency of amplification and purification impact on the results of your work, regardless the input DNA concentration.

#END

