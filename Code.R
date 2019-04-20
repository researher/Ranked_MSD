
### Importing library 

library(caret) # for various machine learning functions
library(dplyr) # for data manipulation
library(class) # # for various classification functions
library(randomForest) # for random forest algorithm
library(pROC) # for Plotting the ROC curves
library(ggplot2) # for impressive plotting
library(tidyverse) # for data manipulation

# load data
Final_Data_133_sub_880_samples_Labeled <- read.csv("Final_Data_133_sub_880_samples_Labeled.csv", header=TRUE, sep=",", check.names = FALSE)

