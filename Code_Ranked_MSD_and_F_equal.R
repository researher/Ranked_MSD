### Importing library
library(caret) # For various machine learning functions
library(dplyr) # For data manipulation
library(class) # # For various classification functions
library(randomForest) # For random forest algorithm
library(pROC) # For plotting the ROC curves
library(ggplot2) # For impressive plotting
library(tidyverse) # For data manipulation

# Loading data
Total_Data <- read.csv("Total_Data.csv", header=TRUE, sep=",", check.names = FALSE)


### Dividing data into training and hold-out test set
set.seed(1234)
# Dividing data set into train (80%) and test (20%) using createDataPartition function of caret package
index_Train <- createDataPartition(y = Total_Data$Label, p = 0.80, list = FALSE)
g_train_data <- Total_Data[index_Train, ]
g_test_data <- Total_Data[-index_Train, ]

# Display the dimensions of training and hold-out test set(rows columns)
(dim(g_train_data))
(dim(g_test_data))

## Using the training data for DEG Identification
### Obtaining healthy reference

# Converting data frame into tibble
Training_Data_tb <- as_tibble(g_train_data)

# Select all samples before inoculation
Healthy_Reference <- Training_Data_tb %>% filter(Label== "negative")


# Display the dimensions of Healthy_Reference matrix (rows columns)
dim(Healthy_Reference)

Healthy_Reference[c(1:7),c(1:7)] # Show first 7 rows


## Here we are calling "Healthy Reference" samples to those samples which are collected from class one (Negative Class).
### Obtaining averaged healthy reference

# Calculating mean value of every gene when subjects are healthy (column name = gene and each row is a different subject)
Averaged_Healthy_Reference <- colMeans(Healthy_Reference[,2:20738])
Averaged_Healthy_Reference <- t(as.matrix(Averaged_Healthy_Reference))

# Display few inicial values of Averaged_Healthy_Reference
Averaged_Healthy_Reference[,1:7]


### Target samples
## Here "Target" is the matrix having all the samples except Healthy Reference samples.

# Select all samples except Healthy Reference samples.
Target <- Training_Data_tb %>% filter(Label== "positive")

# Show Target samples
Target[c(1:7),c(1:7)]

# Show Target dimension
dim(Target)

### Difference between Target and Healthy Reference

# Obtaining difference matrix
TH_Difference <- as.data.frame(Target)
for(i in 1:dim(Target)[1]){
    TH_Difference[i,2:20738] <- as.data.frame(Target[i,2:20738]-Averaged_Healthy_Reference[1,1:20737])
} # end for loop

TH_Difference[1:7,1:7]

write.csv(TH_Difference, file = "Target_Reference_Difference_Matrix.csv")

### Calculating Mean Squared Difference (MSD)
Target_MSD <- colMeans(TH_Difference[,2:20738]^2)

# Show genes with MSD values
head(Target_MSD)

### Rank the Genes based on MSD values

# Ranking in descending order of MSD
Ranked_MSD <- Target_MSD[order(-Target_MSD)]

# Show Ranked_MSD
head(Ranked_MSD)

# Writing/Saving Ranked_MSD
write.csv(Ranked_MSD, file = "Ranked_MSD.csv")

### Prepare the training and test data for feeding into the functions of caret package

# Train and test data
g_train_data_all_genes <- g_train_data
g_test_data_all_genes <- g_test_data

# Converting Label vector into factor as per the requirement of train function of caret package
Labels_g_train_data_all_genes <- as.factor(g_train_data$Label)
Labels_g_test_data_all_genes <- as.factor(g_test_data$Label)

# Converting training and test data into matrix as per the requirement of train function of caret package
g_train_data_all_genes <- as.matrix(g_train_data_all_genes[,-c(1)])
g_test_data_all_genes <- as.matrix(g_test_data_all_genes[,-c(1)])



### Model building using all the Genes

# Enable Parallel Processing
library(doSNOW)
library(doParallel)
cl <- makeCluster(40)
registerDoSNOW(cl)
pt<-proc.time()

set.seed(1234)

# Creating indices for repeated stratified k folds
folds <- createMultiFolds(Labels_g_train_data_all_genes, k = 10, times = 3)

# Creating folds (10 fold cross validation reapeated 3 times)
cross_validation_10_fold <- trainControl(method = "repeatedcv",
                                         number = 10,
                                         repeats = 3,
                                         index = folds)

set.seed(1234)

metric <- "Accuracy"
grid <- expand.grid(k = c(1:30))

# Training the model
train_full_Feature <- train(x= g_train_data_all_genes, # Training Data
                            y = Labels_g_train_data_all_genes,  # Class labels of training data
                            method = "knn", # Train using KNN
                            metric = metric, # Passing "Accuracy" as evaluation matric
                            tuneGrid = grid, # Passing grid for tuning parameters
                            trControl = cross_validation_10_fold) # Passing training control parameters


# Stop Parallel Processing
proc.time()-pt
stopCluster(cl)


# Print trained model
print(train_full_Feature)

# standard deviation for 12023 genes
print((sd(train_full_Feature$resample$Accuracy)))


##################################### test data prediction ####################################


### Test set predicition

# Predicting Test Data
# Passing test data without labels (without fist column which contains labels)
(testPrediction1 <- predict(train_full_Feature, newdata = g_test_data_all_genes))


# Show Test Data
print(Labels_g_test_data_all_genes)


### Performance measure

# Display confusion matrix
(confusionMatrix(testPrediction1, Labels_g_test_data_all_genes))



################### Function to train using desired number of features ####################
##################                                                     ####################

train_any_f_size <- function(f_size){
    
    # Train the Machine Learning classifier using desired number of feature size = f_size
    # print value of feature size
    print(f_size)
    
    # Select top f_size DEGs from Ranked gene vector named Ranked_MSD
    Sel_DEGs <- Ranked_MSD[1:f_size]
    
    # Show number of genes in DEGs
    length(Sel_DEGs)
    
    # Creating a vector containing Differencially Expressed Genes ID as column names
    DEGs_Names_Sel <- names(Sel_DEGs)
    
    
    # Training set that only has DEGs
    g_train_data_DEGs <- g_train_data[ , colnames(g_train_data) %in% DEGs_Names_Sel]
    
    g_train_data_DEGs <- as.matrix(g_train_data_DEGs)
    
    # Assigning column name
    if(f_size==1){
        colnames(g_train_data_DEGs) <- DEGs_Names_Sel
    }
    
    # Show dimension of the training data (having DEGs)
    dim(g_train_data_DEGs)
    
    # Model Building using DEGs
    
    # Enable Parallel Processing
    library(doSNOW)
    library(doParallel)
    cl <- makeCluster(40)
    registerDoSNOW(cl)
    pt<-proc.time()
    
    set.seed(1234)
    
    # creating indices for repeated stratified k folds
    folds <- createMultiFolds(Labels_g_train_data_all_genes, k = 10, times = 3)
    
    # creating folds
    cross_validation_10_fold <- trainControl(method = "repeatedcv",
                                             number = 10,
                                             repeats = 3,
                                             index = folds)
    
    set.seed(1234)
    
    
    metric <- "Accuracy"
    grid <- expand.grid(k = c(1:30))
    
    # training the k nearest neighbour classifier
    train_F_Size <- train(x = g_train_data_DEGs, # Training Data
                          y = Labels_g_train_data_all_genes,  # Class labels of training data
                          method = "knn", # Train using KNN
                          metric = metric, # Passing "Accuracy" as evaluation matric
                          tuneGrid = grid, # Passing grid for tuning parameters
                          trControl = cross_validation_10_fold) # Passing training control parameters
    
    
    # Stop Parallel Processing
    proc.time()-pt
    stopCluster(cl)
    
    # Display the accuracy using selected features
    print(max(train_F_Size$results$Accuracy))
    
    return(train_F_Size)
    
}

###############################  End Function Definition  ########################################



################### Function to return test data using desired number of genes ####################
##################                                                            ####################

test_data_f_size <- function(f_size){
    # evaluate the accuracy using feature size = f_size
    # print value of feature size
    print(f_size)
    
    # Select top f_size DEGs from Ranked gene vector named Ranked_MSD
    Sel_DEGs <- Ranked_MSD[1:f_size]
    
    # Show number of genes in DEGs
    length(Sel_DEGs)
    
    # Creating a vector containing Differencially Expressed Genes ID as column names
    DEGs_Names_Sel <- names(Sel_DEGs)
    
    # Training set that only has DEGs
    g_test_data_DEGs <- g_test_data_all_genes[ , colnames(g_test_data_all_genes) %in% DEGs_Names_Sel]
    
    g_test_data_DEGs <- as.matrix(g_test_data_DEGs)
    
    # assigning column name
    if(f_size==1){
        colnames(g_test_data_DEGs) <- DEGs_Names_Sel
    }
    
    # Show dimension of the training data (having DEGs)
    (dim(g_test_data_DEGs))
    
    return(g_test_data_DEGs)
    
}

###############################  End Function Definition  ########################################


################### Finding F_equal accuracy   ############################################
##################                             ###########################################
# L=left pointer and R is right pointer (pointer denotes feature size)

# inicializtion
First_Sig_equal_Acc <- 0
First_Sig_equal_Model <- 0
First_Sig_equal_Feature_Size <- 0

# L = left index, R = Right index , a = alpha value (t-test significance threshold)

######################## Function definition start ######################################
Search_F_equal <- function(L, R, a){
  
  print(L)
  print(R)
  
  if((L+1)>R){
    print("Left index reaches to right index")
    print("Printing best possible solution")
    return(First_Sig_equal_Model_local)
  }
  
  # Finding mid point
  M = ceiling((L+R)/2)
  
  # Evaluate the accuracy using feature size M-1, M and M+1
  # Show number of DEGs
  print(M)
  
  train_M_Minus_1 <- train_any_f_size(M-1)
  print(train_M_Minus_1)
  
  train_M <- train_any_f_size(M)
  print(train_M)
  
  train_M_Plus_1 <- train_any_f_size(M+1)
  print(train_M_Plus_1)
  
  
  Feature_Acc_Matrix <- matrix(c(M-1, max(train_M_Minus_1$results$Accuracy), M, max(train_M$results$Accuracy), M+1, max(train_M_Plus_1$results$Accuracy)),nrow=1,ncol=6)
  
  print(Feature_Acc_Matrix)
  
  t_test_M <- t.test(train_full_Feature$resample$Accuracy, train_M$resample$Accuracy, paired = TRUE)
  print(t_test_M)
  t_test_M_Plus_1 <- t.test(train_full_Feature$resample$Accuracy, train_M_Plus_1$resample$Accuracy, paired = TRUE)
  print(t_test_M_Plus_1)
  t_test_M_Minus_1 <- t.test(train_full_Feature$resample$Accuracy, train_M_Minus_1$resample$Accuracy, paired = TRUE)
  print(t_test_M_Minus_1)
  
  if(((t_test_M_Minus_1$p.value != "NaN") & (t_test_M_Minus_1$p.value < a) & (max(train_M_Minus_1$results$Accuracy)<max(train_full_Feature$results$Accuracy)) &
      (t_test_M$p.value != "NaN") & (t_test_M_Plus_1$p.value != "NaN") & (t_test_M$p.value > a) & (t_test_M_Plus_1$p.value > a)) |
     ((t_test_M_Minus_1$p.value != "NaN") & (t_test_M_Minus_1$p.value < a) & (max(train_M_Minus_1$results$Accuracy)<max(train_full_Feature$results$Accuracy)) &
      (t_test_M$p.value != "NaN") & (t_test_M$p.value < a) & (max(train_M$results$Accuracy)<max(train_full_Feature$results$Accuracy)) &
      (t_test_M_Plus_1$p.value != "NaN") & (t_test_M_Plus_1$p.value > a))){
    if(t_test_M$p.value > a){
      First_Sig_equal_Acc <<- max(train_M$results$Accuracy)
      First_Sig_equal_Model <<- train_M
      # save the model for First_Sig_equal_Model
      Model_name <- paste("KNN_10CV_First_Sig_equal_Model_",M,".rda", collapse = "",sep="")
      save(train_M, file = Model_name)
      First_Sig_equal_Feature_Size <<- M
      print("First significantly equal accuracy found")
      print("Feature size is: ")
      print(First_Sig_equal_Feature_Size)
      print("First significantly equal accuracy is:")
      print(First_Sig_equal_Acc)
      
      return(train_M)
      
    }else if(t_test_M_Plus_1$p.value > a){
      First_Sig_equal_Acc <<- max(train_M_Plus_1$results$Accuracy)
      First_Sig_equal_Model <<- train_M_Plus_1
      # save the model for First_Sig_equal_Model
      Model_name <- paste("KNN_10CV_First_Sig_equal_Model_",M+1,".rda", collapse = "",sep="")
      save(train_M_Plus_1, file = Model_name)
      First_Sig_equal_Feature_Size <<- M+1
      print("First significantly equal accuracy found")
      print("Feature size is: ")
      print(First_Sig_equal_Feature_Size)
      print("First significantly equal accuracy is:")
      print(First_Sig_equal_Acc)
      
      return(train_M_Plus_1)
    }
    
  }else if(((t_test_M_Minus_1$p.value != "NaN") & (t_test_M_Minus_1$p.value < a) & (max(train_M_Minus_1$results$Accuracy)<max(train_full_Feature$results$Accuracy)) &
            (t_test_M$p.value != "NaN") & (t_test_M$p.value < a) & (max(train_M$results$Accuracy)<max(train_full_Feature$results$Accuracy)) &
            (t_test_M_Plus_1$p.value != "NaN") & (t_test_M_Plus_1$p.value < a) & (max(train_M_Plus_1$results$Accuracy)<max(train_full_Feature$results$Accuracy) ))){
                
    # Update value of L
    L <- M
    # Calling function recursively
    Recall(L=L, R=R, a=a)
    
  }else{
    # Update value of R
    R <- M
    # Calling function recursively
    Recall(L=L, R=R, a=a)
  }
  
}

# Function calling give L, R and alpha in input
# total number of featrues R = lenght of data -1 (all columns - column that has lables)
First_Sig_equal_Model <- Search_F_equal(L=1, R=ncol(g_train_data_all_genes), a=0.05)

# Print the model which gives F_equal accuracy 
print(First_Sig_equal_Model)

# print the First_Sig_equal_Acc
print(First_Sig_equal_Acc)

# print First_Sig_equal_Feature_Size (F_equal)
print(First_Sig_equal_Feature_Size)


### Predicting Test Data using Statistically Equal Model #####################

g_test_data_DEGs <- test_data_f_size(First_Sig_equal_Feature_Size)


# Passing test data without labels (without fist column which contains labels)
testPrediction4 <- predict(First_Sig_equal_Model, newdata = g_test_data_DEGs)

# show test predictions
print(testPrediction4)

# print test data labels
print(Labels_g_test_data_all_genes)


# Performance Measure for DEGs
# Display confusion matrix
(confusionMatrix(testPrediction4, Labels_g_test_data_all_genes))
