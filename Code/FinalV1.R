setwd("C:/Users/Warsec/Desktop/Soil Model")

#Load required libraries
library(ggplot2)#Used for plotting the data
library(reshape2)
library(stringr)
library(grid)
library(prospectr)
library(ptw)
library(RWeka) #Used for PCA, please install rJava for it to work
library(wavelets)#Discrete wavelet transform of spectra variables
library(nnet)
library(e1071)
library(ipred)
library(klaR)
library(caret)#Used for training and predicting models
library(pls)#Used for training model with Partial Least Squares
library(kernlab)#Used for training model using SVM
library(randomForest)#Used for training model using rf


#Additional Functions
#Plotting Spectra
plotSpectra <- function(numberOfSamples, spectralData, subsample, dataDF){
  #based on http://www.kaggle.com/c/afsis-soil-properties/forums/t/10184/first-derivative
  ixs <- sample(which(1:nrow(dataDF) %in% subsample), numberOfSamples)
  trainRawSub <- melt(dataDF[ixs, ], id.vars = "PIDN", measure.vars = spectralData)
  trainRawSub$variable <- as.numeric(str_replace_all(trainRawSub$variable,"m",""))
  ggplot(trainRawSub, aes(x = variable, y = value, colour = PIDN)) + geom_line()
}

# Multiple plot function
#got from : http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))  
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#Load Data
train <- read.csv('training.csv', header = TRUE, stringsAsFactors = FALSE)
test <- read.csv('sorted_test.csv', header = TRUE, stringsAsFactors = FALSE)

#Data Transformation
train <- transform(train, Depth = as.numeric(as.factor(train$Depth)) - 1)
test <- transform(test, Depth = as.numeric(as.factor(test$Depth)) - 1)

#Separating the data to remove the CO2 spectra and predictors
co2spectra <- seq(which(names(train) == 'm2379.76'), which(names(train) == 'm2352.76'))
allspectra <- seq(which(names(train) == 'm7497.96'), which(names(train) == 'm599.76'))
allspectranoco2 <- c(seq(which(names(train) == 'm7497.96'), which(names(train) == 'm2379.76')), 
                       seq(which(names(train) == 'm2352.76'), which(names(train) == 'm599.76')))
spatialpredictors <- seq(which(names(train) == 'BSAN'), which(names(train) == 'TMFI'))
soildepth <- which(names(train) == 'Depth')

#Visualization of data of 10 observations from each sample.
#Subsoil all spectra
samplesubsoil <- plotSpectra(10, spectralData = allspectra, subsample = which(train$Depth == 0), train)
#Topsoil all spectra
sampletopsoil <- plotSpectra(10, spectralData = allspectra, subsample = which(train$Depth == 1), train)
#Subsoil no CO2
samplesubsoilnoco2 <- plotSpectra(10, spectralData = allspectranoco2, subsample = which(train$Depth == 0),train)
#Topsoil no CO2
sampletopsoilnoco2 <- plotSpectra(10, spectralData = allspectranoco2, subsample = which(train$Depth == 1), train)
#Plots the four plots together in 2X2. The top two plots show subsoil sample while the bottom two show topsoil plots 
#The plots in the first column consists of CO2 spectra while the plots in the second column has no Co2 spectra
multiplot(samplesubsoil, sampletopsoil, samplesubsoilnoco2, samplesubsoilnoco2, cols = 2)


#Discrete Wavelet Transforms using Haar Algorithm
#DF1: input matrix for transform
#nTimes: number of iterations
HaarTransform=function(DF1,nTimes=1)
{
  w =function(k)
  {
    s1=dwt(k, filter="haar")
    return (s1@V[[1]])
  }
  Smt=DF1
  for (i in 1:nTimes)
  {
    Smt=t(apply(Smt,1,w))
  }
  return (data.frame(Smt))
}

#Calcule PCA using Weka PrincipalComponents filter
#df: input data frame
#var: variance parameter for PCA
WekPCA=function(df,var)
{
  pc=make_Weka_filter('weka/filters/unsupervised/attribute/PrincipalComponents')
  d1=pc(df[,1]~.,data=df[,-1],control=c('-R',var))
  return (d1[,-ncol(d1)])
}

#Getting Derivatives
#DF1: input matrix for transform
#D: Order
Derivative=function(DF1,D=1)
{
  df1=t(diff(t(DF1), differences = D))
  return(df1)
}


#Seperating the training data into train and validation
iTrain <- createDataPartition(y = train$Depth, p = .80, list = FALSE)
str(iTrain)
training <- train[iTrain,]
testing <- train[-iTrain,]
nrow(training)
nrow(testing)

#Combining the data for preprocessing
trainandtest = rbind(training,testing)

#calcule partial PCs of combined data set
#divide data set into 30 sub frames and then getting PCs
#the combine sub-frames
PC=list()
for (i in 1:30)
{
  j1=(i-1)*119+2
  j2=(i)*119+1
  if (i==30)
  {
    j2=3579
  }
  temp1=trainandtest[,j1:j2]
  flush.console()
  PC[[i]]=WekPCA(cbind(trainandtest['Ca'],temp1),0.999)
}
PComponents=PC[[1]]
for (i in 2:30)
{
  PComponents=cbind(PComponents,PC[[i]])
}

#Multiple Scatter Correction on spectral features(two phases) 
trainandtestreduced = msc(as.matrix(trainandtest[,2:3579]))
trainandtestreduced = msc(as.matrix(trainandtestreduced))

#First Derivatives
trainandtestreduced=Derivative(trainandtestreduced,1)

#Original Dataset 
trainandtestoriginal=trainandtest[2:3600]

#Reduced Dataset
trainandtestreduced=cbind(PComponents, data.frame(HaarTransform(trainandtestreduced,9)),trainandtest[,3580:3600])

#Transformation of P to remove outliers
trainandtestreduced$P = ifelse(trainandtestreduced$P>6,6,trainandtestreduced$P)
trainandtestreduced$P = log(1+trainandtestreduced$P)

#Dividing the Dataset to test and training date frames from reduced Dataset
trainreduced = trainandtestreduced[(1:nrow(training)),]
testreduced = trainandtestreduced[-(1:nrow(training)),]

#Includes the independent variables
xtrain = trainreduced[,1:107]
xtest = testreduced[,1:107]

#Includes the dependent variables (The variables to be predicted)
ytrain = trainreduced[,108:112]
ytest = testreduced[,108:112]

#Setting the cross validation function with 3 fold cross validation 
set.seed(123)
ctrl <- trainControl(method = "repeatedcv",
                     repeats = 3)

#Train and Predict 

#Training using partial least squares (Prediction method p1)                                       
p1Ca = train(x=xtrain, y=ytrain$Ca,
             method="pls",
             tuneLength = 15,
             trControl = ctrl)
p1Capredict = predict(p1Ca, xtest)
p1CaRMSE = RMSE(p1Capredict, ytest$Ca)

p1P = train(x=xtrain, y=ytrain$P,
             method="pls",
             tuneLength = 15,
             trControl = ctrl)    
p1Ppredict = predict(p1P, xtest)
p1PRMSE = RMSE(p1Ppredict, ytest$P)

p1pH = train(x=xtrain, y=ytrain$pH,
             method="pls",
             tuneLength = 15,
             trControl = ctrl)
p1pHpredict = predict(p1pH, xtest)
p1pHRMSE = RMSE(p1pHpredict, ytest$pH)

p1SOC = train(x=xtrain, y=ytrain$SOC,
             method="pls",
             tuneLength = 15,
             trControl = ctrl)
p1SOCpredict = predict(p1SOC, xtest)
p1SOCRMSE = RMSE(p1SOCpredict, ytest$SOC)

p1Sand = train(x=xtrain, y=ytrain$Sand,
             method="pls",
             tuneLength = 15,
             trControl = ctrl)
p1Sandpredict = predict(p1Sand, xtest)
p1SandRMSE = RMSE(p1Sandpredict, ytest$Sand)


#Training using k-Nearest Neighbour (Prediction method p2)  
p2Ca = train(x=xtrain, y=ytrain$Ca,
             method = "knn",
             tuneLength = 20,
             trControl = ctrl)
p2Capredict = predict(p2Ca, xtest)
p2CaRMSE = RMSE(p2Capredict, ytest$Ca)

p2P = train(x=xtrain, y=ytrain$P,
             method = "knn",
             tuneLength = 20,
             trControl = ctrl)
p2Ppredict = predict(p2P, xtest)
p2PRMSE = RMSE(p2Ppredict, ytest$P)

p2pH = train(x=xtrain, y=ytrain$pH,
             method = "knn",
             tuneLength = 20,
             trControl = ctrl)
p2pHpredict = predict(p2pH, xtest)
p2pHRMSE = RMSE(p2pHpredict, ytest$pH)

p2SOC = train(x=xtrain, y=ytrain$SOC,
             method = "knn",
             tuneLength = 20,
             trControl = ctrl)
p2SOCpredict = predict(p2SOC, xtest)
p2SOCRMSE = RMSE(p2SOCpredict, ytest$SOC)

p2Sand = train(x=xtrain, y=ytrain$Sand,
             method = "knn",
             tuneLength = 20,
             trControl = ctrl)
p2Sandpredict = predict(p2Sand, xtest)
p2SandRMSE = RMSE(p2Sandpredict, ytest$Sand)


#Training using Random Forest (Prediction method p3)  
rfGrid = data.frame(.mtry = 2:8)
set.seed(123)
p3Ca = train(x=xtrain, y=ytrain$Ca,
             method = "rf",
             tuneGrid = rfGrid,
             trControl = ctrl)
p3Capredict = predict(p3Ca, xtest)
p3CaRMSE = RMSE(p3Capredict, ytest$Ca)

p3P = train(x=xtrain, y=ytrain$P,
             method = "rf",
            tuneGrid = rfGrid,
            trControl = ctrl)
p3Ppredict = predict(p3P, xtest)
p3PRMSE = RMSE(p3Ppredict, ytest$P)

p3pH = train(x=xtrain, y=ytrain$pH,
             method = "rf",
             tuneGrid = rfGrid,
             trControl = ctrl)
p3pHpredict = predict(p3pH, xtest)
p3pHRMSE = RMSE(p3pHpredict, ytest$pH)

p3SOC = train(x=xtrain, y=ytrain$SOC,
             method = "rf",
             tuneGrid = rfGrid,
             trControl = ctrl)
p3SOCpredict = predict(p3SOC, xtest)
p3SOCRMSE = RMSE(p3SOCpredict, ytest$SOC)

p3Sand = train(x=xtrain, y=ytrain$Sand,
             method = "rf",
             tuneGrid = rfGrid,
             trControl = ctrl)
p3Sandpredict = predict(p3Sand, xtest)
p3SandRMSE = RMSE(p3Sandpredict, ytest$Sand)


#Training using svm (Prediction method p3)  
p4Ca = train(x = xtrain, y=ytrain$Ca,
             method = "svmLinear",
             trControl = ctrl)
p4Capredict = predict(p4Ca, xtest)
p4CaRMSE = RMSE(p4Capredict, ytest$Ca)

p4P = train(x = xtrain, y=ytrain$P,
             method = "svmLinear",
             trControl = ctrl)
p4Ppredict = predict(p4P, xtest)
p4PRMSE = RMSE(p4Ppredict, ytest$P)

p4pH = train(x = xtrain, y=ytrain$pH,
             method = "svmLinear",
             trControl = ctrl)
p4pHpredict = predict(p4pH, xtest)
p4pHRMSE = RMSE(p4pHpredict, ytest$pH)

p4SOC = train(x = xtrain, y=ytrain$SOC,
             method = "svmLinear",
             trControl = ctrl)
p4SOCpredict = predict(p4SOC, xtest)
p4SOCRMSE = RMSE(p4SOCpredict, ytest$SOC)

p4Sand = train(x = xtrain, y=ytrain$Sand,
             method = "svmLinear",
             trControl = ctrl)
p4Sandpredict = predict(p4Sand, xtest)
p4SandRMSE = RMSE(p4Sandpredict, ytest$Sand)


#Training using lm (Prediction method p3)  
p5Ca = train(x = xtrain, y=ytrain$Ca,
             method = "lm",
             trControl = ctrl)
p5Capredict = predict(p5Ca, xtest)
p5CaRMSE = RMSE(p5Capredict, ytest$Ca)

p5P = train(x = xtrain, y=ytrain$P,
             method = "lm",
             trControl = ctrl)
p5Ppredict = predict(p5P, xtest)
p5PRMSE = RMSE(p5Ppredict, ytest$P)

p5pH = train(x = xtrain, y=ytrain$pH,
             method = "lm",
             trControl = ctrl)
p5pHpredict = predict(p5pH, xtest)
p5pHRMSE = RMSE(p5pHpredict, ytest$pH)

p5SOC = train(x = xtrain, y=ytrain$SOC,
             method = "lm",
             trControl = ctrl)
p5SOCpredict = predict(p5SOC, xtest)
p5SOCRMSE = RMSE(p5SOCpredict, ytest$SOC)

p5Sand = train(x = xtrain, y=ytrain$Sand,
             method = "lm",
             trControl = ctrl)
p5Sandpredict = predict(p5Sand, xtest)
p5SandRMSE = RMSE(p5Sandpredict, ytest$Sand)

#Cross Validation for best method using RMSE as metric
bestCamethod = min(p1CaRMSE,p2CaRMSE,p3CaRMSE,p4CaRMSE,p5CaRMSE)
#p3CaRMSE = 0.2600949
bestPmethod = min(p1PRMSE,p2PRMSE,p3PRMSE,p4PRMSE,p5PRMSE)
#p3PEMSE = 0.2724925
bestpHmethod = min(p1pHRMSE,p2pHRMSE,p3pHRMSE,p4pHRMSE,p5pHRMSE)
#p5pHRMSE = 0.3481976
bestSOCmethod = min(p1SOCRMSE,p2SOCRMSE,p3SOCRMSE,p4SOCRMSE,p5SOCRMSE)
#p5SOCRMSE=0.2986313
bestSandmethod = min(p1SandRMSE,p2SandRMSE,p3SandRMSE,p4SandRMSE,p5SandRMSE)
#p5SandRMSE=0.3194639

#RMSE scores of predictor varaibles:
RMSEvalues = cbind(bestCamethod,bestPmethod,bestpHmethod,bestSOCmethod,bestSandmethod)



