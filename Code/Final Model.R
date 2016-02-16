setwd("../Soil Model")

#library 
library(caret)
library(prospectr)

#Import test and training data
trainOriginal=read.csv("training.csv")
testOriginal=read.csv("sorted_test.csv")


#Changing the categories to 0 and 1
trainOriginal$Depth <-  with ( trainOriginal, ifelse ( ( Depth == 'Subsoil' ), 0 , 1 ) )
testOriginal$Depth <-  with ( testOriginal, ifelse ( ( Depth == 'Subsoil' ), 0 , 1 ) ) 

#Delecting the variables with CO2 spectra band
trainOriginal=trainOriginal[,-(2656:2670)]
testOriginal=testOriginal[,-(2656:2670)]

#using only the mid-infrared absorbance measurements
Mytrain=((trainOriginal[,2:3564]))
Mytest=((testOriginal[,2:3564]))

#Smoothing the noise
Mytrain.filtered=savitzkyGolay(X = Mytrain, 0, 3, 11)
Mytest.filtered=savitzkyGolay(X = Mytest, 0, 3, 11)

Mytrain.filtered.haar=HaarTransform(Mytrain,4)
Mytest.filtered.haar=HaarTransform(Mytest,4)

Mytrain.filtered=cbind(Mytrain.filtered, trainOriginal[c(3565:3585)])
Mytest.filtered=cbind(Mytest.filtered, testOriginal[c(3565:3580)])
                                              
set.seed(1234)
# 10 fold cv
indx <- createFolds(trainOriginal[,1], returnTrain = TRUE)
ctrl <- trainControl(method = "cv", index = indx)

#predict Ca
lmTuneCa <- train(x = Xtrainfiltered, y = Ytrain$Ca,
                  method = "lm",
                  trControl = ctrl)
lmTuneCa

write.csv(submission,'FinalSubmission.csv',row.names=F,quote=F,na='')











,