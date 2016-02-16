library(e1071); library(prospectr)

# Data preparation; first derivative of half of spectra via SG (all the 'good stuff' seems to be in that half)

myData <- read.csv("training.csv", header = TRUE)

sg21 <- savitzkyGolay(myData[ , c(1800:3579)], p = 3, w = 21, m = 1)

sg21 <- data.frame(sg21)

# Final models (after much CV exploration of parameters)

svm.Ca <- svm(sg21, myData$Ca, cost = 50, gam = 0.00005, ep = 0.01)

svm.SOC <- svm(sg21, myData$SOC, cost = 40, gam = 0.00005, ep = 0.05)

svm.pH <- svm(sg21, myData$pH, cost = 30, gam = 0.0001, ep = 0.15)

svm.P <- svm(sg21, myData$P, cost = 10, gam = 0.0001, ep = 0.1)

svm.Sand <- svm(sg21, myData$Sand, cost = 5, gam = 0.0001, ep = 0.12)

# Read test data and make predictions

testX <- read.csv("sorted_test.csv", header = TRUE)

test.sg21 <- data.frame(savitzkyGolay(testX[, c(1800:3579)],p = 3, w = 21, m = 1))

sub <- read.csv("sub Oct-04 sg21.csv", header = TRUE)

sub$Ca <- predict(svm.Ca, test.sg21)

sub$SOC <- predict(svm.SOC, test.sg21)

sub$pH <- predict(svm.pH, test.sg21)

sub$P <- predict(svm.P, test.sg21)

sub$Sand <- predict(svm.Sand, test.sg21)

write.csv(sub, file = "sub xxxxx.csv", row.names = FALSE)