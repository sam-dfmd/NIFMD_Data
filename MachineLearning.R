
library(FSelector)
library(e1071)
library(randomForest)

########Set working directory for Data input
#setwd("... /directory file")
#1############### Support Vector Machines
###optimal parameter,  kernel gamma cost
#1 linear   0.5    1
svm.sero <- function(dat.train, dat.test, serotyp.train, serotyp.test){
  model_svm <- svm(x = dat.train, y = serotyp.train, type = "C-classification",kernel = "linear", gamma = 0.5, cost = 1, tolerance = 0.001, epsilon = 0.1)
  #print(model_svm)
  #summary(model_svm)
  ###prediction with trained model
  pred.svm <- predict(model_svm, dat.test)
  #print(pred.svm)
  confus.svm <- table(serotyp.test, pred.svm)
  names(dimnames(confus.svm)) <- c("Actual", "Predicted")
  #accuracy <- function(x) sum(diag(x)/sum(x)) * 100
  #error <- 1 - accuracy
  return(confus.svm)
}
#2##########KNN
library(class)
# best k: 4
knn.sero <- function(dat.train, dat.test, serotyp.train, serotyp.test){
  model.knn <- knn(train = dat.train, test = dat.test, cl = serotyp.train, k = 4, prob=TRUE)
  #print(model.knn)
  #summary(model.knn)
  confus.knn <- table(serotyp.test, model.knn)
  names(dimnames(confus.knn)) <- c("Actual", "Predicted")
  return(confus.knn)
}

#3############ANN
library(neuralnet)
ann.sero <- function(dat.train, dat.test, serotyp.train, serotyp.test){
  form <- as.formula(paste("binary.labels ~", paste(colnames(dat.train), collapse = " + ")))
  model.nn <-  neuralnet(formula = form, data = dat.train, hidden = 3, err.fct = "ce",
                         linear.output = FALSE)
  pred.nn <- compute(model.nn, dat.test, rep = 1)$net.result
  #cbind(pred.nn$net.result, serotyp.test)
  sertype <- levels(serotyp.test)
  colnames(pred.nn) <- sertype
  fx <- function(x) as.factor(sertype[which.max(x)])
  serotyp.pred <- apply(pred.nn, 1, fx)
  out <- cbind(pred.nn, serotyp.pred)
  confus.nn <- table(serotyp.test, serotyp.pred)
  names(dimnames(confus.nn)) <- c("Actual", "Predicted")
  return(confus.nn)
}

#4##############Decision Tree
library(caret)
library(rpart)
DecTree.sero <- function(dat.train, dat.test, serotyp.train, serotyp.test){
  model.tre <- rpart(serotyp.train ~., data = dat.train)
  #rpart.plot(model.tre)
  pred.tre <- predict(model.tre, dat.test)
  fxx <- function(x) as.factor(names(x)[which.max(x)])
  serotyp.pred.tre <- apply(pred.tre, 1, fxx)
  confus.tre <- table(serotyp.test, serotyp.pred.tre)
  names(dimnames(confus.tre)) <- c("Actual", "Predicted")
  return(confus.tre)
}

#5#########Naive Bayes
nb.sero <- function(dat.train, dat.test, serotyp.train, serotyp.test){
  model.nb <- naiveBayes(x = dat.train, y = serotyp.train, laplace = 0)
  pred.nb <- predict(model.nb, dat.test)
  confus.nb <- table (serotyp.test, pred.nb)
  names(dimnames(confus.nb)) <- c("Actual", "Predicted")
  return(confus.nb)
}
#6############Multi-class logistic classification
library(nnet)
mls.sero <- function(dat.train, dat.test, serotyp.train, serotyp.test){
  model.mls <- multinom(serotyp.train ~., data = dat.train)
  #summary(model.mls)
  pred.mls <- predict(model.mls, dat.test)
  confus.mls <- table(serotyp.test, pred.mls)
  names(dimnames(confus.mls)) <- c("Actual", "Predicted")
  return(confus.mls)
}
#7########linear discriminat analysis
library(MASS)
lda.sero <- function(dat.train, dat.test, serotyp.train, serotyp.test){
  model.lda <- lda(serotyp.train ~., data = dat.train)
  #summary(model.lda)
  pred.lda <- predict(model.lda, dat.test)
  confus.lda <- table(serotyp.test, pred.lda$class)
  names(dimnames(confus.lda)) <- c("Actual", "Predicted")
  return(confus.lda)
}