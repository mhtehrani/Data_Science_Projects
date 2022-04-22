################################################################################
#Libraries
################################################################################
library(data.table)
library(stats)
library(caret)
library(Metrics)
library(VIM) 
library(mice)
library(car) 
library(ggplot2)
library(plyr)
library(corrplot)
library(nnet)
library(pls)
library(randomForest)
library(fastmatch)
library(zoo)
library(xgboost)
library(base)

################################################################################
#Functions
################################################################################
#Convert reflectivity (dbz) to mm/hr
marshall_palmer <- function(dbz) {
  ((10**(dbz/10))/200) ** 0.625
}

#user defined metric MAE 
MAE.eval.f<- function(pred,dtrain) {
  labels <- getinfo(dtrain, "label")
  labels = expm1(labels)
  pred = expm1(pred)
  MAE = mean(abs(pred-labels))
  return(list(metric = 'MAE',value=MAE))
}

#User defined MAE.f to calculate overall MAE according to  differnet distance
MAE.f <- function(tr,te,model) {
  cs <- c("ref", "ref1",   "mp", "kdp", "zdr", "rho","rd") # for mdoels based on dist. 
  MAE=rep(0,16)
  dim_te=rep(0,16)
  No.Id=rep(0,16)
  pred<-as.data.frame(te)[,1:2]
  pred[,1:2]<-0
  pred.value=rep(0,dim(te)[1])
  count=1
  for (i in 1:16) {
    tr_temp =tr[tr$rd==i-1,]
    te_temp =te[te$rd==i-1,]
    te_Id_temp=te_temp[,Id]
    te.target_temp=te_temp$target
    mean(te.target_temp) 
    dim=dim(tr_temp)[1]
    val.dim=round(0.2*dim)
    train.dim=dim-val.dim
    ############# "xgb"model
      # prepare traing data
      y_temp<-tr_temp$target  # values to be predict in the training data
      mean(y_temp) 
      y_te_temp<-te_temp$target
      mean(y_te_temp) 
      tr_temp<-as.data.frame(tr_temp)
      tr_temp<-tr_temp[,cs] 
      tr_temp_val=tr_temp[(train.dim+1):(train.dim+val.dim),]
      xgval_temp = xgb.DMatrix(as.matrix(tr_temp[(train.dim+1):(train.dim+val.dim),]), label = y_temp[(train.dim+1):(train.dim+val.dim)], missing = NA)
      # prepare test data
      te_temp<-as.data.frame(te_temp)
      te_temp=te_temp[,cs ]
      xgtest_temp = xgb.DMatrix(as.matrix(te_temp),  missing = NA)
      xgtest_temp2 = xgb.DMatrix(as.matrix(te_temp), label =  te.target_temp, missing = NA)
      gc()
      ##build local model
      watchlist1<-list(val=xgval_temp,train=xgtrain_temp)
      watchlist1<-list(val=xgtest_temp2,train=xgtrain_temp)
      set.seed(1)
      x.mod_temp  <- xgb.train(params = param0, 
                               data = xgtrain_temp ,
                               nrounds =3000,
                               watchlist = watchlist1,
                               maximize = FALSE,
                               feval = MAE.eval.f,
                               verbose = 1,
                               early.stop.round = 50
      ) 
      pr_temp <- predict(x.mod_temp,xgtest_temp)
      mean(pr_temp) 
      pr_temp = expm1(pr_temp)  # transform back log1p
      mean(pr_temp) 
      
      te.target_temp = expm1(te.target_temp)
      mean(te.target_temp)    
      # MAE
      MAE_temp = mean(abs(pr_temp-te.target_temp))

    MAE[i]=MAE_temp
    No.Id[i]=dim(te_temp)[1]
    dim_te[i]=dim(te_temp)[1]
    count_temp=dim(as.data.frame(pr_temp))[1]
    pred[count:(count+count_temp-1),2]=as.data.frame(pr_temp)    #target
    pred[count:(count+count_temp-1),1]=te_Id_temp  # Id
    count=count+count_temp
  }
  # for end
  
  MAE_sum=0
  for (i in 1:16) {
    MAE_sum = MAE_sum + MAE[i]*dim_te[i]
  }
  MAE_sum = MAE_sum/dim(te)[1]
  MAE<-list("MAE"=MAE,"dim_te"=dim_te,"MAE sum"= MAE_sum,"No.Id"=No.Id,"predicted test value"=pred)  # output list
}

#User defined function to calculate duration of each record
duration.f <- function(tr_raw) {
  row_id=1
  k=1
  time_int=rep(0,length(tr_raw))
  duration=rep(0,length(tr_raw))
  time.point.id=rep(0,10)
  
  for (i in 1:length(tr_raw$Id)) {
    minutes_past=rep(0,1)
    
    if (i==1 )  {
      ID=tr_raw$Id[i]
      minutes_past=tr_raw[tr_raw$Id==ID,]$minutes_past
      duration.local=rep(0,length(minutes_past)) # initialize
      
      # calculate time interval
      time_int=rep(0,length(minutes_past)+1)    # initialize
      time_int[length(minutes_past)+1]=60-minutes_past[length(minutes_past)]   #last one
      time_int[1]=minutes_past[1]
      
      if (length(minutes_past)==1) {
        # if observation in id =1   
        duration.local[1]=60
      } else {   
        # if observation in id >=1
        for (j in 2:length(minutes_past)){
          time_int[j]=minutes_past[j] - minutes_past[j-1]
        }
      }  # if end
      # get the local duration time at each time point
      for (j in 1:length(minutes_past)){
        if (j==1){
          if (length(minutes_past)==1){
            duration.local[1]=60
          }else{
            duration.local[j]=time_int[j] + 0.5*time_int[j+1]
          }
        } else if (j==length(minutes_past)) {
          duration.local[j]= 0.5*time_int[j] + time_int[j+1]
          
        } else {
          duration.local[j]= 0.5*time_int[j] + 0.5*time_int[j+1]
        }
      }
      time.point.id[k]=length(duration.local)
      duration[row_id:(row_id+j-1)] = duration.local
      row_id = j+row_id
      k=k+1
      #str(train.1)
      #sum(is.na(train.1))
      
    }  else if (tr_raw$Id[i]>tr_raw$Id[i-1] )  {
      ID=tr_raw$Id[i]
      minutes_past=0
      minutes_past=tr_raw[tr_raw$Id==ID,]$minutes_past
      duration.local=rep(0,length(minutes_past)) # initialize
      
      # calculate time interval
      time_int=rep(0,length(minutes_past)+1)    # initialize
      
      time_int[length(minutes_past)+1]=60-minutes_past[length(minutes_past)]   #last one
      
      #time_int=rep(0,length(minutes_past))
      
      time_int[1]=minutes_past[1]
      
      if (length(minutes_past)==1) {
        # if observation in id only 1
        
        duration.local[1]=60
      } else {   
        
        # if observation in id >=1
        for (j in 2:length(minutes_past)){
          
          time_int[j]=minutes_past[j] - minutes_past[j-1]
        }
      }  # if end
      
      # get the local duration time at each time point
      for (j in 1:length(minutes_past)){
        if (j==1){
          
          if (length(minutes_past)==1){
            
            duration.local[1]=60
          }else{
            duration.local[j]=time_int[j] + 0.5*time_int[j+1]
          }
          
        } else if (j==length(minutes_past)) {
          duration.local[j]= 0.5*time_int[j] + time_int[j+1]
          
        } else {
          duration.local[j]= 0.5*time_int[j] + 0.5*time_int[j+1]
        }
      }# for end
      
      time.point.id[k]=length(minutes_past)
      duration[row_id:(row_id+j-1)] = duration.local
      row_id = j+row_id
      k=k+1
    }        
    
  }
  # for loop END
  
  duration
}

################################################################################
#Data Preparation
################################################################################
Raw.Data <- fread("train.csv")

#Use the complete cases for the rest of process So there would be no missing and there
#is no need of imputation. Still we have more than 2.7 million observation.
Rain.Training <- Raw.Data[complete.cases(Raw.Data)==T,]

#Removing data with expected raining more than 1000mm (It seems unrealistic)
Rain.Training <- Rain.Training[Rain.Training$Expected<1000,]

#Calculate time duration for each observation
Rain.Training$dt <- c(Rain.Training$minutes_past[1], diff(Rain.Training$minutes_past))
Rain.Training$dt[Rain.Training$dt<0] <- Rain.Training$minutes_past[Rain.Training$dt<0]

#Calculate Rain Rate based on Marshal Palmer formula (mm/hr)
Rain.Training$RainRate <- marshall_palmer(Rain.Training$Ref)

Rain.Training <- as.data.frame(Rain.Training)
Rain.Training <- cbind(Rain.Training[,(1:2)], dt=Rain.Training$dt,
                       Rain.Training[,(3:4)], RainRate=Rain.Training$RainRate, 
                       Rain.Training[,(5:24)])

write.csv(Rain.Training, file ="Rain.Training.Full.Data.csv")

################################################################################
#Stochastic Gradient Boosting
################################################################################
#subsetting the data to training and test (Data are already shuffeled)
Training <- Rain.Training[1:2000000,]
Test <- Rain.Training[(2000001:dim(Rain.Training)[1]),]

#Aggregating the data for Training set
trainData <- Training[,.(
  target = log1p(mean(Expected)),
  radardist = mean(radardist_km),
  Ref = mean(Ref),
  Ref10 = mean(Ref_5x5_10th),
  Ref50 = mean(Ref_5x5_50th),
  Ref90 = mean(Ref_5x5_90th),
  RefComposite = mean(RefComposite),
  RefComposite10 = mean(RefComposite_5x5_10th),
  RefComposite50 = mean(RefComposite_5x5_50th),
  RefComposite90 = mean(RefComposite_5x5_90th),
  RhoHV = mean(RhoHV),
  mp = mean(RainRate * dt/60),
  Kdp = mean(Kdp),
  Zdr = mean(Zdr)
), Id]

#Aggregating the data for Test set
testData <- Test[,.(
  target = log1p(mean(Expected)),
  radardist = mean(radardist_km),
  Ref = mean(Ref),
  Ref10 = mean(Ref_5x5_10th),
  Ref50 = mean(Ref_5x5_50th),
  Ref90 = mean(Ref_5x5_90th),
  #sumRef = sum(Ref * dt),
  RefComposite = mean(RefComposite),
  RefComposite10 = mean(RefComposite_5x5_10th),
  RefComposite50 = mean(RefComposite_5x5_50th),
  RefComposite90 = mean(RefComposite_5x5_90th),
  #sumRefComposite = sum(RefComposite * dt),
  RhoHV = mean(RhoHV),
  mp = mean(RainRate * dt/60),
  Kdp = mean(Kdp),
  Zdr = mean(Zdr)
), Id]

testData$Order <- 1:dim(testData)[1]

write.csv(trainData, file ="trainData.csv")
write.csv(testData, file ="testData.csv")

#Predictors we want to use for modeling and test
x <- c(4:15) #Peredictors range

#Matrix containing for each Metrics for training data for each distance
Metric <- data.frame(Distance=rep(0, 16), R2=rep(0, 16), RMSE=rep(0, 16), MAE=rep(0, 16), No=rep(0, 16))

for(i in 1:16){
  trainData <- read.csv("trainData.csv")
  trainData <- trainData[,2:dim(trainData)[2]]
  trainData <- trainData[trainData$radardist==(i-1),]

  Cont <- trainControl(method="repeatedcv", number=5, repeats=1)
  gbmGrid <- expand.grid(interaction.depth=c(2,3,4),  n.trees=(1:40)*100, shrinkage=0.1, n.minobsinnode=20)
  Fit.SGB <- train(x=trainData[,x], y=trainData$target, method='gbm', trControl=Cont, tuneGrid=gbmGrid, verbose=FALSE)
  
  if(i==1){Fit00 <- Fit.SGB}
  if(i==2){Fit01 <- Fit.SGB}
  if(i==3){Fit02 <- Fit.SGB}
  if(i==4){Fit03 <- Fit.SGB}
  if(i==5){Fit04 <- Fit.SGB}
  if(i==6){Fit05 <- Fit.SGB}
  if(i==7){Fit06 <- Fit.SGB}
  if(i==8){Fit07 <- Fit.SGB}
  if(i==9){Fit08 <- Fit.SGB}
  if(i==10){Fit09 <- Fit.SGB}
  if(i==11){Fit10 <- Fit.SGB}
  if(i==12){Fit11 <- Fit.SGB}
  if(i==13){Fit12 <- Fit.SGB}
  if(i==14){Fit13 <- Fit.SGB}
  if(i==15){Fit14 <- Fit.SGB}
  if(i==16){Fit15 <- Fit.SGB}
  
  True.SGB <- exp(trainData$target)-1
  Pred.SGB <- exp(predict(Fit.SGB, newdata=trainData[,x]))-1
  
  Metric[i,1] <- i-1
  Metric[i,2] <- (sum((Pred.SGB-mean(True.SGB)) ** 2))/(sum((True.SGB-mean(True.SGB)) ** 2))
  Metric[i,3] <- rmse(True.SGB, Pred.SGB)
  Metric[i,4] <- mae(True.SGB, Pred.SGB)
  Metric[i,5] <- dim(trainData)[1]
}

#Calculating Predictions on Test Data
Metric.Test <- data.frame(Distance=rep(0, 16), R2=rep(0, 16), RMSE=rep(0, 16), MAE=rep(0, 16), No=rep(0, 16))
Whole.test.Data <- data.frame(Observed=rep(0, 1), Prediction=rep(0, 1), Order=rep(0, 1), Id=rep(0, 1))
for(i in 1:16){
  testData <- read.csv("trainData.csv")
  testData <- testData[,2:dim(testData)[2]]
  testData <- testData[testData$radardist==(i-1),]

  st <- dim(Whole.test.Data)[1]
  
  if(i==1){Fit.SGB <- Fit00
  main="distance=00"
  st <- 0}
  if(i==2){Fit.SGB <- Fit01
  main="distance=01"}
  if(i==3){Fit.SGB <- Fit02
  main="distance=02"}
  if(i==4){Fit.SGB <- Fit03
  main="distance=03"}
  if(i==5){Fit.SGB <- Fit04
  main="distance=04"}
  if(i==6){Fit.SGB <- Fit05
  main="distance=05"}
  if(i==7){Fit.SGB <- Fit06
  main="distance=06"}
  if(i==8){Fit.SGB <- Fit07
  main="distance=07"}
  if(i==9){Fit.SGB <- Fit08
  main="distance=08"}
  if(i==10){Fit.SGB <- Fit09
  main="distance=09"}
  if(i==11){Fit.SGB <- Fit10
  main="distance=10"}
  if(i==12){Fit.SGB <- Fit11
  main="distance=11"}
  if(i==13){Fit.SGB <- Fit12
  main="distance=12"}
  if(i==14){Fit.SGB <- Fit13
  main="distance=13"}
  if(i==15){Fit.SGB <- Fit14
  main="distance=14"}
  if(i==16){Fit.SGB <- Fit15
  main="distance=15"}
  
  plot(Fit.SGB, main=main)
  plot(Fit.SGB, metric="RMSE", plotType="level", scales=list(x=list(rot=90)), main=main)
  
  True.SGB <- exp(testData$target)-1
  Pred.SGB <- exp(predict(Fit.SGB, newdata=testData[,x]))-1
  plot(True.SGB, Pred.SGB, xlim=c(0,50), ylim=c(0,50))
  
  Metric.Test[i,1] <- i-1
  Metric.Test[i,2] <- (sum((Pred.SGB-mean(True.SGB)) ** 2))/(sum((True.SGB-mean(True.SGB)) ** 2))
  Metric.Test[i,3] <- rmse(True.SGB, Pred.SGB)
  Metric.Test[i,4] <- mae(True.SGB, Pred.SGB)
  Metric.Test[i,5] <- dim(testData)[1]
  
  Whole.test.Data[(st+1):(st+dim(testData)[1]),1] <- True.SGB
  Whole.test.Data[(st+1):(st+dim(testData)[1]),2] <- Pred.SGB
  Whole.test.Data[(st+1):(st+dim(testData)[1]),3] <- testData$Order
  Whole.test.Data[(st+1):(st+dim(testData)[1]),4] <- testData$Id
}

Whole.test.Data <- Whole.test.Data[order(Whole.test.Data$Order),] 
MAE <- mae(Whole.test.Data$Observed, Whole.test.Data$Prediction)

#Calculating R2
True.SGB <- log1p(Whole.test.Data$Observed)
Pred.SGB <- log1p(Whole.test.Data$Prediction)
R2 <- (sum((Pred.SGB-mean(True.SGB)) ** 2))/(sum((True.SGB-mean(True.SGB)) ** 2))

plot(Whole.test.Data$Observed, Whole.test.Data$Prediction, xlim=c(0,50), ylim=c(0,50),
     main="Predicted versus Observed values for Test data", xlab="Observed", ylab="Predicted",
     col="dark red")

################################################################################
# Ordinary Least Squares, PLS, Neural Network, Random Forest
################################################################################
# input data set
rawData<-read.table("Rain.Training.Full.Data.csv",header=TRUE,sep=",")
training<-rawData[1:1000000,]
max(training$Expected)
testing<-rawData[2000001:2763588,]
mean(testing$Expected)

# data visualization
# histogram of ref, rhoHV, Zdr, Kdp
hist(training$Ref,main=" ",xlab="Ref")
hist(training$RhoHV,main=" ",xlab="RhoHV")
hist(training$Zdr,main=" ",xlab="Zdr")
hist(training$Kdp,main=" ",xlab="Kdp")

# response transformation
par(mfrow=c(1,2))
#hist(training$Expected)
max(training$Expected)
symbox(training$Expected, data=training$Expected, powers=c(3,2,1,0,-0.5,-1,-2),ylab="Expected")
hist(log1p(training$Expected),main=" ",xlab="log(Expected)")

# Visualize the distribution of predictor
par(mfrow=c(2,6))
#apply(Glass[,1:9],2,hist)
#apply(Glass[,1:9],2,boxplot)
plotTraining<-training[,1:12]
head(training)
head(plotTraining)
for (i in 1:length(plotTraining)) {
  hist(plotTraining[,i], main=names(plotTraining[i]),xlab=" ", type="l")
}
plotTraining<-training[,13:24]
head(training)
head(plotTraining)
for (i in 1:length(plotTraining)) {
  hist(plotTraining[,i], main=names(plotTraining[i]),xlab=" ", type="l")
}

corrplot(cor(training),order="hclust")

# data preparation
training<-training[,-c(1,4,7)]
testing<-testing[,-c(1,4,7)]
selectVar<-c(1,2,3,4,8,12,16,20,24)
training<- training[,selectVar]
testing<-testing[,selectVar]

# training set
training$dt<-c(training$minutes_past[1],diff(training$minutes_past))
training$dt[training$dt<0]=training$minutes_past[training$dt<0]
head(training$dt)
training$RainfallRate <- ((10**(training$Ref/10))/200) ** 0.625
head(training$RainfallRate)
head(training)

trainingNew<- aggregate(cbind(training$radardist_km,training$Ref,
                              training$RefComposite,training$RhoHV,training$Zdr,
                              training$Kdp,training$Expected)~ Id,data=training,mean)
View(trainingNew)
trainingNew <- rename(trainingNew,c("V1"="radarDist","V2"="Ref","V3"="RefComposite",
                                    "V4"="RhoHV","V5"="Zdr","V6"="Kdp","V7"="Expected"))
AccuRain<- aggregate(training$dt*training$RainfallRate~Id,data=training,sum)
trainingNew$Rain<-AccuRain[,2]

# test set
testing$dt<-c(testing$minutes_past[1],diff(testing$minutes_past))
testing$dt[testing$dt<0]=testing$minutes_past[testing$dt<0]
testing$RainfallRate <- ((10**(testing$Ref/10))/200) ** 0.625
testNew<- aggregate(cbind(testing$radardist_km,testing$Ref,
                          testing$RefComposite,testing$RhoHV,testing$Zdr,
                          testing$Kdp,testing$Expected)~ Id,data=testing,mean)
View(testNew)
testNew <- rename(testNew,c("V1"="radarDist","V2"="Ref","V3"="RefComposite",
                            "V4"="RhoHV","V5"="Zdr","V6"="Kdp","V7"="Expected"))
AccuRain<- aggregate(testing$dt*testing$RainfallRate~Id,data=testing,sum)
testNew$Rain<-AccuRain[,2]

trainingNew<-trainingNew[,-1]
head(trainingNew)

testNew<-testNew[,-1]
head(testNew)

# model: OLS
ctrl<-trainControl(method="cv",number=10)
lmFit<-train(log1p(Expected)~.,data=trainingNew,method="lm",preProc=c("center","scale"),trControl=ctrl)
lmFit
plot(lmFit)
lmPred<-predict(lmFit,testNew)
lmValues<-data.frame(obs=log1p(testNew$Expected),pred=lmPred)
#lmPred<-predict(lmFit)
#lmValues<-data.frame(obs=log1p(trainingNew$Expected),pred=lmPred)
MAE<-mean(abs(testNew$Expected-(exp(lmPred)-1)))
MAE
defaultSummary(lmValues)
xyplot(log1p(testNew$Expected)~lmPred,type=c("p","g"),xlab="Predicted",ylab="Observed",xlim=c(0,8))
xyplot(resid(lmFit)~predict(lmFit),type=c("p","g"),xlab="Predicted",ylab="Residual")

# model: PLS
plsFit<-train(log1p(Expected)~.,data=trainingNew,method="pls",tuneLength=20,preProc=c("center","scale"),trControl=ctrl)
plsFit
plot(plsFit)
lmPred2<-predict(plsFit,testNew)
lmValues2<-data.frame(obs=log(testNew$Expected),pred=lmPred2)
defaultSummary(lmValues2)
MAE<-mean(abs(testNew$Expected-(exp(lmPred2)-1)))
MAE
#check model assumption
xyplot(log1p(testNew$Expected)~lmPred2,type=c("p","g"),xlab="Predicted",ylab="Observed")
xyplot(resid(plsFit)~predict(plsFit),type=c("p","g"),xlab="Predicted",ylab="Residual")

# model: Neural network
nnetGrid <- expand.grid(.decay = c(0,0.01, .1),
                        .size = c(12:18), .bag = FALSE)
nnetTune <- train(log(Expected)~.,data=trainingNew,
                  method = "avNNet",
                  tuneGrid = nnetGrid,
                  trControl = ctrl,
                  preProc = c("center", "scale"),
                  linout = TRUE,
                  trace = FALSE,
                  MaxNWts = 18 * ncol(trainingNew)  + 18 + 1)
nnetTune
plot(nnetTune)
lmPred3<-predict(nnetTune)
lmValues3<-data.frame(obs=log1p(testNew$Expected),pred=lmPred3)
MAE<-mean(abs(testNew$Expected-(exp(lmPred3)-1)))
MAE
defaultSummary(lmValues3)

#model: random forest
RFfit<- randomForest(log1p(Expected)~., data = trainingNew, importance = T, ntree=500, mtry=2)
RFfit
plot(RFfit)
lmPred4<-predict(RFfit,testNew)
lmValues4<-data.frame(obs=log1p(testNew$Expected),pred=lmPred4)
defaultSummary(lmValues4)
MAE<-mean(abs(testNew$Expected-(exp(lmPred4)-1)))
MAE
varImpPlot(RFfit)

################################################################################
# xgboost
################################################################################
## read full date set
tr_raw0=fread("Rain.Training.Full.Data.csv")
# assign Train and Test Data
tr_raw=tr_raw0[1:1000000,]
te_raw=tr_raw0[2000001:2763588,]

# Preparation of Training data
# remove the outliers in training data
tr_raw=tr_raw[tr_raw$Expected<200]  
rm(tr_raw0)
gc()

tr_raw$dt <- duration.f(tr_raw)   # self defined duration function
tr_raw$mp <- marshall_palmer(tr_raw$Ref)

#Collapse to one record per Id
tr <- tr_raw[, .(
  target = mean(log1p(Expected), na.rm = T),
  ref = mean(Ref*dt/60, na.rm = T),
  ref1 = mean(RefComposite*dt/60, na.rm = T),
  mp = sum(dt * mp/60, na.rm = T),
  kdp = sum(dt * Kdp/60, na.rm = T),
  zdr = sum(dt * Zdr/60, na.rm = T),
  rho = sum(dt * RhoHV/60, na.rm = T),
  rd = mean(radardist_km, na.rm = T),
  records = .N,
  naCounts = sum(is.na(Ref))
), Id]

write.csv(tr, file = "tr_1M.csv")

# Preparation of Test data #
te_raw$dt <- duration.f(te_raw)
te_raw$mp <- marshall_palmer(te_raw$Ref)
te <- te_raw[, .(
  target = mean(Expected, na.rm = T),  # The original value 
  ref = mean(dt * Ref/60, na.rm = T),
  ref1 = mean(dt * RefComposite/60, na.rm = T),
  mp = sum(dt * mp/60, na.rm = T),
  kdp = sum(dt * Kdp/60, na.rm = T),
  zdr = sum(dt * Zdr/60, na.rm = T),
  rho = sum(dt * RhoHV/60, na.rm = T),
  rd = mean(radardist_km),
  records = .N,
  naCounts = sum(is.na(Ref))
),Id]

write.csv(te, file = "te_1M.csv")

# parameter in xgboost
param0 <- list("objective"  = "reg:linear" 
               #, eval_metric = MAE.eval.f
               #, "eval_metric" = "rmse"
               , "eta" = 0.03            #*optimal
               , "subsample" = 0.7       #*optimal
               , "min_child_weight" =4   #*optimal
               , "max_depth" = 2
               #, "nthreads" = 2
               , "gamma"=0.1
               , "colsample_bytree"=0.55
               , "booster" = "gbtree"
)

# Calculate Training MAE and R2
MAE_outcome = MAE.f(tr,tr,model="xgb")  #
MAE_outcome$"MAE"
MAE_outcome$"MAE sum"
MAE_outcome$"No.Id"
pred.test=MAE_outcome$"predicted test value"

attach(pred.test)
pred.test = pred.test[order(Id),]

plot(tr$target,log1p(pred.test$target),
     xlim=c(0,5),ylim=c(0,5),
     xlab="Observed",
     ylab="Predicted",
     main="Predicted Versus Observed values for Training data",
     col="dark blue")

mean(abs(pred.test$target-expm1(tr$target)))

True.SGB <- tr$target
Pred.SGB <- log1p(pred.test$target)
R2 <- (sum((Pred.SGB-mean(True.SGB)) ** 2))/(sum((True.SGB-mean(True.SGB)) ** 2))

# Calculate Test MAE
# Model based only on distance
MAE_outcome = MAE.f(tr,te,model="xgb")  #
# MAE_outcome
MAE_outcome$"MAE"
MAE_outcome$"MAE sum"
MAE_outcome$"No.Id"
pred.test=MAE_outcome$"predicted test value"

attach(pred.test)
pred.test = pred.test[order(Id),]

plot(te$target,log1p(pred.test$target),
     xlim=c(0,5),ylim=c(0,5),
     xlab="Observed",
     ylab="Predicted",
     main="Predicted Versus Observed values for Test data",
     col="dark red")

mean(abs(pred.test$target-expm1(te$target))) 

True.SGB <- te$target
Pred.SGB <- log1p(pred.test$target)
R2 <- (sum((Pred.SGB-mean(True.SGB)) ** 2))/(sum((True.SGB-mean(True.SGB)) ** 2))