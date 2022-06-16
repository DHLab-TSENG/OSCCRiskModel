library(readxl)
library(nortest)
library(dplyr)
library(MASS)
library(lattice)
library(ggplot2)
library(MLmetrics)
library(randomForest)
library(caret)
library(kernlab)
library(rminer)
library(purrr)
library(iterators)
library(parallel)
library(foreach)
library(doParallel)
library(tidymodels)
library(stacks)
library(DALEXtra)
library(rstatix)
library(tidyr)
library(rstatix)
library(forcats)
Path<-"Output/"
AutoAb <- read_excel("AutoAb.xlsx")

AutoAb$Diagnosis <- 
  factor(ifelse(AutoAb$Diagnosis %in% "Normal", "L",
                ifelse(AutoAb$Diagnosis %in% "Low Risk", "L", "H"))) %>%
  relevel(AutoAb$Diagnosis, ref = "L") # high-risk vs. low-risk

AutoAb$Gender <- factor(AutoAb$Gender)
AutoAb$CIG <- factor(AutoAb$CIG)
AutoAb$ALC <- factor(AutoAb$ALC)
AutoAb$BN <- factor(AutoAb$BN)
AutoAb$AgeOri <- AutoAb$Age
AutoAb$AgeT<- factor(ifelse(AutoAb$Age <= 44, "Y", ifelse(AutoAb$Age >= 65, "O", "M"))) %>%
  relevel(AutoAb$Age, ref = "Y") # Age = ternary
AutoAb$AgeB <- factor(ifelse(AutoAb$Age <= 50, "M", "O")) %>%
  relevel(AutoAb$Age, ref = "M") # Age = binary
MFIBinary<-function(data,train){
  newData<-data
  newData$ANXA2 <- factor(ifelse(newData$ANXA2 >= quantile(train$ANXA2, 0.9), "O", "N"))
  newData$CA2 <- factor(ifelse(newData$CA2 >= quantile(train$CA2, 0.9), "O", "N"))
  newData$HSPA5 <- factor(ifelse(newData$HSPA5 >= quantile(train$HSPA5, 0.9), "O", "N"))
  newData$ISG15 <- factor(ifelse(newData$ISG15 >= quantile(train$ISG15, 0.9), "O", "N"))
  newData$KNG1 <- factor(ifelse(newData$KNG1 >= quantile(train$KNG1, 0.9), "O", "N"))
  newData$MMP1 <- factor(ifelse(newData$MMP1 >= quantile(train$MMP1, 0.9), "O", "N"))
  newData$MMP3 <- factor(ifelse(newData$MMP3 >= quantile(train$MMP3, 0.9), "O", "N"))
  newData$p53 <- factor(ifelse(newData$p53 >= quantile(train$p53, 0.9), "O", "N"))
  newData$PRDX2 <- factor(ifelse(newData$PRDX2 >= quantile(train$PRDX2, 0.9), "O", "N"))
  newData$SPARC <- factor(ifelse(newData$SPARC >= quantile(train$SPARC, 0.9), "O", "N"))
  return (newData)
}

AutoAb$ANXA2Log <- log10(AutoAb$ANXA2)
AutoAb$CA2Log <- log10(AutoAb$CA2)
AutoAb$HSPA5Log <- log10(AutoAb$HSPA5)
AutoAb$ISG15Log <- log10(AutoAb$ISG15)
AutoAb$KNG1Log <- log10(AutoAb$KNG1)
AutoAb$MMP1Log <- log10(AutoAb$MMP1)
AutoAb$MMP3Log <- ifelse(AutoAb$MMP3 == 0, log10((AutoAb$MMP3+1)), log10(AutoAb$MMP3))
AutoAb$p53Log <- log10(AutoAb$p53)
AutoAb$PRDX2Log <- ifelse(AutoAb$PRDX2 == 0, log10((AutoAb$PRDX2+1)), log10(AutoAb$PRDX2))
AutoAb$SPARCLog <- log10(AutoAb$SPARC) #log10

AutoAb$ANXA2Std<- scale(AutoAb$ANXA2)
AutoAb$CA2Std <- scale(AutoAb$CA2)
AutoAb$HSPA5Std <- scale(AutoAb$HSPA5)
AutoAb$ISG15Std <- scale(AutoAb$ISG15)
AutoAb$KNG1Std <- scale(AutoAb$KNG1)
AutoAb$MMP1Std <- scale(AutoAb$MMP1)
AutoAb$MMP3Std <- scale(AutoAb$MMP3)
AutoAb$p53Std <- scale(AutoAb$p53)
AutoAb$PRDX2Std <- scale(AutoAb$PRDX2)
AutoAb$SPARCStd <- scale(AutoAb$SPARC) #scale

# model parameter ------------------------------------------------------------------
Ctrl <- trainControl(method = "repeatedcv", 
                     number = 5, repeats = 5, 
                     classProbs = T,
                     summaryFunction = twoClassSummary,
                     allowParallel = T)
performanceAgg<-function(tst){
  CM <- confusionMatrix(data = tst$raw, 
                        reference = tst$Diagnosis,
                        positive = "H")[[4]]
  SENS <- CM[1]
  SPEC <- CM[2]
  PPV <- CM[3]
  NPV <- CM[4]
  confusionMatrix(data = tst$raw, 
                  reference = tst$Diagnosis,
                  positive = "H",
                  mode = "prec_recall")
  ACC <- postResample(pred = tst$raw, obs = tst$Diagnosis)[1]
  tst$obs <- tst$Diagnosis
  tst$pred <- tst$raw
  tst <- data.frame(tst)
  ROC <- twoClassSummary(tst, lev = levels(tst$obs))[1]
  prSummary(tst, lev = levels(tst$obs))
  return(c(ACC, SENS, SPEC, PPV, NPV, ROC))
}

# model develop and evaluation-xgb ------------------------------------------------------------------

xgb_result <- function(autoAb,seed, mfiBinary) {
  set.seed(seed)
  
  AutoAb_split <- initial_split(autoAb, prop = .80, strata = Diagnosis)
  Trn <- training(AutoAb_split)
  Tst <- testing(AutoAb_split) # split into training and testing sets
  
  #TrnIndex <- createDataPartition(autoAb$Diagnosis, p = .80, 
  #                                list = F, times = 1) # training : test = 80 : 20
  #Trn <- autoAb[TrnIndex,]
  #Tst <- autoAb[-TrnIndex,]
  
  if(mfiBinary==T){
    Tst<- MFIBinary(Tst,Trn)
    Trn<- MFIBinary(Trn,Trn)
  }
  
  set.seed(888) #here
  Trn_xgb <- train(Diagnosis ~ Gender + Age + CIG + ALC + BN +
                     ANXA2 + CA2 + HSPA5 + ISG15 + KNG1 + 
                     MMP1 + MMP3 + p53 + PRDX2 + SPARC,
                   data = Trn,
                   method = "xgbTree", 
                   #trace = F,
                   metric = "ROC",
                   trControl = Ctrl)
  
  bestPar<-Trn_xgb$bestTune
  allPar<-Trn_xgb$results%>% 
    dplyr::select(colnames(Trn_xgb$bestTune)) %>% unique()
  
  Tst$raw <- predict(Trn_xgb, newdata = Tst, type = "raw")
  Tst.prob <- predict(Trn_xgb, newdata = Tst, type = "prob")
  Tst$H <- Tst.prob$H
  Tst$L <- Tst.prob$L
  perf<-performanceAgg(Tst)
  
  # Imp
  VI <- varImp(Trn_xgb)[1]
  
  return(list(bestPar=bestPar,allPar=allPar,perf=perf,Tst=Tst,VI=VI))
}
xgb_result_ponly <- function(autoAb,seed, mfiBinary) {
  set.seed(seed)
  AutoAb_split <- initial_split(autoAb, prop = .80, strata = Diagnosis)
  Trn <- training(AutoAb_split)
  Tst <- testing(AutoAb_split) # split into training and testing sets
  
  #TrnIndex <- createDataPartition(autoAb$Diagnosis, p = .80, 
  #                                list = F, times = 1) # training : test = 80 : 20
  #Trn <- autoAb[TrnIndex,]
  #Tst <- autoAb[-TrnIndex,]
  
  if(mfiBinary==T){
    Tst<- MFIBinary(Tst,Trn)
    Trn<- MFIBinary(Trn,Trn)
  }
  
  set.seed(888) #here
  Trn_xgb <- train(Diagnosis ~ Gender + Age + CIG + ALC + BN ,
                   data = Trn,
                   method = "xgbTree", 
                   #trace = F,
                   metric = "ROC",
                   trControl = Ctrl)
  
  bestPar<-Trn_xgb$bestTune
  allPar<-Trn_xgb$results%>% 
    dplyr::select(colnames(Trn_xgb$bestTune)) %>% unique()
  
  Tst$raw <- predict(Trn_xgb, newdata = Tst, type = "raw")
  Tst.prob <- predict(Trn_xgb, newdata = Tst, type = "prob")
  Tst$H <- Tst.prob$H
  Tst$L <- Tst.prob$L
  perf<-performanceAgg(Tst)
  
  # Imp
  VI <- varImp(Trn_xgb)[1]
  
  return(list(bestPar=bestPar,allPar=allPar,perf=perf,Tst=Tst,VI=VI))
}
# model develop and evaluation-xgb ------------------------------------------------------------------

RF_result <- function(autoAb,seed,mfiBinary) {
  set.seed(seed)
  AutoAb_split <- initial_split(autoAb, prop = .80, strata = Diagnosis)
  Trn <- training(AutoAb_split)
  Tst <- testing(AutoAb_split) # split into training and testing sets
  
  #TrnIndex <- createDataPartition(autoAb$Diagnosis, p = .80, 
  #                                list = F, times = 1) # training : test = 80 : 20
  #Trn <- autoAb[TrnIndex,]
  #Tst <- autoAb[-TrnIndex,]
  if(mfiBinary==T){
    Tst<- MFIBinary(Tst,Trn)
    Trn<- MFIBinary(Trn,Trn)
  }
  set.seed(888)
  Trn_rf <- train(Diagnosis ~ Gender + Age + CIG + ALC + BN +
                    ANXA2 + CA2 + HSPA5 + ISG15 + KNG1 + 
                    MMP1 + MMP3 + p53 + PRDX2 + SPARC,
                  data = Trn,
                  method = "rf", 
                  trace = F,
                  metric = "ROC",
                  trControl = Ctrl)
  bestPar<-Trn_rf$bestTune
  allPar<-Trn_rf$results%>% 
    dplyr::select(colnames(Trn_rf$bestTune)) %>% unique()
  
  Tst$raw <- predict(Trn_rf, newdata = Tst, type = "raw")
  Tst.prob <- predict(Trn_rf, newdata = Tst, type = "prob")
  Tst$H <- Tst.prob$H
  Tst$L <- Tst.prob$L
  
  perf<-performanceAgg(Tst)
  return(list(bestPar=bestPar,allPar=allPar,perf=perf,Tst=Tst))
}
svm_result <- function(autoAb,seed,mfiBinary) {
  set.seed(seed)
  AutoAb_split <- initial_split(autoAb, prop = .80, strata = Diagnosis)
  Trn <- training(AutoAb_split)
  Tst <- testing(AutoAb_split) # split into training and testing sets
  
  #TrnIndex <- createDataPartition(autoAb$Diagnosis, p = .80, 
  #                                list = F, times = 1) # training : test = 80 : 20
  #Trn <- autoAb[TrnIndex,]
  #Tst <- autoAb[-TrnIndex,]
  if(mfiBinary==T){
    Tst<- MFIBinary(Tst,Trn)
    Trn<- MFIBinary(Trn,Trn)
  }
  
  set.seed(888) #here
  Trn_svm <- train(Diagnosis ~ Gender + Age + CIG + ALC + BN +
                     ANXA2 + CA2 + HSPA5 + ISG15 + KNG1 + 
                     MMP1 + MMP3 + p53 + PRDX2 + SPARC,
                   data = Trn,
                   method = "svmRadial", 
                   trace = F,
                   metric = "ROC",
                   trControl = Ctrl)
  bestPar<-Trn_svm$bestTune
  allPar<-Trn_svm$results%>% 
    dplyr::select(colnames(Trn_svm$bestTune)) %>% unique()
  Tst$raw <- predict(Trn_svm, newdata = Tst, type = "raw")
  Tst.prob <- predict(Trn_svm, newdata = Tst, type = "prob")
  Tst$H <- Tst.prob$H
  Tst$L <- Tst.prob$L
  perf<-performanceAgg(Tst)
  return(list(bestPar=bestPar,allPar=allPar,perf=perf,Tst=Tst))
}

LR_result <- function(autoAb,seed,mfiBinary) {
  set.seed(seed)
  AutoAb_split <- initial_split(autoAb, prop = .80, strata = Diagnosis)
  Trn <- training(AutoAb_split)
  Tst <- testing(AutoAb_split) # split into training and testing sets
  
  #TrnIndex <- createDataPartition(autoAb$Diagnosis, p = .80, 
  #                                list = F, times = 1) # training : test = 80 : 20
  #Trn <- autoAb[TrnIndex,]
  #Tst <- autoAb[-TrnIndex,]
  
  if(mfiBinary==T){
    Tst<- MFIBinary(Tst,Trn)
    Trn<- MFIBinary(Trn,Trn)
  }
  set.seed(888)
  Trn_lr <- train(Diagnosis ~ Gender + Age + CIG + ALC + BN +
                    ANXA2 + CA2 + HSPA5 + ISG15 + KNG1 + 
                    MMP1 + MMP3 + p53 + PRDX2 + SPARC,
                  data = Trn,
                  method = "glmStepAIC",
                  family = binomial("logit"),
                  direction = "both",
                  metric = "ROC",
                  trace = F,
                  trControl = Ctrl) # training by glmStepAIC
  bestPar<-Trn_lr$bestTune
  allPar<-Trn_lr$results%>% 
    dplyr::select(colnames(Trn_lr$bestTune)) %>% unique()
  Tst$raw <- predict(Trn_lr, newdata = Tst, type = "raw")
  Tst.prob <- predict(Trn_lr, newdata = Tst, type = "prob")
  Tst$H <- Tst.prob$H
  Tst$L <- Tst.prob$L
  perf<-performanceAgg(Tst)
  return(list(bestPar=bestPar,allPar=allPar,perf=perf,Tst=Tst))
}


stk_result <- function(autoAb,seed,mfiBinary) {
  set.seed(seed)
  AutoAb_split <- initial_split(autoAb, prop = .80, strata = Diagnosis)
  Trn <- training(AutoAb_split)
  Tst <- testing(AutoAb_split) # split into training and testing sets
  
  #TrnIndex <- createDataPartition(autoAb$Diagnosis, p = .80, 
  #                                list = F, times = 1) # training : test = 80 : 20
  #Trn <- autoAb[TrnIndex,]
  #Tst <- autoAb[-TrnIndex,]
  
  if(mfiBinary==T){
    Tst<- MFIBinary(Tst,Trn)
    Trn<- MFIBinary(Trn,Trn)
  }
  
  set.seed(444)
  
  rec <- recipe(Diagnosis ~Gender + Age + CIG + ALC + BN +
                  ANXA2 + CA2 + HSPA5 + ISG15 + KNG1 + 
                  MMP1 + MMP3 + p53 + PRDX2 + SPARC, data = Trn) %>%
    step_dummy(all_nominal(), -all_outcomes()) %>% 
    step_zv(all_predictors()) %>% 
    step_normalize(all_predictors()) #set up a basic recipe
  #20220613 start
  TrnOri<-Trn
  Trn<-rec %>% prep(TrnOri) %>%bake(new_data = Trn) 
  Tst<-rec %>% prep(TrnOri) %>%bake(new_data = Tst) 
  
  new_rec <- recipe(Diagnosis ~., data = Trn) 
  #20220613 end
  folds <- rsample::vfold_cv(Trn, v = 5, repeats = 5, strata = Diagnosis) #use 5-fold cv repeated 5 times
  metric <- metric_set(accuracy, sens, spec, ppv, npv, roc_auc) #define the metrics
  ctrl_grid <- control_stack_grid() ##tuning and fitting results
  ctrl_res <- control_stack_resamples() ##tuning and fitting results
  
  # logistic regression model
  lr_spec <- logistic_reg() %>%
    set_engine("glm") %>%
    set_mode("classification") # create a logistic regression model definition
  
  lr_wflow <- workflow() %>%
    add_recipe(new_rec) %>%
    add_model(lr_spec) # add to a workflow
  
  set.seed(888)
  lr_model <- fit_resamples(lr_wflow,
                            resamples = folds,
                            metrics = metric,
                            control = ctrl_res) # fit to the 5-fold cv repeated 5 times
  
  # random forest model
  rf_spec <- rand_forest(mtry = tune(),
                         trees = 1000,
                         min_n = tune()) %>% #tune parameters
    set_engine("randomForest") %>%
    set_mode("classification") # create a random forest model definition
  
  rf_wflow <- workflow() %>%
    add_recipe(new_rec) %>%
    add_model(rf_spec) # add to a workflow
  
  set.seed(666)
  rf_res <- tune_grid(rf_wflow,
                      resamples = folds,
                      grid = 50,#50
                      metrics = metric,
                      control = ctrl_grid) # fit to the 5-fold cv repeated 5 times
  
  best_rf <- rf_res %>%
    select_best("roc_auc") #select the best parameter
  
  rf_bestwf <- rf_wflow %>%
    finalize_workflow(best_rf) #add the best parameter to a workflow
  
  set.seed(888)
  best_rf_model <- fit_resamples(rf_bestwf,
                                 resamples = folds,
                                 metrics = metric,
                                 control = ctrl_res) # fit the best random forest model to the 5-fold cv repeated 5 times
  
  # svm with rbf model
  svm_spec <- svm_rbf(cost = tune(), 
                      rbf_sigma = tune()) %>% #tune parameters
    set_engine("kernlab") %>%
    set_mode("classification") # create a svm with rbf model definition
  
  svm_wflow <- workflow() %>%
    add_recipe(new_rec) %>%
    add_model(svm_spec) # add to a workflow
  
  set.seed(666)
  svm_res <- tune_grid(svm_wflow,
                       resamples = folds,
                       grid = 50,
                       metrics = metric,
                       control = ctrl_grid) #fit to the 5-fold cv repeated 5 times
  
  best_svm <- svm_res %>%
    select_best("roc_auc") #select the best parameter
  
  svm_bestwf <- svm_wflow %>%
    finalize_workflow(best_svm) #add the best parameter to a workflow
  
  set.seed(888)
  best_svm_model <- fit_resamples(svm_bestwf,
                                  resamples = folds,
                                  metrics = metric,
                                  control = ctrl_res) #fit the best svm model to the 5-fold cv repeated 5 times
  
  #XGBoost model
  xgb_spec <- boost_tree(mtry = tune(), 
                         trees = 1000,
                         min_n = tune(),
                         tree_depth = tune(),
                         learn_rate = tune(),
                         loss_reduction = tune(),
                         sample_size = tune()) %>% #tune parameters
    set_engine("xgboost") %>%
    set_mode("classification") #create a xgboost model definition
  
  xgb_wflow <- workflow() %>% 
    add_recipe(new_rec) %>%
    add_model(xgb_spec) #add to a workflow
  
  set.seed(666)
  xgb_res <- tune_grid(xgb_wflow,
                       resamples = folds,
                       grid = 50,
                       metrics = metric,
                       control = ctrl_grid) #fit to the 5-fold cv repeated 5 times
  
  best_xgb <- xgb_res %>%
    select_best("roc_auc") #select the best parameter
  
  xgb_bestwf <- xgb_wflow %>%
    finalize_workflow(best_xgb) #add the best parameter to a workflow
  
  set.seed(888)
  best_xgb_model <- fit_resamples(xgb_bestwf,
                                  resamples = folds,
                                  metrics = metric,
                                  control = ctrl_res) #fit the best xgboost model to the 5-fold cv repeated 5 times
  
  #initialize the stack
  AutoAb_stk <- stacks() %>% 
    add_candidates(lr_model) %>% 
    add_candidates(best_rf_model) %>%
    add_candidates(best_svm_model) %>%
    add_candidates(best_xgb_model) #add candidate members
  
  AutoAb_model_stk <- AutoAb_stk %>% # determine how to combine their predictions
    blend_predictions(penalty = 10^(-6:-1)) %>% #tune penalty from 10^(-6:-1) ## here
    fit_members() #fit the candidates with nonzero stacking coefficients
  
  bestPar<-data.frame(penalty=AutoAb_model_stk$penalty$penalty) #?
  allPar<-data.frame(penalty = 10^(-6:-1))
  
  #predict testing set
  Tst_pred <- Tst %>%
    bind_cols(predict(AutoAb_model_stk, .))
  Tst_prob <- Tst %>%
    bind_cols(predict(AutoAb_model_stk, ., type = "prob"))
  
 
  
  explainer_stk <- 
    explain_tidymodels(
      AutoAb_model_stk, 
      data =  Trn %>% select(-Diagnosis) , 
      y = as.numeric(Trn$Diagnosis),
      label = "Stacking",
      verbose = FALSE
    )
  set.seed(1804)
  vip_stk <- model_parts(explainer_stk)
  
  #extract the output
  ACC <- yardstick::accuracy(Tst_pred, truth = Diagnosis, .pred_class)[[3]]
  SENS <- yardstick::sens(Tst_pred, truth = Diagnosis, .pred_class,event_level="second")[[3]]
  SPEC <- yardstick::spec(Tst_pred, truth = Diagnosis, .pred_class,event_level="second")[[3]]
  PPV <- yardstick::ppv(Tst_pred, truth = Diagnosis, .pred_class,event_level="second")[[3]]
  NPV <- yardstick::npv(Tst_pred, truth = Diagnosis, .pred_class,event_level="second")[[3]]
  ROC <- yardstick::roc_auc(Tst_prob, truth = Diagnosis, .pred_H,event_level="second")[[3]]
  perf<-c(Accuracy=ACC, Sensitivity=SENS, Specificity=SPEC, 
          `Pos Pred Value`=PPV, `Neg Pred Value`=NPV, ROC=ROC)
  
  return(list(bestPar=bestPar,allPar=allPar,perf=perf,
              Tst=Tst_prob[,c("Diagnosis",".pred_L",".pred_H")],
              VI=vip_stk))
  #return(c(ACC, SENS, SPEC, PPV, NPV, ROC))
}

stk_result_ponly <- function(autoAb,seed,mfiBinary) {
  set.seed(seed)
  AutoAb_split <- initial_split(autoAb, prop = .80, strata = Diagnosis)
  Trn <- training(AutoAb_split)
  Tst <- testing(AutoAb_split) # split into training and testing sets
  
  #TrnIndex <- createDataPartition(autoAb$Diagnosis, p = .80, 
  #                                list = F, times = 1) # training : test = 80 : 20
  #Trn <- autoAb[TrnIndex,]
  #Tst <- autoAb[-TrnIndex,]
  
  if(mfiBinary==T){
    Tst<- MFIBinary(Tst,Trn)
    Trn<- MFIBinary(Trn,Trn)
  }
  
  set.seed(444)
  
  rec <- recipe(Diagnosis ~Gender + Age + CIG + ALC + BN , data = Trn) %>%
    step_dummy(all_nominal(), -all_outcomes()) %>% 
    step_zv(all_predictors()) %>% 
    step_normalize(all_predictors()) #set up a basic recipe
  #20220613 start
  TrnOri<-Trn
  Trn<-rec %>% prep(TrnOri) %>%bake(new_data = Trn) 
  Tst<-rec %>% prep(TrnOri) %>%bake(new_data = Tst) 
  new_rec <- recipe(Diagnosis ~., data = Trn) 
  #20220613 end
  folds <- rsample::vfold_cv(Trn, v = 5, repeats = 5, strata = Diagnosis) #use 5-fold cv repeated 5 times
  metric <- metric_set(accuracy, sens, spec, ppv, npv, roc_auc) #define the metrics
  ctrl_grid <- control_stack_grid() ##tuning and fitting results
  ctrl_res <- control_stack_resamples() ##tuning and fitting results
  
  # logistic regression model
  lr_spec <- logistic_reg() %>%
    set_engine("glm") %>%
    set_mode("classification") # create a logistic regression model definition
  
  lr_wflow <- workflow() %>%
    add_recipe(new_rec) %>%
    add_model(lr_spec) # add to a workflow
  
  set.seed(888)
  lr_model <- fit_resamples(lr_wflow,
                            resamples = folds,
                            metrics = metric,
                            control = ctrl_res) # fit to the 5-fold cv repeated 5 times
  
  # random forest model
  rf_spec <- rand_forest(mtry = tune(),
                         trees = 1000,
                         min_n = tune()) %>% #tune parameters
    set_engine("randomForest") %>%
    set_mode("classification") # create a random forest model definition
  
  rf_wflow <- workflow() %>%
    add_recipe(new_rec) %>%
    add_model(rf_spec) # add to a workflow
  
  set.seed(666)
  rf_res <- tune_grid(rf_wflow,
                      resamples = folds,
                      grid = 50,#50
                      metrics = metric,
                      control = ctrl_grid) # fit to the 5-fold cv repeated 5 times
  
  best_rf <- rf_res %>%
    select_best("roc_auc") #select the best parameter
  
  rf_bestwf <- rf_wflow %>%
    finalize_workflow(best_rf) #add the best parameter to a workflow
  
  set.seed(888)
  best_rf_model <- fit_resamples(rf_bestwf,
                                 resamples = folds,
                                 metrics = metric,
                                 control = ctrl_res) # fit the best random forest model to the 5-fold cv repeated 5 times
  
  # svm with rbf model
  svm_spec <- svm_rbf(cost = tune(), 
                      rbf_sigma = tune()) %>% #tune parameters
    set_engine("kernlab") %>%
    set_mode("classification") # create a svm with rbf model definition
  
  svm_wflow <- workflow() %>%
    add_recipe(new_rec) %>%
    add_model(svm_spec) # add to a workflow
  
  set.seed(666)
  svm_res <- tune_grid(svm_wflow,
                       resamples = folds,
                       grid = 50,
                       metrics = metric,
                       control = ctrl_grid) #fit to the 5-fold cv repeated 5 times
  
  best_svm <- svm_res %>%
    select_best("roc_auc") #select the best parameter
  
  svm_bestwf <- svm_wflow %>%
    finalize_workflow(best_svm) #add the best parameter to a workflow
  
  set.seed(888)
  best_svm_model <- fit_resamples(svm_bestwf,
                                  resamples = folds,
                                  metrics = metric,
                                  control = ctrl_res) #fit the best svm model to the 5-fold cv repeated 5 times
  
  #XGBoost model
  xgb_spec <- boost_tree(mtry = tune(), 
                         trees = 1000,
                         min_n = tune(),
                         tree_depth = tune(),
                         learn_rate = tune(),
                         loss_reduction = tune(),
                         sample_size = tune()) %>% #tune parameters
    set_engine("xgboost") %>%
    set_mode("classification") #create a xgboost model definition
  
  xgb_wflow <- workflow() %>% 
    add_recipe(new_rec) %>%
    add_model(xgb_spec) #add to a workflow
  
  set.seed(666)
  xgb_res <- tune_grid(xgb_wflow,
                       resamples = folds,
                       grid = 50,
                       metrics = metric,
                       control = ctrl_grid) #fit to the 5-fold cv repeated 5 times
  
  best_xgb <- xgb_res %>%
    select_best("roc_auc") #select the best parameter
  
  xgb_bestwf <- xgb_wflow %>%
    finalize_workflow(best_xgb) #add the best parameter to a workflow
  
  set.seed(888)
  best_xgb_model <- fit_resamples(xgb_bestwf,
                                  resamples = folds,
                                  metrics = metric,
                                  control = ctrl_res) #fit the best xgboost model to the 5-fold cv repeated 5 times
  
  #initialize the stack
  AutoAb_stk <- stacks() %>% 
    add_candidates(lr_model) %>% 
    add_candidates(best_rf_model) %>%
    add_candidates(best_svm_model) %>%
    add_candidates(best_xgb_model) #add candidate members
  
  AutoAb_model_stk <- AutoAb_stk %>% # determine how to combine their predictions
    blend_predictions(penalty = 10^(-6:-1)) %>% #tune penalty from 10^(-6:-1) ## here
    fit_members() #fit the candidates with nonzero stacking coefficients
  
  bestPar<-data.frame(penalty=AutoAb_model_stk$penalty$penalty) #?
  allPar<-data.frame(penalty = 10^(-6:-1))
  
  #predict testing set
  Tst_pred <- Tst %>%
    bind_cols(predict(AutoAb_model_stk, .))
  Tst_prob <- Tst %>%
    bind_cols(predict(AutoAb_model_stk, ., type = "prob"))
  
  
  
  explainer_stk <- 
    explain_tidymodels(
      AutoAb_model_stk, 
      data =  Trn%>% select(-Diagnosis) , 
      y = as.numeric(Trn$Diagnosis),
      label = "Stacking",
      verbose = FALSE
    )
  set.seed(1804)
  vip_stk <- model_parts(explainer_stk)
  
  #extract the output
  
  ACC <- yardstick::accuracy(Tst_pred, truth = Diagnosis, .pred_class)[[3]]
  SENS <- yardstick::sens(Tst_pred, truth = Diagnosis, .pred_class,event_level="second")[[3]]
  SPEC <- yardstick::spec(Tst_pred, truth = Diagnosis, .pred_class,event_level="second")[[3]]
  PPV <- yardstick::ppv(Tst_pred, truth = Diagnosis, .pred_class,event_level="second")[[3]]
  NPV <- yardstick::npv(Tst_pred, truth = Diagnosis, .pred_class,event_level="second")[[3]]
  ROC <- yardstick::roc_auc(Tst_prob, truth = Diagnosis, .pred_H,event_level="second")[[3]]
  perf<-c(Accuracy=ACC, Sensitivity=SENS, Specificity=SPEC, 
          `Pos Pred Value`=PPV, `Neg Pred Value`=NPV, ROC=ROC)
  
  return(list(bestPar=bestPar,allPar=allPar,perf=perf,
              Tst=Tst_prob[,c("Diagnosis",".pred_L",".pred_H")],
              VI=vip_stk))
}

resultOutput<-function(modelObj,model,class){
  model_results<-map_dfr(modelObj,3)
  model_bestPar<-map_dfr(modelObj,1)
  model_allPar<-map_dfr(modelObj,2)
  model_testRaw<-map_dfr(modelObj,4)
  colnames(model_results) <- c("Accuracy", "Sensitivity", "Specificity", "PPV", "NPV", "ROC")
  model_results$Class<-class
  model_results$Algorithm<-model
  model_results$Seed<-SeedList
  write.csv(model_results, file = paste0(Path,model,"5f_",class,"_new.csv"), row.names = F)
  model_bestPar$Class<-class
  model_bestPar$Algorithm<-model
  write.csv(model_bestPar, file = paste0(Path,model,"bestPar_",class,"_new.csv"), row.names = F)
  model_allPar$Class<-class
  model_allPar$Algorithm<-model
  write.csv(model_allPar, file = paste0(Path,model,"allPar_",class,"_new.csv"), row.names = F)
  model_testRaw$Class<-class
  model_testRaw$Algorithm<-model
  write.csv(model_testRaw, file = paste0(Path,model,"testRaw_",class,"_new.csv"), row.names = F)
}
se <- function(x) sqrt(var(x) / length(x))

ggplot_imp <- function(...) {
  obj <- list(...)
  metric_name <- attr(obj[[1]], "loss_name")
  metric_lab <- paste(metric_name, 
                      "after permutations\n(higher indicates more important)")
  
  full_vip <- bind_rows(obj) %>%
    filter(variable != "_baseline_")
  
  perm_vals <- full_vip %>% 
    filter(variable == "_full_model_") %>% 
    group_by(label) %>% 
    summarise(dropout_loss = mean(dropout_loss))
  
  p <- full_vip %>%
    filter(variable != "_full_model_") %>% 
    mutate(variable = forcats::fct_reorder(variable, dropout_loss)) %>%
    ggplot(aes(dropout_loss, variable)) 
  if(length(obj) > 1) {
    p <- p + 
      facet_wrap(vars(label)) +
      geom_vline(data = perm_vals, aes(xintercept = dropout_loss, color = label),
                 size = 1.4, lty = 2, alpha = 0.7) +
      geom_boxplot(aes(color = label, fill = label), alpha = 0.2)
  } else {
    p <- p + 
      geom_vline(data = perm_vals, aes(xintercept = dropout_loss),
                 size = 1.4, lty = 2, alpha = 0.7) +
      geom_boxplot(fill = "#91CBD765", alpha = 0.4)
    
  }
  p +
    theme(legend.position = "none") +
    labs(x = metric_lab, 
         y = NULL,  fill = NULL,  color = NULL)
}
