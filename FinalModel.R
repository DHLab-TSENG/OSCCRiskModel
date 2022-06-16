
stk_output_shiny <- function(autoAb,seed,mfiBinary) {
  set.seed(seed)
  # TrnIndex <- createDataPartition(autoAb$Diagnosis, p = .80, 
  #                                 list = F, times = 1)
  # Trn <- autoAb[TrnIndex,]
  # Tst <- autoAb[-TrnIndex,]
  
  # new split
  AutoAb_split <- initial_split(AutoAb, prop = .80, strata = Diagnosis)
  Trn <- training(AutoAb_split)
  Tst <- testing(AutoAb_split) # split into training and testing sets
  # end- new split
  
  
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
  #TrnOri<-Trn
  #Trn<-rec %>% prep(TrnOri) %>%bake(new_data = Trn) 
  #Tst<-rec %>% prep(TrnOri) %>%bake(new_data = Tst) 
  
  #new_rec <- recipe(Diagnosis ~., data = Trn) 
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
    add_recipe(rec) %>%
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
    add_recipe(rec) %>%
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
    add_recipe(rec) %>%
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
    add_recipe(rec)%>%
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
  
  #initialize the stack here
  AutoAb_stk <- stacks() %>% 
    add_candidates(lr_model) %>% 
    add_candidates(best_rf_model) %>%
    add_candidates(best_svm_model) %>%
    add_candidates(best_xgb_model) #add candidate members
  
  AutoAb_model_stk <- AutoAb_stk %>% # determine how to combine their predictions
    blend_predictions(penalty = 10^(-6:-1)) %>% #tune penalty from 10^(-6:-1) ## here
    fit_members() #fit the candidates with nonzero stacking coefficients
  
  return(AutoAb_model_stk)
}
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
    mutate(variable = fct_reorder(variable, dropout_loss)) %>%
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
source("libAndData.R")
bestSeed<-4388 ##here


#LO
AutoAb$Age<- AutoAb$AgeOri 

AutoAb$ANXA2 <- AutoAb$ANXA2Log
AutoAb$CA2 <- AutoAb$CA2Log
AutoAb$HSPA5 <- AutoAb$HSPA5Log
AutoAb$ISG15 <- AutoAb$ISG15Log
AutoAb$KNG1 <- AutoAb$KNG1Log
AutoAb$MMP1 <- AutoAb$MMP1Log
AutoAb$MMP3 <- AutoAb$MMP3Log
AutoAb$p53 <- AutoAb$p53Log
AutoAb$PRDX2 <- AutoAb$PRDX2Log
AutoAb$SPARC <- AutoAb$SPARCLog #log10
if(file.exists(paste0("stk_model",bestSeed,".rds"))){
  stk_model<-readRDS(paste0("stk_model",bestSeed,".rds"))
}else{
  cpu.core<-parallel::detectCores()
  cl <- makeCluster(cpu.core)
  registerDoParallel(cl)
  stk_model<-stk_output_shiny(AutoAb,bestSeed,FALSE)
  stopCluster(cl)
  stk_model$metrics<-NULL
  stk_model$train<-NULL
  stk_model$splits<-NULL
  stk_model$model_metrics<-NULL 
  stk_model$model_defs$lr_model$pre$actions$recipe$recipe$template<-NULL
  stk_model$model_defs$best_rf_model$pre$actions$recipe$recipe$template<-NULL
  stk_model$model_defs$best_svm_model$pre$actions$recipe$recipe$template<-NULL
  stk_model$model_defs$best_xgb_model$pre$actions$recipe$recipe$template<-NULL
  stk_model$member_fits$lr_model_1_1$pre$actions$recipe$recipe$template<-NULL
  stk_model$member_fits$best_rf_model_1_1$pre$actions$recipe$recipe$template<-NULL
  stk_model$member_fits$best_xgb_model_1_1$pre$actions$recipe$recipe$template<-NULL
  saveRDS(stk_model,paste0("stk_model",bestSeed,".rds"))
}
explainer_stk <- 
  explain_tidymodels(
    stk_model, 
    data = AutoAb %>% dplyr::select(Gender,Age,CIG,ALC,BN ,
                                    ANXA2,CA2,HSPA5 , ISG15 ,KNG1 , 
                                    MMP1,MMP3 , p53 , PRDX2 , SPARC), 
    y = as.numeric(AutoAb$Diagnosis)-1,
    label = "Stacking",
    verbose = FALSE
  )
set.seed(bestSeed)
vip_stk <- model_parts(explainer_stk)
vip_stk$variable<-gsub("^BN","Betel nut chewing",vip_stk$variable)
vip_stk$variable<-gsub("^ALC","Alcohol consumption",vip_stk$variable)
vip_stk$variable<-gsub("^CIG","Smoking",vip_stk$variable)
ggplot_imp(vip_stk)+theme_bw()
ggsave("FinalVarImpLO.pdf",device="pdf",width = 6,height = 6)
ggsave("FinalVarImpLO.png",device="png",width = 6,height = 6)
AutoAb_prob <- AutoAb %>%
  bind_cols(predict(stk_model, ., type = "prob"))
AutoAb_pred <- AutoAb_prob %>%
  bind_cols(predict(stk_model, .))
yardstick::roc_auc(AutoAb_prob, truth = Diagnosis, .pred_H,event_level="second")[[3]]
yardstick::sens(AutoAb_pred, truth = Diagnosis, .pred_class,event_level="second")[[3]]
yardstick::spec(AutoAb_pred, truth = Diagnosis, .pred_class,event_level="second")[[3]]

AutoAb_prob_out<-AutoAb_pred %>% dplyr::select(Diagnosis,.pred_H,.pred_class)

saveRDS(AutoAb_prob_out,"AutoAb_prob_out.rds")
