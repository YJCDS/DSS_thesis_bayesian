##################################
### Thesis Y.J.C. Bragt        ###
### Data Science & Society MSc ###
### Tilburg University         ###
### 16 December 2024           ###
##########################################
###                                    ###
### Analysis 1: Stratified predictions ###
### Analysis 2: LOPOCV                 ###
###                                    ###
##########################################

# Define output location
# outputloc <- XXXX

set.seed(1995)
setwd(outputloc)

# Packages
library(dplyr)
library(rstan)         # For Bayesian models
library(StanHeaders)   # For Bayesian models
library(cmdstanr)      # For Bayesian models
library(parallel)      # For Bayesian models
library(brms)          # For Bayesian models
library(glmnet)        # For Ridge regression
library(Metrics)
library(ranger)        # Fastest random forests
library(lme4)          # For multilevel linear regression

# Bayesian model options
options(mc.cores = parallel::detectCores(), brms.backend = "cmdstanr")
rstan_options(auto_write = TRUE)

# Which target vars, models to run?
analysis1  <- T
analysis2  <- T
targetvars <- c("pam","phq4_score")
methods    <- c("linear regression","ridge regression","random forest","multilevel linear regression","multilevel random forest","bayesian multilevel regression")

for(loop_targetvar in targetvars){
  print(paste0("target var = ",loop_targetvar," | Time: ",Sys.time()))
  print(paste0("methods to be evaluated: ",list(methods)))
  print(paste0("Analysis 1: ",ifelse(analysis1,"yes","no"),". Analysis 2: ",ifelse(analysis2,"yes","no"),"."))
  
  # Customize the folliwng vars to tweak the (stratified) CV (analysis 1)
  kfolds          <- 10
  dataset         <- "df_android"
  loop_data       <- get(dataset)
  predictors      <- intersect(names(features_all), names(get(dataset)))
  targetvar       <- loop_targetvar
  bol_reduce_data <- T       # Whole dataset or limited number of rows to use for CV?
  nr_of_rows      <- 2000    # nrows if data is reduced
  bol_scaling     <- T       # Mean-scaling and centering
  bol_stratify    <- T       # Stratify on target var if needed
  stratifying_on  <- c("pid")
  bol_write       <- T       #write results as csv
  
  # My models are:
  # Linear (frequentist) regression with selected features and with all features
  # Ridge regression with selected features and (with all feautures (= rr_full)
  # Random forest with selected features and all features (= rf_full)
  # Multilevel models: LR and RF
  # Bayesian multilevel model with selected features only. Full is computationally to expensive
  
  data_cv_all     <- loop_data %>% drop_na() %>% tibble::rowid_to_column("case_id") %>% # to check later if cases are leaking
    mutate(across(all_of(stratifying_on), as.factor))
  dummys          <- tibble(pid = unique(data_cv_all$pid)) %>% mutate(pid_pivot = pid, one = 1) %>%
    pivot_wider(names_from = pid_pivot, names_prefix = "dummy_pid",values_from = one, values_fill = 0)
  data_cv_all     <- left_join(data_cv_all,dummys, by = "pid")
  data_cv         <- data_cv_all
  formula_full    <- as.formula(paste(targetvar, "~", paste(predictors, collapse = "+")))
  
  
  # Feature selection: fit a full frequentist regression, stepwise selection of vars, keep vals with t > 3
  full_model           <- lm(formula_full, data = data_cv_all)
  stepwise_model       <- step(full_model, direction = "both", trace = 0)
  tvals                <- abs(summary(stepwise_model)$coefficients[-1,"t value"])
  assign(paste0("tvals_",loop_targetvar),tvals)
  predictors_selection <- unique(c(names(tvals[tvals > 3])))
  assign(paste0("predictors_selection_",loop_targetvar),predictors_selection)
  predictors_all       <- predictors
  
  formula_selection    <- as.formula(paste(targetvar, "~", paste(predictors_selection, collapse = "+")))
  formula_multilevel   <- as.formula(paste(targetvar, "~", paste(predictors_selection, collapse = "+"),"+ (1 +",paste(predictors_selection, collapse = "+"),"| pid)"))
  
  names_dummys <- names(select(data_cv, starts_with("dummy")))
  
  rm(interactions)
  for (pred in predictors_selection){
    for(dummy in names_dummys){
      interaction <- as.tibble(data_cv[[pred]] * data_cv[[dummy]])
      names(interaction) <- paste0("interaction",pred,"_x_",dummy)
      if(exists("interactions")){interactions <- bind_cols(interactions, interaction)}else{interactions <- tibble(interaction)}
    }
  }
  names_interactions <- names(interactions)
  interactions$case_id <- data_cv_all$case_id
  
  formula_interactions    <- as.formula(paste(targetvar, "~",
                                              paste(predictors_selection, collapse = "+"), "+",
                                              paste(names_dummys, collapse = "+"), "+",
                                              paste(names_interactions, collapse = "+")))
  
  # My formulas are:
  # (1) full formula: intercept + features
  # (2) selected formula: intercept + selected features
  # (3) Multilevel formula: intercept + selected features + (1 + selected features | pid)
  
  ###################################
  ### Hyperparams for grid-search ###
  ###################################
  hyperparams_rr <- bind_cols(iter = c(1:100), lambda = 10^seq(10, -6, length = 100))
  hyperparams_rr2 <- bind_cols(iter = c(1:100), lambda = 10^seq(4, -9, length = 100))
  
  hyperparams_rf <- expand.grid(mtry_factor = c(0.3,0.4,0.5),ntree = c(500,700,900,1100),minnodesize = c(17,23,30), maxdepth = c(10,20,30)) %>% rowid_to_column("iter")
  hyperparams_rf2 <- expand.grid(mtry_factor = c(0.3,0.4),ntree = c(500,700,900),minnodesize = c(17,23), maxdepth = c(10,20)) %>% rowid_to_column("iter")
  
  hyperparams_rfml <- expand.grid(mtry_factor = c(0.4,0.5,0.6),ntree = c(500,700,900,1100),minnodesize = c(17,23,30), maxdepth = c(10,20,30)) %>% rowid_to_column("iter")
  hyperparams_rfml2 <- expand.grid(mtry_factor = c(0.4,0.5),ntree = c(700,900,1100),minnodesize = c(23,30), maxdepth = c(10)) %>% rowid_to_column("iter")
  
  ######################################
  ### Generating train-val/test data ###
  ######################################
  
  # 1. shuffle (+ reduce) data
  data_cv <- data_cv[sample(1:nrow(data_cv)),]
  if (bol_stratify)   {data_cv <- data_cv %>% group_by_at(stratifying_on) %>% mutate(rn = row_number()) %>% 
    ungroup() %>% arrange(rn) %>% select(!rn) %>%
    select(pid,gender, everything())}
  if (bol_reduce_data){data_cv <- data_cv[c(1:nr_of_rows),]}
  
  # 2. CV schema
  schema <- data_cv %>% select(case_id, all_of(stratifying_on)) %>%
    mutate(fold = ceiling(row_number()/nrow(.)*kfolds)) #might be convoluted but makes sure that no case gets in more than one fold
  
  # 3. Prediction table
  preds <- tibble()
  
  
if(analysis1){  
  #######################################
  ### Analysis 1: K-fold CV of models ###
  #######################################

  for (k in c(1:(kfolds+1))){
    
    #####################################
    ### Selecting train-val/test data ###
    #####################################
    
    loop_id <- paste0("k",k, "_",targetvar)
    
    if (k %in% (1:kfolds)){
      print(paste0("running fold ",k," / ",kfolds," | Loop ID = ",loop_id," | Time: ",Sys.time()))
      train <- data_cv %>% filter(case_id %in% unlist(schema[schema$fold != k,"case_id"]))
      val   <- data_cv %>% filter(case_id %in% unlist(schema[schema$fold == k,"case_id"]))
    }else{
      print(paste0("running test fold | Var = ",targetvar," | Time: ",Sys.time()))
      train <- data_cv
      val <- anti_join(data_cv_all, data_cv, by = "case_id")
    }
    
    if(bol_scaling){
      train <- train %>% mutate(across(all_of(predictors_all), ~ (. - mean(.)) / sd(.) ^ as.logical(sd(.))))
      val <- val %>% mutate(across(all_of(predictors_all), ~ (. - mean(.)) / sd(.) ^ as.logical(sd(.))))
    }
    
    ##################################
    ### Model 1: Linear regression ###
    ##################################
    
    if("linear regression" %in% methods){
      print("linear regression fitting")
      
      model_lr_selection  <- lm(formula_selection, data = train)
      model_lr_full       <- lm(formula_full, data = train)
      
      pred  <- bind_cols(fold = k, model = "lr",      pred = c(predict(model_lr_selection, newdata = val)),actual = val[[targetvar]], iter = 1)
      pred2 <- bind_cols(fold = k, model = "lr_full", pred = c(predict(model_lr_full, newdata = val)),actual = val[[targetvar]], iter = 1)
      preds <- bind_rows(preds, pred, pred2)
    }
    #################################
    ### Model 2: Ridge regression ###
    #################################
    
    if("ridge regression" %in% methods){
      print("ridge regression fitting")

      for(i in c(1:nrow(hyperparams_rr))){
        if (k %in% (1:kfolds)){
          j <- i
        }else{
          i <- unlist(best_models[best_models$model == "rr","iter"])
          j <- unlist(best_models[best_models$model == "rr_full","iter"])
          print(paste0("hyperparam rr selection model (test set): ",list(hyperparams_rr[i,])))  
          print(paste0("hyperparam rr full model (test set): ",list(hyperparams_rr[j,])))
        }
        
        model_rr_sel  <- glmnet(as.matrix(train[predictors_selection]), train[[targetvar]], alpha = 0, lambda = hyperparams_rr[[i,"lambda"]], family = "gaussian")
        model_rr_full <- glmnet(as.matrix(train[predictors_all]), train[[targetvar]], alpha = 0, lambda = hyperparams_rr[[j,"lambda"]], family = "gaussian")
        
        pred  <- bind_cols(fold = k, model = "rr",
                           pred = c(predict(model_rr_sel, s = "lambda.min", newx = as.matrix(val[predictors_selection]))),
                           actual = val[[targetvar]], iter = i)
        pred2 <- bind_cols(fold = k, model = "rr_full",
                           pred = c(predict(model_rr_full, s = "lambda.min", newx = as.matrix(val[predictors_all]))),
                           actual = val[[targetvar]], iter = j)
        preds <- bind_rows(preds, pred, pred2)
        if (k == kfolds+1){
          print("Test set fitted. Grid search skipped")
          break
        }
      }
    }
    
    ##############################
    ### Model 3: Random Forest ###
    ##############################
    
    if("random forest" %in% methods){
      print("random forest fitting")
      t0 <- as.numeric(Sys.time())
      
      for(i in c(1:nrow(hyperparams_rf))){  
        if (k %in% (1:kfolds)){
          j <- i
        }else{
          i <- unlist(best_models[best_models$model == "rf","iter"])
          j <- unlist(best_models[best_models$model == "rf_full","iter"])
          print(paste0("hyperparam rf selection model (test set): ",list(hyperparams_rf[i,])))  
          print(paste0("hyperparam rf full model (test set): ",list(hyperparams_rf[j,])))
        }
        if(i%%10 == 0){
          print(paste0("progress: ",i," / ",nrow(hyperparams_rf),". Elapsed time: ",round(as.numeric(Sys.time())-t0), " seconds. Estimated remaining time: ",
                       round(round(as.numeric(Sys.time())-t0)/i*(nrow(hyperparams_rf)-i)*(2-i/nrow(hyperparams_rf))^1), " seconds"))
        }
        
        model_rf_sel <- ranger(
          formula_selection,
          data = train,
          mtry = floor(length(predictors_selection)*hyperparams_rf[[i,"mtry_factor"]]),
          num.trees = hyperparams_rf[[i,"ntree"]],
          min.node.size = hyperparams_rf[[i,"minnodesize"]],
          max.depth = hyperparams_rf[[i,"maxdepth"]])
        
        model_rf_full <- ranger(
          formula_full,
          data = train,
          mtry = floor(length(predictors_all)*hyperparams_rf[[j,"mtry_factor"]]),
          num.trees = hyperparams_rf[[j,"ntree"]],
          min.node.size = hyperparams_rf[[j,"minnodesize"]],
          max.depth = hyperparams_rf[[j,"maxdepth"]])
        
        pred  <- bind_cols(fold = k, model = "rf",
                           pred = predict(model_rf_sel, data = val)$predictions,
                           actual = val[[targetvar]], iter = i)
        pred2 <- bind_cols(fold = k, model = "rf_full",
                           pred = predict(model_rf_full, data = val)$predictions,
                           actual = val[[targetvar]], iter = j)
        preds <- bind_rows(preds, pred, pred2)
        if (k == kfolds+1){
          print("Test set fitted. Grid search skipped")
          break
        }
      }
    }
    #############################################
    ### Model 4: Multilevel linear regression ###
    #############################################
    
    if("multilevel linear regression" %in% methods){
      print("multilevel linear regression fitting")
      model_lrml <- lmer(formula_multilevel, data = train)
      
      pred  <- bind_cols(fold = k, model = "lr_ml", pred = predict(model_lrml, newdata = val),actual = val[[targetvar]], iter = 1)
      preds <- bind_rows(preds, pred)
    }
    
    #########################################
    ### Model 5: Multilevel Random Forest ###
    #########################################
    
    if("multilevel random forest" %in% methods){
      print("multilevel random forest fitting")
      t0 <- as.numeric(Sys.time())
      
      train_interactions <- train %>% left_join(interactions, by = "case_id")
      val_interactions <- val %>% left_join(interactions, by = "case_id")
      
      for(i in c(1:nrow(hyperparams_rfml))){  
        if (k == kfolds+1){
          i <- unlist(best_models[best_models$model == "rf_ml","iter"])
          print(paste0("hyperparam rf multilevel model (test set): ",list(hyperparams_rfml[i,])))  
        }
        if(i%%5 == 0){
          print(paste0("progress: ",i," / ",nrow(hyperparams_rfml),". Elapsed time: ",round(as.numeric(Sys.time())-t0), " seconds. Estimated remaining time: ",
                       round(round(as.numeric(Sys.time())-t0)/i*(nrow(hyperparams_rfml)-i)*(2-i/nrow(hyperparams_rfml))^1.5), " seconds"))
        }
        
        model_rfml <- ranger(
          formula_interactions,
          data = train_interactions,
          mtry = floor(length(train_interactions)*hyperparams_rfml[[i,"mtry_factor"]]),
          num.trees = hyperparams_rfml[[i,"ntree"]],
          min.node.size = hyperparams_rfml[[i,"minnodesize"]],
          max.depth = hyperparams_rfml[[i,"maxdepth"]])
        
        pred  <- bind_cols(fold = k, model = "rf_ml",
                           pred = predict(model_rfml, data = val_interactions)$predictions,
                           actual = val_interactions[[targetvar]], iter = i)
        preds <- bind_rows(preds, pred)
        if (k == kfolds+1){
          print("Test set fitted. Grid search skipped")
          break
        }
      }
    }
    
    ###############################################
    ### MODEL 6: Bayesian Multilevel regression ###
    ###############################################
    
    if("bayesian multilevel regression" %in% methods){
      print("bayesian model fitting")
      priors <- c(prior("normal(0, 1.5)", class = "b"),
                  prior("lkj(2)", class = "cor"))
      
      setwd("E:/Lesstof/DSS MSc/Thesis/output")
      model_br <- brm(
        formula = formula_multilevel,
        data = train,
        prior = priors,
        family = gaussian(),
        chains = 12,
        cores = 12,
        iter = 4000,
        warmup = 1000,
        seed = 1995,
        silent = 2,
        file = paste0("fits/",loop_id),
        control = list(adapt_delta = 0.98, max_treedepth = 15)
      )
      
      pred  <- bind_cols(fold = k, model = "br", pred = fitted(model_br, newdata = val)[,1],actual = val[[targetvar]], iter = 1)
      val_m <- val[val$gender == "M",]
      pred2  <- bind_cols(fold = k, model = "br_male", pred = fitted(model_br, newdata = val_m)[,1],actual = val_m[[targetvar]], iter = 1)
      val_f <- val[val$gender == "F",]
      pred3  <- bind_cols(fold = k, model = "br_female", pred = fitted(model_br, newdata = val_f)[,1],actual = val_f[[targetvar]], iter = 1)
      preds <- bind_rows(preds, pred,pred2,pred3)
    }
    
    ##########################################
    ### Computing best models for test set ###
    ##########################################
    
    if(k == kfolds){
      print("CV finished. Best model determining for test set...")
      best_models <- preds %>%
        mutate(sqerror = (pred-actual)^2,
               aerror = abs(pred-actual)) %>%
        group_by(model,iter) %>%
        summarize(avg_error = mean(sqerror)) %>%
        group_by(model) %>%
        arrange(avg_error) %>% slice(1) %>%
        select(model, iter) %>% ungroup() %>%
        mutate(keep = T)
    }
    
    ############################################
    ### Saving predictions under unique name ###
    ############################################
    
    assign(paste0("preds_",loop_targetvar,"_",dataset),preds) 
    
    if (bol_write){
      setwd(outputloc)
      write_csv2(preds,paste0("preds_",loop_targetvar,"_",dataset,".csv"))
    }
  }
}else{
  print("skipping Analysis 1: stratified CV")
  Sys.sleep(5)
}
if(analysis2){
  ##########################################
  ### Analysis 2: Leave-on-person-out-CV ###
  ##########################################
  rds_name <- paste0("preds2_",loop_targetvar,"_",dataset,".rds")
   #tibble to save predictions in
  
  data_lopocv <- data_cv_all %>% mutate(across(all_of(predictors_all), ~ (. - mean(.)) / sd(.) ^ as.logical(sd(.))))
  pid <- unique(data_lopocv$pid)
  w <- c(1:7)
  loop_grid2 <- expand_grid(pid,w) %>% rowid_to_column("iter")
  
  # Check if analysis was already (partially) run, so analysis will not be performed twice
  setwd(outputloc)
  if (file.exists(rds_name)){
    print("preds loaded from file")
    preds2 <- readRDS(rds_name)
  }else{
    print("empty preds-tibble created")  
    preds2 <- tibble()
  }
  
  # Starts fitting models for each pid, window
  for (iter in loop_grid2$iter){
    print(paste0("iter = ",loop_grid2$iter[iter]," | pid = ",loop_grid2$pid[iter]," | w = ",loop_grid2$w[iter]))
    w <- loop_grid2[[iter,"w"]]
    loop_pid <- loop_grid2[[iter,"pid"]]
    
    #skip iters that have already been fitted
    if(iter %in% unique(preds2$fold)){
      print("skip. Iter present.")
      next
    }
    data <- data_lopocv[sample(1:nrow(data_lopocv)),] %>%
      group_by(pid) %>%
      filter(row_number() <= 25) %>%
      ungroup()
    data_pid       <- filter(data, pid == loop_pid)
    train          <- suppressMessages(anti_join(data,data_pid))
    train_pid      <- slice_sample(data_pid, n = w)
    train          <- bind_rows(train,train_pid)
    val            <- suppressMessages(anti_join(filter(data_lopocv, pid == loop_pid),train_pid))
    
    ##################################
    ### Model 1: Linear regression ###
    ##################################
    
    if("linear regression" %in% methods){
      print("linear regression fitting")
      
      model_lr_selection  <- lm(formula_selection, data = train)
      model_lr_full       <- lm(formula_full, data = train)
      
      pred  <- bind_cols(fold = iter, model = "lr",      pred = c(predict(model_lr_selection, newdata = val)),actual = val[[targetvar]], iter = 1)
      pred2 <- bind_cols(fold = iter, model = "lr_full", pred = c(predict(model_lr_full, newdata = val)),actual = val[[targetvar]], iter = 1)
      preds <- bind_rows(preds2, pred, pred2)
    }
    #################################
    ### Model 2: Ridge regression ###
    #################################
    
    if("ridge regression" %in% methods){
      print("ridge regression fitting")
      
      for(i in c(1:nrow(hyperparams_rr2))){
        
        model_rr_sel  <- glmnet(as.matrix(train[predictors_selection]), train[[targetvar]], alpha = 0, lambda = hyperparams_rr2[[i,"lambda"]], family = "gaussian")
        model_rr_full <- glmnet(as.matrix(train[predictors_all]), train[[targetvar]], alpha = 0, lambda = hyperparams_rr2[[i,"lambda"]], family = "gaussian")
        
        pred  <- bind_cols(fold = iter, model = "rr",
                           pred = c(predict(model_rr_sel, s = "lambda.min", newx = as.matrix(val[predictors_selection]))),
                           actual = val[[targetvar]], iter = i)
        pred2 <- bind_cols(fold = iter, model = "rr_full",
                           pred = c(predict(model_rr_full, s = "lambda.min", newx = as.matrix(val[predictors_all]))),
                           actual = val[[targetvar]], iter = i)
        preds2 <- bind_rows(preds2, pred, pred2)
      }
    }
    
    ##############################
    ### Model 3: Random Forest ###
    ##############################
    
    if("random forest" %in% methods){
      print("random forest fitting")
      t0 <- as.numeric(Sys.time())
      
      for(i in c(1:nrow(hyperparams_rf2))){  
        if(i%%10 == 0){
          print(paste0("progress: ",i," / ",nrow(hyperparams_rf2),". Elapsed time: ",round(as.numeric(Sys.time())-t0), " seconds. Estimated remaining time: ",
                       round(round(as.numeric(Sys.time())-t0)/i*(nrow(hyperparams_rf2)-i)*(2-i/nrow(hyperparams_rf2))^1), " seconds"))
        }
        
        model_rf_sel <- ranger(
          formula_selection,
          data = train,
          mtry = floor(length(predictors_selection)*hyperparams_rf2[[i,"mtry_factor"]]),
          num.trees = hyperparams_rf2[[i,"ntree"]],
          min.node.size = hyperparams_rf2[[i,"minnodesize"]],
          max.depth = hyperparams_rf2[[i,"maxdepth"]])
        
        model_rf_full <- ranger(
          formula_full,
          data = train,
          mtry = floor(length(predictors_all)*hyperparams_rf2[[i,"mtry_factor"]]),
          num.trees = hyperparams_rf2[[i,"ntree"]],
          min.node.size = hyperparams_rf2[[i,"minnodesize"]],
          max.depth = hyperparams_rf2[[i,"maxdepth"]])
        
        pred  <- bind_cols(fold = iter, model = "rf",
                           pred = predict(model_rf_sel, data = val)$predictions,
                           actual = val[[targetvar]], iter = i)
        pred2 <- bind_cols(fold = iter, model = "rf_full",
                           pred = predict(model_rf_full, data = val)$predictions,
                           actual = val[[targetvar]], iter = i)
        preds2 <- bind_rows(preds2, pred, pred2)
      }
    }
    #############################################
    ### Model 4: Multilevel linear regression ###
    #############################################
    
    if("multilevel linear regression" %in% methods){
      print("multilevel linear regression fitting")
      model_lrml <- lmer(formula_multilevel, data = train)
      
      pred  <- bind_cols(fold = iter, model = "lr_ml", pred = predict(model_lrml, newdata = val),actual = val[[targetvar]], iter = 1)
      preds2 <- bind_rows(preds2, pred)
    }
    
    #########################################
    ### Model 5: Multilevel Random Forest ###
    #########################################
    
    if("multilevel random forest" %in% methods){
      print("multilevel random forest fitting")
      t0 <- as.numeric(Sys.time())
      
      train_interactions <- train %>% left_join(interactions, by = "case_id")
      val_interactions <- val %>% left_join(interactions, by = "case_id")
      
      for(i in c(1:nrow(hyperparams_rfml2))){  
        if(i%%5 == 0){
          print(paste0("progress: ",i," / ",nrow(hyperparams_rfml2),". Elapsed time: ",round(as.numeric(Sys.time())-t0), " seconds. Estimated remaining time: ",
                       round(round(as.numeric(Sys.time())-t0)/i*(nrow(hyperparams_rfml2)-i)*(2-i/nrow(hyperparams_rfml2))^1.5), " seconds"))
        }
        
        model_rfml <- ranger(
          formula_interactions,
          data = train_interactions,
          mtry = floor(length(train_interactions)*hyperparams_rfml2[[i,"mtry_factor"]]),
          num.trees = hyperparams_rfml2[[i,"ntree"]],
          min.node.size = hyperparams_rfml2[[i,"minnodesize"]],
          max.depth = hyperparams_rfml2[[i,"maxdepth"]])
        
        pred  <- bind_cols(fold = iter, model = "rf_ml",
                           pred = predict(model_rfml, data = val_interactions)$predictions,
                           actual = val_interactions[[targetvar]], iter = i)
        preds2 <- bind_rows(preds2, pred)
      }
    }
    
    ###############################################
    ### MODEL 6: Bayesian Multilevel regression ###
    ###############################################
    
    if("bayesian multilevel regression" %in% methods){
      print("bayesian model fitting")
      priors <- c(prior("normal(0, 1.5)", class = "b"),
                  prior("lkj(2)", class = "cor"))
      
      setwd("E:/Lesstof/DSS MSc/Thesis/output")
      model_br <- brm(
        formula = formula_multilevel,
        data = train,
        prior = priors,
        family = gaussian(),
        chains = 8,
        cores = 8,
        iter = 2000,
        warmup = 500,
        seed = 1995,
        silent = 2,
        file = paste0("fits/lopocv",iter,loop_targetvar),
        control = list(adapt_delta = 0.97, max_treedepth = 12)
      )
      
      pred  <- bind_cols(fold = iter, model = "br", pred = fitted(model_br, newdata = val)[,1],actual = val[[targetvar]], iter = 1)
      preds2 <- bind_rows(preds2, pred)
    }
    
    ############################################
    ### Saving predictions under unique name ###
    ############################################
    
    assign(paste0("preds2_",loop_targetvar,"_",dataset),preds) 
    
    if (bol_write & iter%%20 == 0 | bol_write & iter == max(loop_grid2$iter)){
      print(paste0("saving backup of preds2 as iter = ",iter," ..."))
      setwd(outputloc)
      saveRDS(preds2,rds_name)
      saveRDS(preds2,paste0("backup_",format(Sys.time(),"%Y-%m-%d, %H-%M"),"_",rds_name))
    }
  }
}else{
  print("Skipping analysis 2: LOPOCV")}
}
