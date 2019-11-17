###### Evaluation Funciton


eval_data = function(x, y, seed, prop, repeats, folds = 5){
  
  trainIndex = createDataPartition(y,
                                   p = .67,
                                   list = FALSE,
                                   times = repeats)
  
  for(b in 1:repeats){
    x_train = x[trainIndex[,i],]
    y_train = x[trainIndex[,i]]
    
    x_test = x[-trainIndex[,i],]
    y_test = x[-trainIndex[,i]]
    
    preProcValues = preProcess(x_train, method = c("center", "scale"))
    
    x_train = predict(preProcValues, x_train)
    x_test = predict(preProcValues, x_test)
    
    # Vivid
    
    vivid_start = Sys.time()
    vivid_fit = vivid(x = x_train,
                      y = y_train,
                      bootstraps = 100,
                      cores = detectCores() - 1,
                      seed = 1234567,
                      compareMethod = "BIC")
    end_time = Sys.time()
    
    new_x = x_train[,unlist(vivid_fit$optModel) == 1]
    
    glm_fit = cv.glmnet(new_x, y_train, alpha = 0, family = "binomial")
    
    new_x_test = x_test[,unlist(vivid_fit$optModel) == 1]
    
    vivid_pred = predict(glm.fit, s = "lambda.1se", newx = new_x_test, type = "response")
    
    model_pred = 1*(pred > 0.5)
    gender_binary = 1*(y_test == "F")
    
    # Random Forest
    
    rf_start = Sys.time()
    
    rf_fit = train(x_train,
                   y_train,
                   method = "rf")
    
    rf_finish = Sys.time()
    
    
    
    # LDA
    
    lda_start = Sys.time()
    
    lda_finish = Sys.time()
    
    # Boosting
    
    
    
  }
  
}
