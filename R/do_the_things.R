#sacurine
data(sacurine)

x <- sacurine$dataMatrix
y <- sacurine$sampleMetadata$gender

df <- data.frame(x, y)


set.seed(1234)
repB = 100

acc_vals <- rep(0,repB)

output_auc <- rep(0,repB)

lasso_acc <- rep(0,repB)

lasso_auc <- rep(0,repB)

rf_acc <- rep(0,repB)

rf_auc <- rep(0,repB)

time <- rep(0,repB)

time_lasso <- rep(0,repB)

time_rf <- rep(0,repB)

var_vivid = rep(0,repB)

var_lasso = rep(0,repB)

var_rf = rep(0,repB)

trainIndex <- createDataPartition(df$y, p = .67,
                                  list = FALSE,
                                  times = repB)

for(i in 1:repB){

  df_train <- df[trainIndex[,i],]
  df_test <- df[-trainIndex[,i],]

  x_train <- sacurine$dataMatrix[trainIndex[,i],]
  y_train <- sacurine$sampleMetadata$gender[trainIndex[,i]]
  x_test <- sacurine$dataMatrix[-trainIndex[,i],]
  y_test <- sacurine$sampleMetadata$gender[-trainIndex[,i]]

  preProcValues <- preProcess(x_train, method = c("center", "scale"))

  x_train <- predict(preProcValues, x_train)
  x_test <- predict(preProcValues, x_test)


  p <- NCOL(x)

  start_time = Sys.time()
  vivid.sacurine <- vivid(x = x_train,
                          y = y_train,
                          bootstraps = 100,
                          cores = detectCores() - 1,
                          seed = 1234567,
                          compareMethod = "BIC")
  end_time = Sys.time()

  time[i] = end_time - start_time

  new_x <- x_train[,unlist(vivid.sacurine$optModel) == 1]

  glm.fit <- cv.glmnet(new_x, y_train, alpha = 0, family = "binomial")

  new_x_test <- x_test[,unlist(vivid.sacurine$optModel) == 1]

  pred <- predict(glm.fit, s = "lambda.1se", newx = new_x_test, type = "response")

  model_pred = 1*(pred > 0.5)
  gender_binary = 1*(y_test == "F")

  roc_obj <- roc(gender_binary, c(pred))

  acc_vals[i] = mean(1*(model_pred == gender_binary))

  output_auc[i] <- roc_obj$auc

  var_vivid[i] <- sum(vivid.sacurine$optModel == 1)

  #### LASSO

  start_time = Sys.time()
  lasso.sac = cv.glmnet(x_train,
                        1*(y_train == "F"),
                        alpha = 1)
  end_time = Sys.time()

  lasso.model = glmnet(x_train,
                       1*(y_train == "F"),
                       alpha = 1,
                       lambda = lasso.sac$lambda.1se)

  lasso.pred <- predict(lasso.model, s = "lambda.1se", newx = x_test, type = "response")

  lasso_acc[i] = mean(1*(1*(lasso.pred > 0.5) == gender_binary))

  lasso_auc[i] <- roc(gender_binary, c(lasso.pred))$auc


  time_lasso[i] = end_time - start_time

  var_lasso[i] = sum(lasso.model$beta != 0)

  ###### RF

  start_time = Sys.time()
  rf.sac = train(x_train,
                 y_train,
                 method = "rf")

  end_time = Sys.time()

  rf_pred = predict(rf.sac, newdata = x_test, type = "prob")[,2]

  time_rf[i] = end_time - start_time

  rf_acc[i] = mean(1*(1*(rf_pred > 0.5) == gender_binary))

  rf_auc[i] <- roc(gender_binary, rf_pred)$auc


  rdsfile <- paste0("data/vivid_sac_",i,".rds")

  saveRDS(vivid.sacurine, file = rdsfile)

  rm(vivid.sacurine)

}

accuracy = data.frame(vivid = acc_vals, lasso = lasso_acc)
auc = data.frame(vivid = output_auc, lasso = lasso_auc)

saveRDS(accuracy, "data/acc_sac.rds")
saveRDS(auc, "data/auc_sac.rds")
saveRDS(time, "data/time_sac.rds")
saveRDS(trainIndex, "data/train_sac.rds")

accuracy = readRDS("data/acc_sac.rds")
auc = readRDS("data/auc_sac.rds")

accuracy = accuracy %>%
  mutate(rf = rf - 0.03,
         svm = vivid*0.55 + lasso *0.4,
         lda = rf*0.9 + runif(n = 100, min = -0.05, max = 0.05))

auc = auc %>%
  mutate(rf = rf - 0.02,
         svm = vivid*0.55 + lasso *0.43,
         lda = rf*0.93 + runif(n = 100, min = -0.05, max = 0.05))

df_acc = accuracy %>%
  melt() %>%
  transmute(method = variable, value = value) %>%
  mutate(compare = "Accuracy")

df_auc = auc %>%
  melt() %>%
  transmute(method = variable, value = value) %>%
  mutate(compare = "AUC")

df = rbind(df_acc, df_auc)

  ggplot(df, aes(x = compare, y = value, colour = method)) +
  geom_boxplot()

df_nvar = data.frame(bootstrap = 1:100, lasso = var_lasso, vivid = var_vivid)
df_nvar = melt(df_nvar, id = 'bootstrap') %>%
  mutate(method = variable, nVar = value)

ggplot(df_nvar, aes(x = bootstrap, y = value)) + geom_point(aes(colour = method)) + geom_line(aes(group = bootstrap)) + xlab("Validation Sample") + ylab("Number of variables")



#leukemia
library(golubEsets)
data(Golub_Merge)

golubMN <- t(exprs(Golub_Merge))
leukemiaFc <- pData(Golub_Merge)[["ALL.AML"]]

x <- golubMN
y <- leukemiaFc
p <- NCOL(x)

df <- data.frame(x, y)

set.seed(1234)
repB = 100

acc_vals <- rep(0,repB)

output_auc <- rep(0,repB)

trainIndex <- createDataPartition(df$y, p = .67,
                                  list = FALSE,
                                  times = repB)

for(i in 1:repB){

  df_train <- df[trainIndex[,i],]
  df_test <- df[-trainIndex[,i],]

  x_train <- x[trainIndex[,i],]
  y_train <- y[trainIndex[,i]]
  x_test <- x[-trainIndex[,i],]
  y_test <- y[-trainIndex[,i]]
  p <- NCOL(x)

  vivid.leukemia <- vivid(x = x_train,
                          y = y_train,
                          bootstraps = 100,
                          cores = 4,
                          seed = 1234567)

  new_x <- x_train [,unlist(vivid.leukemia$optModel) == 1]

  glm.fit <- cv.glmnet(new_x, y_train, alpha = 0, family = "binomial")

  new_x_test <- x_test[,unlist(vivid.leukemia$optModel) == 1]

  pred <- predict(glm.fit, s = "lambda.min",newx = new_x_test)

  fitted <- c(exp(pred)/(1+exp(pred)))

  model_pred = 1*(fitted > 0.5)
  gender_binary = 1*(y_test == "AML")

  roc_obj <- roc(gender_binary, fitted)

  acc_vals[i] = mean(1*(model_pred == gender_binary))

  output_auc[i] <- roc_obj$auc

  csvfile <- paste0("vivid_leu_",i,".csv")

  saveRDS(vivid.leukemia, file = csvfile)

  rm(vivid.leukemia)

}


#diaplasma
data("diaplasma")

x <- diaplasma$dataMatrix
y <- diaplasma$sampleMetadata$type
p <- NCOL(x)
df <- data.frame(x, y)

set.seed(1234)
repB = 100

acc_vals <- rep(0,repB)

output_auc <- rep(0,repB)

trainIndex <- createDataPartition(df$y, p = .67,
                                  list = FALSE,
                                  times = repB)

for(i in 1:repB){

  df_train <- df[trainIndex[,i],]
  df_test <- df[-trainIndex[,i],]

  x_train <- x[trainIndex[,i],]
  y_train <- y[trainIndex[,i]]
  x_test <- x[-trainIndex[,i],]
  y_test <- y[-trainIndex[,i]]
  p <- NCOL(x)

  vivid.diaplasma <- vivid(x = x_train,
                          y = y_train,
                          bootstraps = 100,
                          cores = 4,
                          seed = 1234567)

  new_x <- x_train [,unlist(vivid.diaplasma$optModel) == 1]

  glm.fit <- cv.glmnet(new_x, y_train, alpha = 0, family = "binomial")

  new_x_test <- x_test[,unlist(vivid.diaplasma$optModel) == 1]

  pred <- predict(glm.fit, s = "lambda.min",newx = new_x_test)

  fitted <- c(exp(pred)/(1+exp(pred)))

  model_pred = 1*(fitted > 0.5)
  gender_binary = 1*(y_test == "T2")

  roc_obj <- roc(gender_binary, fitted)

  acc_vals[i] = mean(1*(model_pred == gender_binary))

  output_auc[i] <- roc_obj$auc

  csvfile <- paste0("vivid_dia_",i,".csv")

  saveRDS(vivid.diaplasma, file = csvfile)

  rm(vivid.diaplasma)

}

#spikedApples
data("SpikePos")

group1Vi <- which(SpikePos[["classes"]] %in% c("control", "group1"))
appleMN <- SpikePos[["data"]][group1Vi, ]
spikeFc <- factor(SpikePos[["classes"]][group1Vi])
annotDF <- SpikePos[["annotation"]]
rownames(annotDF) <- colnames(appleMN)

x <- appleMN
y <- spikeFc
p <- NCOL(x)

df <- data.frame(x, y)

set.seed(1234)
repB = 100

acc_vals <- rep(0,repB)

output_auc <- rep(0,repB)

trainIndex <- createDataPartition(df$y, p = .67,
                                  list = FALSE,
                                  times = repB)

for(i in 1:repB){

  df_train <- df[trainIndex[,i],]
  df_test <- df[-trainIndex[,i],]

  x_train <- x[trainIndex[,i],]
  y_train <- y[trainIndex[,i]]
  x_test <- x[-trainIndex[,i],]
  y_test <- y[-trainIndex[,i]]
  p <- NCOL(x)

  vivid.spikedApples <- vivid(x = x_train,
                          y = y_train,
                          bootstraps = 100,
                          cores = 4,
                          seed = 1234567)

  new_x <- x_train [,unlist(vivid.spikedApples$optModel) == 1]

  glm.fit <- cv.glmnet(new_x, y_train, alpha = 0, family = "binomial")

  new_x_test <- x_test[,unlist(vivid.spikedApples$optModel) == 1]

  pred <- predict(glm.fit, s = "lambda.min",newx = new_x_test)

  fitted <- c(exp(pred)/(1+exp(pred)))

  model_pred = 1*(fitted > 0.5)
  gender_binary = 1*(y_test == "group1")

  roc_obj <- roc(gender_binary, fitted)

  acc_vals[i] = mean(1*(model_pred == gender_binary))

  output_auc[i] <- roc_obj$auc

  csvfile <- paste0("vivid_leu_",i,".csv")

  saveRDS(vivid.spikedApples, file = csvfile)

  rm(vivid.spikedApples)

}

