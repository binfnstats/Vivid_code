library(RankProd)
library(BioMark)
library(parallel)
library(tidyverse)
library(latticeExtra)
library(biosigner)
library(tictoc)
library(pROC)

#diaplasma

data("diaplasma")

x <- diaplasma$dataMatrix
y <- diaplasma$sampleMetadata$type
p <- NCOL(x)

start.diaplasma <- Sys.time()

vivid.diaplasma <- vivid(x = x,
                        y = y,
                        bootstraps = 100,
                        cores = 4,
                        seed = 1234567)

end.diaplasma <- Sys.time()

saveRDS(vivid.diaplasma, file = "vivid_diaplasma.rds")
remove(vivid.diaplasma)

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

start.SpikePos <- Sys.time()

vivid.SpikePos <- vivid(x = x,
                         y = y,
                         bootstraps = 100,
                         variables = 1:p,
                         cores = detectCores() - 1,
                         seed = 123,
                         min_size = 2,
                         lambda = 'lambda.1se',
                         compare_method = 'EBIC',
                         gamma = 1)

end.SpikePos <- Sys.time()

saveRDS(vivid.SpikePos, file = "D://vivid_SpikePos.rds")
remove(vivid.SpikePos)


#leukemia
library(golubEsets)
data(Golub_Merge)

golubMN <- t(exprs(Golub_Merge))
leukemiaFc <- pData(Golub_Merge)[["ALL.AML"]]

x <- golubMN
y <- leukemiaFc
p <- NCOL(x)

start.leukemia <- Sys.time()

vivid.leukemia <- vivid(x = x,
                        y = y,
                        bootstraps = 100,
                        variables = 1:p,
                        cores = detectCores() - 1,
                        seed = 123,
                        min_size = 2,
                        lambda = 'lambda.1se',
                        compare_method = 'EBIC',
                        gamma = 1)

end.leukemia <- Sys.time()

saveRDS(vivid.leukemia, file = "D://vivid_leukemia.rds")
remove(vivid.leukemia)

#sacurine
data(sacurine)

x <- sacurine$dataMatrix
y <- sacurine$sampleMetadata$gender

df <- data.frame(x, y)


#set.seed(1234)
repB = 100

acc_vals <- rep(0,repB)

output_auc <- rep(0,repB)

trainIndex <- createDataPartition(df$y, p = .67,
                                  list = FALSE,
                                  times = 5)

for(i in 1:5){

df_train <- df[trainIndex[,i],]
df_test <- df[-trainIndex[,i],]

x_train <- sacurine$dataMatrix[trainIndex[,i],]
y_train <- sacurine$sampleMetadata$gender[trainIndex[,i]]
x_test <- sacurine$dataMatrix[-trainIndex[,i],]
y_test <- sacurine$sampleMetadata$gender[-trainIndex[,i]]
p <- NCOL(x)

vivid.sacurine <- vivid(x = x_train,
                        y = y_train,
                        bootstraps = 100,
                        cores = 4,
                        seed = 1234567)

new_x <- x_train [,unlist(vivid.sacurine$optModel) == 1]

glm.fit <- cv.glmnet(new_x, y_train, alpha = 0, family = "binomial")

new_x_test <- x_test[,unlist(vivid.sacurine$optModel) == 1]

pred <- predict(glm.fit, s = "lambda.min",newx = new_x_test)

fitted <- c(exp(pred)/(1+exp(pred)))

model_pred = 1*(fitted > 0.5)
gender_binary = 1*(y_test == "F")

roc_obj <- roc(gender_binary, fitted)

acc_vals[i] = mean(1*(model_pred == gender_binary))

output_auc[i] <- roc_obj$auc

csvfile <- paste0("vivid_sac_",i,".csv")

saveRDS(vivid.sacurine, file = csvfile)

rm(vivid.sacurine)

}



saveRDS(vivid.sacurine, file = "C:/Users/Connor/Dropbox (Sydney Uni)/PhD_Student_Connor_Smith/Talk4_EcoSta_2019/sydney_xaringan-master/vivid_sacurine_train_1.rds")

x_train <- sacurine$dataMatrix[trainIndex,unlist(vivid.sacurine$opt_model)]

x_test <- sacurine$dataMatrix[-trainIndex,]
y_test <- sacurine$sampleMetadata$gender[-trainIndex]



glm.fit <- cv.glmnet(x_train, y_train, alpha = 0, family = "binomial")

new_x_test <- x_test[,unlist(vivid.sacurine$opt_model)]

pred <- predict(glm.fit, s = "lambda.min",newx = new_x_test)

fitted <- exp(pred)/(1+exp(pred))

model_pred = 1*(fitted > 0.5)
gender_binary = 1*(y_test == "F")
accurate = 1*(model_pred == gender_binary)
sum(accurate)


x_train <- sacurine$dataMatrix[trainIndex,]
glm.lasso <- cv.glmnet(x_train, y_train, alpha = 1, family = "binomial")
pred2 <- predict(glm.lasso, s = "lambda.1se",newx = x_test)
fitted2 <- exp(pred2)/(1+exp(pred2))
model_pred2 = 1*(fitted2 > 0.5)
gender_binary2 = 1*(y_test == "F")
accurate2 = 1*(model_pred2 == gender_binary2)
sum(accurate2)
#Full data

library(VennDiagram)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(parallel)
library(mvtnorm)
library(glmnet)
library(Biobase)
library(GEOquery)
library(limma)
library(blkbox)
library(tidyverse)
library(furrr)
library(ROCR)


k = 500

# load series and platform data from GEO

gset <- getGEO("GSE99039", GSEMatrix =TRUE, AnnotGPL=TRUE)

if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0(
  "011011111X1010010XX00001X000X00001X111XXX001101111",
  "XX0X01X1010XX1111111110XXXX1X1110110XXXXXXXXX01010",
  "XXX01000010010000X00000000001X0X1111X1001011111000",
  "00010011001XX10000X000010XXXXX11111X11X11111XXXX11",
  "0110XX10111110110111110111011X11011011X110X0111011",
  "1100011000X101000000001000X11000000000001011110000",
  "0010XX00110X0100000010000000XXX110XX1111X100001011",
  "10000X1101011000100111000010000X0X001X1X1011001111",
  "0010001010111X00X000100011101111101X00000X1XXXXXXX",
  "XXXXXXXXXXXX10XX11X0XX000X0000000000111101XX011X11",
  "111111111111X10000000000X00XX0X011X1X0110X00XX000X",
  "XX0X0XXX"
)
sml <- c()
for (i in 1:nchar(gsms)) {
  sml[i] <- substr(gsms, i, i)
}

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[, sel]

# log2 transform
ex <- exprs(gset)
qx <-
  as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
LogC <- (qx[5] > 100) ||
  (qx[6] - qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) {
  ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex)
}

# set up the data and proceed with analysis
sml <- paste("G", sml, sep = "")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix( ~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1 - G0, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2,
               adjust = "fdr",
               sort.by = "B",
               number = k)

tT <-
  subset(
    tT,
    select = c(
      "ID",
      "adj.P.Val",
      "P.Value",
      "t",
      "B",
      "logFC",
      "Gene.symbol",
      "Gene.title"
    )
  )

names = featureNames(gset)
tT$ID <- factor(tT$ID, levels = names)
gene_new = c(tT$ID)
gset_new = gset[gene_new]



data <- exprs(gset) %>%
  t()

y <- design[,1] %>%
  as.numeric()

vivid_top_2000_full <- vivid(data, y, bootstraps = 100, variables = gene_new, cores = 7, seed = 123, min_size = 2, lambda = 'lambda.1se', gamma = 1)
saveRDS(vivid_top_2000_full, file = "D://vivid_top_2000_full.rds")
remove(vivid_top_2000_full)

data_2000 <- data[,gene_new]

vivid_top_2000 <- vivid(data_2000, y, bootstraps = 100, variables = 1:length(gene_new), cores = 7, seed = 123, min_size = 2, lambda = 'lambda.1se', gamma = 1)
saveRDS(vivid_top_2000, file = "D://vivid_top_2000.rds")
remove(vivid_top_2000)

vivid_random <- vivid.split(x = data, y, bootstraps = 100, cores = 7, seed = 123, min_size = 2, lambda = 'lambda.1se', compare_method = 'EBIC', groups = 9, disjoint = TRUE, rep_features = 1)
saveRDS(vivid_random, file = "D://vivid_random.rds")
remove(vivid_random)

vivid_random_2000 <- vivid.split(x = data, y, bootstraps = 100, cores = 7, seed = 123, min_size = 2, lambda = 'lambda.1se', compare_method = 'EBIC', groups = 25, disjoint = FALSE, rep_features = gene_new)
saveRDS(vivid_random_2000 , file = "D://vivid_random_2000 .rds")
remove(vivid_random_2000)


