
#####################################################################################################
# Determine broad and narrow sense heritability of FDKL, FDKLhat, and DON in training and test sets #
#####################################################################################################

setwd("/Users/Dir/")

library(asreml)
library(data.table)
library(magrittr)
library(rrBLUP)
library(tidyverse)

# Upload geno and pheno data ----------------------------------------------

pheno <- read_tsv("FHBpheno_2020-21.tsv") %>% 
  mutate(across(c(1:15), factor))

geno <- as.matrix(fread("FHBgeno_2020-21.csv"), rownames = TRUE) 

# Narrow-Sense with only entries with GBS data -----------------------------
 # Seperate pheno file into training and validation

train_h2 <- pheno %>%  filter(GS_trainingset == 1) %>%
  print(width = Inf) 

val_h2 <- pheno %>%  filter(GS_validset == 1) %>% 
  print(width = Inf) 

# Create kinship matrix with gmat --------------------------------------------------------------------

genoTrain_h2 <- geno[rownames(geno) %in% train_h2$name2, ] # Filter gbs file for phenotyped lines
gmatTrain_h2 <- A.mat(genoTrain_h2)
diag(gmatTrain_h2)<- diag(gmatTrain_h2)+0.00001
diagTrain <- mean(diag(gmatTrain_h2))

trainGBS_h2 <- train_h2 %>% filter(name2 %in% row.names(genoTrain_h2)) # Filter pheno file for lines with gbs data

genoVal_h2 <- geno[rownames(geno) %in% val_h2$name2, ] # Filter gbs file for phenotyped lines
gmatVal_h2 <- A.mat(genoVal_h2)
diag(gmatVal_h2)<- diag(gmatVal_h2)+0.00001
diagVal <- mean(diag(gmatVal_h2))

valGBS_h2 <- val_h2 %>% filter(name2 %in% row.names(genoVal_h2)) # Filter pheno file for lines with gbs data

# Fit mixed model ---------------------------------------------------------
  # Fit model for training set
    # FDKV

amod_Train_FDKV_h2 <- asreml(fixed = FDK_V ~ studyName:blockNumber,
                             random= ~  vm(name2, gmatTrain_h2),
                             residual = ~ id(units),
                             workspace = 64e6, data = trainGBS_h2,  
                             na.action = na.method(y='omit', x='omit'))

summary(amod_Train_FDKV_h2)$varcomp
V1 <- summary(amod_Train_FDKV_h2)$varcomp[1,1]*diagTrain
V2 <- summary(amod_Train_FDKV_h2)$varcomp[2,1]
train_FDKV_h2 <- V1 / (V1+V2)

    # FDKL

amod_Train_FDKL_h2 <- asreml(fixed = FDK_L ~ studyName:blockNumber,
                             random= ~ vm(name2, gmatTrain_h2),
                             residual = ~ id(units),
                             workspace = 64e6, data = trainGBS_h2,  
                             na.action = na.method(y='omit', x='omit'))

summary(amod_Train_FDKL_h2)$varcomp
V1 <- summary(amod_Train_FDKL_h2)$varcomp[1,1]*diagTrain
V2 <- summary(amod_Train_FDKL_h2)$varcomp[2,1]
train_FDKL_h2 <- V1 / (V1+V2)

    # DON

amod_Train_DON_h2 <- asreml(fixed = DON ~ studyName:blockNumber,
                            random= ~ vm(name2, gmatTrain_h2),
                            residual = ~ id(units),
                            workspace = 64e6, data = trainGBS_h2,  
                            na.action = na.method(y='omit', x='omit'))

summary(amod_Train_DON_h2)$varcomp
V1 <- summary(amod_Train_DON_h2)$varcomp[1,1]*diagTrain
V2 <- summary(amod_Train_DON_h2)$varcomp[2,1]
train_DON_h2 <- V1 / (V1+V2)

# Fit model for validation set
    # FDKv

amod_Val_FDKV_h2 <- asreml(fixed = FDK_V ~ blockNumber,
               random= ~ vm(name2, gmatVal_h2),
               residual = ~ id(units),
               workspace = "8gb", data=val_h2,  
               na.action = na.method(y='omit', x='omit'))

summary(amod_Val_FDKV_h2)$varcomp
V1 <- summary(amod_Val_FDKV_h2)$varcomp[1,1]*diagVal
V2 <- summary(amod_Val_FDKV_h2)$varcomp[2,1]
val_FDKV_h2 <- V1 / (V1+V2)

    # FDKL-hat

amod_Val_FDKLhat_h2 <- asreml(fixed = FDK_Lhat ~ blockNumber,
                    random= ~ vm(name2, gmatVal_h2),
                    residual = ~ id(units),
                    workspace = "8gb", data=val_h2,  
                    na.action = na.method(y='omit', x='omit'))

summary(amod_Val_FDKLhat_h2)$varcomp
V1 <- summary(amod_Val_FDKLhat_h2)$varcomp[1,1]*diagVal
V2 <- summary(amod_Val_FDKLhat_h2)$varcomp[2,1]
val_FDKLhat_h2 <- V1 / (V1 + V2)
   
    # DON

amod_Val_DON_h2 <- asreml(fixed = DON ~ blockNumber,
                     random= ~ vm(name2, gmatVal_h2),
                     residual = ~ id(units),
                     workspace = 64e6, data=valGBS_h2,  
                     na.action = na.method(y='omit', x='omit'))

summary(amod_Val_DON_h2)$varcomp
V1 <- summary(amod_Val_DON_h2)$varcomp[1,1]*diagVal
V2 <- summary(amod_Val_DON_h2)$varcomp[1,2]
val_DON_h2 <- V1 / (V1 + V2)

# Save results for narrow sense heritability --------------------------------

c1 <- c("train_FDKV_h2", "train_FDKL_h2", "train_DON_h2",
        "val_FDKV_h2", "val_FDKLhat_h2", "val_DON_h2")
c2 <- c(train_FDKV_h2, train_FDKL_h2, train_DON_h2,
       val_FDKV_h2, val_FDKLhat_h2, val_DON_h2)
h2 <- tibble(populationMethod = c1, Estimate = c2)

print(h2)

#### END NARROW SENSE SEGMENT ####

# Broad sense with all entries with phenotype data for FDKV, FDKL/FDKLhat, and DON --------

train_H2 <- pheno %>%  filter(GS_trainingset == 1) %>% #filter for studyName
  print(width = Inf) 

val_H2 <- pheno %>%  filter(GS_validset == 1) %>% 
  print(width = Inf) 

# Fit mixed model ---------------------------------------------------------
  # Fit model for training set
    # FDKV

amod_Train_FDKV_H2 <- asreml(fixed = FDK_V ~ studyName:blockNumber,
                       random = ~  name2,
                       residual = ~ id(units), workspace=64e6,
                       data=train_H2,  na.action = na.method(y='omit', x='omit'))

summary(amod_Train_FDKV_H2)$varcomp

train_FDKV_H2 <- vpredict(amod_Train_FDKV_H2, train_FDKV_H2 ~ V1 / (V1 + V2))

    # FDKL

amod_Train_FDKL_H2 <- asreml(fixed = FDK_L ~ studyName:blockNumber,
                         random = ~ name2,
                         residual = ~ id(units), workspace=64e6,
                         data=train_H2,  na.action = na.method(y='omit', x='omit'))

summary(amod_Train_FDKL_H2)$varcomp

train_FDKL_H2 <- vpredict(amod_Train_FDKL_H2, train_FDKL_H2 ~ V1 / (V1 + V2))

    # DON 

amod_Train_DON_H2 <- asreml(fixed = DON ~ studyName:blockNumber,
                       random = ~ name2,
                       residual = ~ id(units), workspace=64e6,
                       data=train_H2,  na.action = na.method(y='omit', x='omit'))

summary(amod_Train_DON_H2)$varcomp

train_DON_H2 <- vpredict(amod_Train_DON_H2, train_DON_H2 ~ V1 / (V1 + V2))

# Fit model for validation set
  # FDKV

amod_Val_FDKV_H2 <- asreml(fixed = FDK_V ~ blockNumber,
                      random = ~ name2,
                      residual = ~ id(units),workspace=64e6,
                      data=val_H2,  na.action = na.method(y='omit', x='omit'))

summary(amod_Val_FDKV_H2)$varcomp

val_FDKV_H2 <- vpredict(amod_Val_FDKV_H2, val_Vis_H2 ~ V1 / (V1 + V2))

  # FDKLhat

amod_Val_FDKLhat_H2 <- asreml(fixed = FDK_Lhat ~ blockNumber,
                        random = ~ name2,
                        residual = ~ id(units), workspace=64e6,
                        data=val_H2,  na.action = na.method(y='omit', x='omit'))

summary(amod_Val_FDKLhat_H2)$varcomp

val_FDKLhat_H2 <- vpredict(amod_Val_FDKLhat_H2, val_FDKLhat_H2 ~ V1 / (V1 + V2))

  # DON

amod_Val_DON_H2 <- asreml(fixed = DON ~ studyName:blockNumber,
                     random = ~ name2,
                     residual = ~ id(units), workspace=64e6,
                     data=val_H2,  na.action = na.method(y='omit', x='omit'))

summary(amod_Val_DON_H2)$varcomp

val_DON_H2 <- vpredict(amod_Val_DON_H2, val_DON_H2 ~ V1 / (V1 + V2))

#### END BROAD SENSE SEGMENT ####

# Save results for broad sense heritabilities----------------------------

H2 <- tibble(rownames_to_column(bind_rows(train_FDKV_H2, train_FDKL_H2, train_DON_H2, 
                                          val_FDKV_H2, val_FDKLhat_H2, val_DON_H2), 
                                          var = "populationMethod"))

print(H2)

# Combine broad sense and narrow sense and export -------------------------

FHB_heritabilities <- bind_rows(H2, h2) %>% 
  separate(col = populationMethod, into = c("set", "trait", "heritability" ), sep = "_") %>% 
  write_csv("FHB_heritabilities.csv")

#### END ####
