setwd("~/Documents/Scab phenotyping project")
library(rrBLUP) #used for the relationship matrix
library(asreml) #used to fit mixed models
library(readxl) #used to read phenotypic data
asreml.options(maxit=50)

#read in data
pheno <- as.data.frame(read_xls("FHBpheno_2020-21.xls"))#phenotypic data
geno<- read.csv('FHBgeno_2020-21.csv', row.names=1) #genotypic data

#select data for model training and validation
pheno<- droplevels.data.frame(pheno[c(which(pheno$GS_trainingset==1), 
                                      which(pheno$GS_validset==1)),])
geno<- geno[which(row.names(geno) %in% pheno$name2),]
 
#make relationship matrix
Gmat<- A.mat(as.matrix(geno)) #make relationship matrix
diag(Gmat)<- diag(Gmat)+0.00001 #bend matrix

#Change characters to factors
pheno$studyName<- as.factor(pheno$studyName)
pheno$blockNumber<- as.factor(pheno$blockNumber)
pheno$name2<- as.factor(pheno$name2)

#Change phenotypic data to numeric
pheno$DON<- as.numeric(pheno$DON)
pheno$FDK_V<- as.numeric(pheno$FDK_V)
pheno$FDK_L<- as.numeric(pheno$FDK_L)
pheno$FDK_Lhat<- as.numeric(pheno$FDK_Lhat)

#Consider FDK_L and FDK_Lhat as one trait
pheno$FDK_Lall<- pheno$FDK_L
pheno[which(pheno$GS_validset==1),'FDK_Lall']<- 
  pheno[which(pheno$GS_validset==1),'FDK_Lhat']

#Mask DON data on the validation set
pheno0<- pheno
pheno[which(pheno$GS_validset==1),'DON']<- NA


##########################
##     Fit GS models    ##
##########################

#Multi-trait model with both FDK_V, FDK_Lall, and DON
amod1<- asreml(fixed=cbind(DON, FDK_V, FDK_Lall)~trait+trait:studyName,
               random= ~trait:studyName:blockNumber+us(trait):vm(name2, Gmat),
               residual = ~id(units):us(trait),workspace=16e6, 
               data=pheno,  na.action = na.method(y='include', x='include'))

#Multi-trait model with FDK_Lall and DON
amod2<- asreml(fixed=cbind(DON, FDK_Lall)~trait+trait:studyName,
               random= ~trait:studyName:blockNumber+us(trait):vm(name2, Gmat),
               residual = ~id(units):us(trait),workspace=16e6,
               data=pheno,  na.action = na.method(y='include', x='include'))

#Multi-trait model FDK_V and DON
amod3<- asreml(fixed=cbind(DON, FDK_V)~trait+trait:studyName,
               random= ~trait:studyName:blockNumber+us(trait):vm(name2, Gmat),
               residual = ~id(units):us(trait),workspace=16e6, 
               data=pheno,  na.action = na.method(y='include', x='include'))

#Single trait model with DON
amod4<- asreml(fixed=DON~1+studyName,
               random= ~studyName:blockNumber+vm(name2, Gmat),workspace=16e6, 
               data=pheno,  na.action = na.method(y='include', x='include'))

#Single trait model for DON with no missing DON data for validation
amod0<- asreml(fixed=DON~1+studyName,
               random= ~studyName:blockNumber+vm(name2, Gmat),workspace=16e6, 
               data=pheno0,  na.action = na.method(y='include', x='include'))

#get genomic BLUPs from each model
mod1blup<- predict(amod1, classify='name2:trait', ignore=c('trait', 'trait:studyName'))$pvals
mod2blup<- predict(amod2, classify='name2:trait', ignore=c('trait', 'trait:studyName'))$pvals
mod3blup<- predict(amod3, classify='name2:trait', ignore=c('trait', 'trait:studyName'))$pvals
mod4blup<- predict(amod4, classify='name2', ignore=c('(Intercept)','studyName'))$pvals
donBLUP<- predict(amod0, classify='name2', ignore=c('(Intercept)','studyName'))$pvals

#names of the validation set
valid_names<- unique(pheno[which(pheno$GS_validset==1),'name2'])

#####################################
## Predictive abilities of models  ##
#####################################

##FDK_Lall and DON
#merge predicted DON with 'true' DON
rslt<- merge(mod2blup[which(mod2blup$trait=='DON'),], 
             donBLUP, by='name2') 
a<- cor.test(rslt[which(rslt$name2 %in%valid_names),'predicted.value.x'],
             rslt[which(rslt$name2 %in%valid_names),'predicted.value.y']) #correlation between values
rslt$model<- 'FDK_Lall' #add secondary trait information to the data.frame
rsltAll<- rslt[which(rslt$name2 %in%valid_names),] #save predicted and actual values to results data.frame

##FDK_V and DON
rslt<- merge(mod3blup[which(mod3blup$trait=='DON'),], 
             donBLUP, by='name2') #merge predicted DON with 'true' DON
b<- cor.test(rslt[which(rslt$name2 %in%valid_names),'predicted.value.x'],
             rslt[which(rslt$name2 %in%valid_names),'predicted.value.y']) #correlation between values
rslt$model<- 'FDK_V' #add secondary trait information to the data.frame
rsltAll<- rbind(rsltAll, 
                rslt[which(rslt$name2 %in%valid_names),]) #save predicted and actual values to results data.frame

##DON alone
rslt<- merge(mod4blup, donBLUP, by='name2') #merge predicted DON with 'true' DON
d<- cor.test(rslt[which(rslt$name2 %in%valid_names),'predicted.value.x'],
             rslt[which(rslt$name2 %in%valid_names),'predicted.value.y']) #correlation between values
rslt$model<- 'None' #add secondary trait information to the data.frame
rsltAll<- plyr::rbind.fill(rsltAll, 
                     rslt[which(rslt$name2 %in%valid_names),]) #save predicted and actual values to results data.frame

##both FDK_V and FDK_Lall 
rslt<- merge(mod1blup[which(mod1blup$trait=='DON'),], donBLUP, by='name2') #correlation between values
e<- cor.test(rslt[which(rslt$name2 %in%valid_names),'predicted.value.x'],
             rslt[which(rslt$name2 %in%valid_names),'predicted.value.y']) #correlation between values
rslt$model<- 'FDK_V and FDK_Lall' #add secondary trait information to the data.frame
rsltAll<- rbind(rsltAll, 
                rslt[which(rslt$name2 %in%valid_names),]) #save predicted and actual values to results data.frame

#Summary of accuracy results
df<- data.frame(Secondary_trait= c('FDK_Lall', 'FDK_V', 'None'),
                Accuracy= c(a$estimate, b$estimate, d$estimate))


#####################################################
## Correlations between traits and narrow sense h2 ##
##                   -Validation set-              ##
#####################################################
pheno0v<- droplevels.data.frame(pheno0[which(pheno0$GS_validset == 1),]) #subset data
#fit model
amodFull<- asreml(fixed=cbind(DON, FDK_V, FDK_Lall)~trait,
                  random= ~trait:blockNumber+us(trait):vm(name2, Gmat),
                  residual = ~id(units):us(trait),workspace=16e6, 
                  data=pheno0v,  na.action = na.method(y='include', x='include'))

#get genetic correlations from model results
traits_sub<- c('DON', 'FDK_V', 'FDK_Lall')
gencor<- matrix(nrow=length(traits_sub), ncol=length(traits_sub))
for(i in 1:length(traits_sub)){
  for(j in 1:length(traits_sub)){
    num1<-summary(amodFull)$varcomp[paste('trait:vm(name2, Gmat)!trait_', traits_sub[i], ":", traits_sub[j], sep=""),'component']
    num2<-summary(amodFull)$varcomp[paste('trait:vm(name2, Gmat)!trait_', traits_sub[j], ":", traits_sub[i], sep=""),'component']
    num<- unique(na.omit(c(num1,num2)))
    denom_a<- summary(amodFull)$varcomp[paste('trait:vm(name2, Gmat)!trait_', traits_sub[i], ":", traits_sub[i], sep=""),'component']
    denom_b<- summary(amodFull)$varcomp[paste('trait:vm(name2, Gmat)!trait_', traits_sub[j], ":", traits_sub[j], sep=""),'component']
    cr<- num/c(sqrt(denom_a)*sqrt(denom_b))
    gencor[i,j]<- cr
  }
}
colnames(gencor)<- traits_sub
rownames(gencor)<- traits_sub

#matrix of genetic correlations
gencor

#calculate narrow sense heritability for each trait
ixG<- match(pheno0v$name2, row.names(Gmat))
meanD<- mean(diag(Gmat[ixG, ixG]))
Vg<- summary(amodFull)$varcomp[paste('trait:vm(name2, Gmat)!trait_', traits_sub, ":", traits_sub, sep=""),'component']
Ve<- summary(amodFull)$varcomp[paste('units:trait!trait_', traits_sub, ":", traits_sub, sep=""),'component']
h2<- (meanD*Vg)/((meanD*Vg)+Ve)
names(h2)<- traits_sub
#vector of heritabilities
h2

#####################################################
## Correlations between traits and narrow sense h2 ##
##                   -Training set-                ##
#####################################################

#Fit model
pheno0t<- droplevels.data.frame(pheno0[which(pheno0$GS_trainingset == 1),]) #subset data
amodFull<- asreml(fixed=cbind(DON, FDK_V, FDK_Lall)~trait+trait:studyName,
                  random= ~trait:studyName:blockNumber+us(trait):vm(name2, Gmat),
                  residual = ~id(units):us(trait),workspace=16e6, 
                  data=pheno0t,  na.action = na.method(y='include', x='include'))
amodFull<- update(amodFull)

#get genetic correlations from model results
traits_sub<- c('DON', 'FDK_V', 'FDK_Lall')
gencor<- matrix(nrow=length(traits_sub), ncol=length(traits_sub))
for(i in 1:length(traits_sub)){
  for(j in 1:length(traits_sub)){
    num1<-summary(amodFull)$varcomp[paste('trait:vm(name2, Gmat)!trait_', traits_sub[i], ":", traits_sub[j], sep=""),'component']
    num2<-summary(amodFull)$varcomp[paste('trait:vm(name2, Gmat)!trait_', traits_sub[j], ":", traits_sub[i], sep=""),'component']
    num<- unique(na.omit(c(num1,num2)))
    denom_a<- summary(amodFull)$varcomp[paste('trait:vm(name2, Gmat)!trait_', traits_sub[i], ":", traits_sub[i], sep=""),'component']
    denom_b<- summary(amodFull)$varcomp[paste('trait:vm(name2, Gmat)!trait_', traits_sub[j], ":", traits_sub[j], sep=""),'component']
    cr<- num/c(sqrt(denom_a)*sqrt(denom_b))
    gencor[i,j]<- cr
  }
}
colnames(gencor)<- traits_sub
rownames(gencor)<- traits_sub
#matrix of genetic correlations
gencor

#calculate narrow sense heritability for each trait
ixG<- match(pheno0t$name2, row.names(Gmat))
meanD<- mean(diag(Gmat[ixG, ixG]))
Vg<- summary(amodFull)$varcomp[paste('trait:vm(name2, Gmat)!trait_', traits_sub, ":", traits_sub, sep=""),'component']
Ve<- summary(amodFull)$varcomp[paste('units:trait!trait_', traits_sub, ":", traits_sub, sep=""),'component']
h2<- (meanD*Vg)/((meanD*Vg)+Ve)
names(h2)<- traits_sub
#vector of heritabilities
h2
