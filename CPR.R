#######################################################
# Code for SARS-CoV-2 clinical prediction rule (Jody Reimer first author)
#######################################################
rm(list=ls())
graphics.off();

library(readr)
library(tidyverse)
library(ranger)
library(pROC)
library(cvAUC)
setwd()

####################
#import and merging data
###################
data <- read_csv("PatientData.csv")

################### merge in ADI - not needed for fake data ####
# ADI <- read.csv("C://UT_2015_ADI_9 Digit Zip Code_v2.0.txt")
# 
# #drop the extra G at the beginning of each ZIPID
# ADI$zip <- as.numeric(substr(ADI$ZIPID,2,10))
# 
# #remove the hyphen from data pat_zipcode
# data$zip <- as.numeric(gsub("-", "", data$pat_zipcode_x ))
# 
# data$zip_length <- nchar(sub('^0+','',sub('\\.','',data$zip)))
# 
# ADI_5digit <- read.csv("C://Zip Code Agg Data.csv")
# #from student  manually averaged ADI's into 5-digit zips
# 
# #using 9-digit zip when have it, using 5-digit zip when don't
# data <- data %>% left_join(ADI, by="zip")
# data <- data %>% left_join(ADI_5digit, by=c("zip"="ZIP5"))
# #if have 9-digit zip in our data but no matching 9-digit ADI zip,
# #truncate to 5-digit ADI
# data <- data %>% mutate(zip_truncated=(case_when(((is.na(ADI_NATRANK))&(is.na(AVG.ADI_NATRANK.)))~as.numeric(substr(data$zip,1,5)), TRUE~zip)))
# data <- data %>% left_join(ADI_5digit, by=c("zip_truncated"="ZIP5"))
# #if had 9-digit zip, want ADI_NATRANK
# #if originally had 5-digit zip, want AVG.ADI_NATRANK..x
# #if originally had neither, truncated to 5-digit zip and tried again, want AVG.ADI_NATRANK..y
# data <- data %>% mutate(ADI_natl=(case_when((!is.na(ADI_NATRANK))~as.numeric(ADI_NATRANK),
#                                             (!is.na(AVG.ADI_NATRANK..x))~as.numeric(AVG.ADI_NATRANK..x),
#                                             (!is.na(AVG.ADI_NATRANK..y))~as.numeric(AVG.ADI_NATRANK..y))))
# data <- data %>% mutate(ADI_UT=(case_when((!is.na(ADI_STATERNK))~as.numeric(ADI_STATERNK),
#                                           (!is.na(AVG.ADI_STATERNK..x))~as.numeric(AVG.ADI_STATERNK..x),
#                                           (!is.na(AVG.ADI_STATERNK..y))~as.numeric(AVG.ADI_STATERNK..y))))
# 
# data <- data %>% mutate(zip_5digit=(case_when((zip_length==9)~as.numeric(substr(data$zip,1,5)), TRUE~zip)))

####################
#data cleaning
################### define variables - unneeded for fake data ####
# data=data %>% mutate(COVID=as.numeric(interpreted_result=="POSITIVE"))
# 
# data=data %>% mutate(gender=(case_when(is.na(pat_gender_desc)~2,pat_gender_desc=="Male"~1,pat_gender_desc=="Female"~0,pat_gender_desc=="Unknown"~2)))
# #0=Female, 1=Male, 2=Unknown/missing
# 
# yesnocols<-c("prior_exp","travel","cough","fever","nasal","sore_throat",         
#              "headache","sob","hc_worker",
#              "diarrhea","myalgia","lethargy","nausea_vomit")
# #for symptoms, 0=no, 1=yes, 2=not recorded or missing
# #not recorded/missin(2) to be no(0) for symptoms, 
# #if no symptom response (if "NA") drop from analysis
# data<-data %>% mutate_at(vars(yesnocols),~case_when(is.na(.)~.,.==1~1,.==0~0,.==2~0))
# 
# #smoker is: never(0), former(1), current(2), missing(3)
# data=data %>% mutate(smoker=(case_when(is.na(smoker) ~ 3, TRUE~smoker)))
# 
# #a single variable for any comorbidities, not including obesity
# comorbids<-c("CHF_ind","COAG_ind","HTNCX_ind","HTNWCHF_ind","LIVER_ind",
#              "LYTES_ind","NEURO_ind","PULMCIRC_ind","RENLFAIL_ind","ULCER_ind",
#              "ALCOHOL_ind","ANEMDEF_ind","ARTH_ind","BLDLOSS_ind","CHRNLUNG_ind",
#              "DEPRESS_ind","DM_ind","DMCX_ind","DRUG_ind","HHRWCHF_ind",
#              "HHRWHRF_ind","HHRWOHRF_ind","HRENWORF_ind","HRENWRF_ind","HTN_ind",
#              "HTNPREG_ind","HTNWOCHF_ind","HYPOTHY_ind","LYMPH_ind","METS_ind",
#              "OHTNPREG_ind","PARA_ind","PERIVASC_ind","PSYCH_ind","TUMOR_ind",
#              "VALVE_ind")
# data$comorbid <- ((rowSums(data[comorbids] == 1, na.rm=T) > 0) * 1)
# #comorbid==0 if no comorbidities, ==1 if have at least one comorbidity
# 
# #not included: OBESE_ind, WGHTLOSS_ind
# 
# #adding a Hispanic category to race for anyone who said race==other or unknown AND said ethnicity==Hispanic
# data <- data %>% mutate(race=(case_when(((pat_race_desc=="Other" | pat_race_desc=="Unknown") & 
#                                            pat_ethnicity_desc=="Hispanic/Latino")~"Hispanic", TRUE~pat_race_desc)))
# data <- data %>% mutate(race=(case_when(is.na(race) ~ "Unknown", TRUE~race)))
# 
# #making white reference
# data$race.factor <- as.factor(data$race)
# data <- within(data, race.factor<-relevel(race.factor,ref="White or Caucasian"))
# 
# #adding Native into Other since small cell sizes
# data <- data %>% mutate(race_temp=(case_when((race=="American Indian and Alaska Native") ~ "Other", TRUE~race)))

################### drop missing - unneeded for fake data ####
# data <- data %>% filter(!is.na(COVID))
# 
# #dropping non-binary gender for now since too few ppl
# data <- data %>% filter(gender !=2)
# 
# #drop those missing symptom data (were originally no response, keep those ==2)
# data <- data %>% filter(!is.na(prior_exp)&!is.na(travel)&!is.na(nasal)&!is.na(sore_throat)&
#                           !is.na(headache)&!is.na(hc_worker)&!is.na(diarrhea)&!is.na(myalgia)&
#                           !is.na(lethargy)&!is.na(nausea_vomit))

################### subset, list of vars considering ####
data <- data %>% filter(!is.na(ADI_UT))

#easier to just choose which variables to explore rather than drop some
names.main <- c("indexage","gender",
                     "prior_exp","travel","cough","fever","nasal",
                      "sore_throat","headache","sob","hc_worker","smoker",
                      "diarrhea","myalgia","lethargy","nausea_vomit", "comorbid",
                      "ADI_UT")
names.supplA <- c("indexage","gender",
                     "prior_exp","travel","hc_worker","smoker",
                     "comorbid","ADI_UT")
names.supplB <- c("indexage","gender",
                     "prior_exp","travel","cough","fever","nasal",
                     "sore_throat","headache","sob","hc_worker","smoker",
                     "diarrhea","myalgia","lethargy","nausea_vomit", "comorbid",
                     "ADI_UT","spo2_value","pulse_value")
names.supplC <- c("indexage","gender",
                     "prior_exp","travel","cough","fever","nasal",
                     "sore_throat","headache","sob","hc_worker","smoker",
                     "diarrhea","myalgia","lethargy","nausea_vomit", "comorbid",
                     "ADI_UT","race_temp")

####################
#Random Forest (RF)
################### main model ####
#need complete cases, are a handful of observations missing for these predictors
data <- data[which(!is.na(data$fever)&!is.na(data$sob)),]

COVID <- "COVID"
out=ranger(as.formula(paste(COVID, '~' ,paste(names.main,collapse="+"),sep="")),data=data,importance="impurity",num.trees=1000)
imps=importance(out)

df_imps=data.frame(names=names(imps),var_red=as.numeric(imps)) %>% arrange(desc(var_red))
colnames(df_imps) <- c(paste("n=",dim(data)[1],"    ",colnames(df_imps)[1]),
                       colnames(df_imps)[2])
df_imps
# n= 1928      names   var_red
# 1            indexage 23.725018
# 2              ADI_UT 17.880405
# 3           prior_exp  5.005538
# 4              smoker  4.748546
# 5              travel  4.120636
# 6                 sob  3.595011
# 7            headache  3.374923
# 8              gender  3.315168
# 9             myalgia  3.273567
# 10              nasal  3.100431
# 11           lethargy  2.879487
# 12        sore_throat  2.737227
# 13              fever  2.640746
# 14           comorbid  2.331370
# 15          hc_worker  2.245508
# 16           diarrhea  1.979109
# 17              cough  1.923026
# 18       nausea_vomit  1.367887

nvars_opts=c(1:dim(df_imps)[1])
result=data.frame(iter=NA,nvar=NA,true=NA,pred_glm=NA,pred_RF=NA)
data$index=1:dim(data)[1]

test_record <- NA
train_record <- NA

for (each in 1:100){
  print(each)
  train=data %>% sample_frac(.80,replace=F)
  
  test=data[-which(data$index %in% train$index),]
  
  train_record <- c(train_record,table(train$COVID)[["1"]])
  test_record <- c(test_record,table(test$COVID)[["1"]])
  
  out=ranger(as.formula(paste(COVID,'~',paste(names.main,collapse="+"),sep="")),data=train,num.trees=1000,importance="impurity")
  df_imps=data.frame(names=names(ranger::importance(out)),imps=ranger::importance(out)) %>% arrange(desc(imps))
  for (nvars in nvars_opts){
    
    print(nvars)
    out1=glm(as.formula(paste(COVID,'~',paste(df_imps$names[1:nvars],collapse="+"),sep="")),data=train,family="binomial")
    out2=ranger(as.formula(paste(COVID,'~',paste(df_imps$names[1:nvars],collapse="+"),sep="")),data=train,num.trees=1000)
    
    df=data.frame(iter=each,nvar=nvars,true=test[[COVID]],pred_glm=as.numeric(predict(out1,newdata=test,type="response")),pred_RF=as.numeric(predict(out2,data=test,type="response")$predictions))
    result=rbind(result,df)
  }
}
result=result[-1,]
#write.csv(result,"C:/.csv")


#check number of cases in the test/training dataset
train_record=train_record[-1]
test_record=test_record[-1]
train_record
test_record

AUCs=result %>% split(.$nvar) %>% purrr::map(~ci.cvAUC(.$pred_glm,.$true,folds=.$iter))
AUCs2=result %>% split(.$nvar) %>% purrr::map(~ci.cvAUC(.$pred_RF,.$true,folds=.$iter))

AUC_df=rbind(bind_rows(AUCs %>% purrr::map(~data.frame(AUC=.$cvAUC,SE=.$se,lower=.$ci[1],upper=.$ci[2],level=.$confidence,Model="LR"))),
             bind_rows(AUCs2 %>% purrr::map(~data.frame(AUC=.$cvAUC,SE=.$se,lower=.$ci[1],upper=.$ci[2],level=.$confidence,Model="RF"))))
AUC_df$nvar=rep(nvars_opts,2)
AUC_df
#          AUC          SE     lower     upper level Model nvar
# 1  0.4940545 0.005823689 0.4826403 0.5054688  0.95    LR    1
# 2  0.5940151 0.005851717 0.5825459 0.6054842  0.95    LR    2
# 3  0.6239789 0.005694941 0.6128170 0.6351407  0.95    LR    3
# 4  0.6511177 0.005578162 0.6401847 0.6620507  0.95    LR    4
# 5  0.6744313 0.005304758 0.6640341 0.6848284  0.95    LR    5
# 6  0.6841815 0.005200674 0.6739884 0.6943746  0.95    LR    6
# 7  0.6932918 0.005059162 0.6833760 0.7032076  0.95    LR    7
# 8  0.6981037 0.005038500 0.6882285 0.7079790  0.95    LR    8
# 9  0.7007493 0.005022239 0.6909059 0.7105927  0.95    LR    9
# 10 0.7015145 0.005024585 0.6916665 0.7113625  0.95    LR   10
# 11 0.7018568 0.004996626 0.6920636 0.7116501  0.95    LR   11
# 12 0.7044049 0.004957406 0.6946886 0.7141212  0.95    LR   12
# 13 0.7097677 0.004830257 0.7003005 0.7192348  0.95    LR   13
# 14 0.7095464 0.004851231 0.7000381 0.7190546  0.95    LR   14
# 15 0.7075866 0.004913274 0.6979568 0.7172165  0.95    LR   15
# 16 0.7056223 0.004955989 0.6959087 0.7153358  0.95    LR   16
# 17 0.7022946 0.005018493 0.6924585 0.7121306  0.95    LR   17
# 18 0.7003548 0.005073537 0.6904109 0.7102988  0.95    LR   18
# 19 0.5296588 0.007474050 0.5150099 0.5443077  0.95    RF    1
# 20 0.5949745 0.005442716 0.5843070 0.6056421  0.95    RF    2
# 21 0.6247781 0.005772281 0.6134646 0.6360915  0.95    RF    3
# 22 0.6171738 0.005772828 0.6058592 0.6284883  0.95    RF    4
# 23 0.6542839 0.005560993 0.6433846 0.6651832  0.95    RF    5
# 24 0.6606023 0.005621922 0.6495835 0.6716210  0.95    RF    6
# 25 0.6728117 0.005486344 0.6620587 0.6835647  0.95    RF    7
# 26 0.6779033 0.005420304 0.6672797 0.6885269  0.95    RF    8
# 27 0.6795914 0.005320259 0.6691638 0.6900189  0.95    RF    9
# 28 0.6803568 0.005340119 0.6698903 0.6908232  0.95    RF   10
# 29 0.6787141 0.005332701 0.6682622 0.6891660  0.95    RF   11
# 30 0.6819800 0.005269548 0.6716518 0.6923081  0.95    RF   12
# 31 0.6943995 0.005127594 0.6843496 0.7044494  0.95    RF   13
# 32 0.6931133 0.005179461 0.6829618 0.7032649  0.95    RF   14
# 33 0.6920084 0.005230663 0.6817565 0.7022603  0.95    RF   15
# 34 0.6890994 0.005233351 0.6788422 0.6993566  0.95    RF   16
# 35 0.6924501 0.005163638 0.6823296 0.7025707  0.95    RF   17
# 36 0.6907261 0.005202257 0.6805299 0.7009224  0.95    RF   18
#write.csv(AUC_df,"C:/.csv")

#pdf(file="C:/AUCnumVar.pdf",width=6, height=4)
plot(AUC_df$nvar[1:dim(df_imps)[1]],AUC_df$AUC[1:dim(df_imps)[1]],
     xlab="number of variables",ylab="AUC",
#     main="Predicting COVID+ Test in Utah",
     ylim=c(0.4,0.8),
     pch=1,col="red")
points(AUC_df$nvar[1:dim(df_imps)[1]],AUC_df$AUC[(dim(df_imps)[1]+1):dim(AUC_df)[1]],
       pch=2,col="blue")
legend("topleft",c("logistic regression","random forest"),col=c("red","blue"),pch=c(1,2))
dev.off()



top5 <- c("indexage","ADI_UT","prior_exp","smoker","travel")

COVID <- "COVID"
out=ranger(as.formula(paste(COVID, '~' ,paste(top5,collapse="+"),sep="")),data=data,importance="impurity",num.trees=1000)
imps=importance(out)

df_imps=data.frame(names=names(imps),var_red=as.numeric(imps)) %>% arrange(desc(var_red))
colnames(df_imps) <- c(paste("n=",dim(data)[1],"    ",colnames(df_imps)[1]),
                       colnames(df_imps)[2])
df_imps
# n= 1928      names   var_red
# 1           indexage 27.988744
# 2             ADI_UT 14.624491
# 3          prior_exp  4.302415
# 4             smoker  4.034006
# 5             travel  3.103529

result=data.frame(iter=NA,nvar=NA,true=NA,pred_glm=NA,pred_RF=NA)

data$index=1:dim(data)[1]

for (each in 1:100){
  print(each)
  train=data %>% sample_frac(.80,replace=F)
  
  test=data[-which(data$index %in% train$index),]
  
  for (nvars in length(top5)){
    
    print(nvars)
    out1=glm(as.formula(paste(COVID,'~',paste(top5,collapse="+"),sep="")),data=train,family="binomial")
    out2=ranger(as.formula(paste(COVID,'~',paste(top5,collapse="+"),sep="")),data=train,num.trees=1000)
    
    df=data.frame(iter=each,nvar=nvars,true=test[[COVID]],pred_glm=as.numeric(predict(out1,newdata=test,type="response")),pred_RF=as.numeric(predict(out2,data=test,type="response")$predictions))
    result=rbind(result,df)
  }
}
result=result[-1,]
#write.csv(result,"C:/.csv")

AUCs=result %>% split(.$nvar) %>% purrr::map(~ci.cvAUC(.$pred_glm,.$true,folds=.$iter))
AUCs2=result %>% split(.$nvar) %>% purrr::map(~ci.cvAUC(.$pred_RF,.$true,folds=.$iter))

AUC_df=rbind(bind_rows(AUCs %>% purrr::map(~data.frame(AUC=.$cvAUC,SE=.$se,lower=.$ci[1],upper=.$ci[2],level=.$confidence,Model="LR"))),
             bind_rows(AUCs2 %>% purrr::map(~data.frame(AUC=.$cvAUC,SE=.$se,lower=.$ci[1],upper=.$ci[2],level=.$confidence,Model="RF"))))
AUC_df
#         AUC          SE     lower     upper level Model
# 1 0.6939016 0.005269242 0.6835740 0.7042291  0.95    LR
# 2 0.6694827 0.005488610 0.6587252 0.6802402  0.95    RF




################### SupplA ####
COVID <- "COVID"
out=ranger(as.formula(paste(COVID, '~' ,paste(names.supplA,collapse="+"),sep="")),data=data,importance="impurity",num.trees=1000)
imps=importance(out)

df_imps=data.frame(names=names(imps),var_red=as.numeric(imps)) %>% arrange(desc(var_red))
colnames(df_imps) <- c(paste("n=",dim(data)[1],"    ",colnames(df_imps)[1]),
                       colnames(df_imps)[2])
df_imps
# n= 1930      names   var_red
# 1           indexage 19.031252
# 2             ADI_UT 12.338648
# 3          prior_exp  4.302280
# 4             smoker  3.406510
# 5             travel  3.100879
# 6             gender  1.991362
# 7           comorbid  1.856905
# 8          hc_worker  1.603600

nvars_opts=c(1:dim(df_imps)[1])
result=data.frame(iter=NA,nvar=NA,true=NA,pred_glm=NA,pred_RF=NA)
data$index=1:dim(data)[1]

test_record <- NA
train_record <- NA

result=data.frame(iter=NA,nvar=NA,true=NA,pred_glm=NA,pred_RF=NA)

data$index=1:dim(data)[1]

for (each in 1:100){
  print(each)
  train=data %>% sample_frac(.80,replace=F)
  
  test=data[-which(data$index %in% train$index),]
  
  for (nvars in length(names.supplA)){
    
    print(nvars)
    out1=glm(as.formula(paste(COVID,'~',paste(names.supplA,collapse="+"),sep="")),data=train,family="binomial")
    out2=ranger(as.formula(paste(COVID,'~',paste(names.supplA,collapse="+"),sep="")),data=train,num.trees=1000)
    
    df=data.frame(iter=each,nvar=nvars,true=test[[COVID]],pred_glm=as.numeric(predict(out1,newdata=test,type="response")),pred_RF=as.numeric(predict(out2,data=test,type="response")$predictions))
    result=rbind(result,df)
  } 
} 
result=result[-1,]
#write.csv(result,"C:/.csv")

AUCs=result %>% split(.$nvar) %>% purrr::map(~ci.cvAUC(.$pred_glm,.$true,folds=.$iter))
AUCs2=result %>% split(.$nvar) %>% purrr::map(~ci.cvAUC(.$pred_RF,.$true,folds=.$iter))

AUC_df=rbind(bind_rows(AUCs %>% purrr::map(~data.frame(AUC=.$cvAUC,SE=.$se,lower=.$ci[1],upper=.$ci[2],level=.$confidence,Model="LR"))),
             bind_rows(AUCs2 %>% purrr::map(~data.frame(AUC=.$cvAUC,SE=.$se,lower=.$ci[1],upper=.$ci[2],level=.$confidence,Model="RF"))))
AUC_df
#         AUC          SE     lower     upper level Model
# 1 0.6825466 0.005479498 0.6718070 0.6932862  0.95    LR
# 2 0.6684545 0.005717841 0.6572477 0.6796612  0.95    RF

################### SupplB ####
#need complete cases, are a handful of observations missing for these predictors
data <- data[which(!is.na(data$fever)&!is.na(data$sob)&!is.na(data$spo2_value)&!is.na(data$pulse_value)),]

COVID <- "COVID"
out=ranger(as.formula(paste(COVID, '~' ,paste(names.supplB,collapse="+"),sep="")),data=data,importance="impurity",num.trees=1000)
imps=importance(out)

df_imps=data.frame(names=names(imps),var_red=as.numeric(imps)) %>% arrange(desc(var_red))
colnames(df_imps) <- c(paste("n=",dim(data)[1],"    ",colnames(df_imps)[1]),
                       colnames(df_imps)[2])
df_imps
# n= 1641      names    var_red
# 1            indexage 15.1685603
# 2         pulse_value 12.0566346
# 3              ADI_UT 11.2662175
# 4          spo2_value  7.8497088
# 5              smoker  3.2328394
# 6              travel  3.0353385
# 7           prior_exp  2.9846529
# 8                 sob  2.7231327
# 9            headache  2.5187449
# 10              nasal  2.3065771
# 11             gender  2.2464674
# 12            myalgia  2.2260639
# 13           lethargy  1.9239307
# 14              fever  1.8890387
# 15        sore_throat  1.7996984
# 16           diarrhea  1.4691262
# 17          hc_worker  1.3547858
# 18           comorbid  1.2049258
# 19              cough  1.1155810
# 20       nausea_vomit  0.9705187

nvars_opts=c(1:dim(df_imps)[1])
result=data.frame(iter=NA,nvar=NA,true=NA,pred_glm=NA,pred_RF=NA)
data$index=1:dim(data)[1]

test_record <- NA
train_record <- NA

for (each in 1:100){
  print(each)
  train=data %>% sample_frac(.80,replace=F)
  
  test=data[-which(data$index %in% train$index),]
  
  for (nvars in length(names.supplB)){
    
    print(nvars)
    out1=glm(as.formula(paste(COVID,'~',paste(names.supplB,collapse="+"),sep="")),data=train,family="binomial")
    out2=ranger(as.formula(paste(COVID,'~',paste(names.supplB,collapse="+"),sep="")),data=train,num.trees=1000)
    
    df=data.frame(iter=each,nvar=nvars,true=test[[COVID]],pred_glm=as.numeric(predict(out1,newdata=test,type="response")),pred_RF=as.numeric(predict(out2,data=test,type="response")$predictions))
    result=rbind(result,df)
  }
}
result=result[-1,]
#write.csv(result,"C:/.csv")

AUCs=result %>% split(.$nvar) %>% purrr::map(~ci.cvAUC(.$pred_glm,.$true,folds=.$iter))
AUCs2=result %>% split(.$nvar) %>% purrr::map(~ci.cvAUC(.$pred_RF,.$true,folds=.$iter))

AUC_df=rbind(bind_rows(AUCs %>% purrr::map(~data.frame(AUC=.$cvAUC,SE=.$se,lower=.$ci[1],upper=.$ci[2],level=.$confidence,Model="LR"))),
             bind_rows(AUCs2 %>% purrr::map(~data.frame(AUC=.$cvAUC,SE=.$se,lower=.$ci[1],upper=.$ci[2],level=.$confidence,Model="RF"))))
AUC_df
#         AUC          SE     lower     upper level Model
# 1 0.7054032 0.005327196 0.6949621 0.7158443  0.95    LR
# 2 0.6726852 0.005681347 0.6615499 0.6838204  0.95    RF

################### SupplC ####
#need complete cases, are a handful of observations missing for these predictors
data <- data[which(!is.na(data$fever)&!is.na(data$sob)),]

COVID <- "COVID"
out=ranger(as.formula(paste(COVID, '~' ,paste(names.supplC,collapse="+"),sep="")),data=data,importance="impurity",num.trees=1000)
imps=importance(out)

df_imps=data.frame(names=names(imps),var_red=as.numeric(imps)) %>% arrange(desc(var_red))
colnames(df_imps) <- c(paste("n=",dim(data)[1],"    ",colnames(df_imps)[1]),
                       colnames(df_imps)[2])
df_imps
# n= 1928      names   var_red
# 1            indexage 22.270487
# 2              ADI_UT 16.774334
# 3           race_temp  5.319935
# 4           prior_exp  4.730102
# 5              smoker  4.707467
# 6              travel  3.856690
# 7                 sob  3.424887
# 8            headache  3.307009
# 9             myalgia  3.238331
# 10             gender  3.202329
# 11              nasal  3.027810
# 12           lethargy  2.723203
# 13        sore_throat  2.711674
# 14              fever  2.439772
# 15           comorbid  2.309846
# 16          hc_worker  2.091206
# 17           diarrhea  1.855407
# 18              cough  1.841688
# 19       nausea_vomit  1.288514

nvars_opts=c(1:dim(df_imps)[1])
result=data.frame(iter=NA,nvar=NA,true=NA,pred_glm=NA,pred_RF=NA)
data$index=1:dim(data)[1]

test_record <- NA
train_record <- NA

for (each in 1:100){
  print(each)
  train=data %>% sample_frac(.80,replace=F)
  
  test=data[-which(data$index %in% train$index),]
  
  for (nvars in length(names.supplC)){
    
    print(nvars)
    out1=glm(as.formula(paste(COVID,'~',paste(names.supplC,collapse="+"),sep="")),data=train,family="binomial")
    out2=ranger(as.formula(paste(COVID,'~',paste(names.supplC,collapse="+"),sep="")),data=train,num.trees=1000)
    
    df=data.frame(iter=each,nvar=nvars,true=test[[COVID]],pred_glm=as.numeric(predict(out1,newdata=test,type="response")),pred_RF=as.numeric(predict(out2,data=test,type="response")$predictions))
    result=rbind(result,df)
  }
}
result=result[-1,]
#write.csv(result,"C:/.csv")

AUCs=result %>% split(.$nvar) %>% purrr::map(~ci.cvAUC(.$pred_glm,.$true,folds=.$iter))
AUCs2=result %>% split(.$nvar) %>% purrr::map(~ci.cvAUC(.$pred_RF,.$true,folds=.$iter))

AUC_df=rbind(bind_rows(AUCs %>% purrr::map(~data.frame(AUC=.$cvAUC,SE=.$se,lower=.$ci[1],upper=.$ci[2],level=.$confidence,Model="LR"))),
             bind_rows(AUCs2 %>% purrr::map(~data.frame(AUC=.$cvAUC,SE=.$se,lower=.$ci[1],upper=.$ci[2],level=.$confidence,Model="RF"))))
AUC_df
#         AUC          SE     lower     upper level Model
# 1 0.7224096 0.004776301 0.7130482 0.7317709  0.95    LR
# 2 0.6959231 0.005059139 0.6860074 0.7058389  0.95    RF

