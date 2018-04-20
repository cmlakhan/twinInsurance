#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)



library(broom)
library(tidyr)
library(dplyr)
library(magrittr)
library(boot)
library(readr)

#library(glm)
library(binomTools)

library(gtools)

library(polycor)


library(lme4)
library(HLMdiag)
library(reshape)

library(psych)
library(mvtnorm)
library(fMultivar)


library(foreach)
library(doParallel)


WD <- getwd()
setwd(WD)

library(SqlServerJtds)

cn <- connect.sql.server()



options(dplyr.width = Inf)


library(psych)

column_selection <- c('phewas_code','genderPairType','individual_idT1','individual_idT2','GenderT1','GenderT2','ageAverageT1','ageAverageT2', 
                      'PheWAS_Indicator_T1','PheWAS_Indicator_T2','monthEnrollmentT1','monthEnrollmentT2',
                      'avgMonthlyVisitsT1','avgMonthlyVisitsT2','numComorbidT1', 'numComorbidT2', 'ZipCodeT1', 'ZipCodeT2')

T1_column <- c('phewas_code','individual_idT1','GenderT1','SubscriberIdT1','ageAverageT1',
               'PheWAS_Indicator_T1', 'genderConfig','extraSSIDT1','monthEnrollmentT1','avgMonthlyVisitsT1','numComorbidT1','aqi_group','temp_group','income_group' )
T2_column <- c('phewas_code','individual_idT2','GenderT2','SubscriberIdT2','ageAverageT2', 
               'PheWAS_Indicator_T2', 'genderConfig','extraSSIDT2','monthEnrollmentT2','avgMonthlyVisitsT2','numComorbidT2','aqi_group','temp_group','income_group')

names_RF <- c('phewas_code','MemberId','Gender', 'SubscriberId', 'MemberBirthYear',  
              'Indicator', 'genderConfig', 'extraSSID', 'monthEnrollment', 'avgMonthlyVisits', 'numComorbid','aqi_group','temp_group','income_group')


regression_column <- c('phewas_code', 'Indicator', 'MemberBirthYear', 'Gender',  
                       'SubscriberId','genderConfig', 'extraSSID','monthEnrollment', 'avgMonthlyVisits', 'numComorbid','aqi_group','temp_group','income_group')


lmer_formula <- as.formula(Indicator ~  (1|SubscriberId) + (1|extraSSID) + MemberBirthYear + as.factor(Gender) + 
                             monthEnrollment + monthEnrollment^2 + (1|aqi_group) + (1|temp_group) + (1|income_group) )


twinOS_tetrachoric <- function(n11,n10,n01,n00) {
  
  #  the input to compute polychoric correlation for twins in a contingency table:
  #             Disease | !Disease
  #------------------------
  #  Disease |   ConcordantD      |     DisconcordantD/2
  # !Disease |  DisconcordantD/2       |     ConcordantND
  #
  if(n11 == 0) {n11 <- 0.5; } ## add in a constant for zero cells
  if(n10 == 0) {n10 <- 0.5; } ## this will be 0.5 when input to the function
  if(n01 == 0) {n01 <- 0.5; }
  
  #n01 <- (n01+n10)/2
  #n10 <- n01
  
  X <- matrix(c(n11,n01,n10,n00), nrow=2, byrow=T)
  tetraChorRho <- polychor(X, ML=T)
  return(tetraChorRho)
}

twinOS_tetrachoric_SE_kendall_stewart <- function(n11,n10,n01,n00) {
  if(n11 == 0) {n11 <- 0.5; } ## add in a constant for zero cells
  if(n10 == 0) {n10 <- 1; }
  if(n01 == 0) {n01 <- 1; }
  n <- n00 + n11 + n10 + n01
  

  p0. <- (n00 + n01)/n
  p.0 <- (n00 + n10)/n
  #Threshold t1 and t2
  t1 <- qnorm(p0.)
  t2 <- qnorm(p.0)
  
  r_ML <- twinOS_tetrachoric(n11,n10,n01,n00)
  Vr <- as.numeric(n^2*(1/n00 + 1/n01 + 1/n10 + 1/n11)*(dnorm2d(t1, t2, rho = r_ML))^2)
  se= sqrt(1/Vr)
  return(se)
}

heritability_from_tetrachor <- function(tetraChorRhoSameSex, tetraChorRhoOppositeSex,p) {
  # p ~ 0.5; this is from the estimation of # of MZ in same sex twins
  # h2 = 2(rho_ss - rho_os) / p
  # c2 = (rho_os * (1 + p) - rho_ss) / p
  #hsq <- 4 * (tetraChorRhoSameSex - tetraChorRhoOppositeSex)
  
  hsq <- 2*(tetraChorRhoSameSex - tetraChorRhoOppositeSex) / p
  
  return(hsq)
}

shared_environment <- function(tetraChorRhoSameSex, tetraChorRhoOppositeSex,p) {
  

  csq = (tetraChorRhoOppositeSex * (1 + p) - tetraChorRhoSameSex) / p
  return(csq)
}

twin_heritability_SE <- function(se.rOS, se.rSS,p) {
  ## compute the SE of the H2
  #SE(h2) = 4 * sqrt( SE(rSS)^2 + SE(rOS)^2 )
  se_h <- (2/p) * sqrt(se.rOS^2 + se.rSS^2)
  return(se_h)
}

twin_environment_SE <- function(se.rOS, se.rSS,p) {
  ## compute the SE of the E
  # E = 2 * ((1.5 * rOS) - rSS)
  # SE(E) = 2 * sqrt( (1.5^2) * SE(rOS)^2 + SE(rOS)^2 ))
  #se_e <- 2 * sqrt(2.25 * (se.rOS^2) + (se.rSS^2))
  
  se_e <- (1/p)*sqrt((1+p)^2*(se.rOS^2) + (se.rSS^2))
  return(se_e)
}

getPheWAS_data <- function(phewascode,cn,column_selection){
  phewas_code <- paste0("select * from cl313.dbo.twinPair_PheWAS_binary_indicator_single_pair a where a.phewas_code = '", phewascode, "'" )
  sqlResults <- dbGetQuery(cn, phewas_code)
  sqlResults <- sqlResults[,column_selection]
  sqlResults <- sqlResults[complete.cases(sqlResults),]
  sqlResults$genderConfig <- ifelse(sqlResults$genderPairType == 'MF', 'OS', 'SS')
  
  return(sqlResults)
}

mixData <- function(sqlResults_P1){
  twinPairsSS <- sqlResults_P1[sqlResults_P1$genderPairType != 'MF',]
  twinPairsOS <- sqlResults_P1[sqlResults_P1$genderPairType == 'MF',]
  index <- sample(nrow(twinPairsOS),size=nrow(twinPairsOS)/2,replace=FALSE)
  twinPairsOS1 <- twinPairsOS[index,]
  twinPairsOS2 <- twinPairsOS[-index,]
  twinPairsOS2_copies <- twinPairsOS2
  twinPairsOS2_copies$GenderT1 <- twinPairsOS2$GenderT2
  twinPairsOS2_copies$GenderT2 <- twinPairsOS2$GenderT1
  twinPairsOS2_copies$PheWAS_Indicator_T1 <- twinPairsOS2$PheWAS_Indicator_T2
  twinPairsOS2_copies$PheWAS_Indicator_T2 <- twinPairsOS2$PheWAS_Indicator_T1
  twinPairsOS2_copies$MemberIdT1 <- twinPairsOS2$MemberIdT2
  twinPairsOS2_copies$MemberIdT2 <- twinPairsOS2$MemberIdT1
  twinPairsOS2_copies$monthEnrollmentT1 <- twinPairsOS2$monthEnrollmentT2
  twinPairsOS2_copies$monthEnrollmentT2 <- twinPairsOS2$monthEnrollmentT1
  twinPairs <-rbind(twinPairsSS,twinPairsOS1,twinPairsOS2_copies)
  return(twinPairs)
}

prepareDatasets_RF <- function(sqlData, T1_column, T2_column, names_RF, regression_column){
  mixData <- sqlData
  mixData <- mixData %>% mutate(SubscriberIdT1=1:nrow(.)) %>% mutate(SubscriberIdT2=1:nrow(.)) %>%
    mutate(individual_idT1 = paste0(SubscriberIdT1,'_1')) %>%
    mutate(individual_idT2 = paste0(SubscriberIdT1,'_2'))
    
    
  mixData$genderConfig <- ifelse(mixData$genderPairType == 'MF', 'OS', 'SS')
  randomID <- mixData %>%  select(individual_idT1,individual_idT2, genderConfig)
  SSData <- subset(randomID, genderConfig == 'SS')
  SSData$extraSSIDT1 <- 1:nrow(SSData)
  SSData$extraSSIDT2 <- 1:nrow(SSData)
  OSData <- subset(randomID, genderConfig == 'OS')
  OSData$extraSSIDT1 <-1:nrow(OSData)
  OSData$extraSSIDT2 <- 1:nrow(OSData)
  OSData$extraSSIDT1 <- OSData$extraSSIDT1 + 800000
  OSData$extraSSIDT2 <- OSData$extraSSIDT2 + 1600000
  bothData <- bind_rows(SSData,OSData) %>% select(individual_idT1,individual_idT2,extraSSIDT1,extraSSIDT2)
  mixData <- mixData %>% inner_join(.,bothData, by=c("individual_idT1","individual_idT2"))
  twinPair1 <- mixData %>% select_(.,.dots=T1_column)
  twinPair2 <- mixData %>% select_(.,.dots=T2_column)
  names(twinPair1) <- names_RF
  names(twinPair2) <- names_RF
  bothPairs2 <- bind_rows(twinPair1,twinPair2)
  regression_all_b1_b2_rf_ssextra <- bothPairs2[,regression_column]
  regression_all_b1_b2_rf_ssextra <- regression_all_b1_b2_rf_ssextra %>% mutate(SubscriberId = as.factor(SubscriberId)) %>% mutate(extraSSID = as.factor(extraSSID))
  return(regression_all_b1_b2_rf_ssextra)
}

mixDataCount <- function(sqlData,phewas_code){
  
  

  mixData <- sqlData
  mixData <- mixData %>% mutate(SubscriberIdT1=1:nrow(.)) %>% mutate(SubscriberIdT2=1:nrow(.))
  mixData$genderConfig <- ifelse(mixData$genderPairType == 'MF', 'OS', 'SS')
  phewas <- mixData$phewas_code[1]
  mixData_T1_0_count <- length(which(mixData$PheWAS_Indicator_T1 =='0'))
  mixData_T1_0 <- data.frame(phewas, mixData_T1_0_count)
  names(mixData_T1_0) <- c('phewas_code','T1_0')
  mixData_T1_1_count <- length(which(mixData$PheWAS_Indicator_T1 =='1'))
  mixData_T1_1 <- data.frame(phewas, mixData_T1_1_count)
  names(mixData_T1_1) <- c('phewas_code','T1_1')
  mixData_T2_0_count <- length(which(mixData$PheWAS_Indicator_T2 =='0'))
  mixData_T2_0 <- data.frame(phewas, mixData_T2_0_count)
  names(mixData_T2_0) <- c('phewas_code','T2_0')
  mixData_T2_1_count <- length(which(mixData$PheWAS_Indicator_T2 =='1'))
  mixData_T2_1 <- data.frame(phewas, mixData_T2_1_count)
  names(mixData_T2_1) <- c('phewas_code','T2_1')
  mixData_count <- mixData_T2_0 %>%  left_join(.,mixData_T2_1, by="phewas_code") %>% left_join(.,mixData_T1_1, by="phewas_code") %>% left_join(.,mixData_T1_0, by="phewas_code")
  mixData_count[is.na(mixData_count)] <- 0
  mixData_count$total <- mixData_count$T2_1 + mixData_count$T2_0 + mixData_count$T1_1 + mixData_count$T1_0
  mixData_count$prevalence <-  (mixData_count$T2_1 + mixData_count$T1_1 ) / mixData_count$total
  mixData_count$Evary01 <- mixData_count$prevalence * (1-mixData_count$prevalence)
  mixData_count$T <- qnorm(1-mixData_count$prevalence)
  mixData_count$z_height <- dnorm(mixData_count$T)
  mixData_count$i <- mixData_count$z_height/mixData_count$prevalence
  return(mixData_count)
  
  
  
  
  
  
}

createMixData <- function(sqlData){
  mixData <- sqlData
  mixData <- mixData %>% mutate(SubscriberIdT1=1:nrow(.)) %>% mutate(SubscriberIdT2=1:nrow(.))
  mixData$genderConfig <- ifelse(mixData$genderPairType == 'MF', 'OS', 'SS')
  phewas <- mixData$phewas_code[1]
  
mixData_T1_0_count <- length(which(mixData$PheWAS_Indicator_T1 =='0'))
mixData_T1_0 <- data.frame(phewas, mixData_T1_0_count)
names(mixData_T1_0) <- c('phewas_code','T1_0')

mixData_T1_1_count <- length(which(mixData$PheWAS_Indicator_T1 =='1'))
mixData_T1_1 <- data.frame(phewas, mixData_T1_1_count)
names(mixData_T1_1) <- c('phewas_code','T1_1')

mixData_T2_0_count <- length(which(mixData$PheWAS_Indicator_T2 =='0'))
mixData_T2_0 <- data.frame(phewas, mixData_T2_0_count)
names(mixData_T2_0) <- c('phewas_code','T2_0')

mixData_T2_1_count <- length(which(mixData$PheWAS_Indicator_T2 =='1'))
mixData_T2_1 <- data.frame(phewas, mixData_T2_1_count)
names(mixData_T2_1) <- c('phewas_code','T2_1')

mixData_count <- mixData_T2_0 %>%  left_join(.,mixData_T2_1, by="phewas_code") %>% left_join(.,mixData_T1_1, by="phewas_code") %>% left_join(.,mixData_T1_0, by="phewas_code")

mixData_count[is.na(mixData_count)] <- 0

mixData_count$total <- mixData_count$T2_1 + mixData_count$T2_0 + mixData_count$T1_1 + mixData_count$T2_0

mixData_count$prevalence <- (mixData_count$T2_1 + mixData_count$T1_1) / mixData_count$total

mixData_count$Evary01 <- mixData_count$prevalence * (1-mixData_count$prevalence)

mixData_count$T <- qnorm(1-mixData_count$prevalence)

mixData_count$z_height <- dnorm(mixData_count$T)

mixData_count$i <- mixData_count$z_height/ifelse(mixData_count$prevalence == 0, .00000001,mixData_count$prevalence)

randomID <- mixData %>%  select(individual_idT1,individual_idT2, genderConfig)

SSData <- subset(randomID, genderConfig == 'SS')
SSData$extraSSIDT1 <- 1:nrow(SSData)
SSData$extraSSIDT2 <- 1:nrow(SSData)

OSData <- subset(randomID, genderConfig == 'OS')
OSData$extraSSIDT1 <-1:nrow(OSData)
OSData$extraSSIDT2 <- 1:nrow(OSData)
OSData$extraSSIDT1 <- OSData$extraSSIDT1 + 80000
OSData$extraSSIDT2 <- OSData$extraSSIDT2 + 160000

bothData <- bind_rows(SSData,OSData) %>% select(individual_idT1,individual_idT2,extraSSIDT1,extraSSIDT2)
mixData <- mixData %>% inner_join(.,bothData, by=c("individual_idT1","individual_idT2"))
return(mixData)
}

h2_c2_calc <- function(mixData){
  
  phewas_code <- mixData %>% select(phewas_code) %>% unique()
  
  
  
  SS_0_0 <- length(which(mixData$genderConfig == 'SS' & mixData$PheWAS_Indicator_T1 == '0' & mixData$PheWAS_Indicator_T2 == '0'))
  SS_0_1 <- length(which(mixData$genderConfig == 'SS' & mixData$PheWAS_Indicator_T1 == '0' & mixData$PheWAS_Indicator_T2 == '1'))
  SS_1_0 <- length(which(mixData$genderConfig == 'SS' & mixData$PheWAS_Indicator_T1 == '1' & mixData$PheWAS_Indicator_T2 == '0'))
  SS_1_1 <- length(which(mixData$genderConfig == 'SS' & mixData$PheWAS_Indicator_T1 == '1' & mixData$PheWAS_Indicator_T2 == '1'))
  

  OS_0_0 <- length(which(mixData$genderConfig == 'OS' & mixData$PheWAS_Indicator_T1 == '0' & mixData$PheWAS_Indicator_T2 == '0'))
  OS_0_1 <- length(which(mixData$genderConfig == 'OS' & mixData$PheWAS_Indicator_T1 == '0' & mixData$PheWAS_Indicator_T2 == '1'))
  OS_1_0 <- length(which(mixData$genderConfig == 'OS' & mixData$PheWAS_Indicator_T1 == '1' & mixData$PheWAS_Indicator_T2 == '0'))
  OS_1_1 <- length(which(mixData$genderConfig == 'OS' & mixData$PheWAS_Indicator_T1 == '1' & mixData$PheWAS_Indicator_T2 == '1'))
  
  
    

    
  data <- cbind(phewas_code,SS_0_0,SS_0_1,SS_1_0,SS_1_1,OS_0_0,OS_0_1,OS_1_0,OS_1_1) %>%
    mutate(r_SS=twinOS_tetrachoric(SS_1_1, SS_1_0, SS_0_1, SS_0_0), 
           r_OS=twinOS_tetrachoric(OS_1_1, OS_1_0, OS_0_1, OS_0_0)) %>%
    rowwise() %>%
    mutate(r_SS_SE=twinOS_tetrachoric_SE_kendall_stewart(SS_1_1, SS_1_0, SS_0_1, SS_0_0), r_OS_SE=twinOS_tetrachoric_SE_kendall_stewart(OS_1_1, OS_1_0, OS_0_1, OS_0_0)) %>%
    rowwise() %>%
    mutate(h2=heritability_from_tetrachor(r_SS, r_OS,p), c2=shared_environment(r_SS, r_OS,p)) %>%
    rowwise() %>%
    mutate(h2_SE=twin_heritability_SE(r_OS_SE, r_SS_SE,p), c2_SE=twin_environment_SE(r_OS_SE, r_SS_SE,p))
  
  
  
  return(data)  
  
}

randomEffects <- function(dataAnalysis, lmer_formula){
      fit <- lmer(lmer_formula, data=dataAnalysis)
      return(fit)
}

fit_h2_c2 <- function(d,K,p, i, formula_lmer){
  

  fit <- lmer(formula_lmer, data = d)
  
  
  tidyDF_boot <- tidy(fit)
  
  analysisDF <- tidyDF_boot[,c('term','estimate')] %>% spread(term, estimate) %>% mutate(varpair=`sd_(Intercept).SubscriberId` * `sd_(Intercept).SubscriberId`) %>% 
    mutate(varpair=`sd_(Intercept).SubscriberId` * `sd_(Intercept).SubscriberId`) %>%
    mutate(varextraSS=`sd_(Intercept).extraSSID` * `sd_(Intercept).extraSSID`) %>%
    mutate(varincome=`sd_(Intercept).income_group` * `sd_(Intercept).income_group`) %>%
    mutate(varaqi=`sd_(Intercept).aqi_group` * `sd_(Intercept).aqi_group`) %>%
    mutate(vartemp=`sd_(Intercept).temp_group` * `sd_(Intercept).temp_group`) %>%
    mutate(varres=`sd_Observation.Residual` * `sd_Observation.Residual`) %>%
    mutate(vartot=varpair + varextraSS + varres + varincome + varaqi + vartemp) %>%
    mutate(varos=varpair) %>%
    mutate(varss=varpair + varextraSS) %>%                      
    mutate(rss01=varss/vartot) %>%
    mutate(ros01=varos/vartot) %>%
    mutate(rincome01=varincome/vartot) %>%
    mutate(raqi01=varaqi/vartot) %>%
    mutate(rtemp01=vartemp/vartot) %>%
    rowwise() %>%
    mutate(Eb1= K + (varss/K)) %>%
    rowwise() %>%
    mutate(Eb2= K + (varos/K)) %>%
    rowwise() %>%
    mutate(Ebincome= K + (varincome/K)) %>%
    rowwise() %>%
    mutate(Ebaqi= K + (varaqi/K)) %>%
    rowwise() %>%
    mutate(Ebtemp= K + (vartemp/K)) %>%
    rowwise() %>%
    mutate(K=K) %>%
    mutate(Tss= qnorm(1 - Eb1)) %>%
    mutate(Tos= qnorm(1 - Eb2)) %>% 
    mutate(Tincome= qnorm(1 - Ebincome)) %>% 
    mutate(Taqi= qnorm(1 - Ebaqi)) %>% 
    mutate(Ttemp= qnorm(1 - Ebtemp)) %>% 
    mutate(rliabss= ((T - Tss)*sqrt(1 - (T^2 - Tss^2)*(1 - T/i))) / (i + Tss^2*(i - T))) %>% 
    mutate(rliabos= ((T - Tos)*sqrt(1 - (T^2 - Tos^2)*(1 - T/i))) / (i + Tos^2*(i - T))) %>%
    mutate(rliabincome= ((T - Tincome)*sqrt(1 - (T^2 - Tincome^2)*(1 - T/i))) / (i + Tincome^2*(i - T))) %>%
    mutate(rliabaqi= ((T - Taqi)*sqrt(1 - (T^2 - Taqi^2)*(1 - T/i))) / (i + Taqi^2*(i - T))) %>%
    mutate(rliabtemp= ((T - Ttemp)*sqrt(1 - (T^2 - Ttemp^2)*(1 - T/i))) / (i + Ttemp^2*(i - T))) %>%
    rowwise() %>%
    mutate(h2_liab=heritability_from_tetrachor(rliabss,rliabos,p)) %>%
    rowwise() %>%
    mutate(c2_liab=shared_environment(rliabss,rliabos,p)) %>%
    mutate(rres01=varres/vartot) %>%
    mutate(Ebres=K + (varres/K)) %>%
    mutate(Tres= qnorm(1 - Ebres)) %>%
    mutate(rliabres= ((T - Tres)*sqrt(1 - (T^2 - Tres^2)*(1 - T/i))) / (i + Tres^2*(i - T)))   
  
  return(analysisDF)
  
}

boostrap_SE <- function(dataAnalysis,iterations,K,p, i,formula_lmer, T1_column, T2_column, names_RF, regression_column ){
  bootIter <- 0
  
  
  while(bootIter < iterations)
  {

    dataAnalysis_1 <- dataAnalysis %>% sample_n(., nrow(.),replace=TRUE)
    dataAnalysis2 <- prepareDatasets_RF(dataAnalysis_1,T1_column, T2_column, names_RF, regression_column)
    d <- dataAnalysis2
    

    num_1 <- nrow(d[d$Indicator==1,])
      data_h2_c2_row <- fit_h2_c2(d,K,p, i, formula_lmer)
      data_h2_c2_row <- data_h2_c2_row %>% dplyr::select(rliabss_SE=rliabss,rliabos_SE=rliabos,
                                                  h2_liab_SE=h2_liab,c2_liab_SE=c2_liab,
                                                  lambda_OS_SE=lambda_OS,lambda_SS_SE=lambda_SS,
                                                  rliabincome_SE=rliabincome,rliabaqi_SE=rliabaqi,
                                                  rliabtemp_SE=rliabtemp, rliabres_SE=rliabres)
      if(bootIter == 0)
      {
        data_h2_c2_all <- data_h2_c2_row
      }
      data_h2_c2_all <- data_h2_c2_all %>% bind_rows(.,data_h2_c2_row)
      bootIter <- bootIter+1
     
  }
  return(data_h2_c2_all)
}




flattenTidyFile <- function(tidyFrame,K,i,T){
  tidyFrame
  
  analysisDF <- tidyFrame[,c('term','estimate')] %>% spread(term, estimate) %>% mutate(varpair=`sd_(Intercept).SubscriberId` * `sd_(Intercept).SubscriberId`) %>% 
    mutate(varpair=`sd_(Intercept).SubscriberId` * `sd_(Intercept).SubscriberId`) %>%
    mutate(varextraSS=`sd_(Intercept).extraSSID` * `sd_(Intercept).extraSSID`) %>%
    mutate(varincome=`sd_(Intercept).income_group` * `sd_(Intercept).income_group`) %>%
    mutate(varaqi=`sd_(Intercept).aqi_group` * `sd_(Intercept).aqi_group`) %>%
    mutate(vartemp=`sd_(Intercept).temp_group` * `sd_(Intercept).temp_group`) %>%
    mutate(varres=`sd_Observation.Residual` * `sd_Observation.Residual`) %>%
    mutate(vartot=varpair + varextraSS + varres + varincome + varaqi + vartemp) %>%
    mutate(varos=varpair) %>%
    mutate(varss=varpair + varextraSS) %>%                      
    mutate(rss01=varss/vartot) %>%
    mutate(ros01=varos/vartot) %>%
    mutate(rincome01=varincome/vartot) %>%
    mutate(raqi01=varaqi/vartot) %>%
    mutate(rtemp01=vartemp/vartot) %>%
    rowwise() %>%
    mutate(Eb1= K + (varss/K)) %>%
    rowwise() %>%
    mutate(Eb2= K + (varos/K)) %>%
    rowwise() %>%
    mutate(Ebincome= K + (varincome/K)) %>%
    rowwise() %>%
    mutate(Ebaqi= K + (varaqi/K)) %>%
    rowwise() %>%
    mutate(Ebtemp= K + (vartemp/K)) %>%
    rowwise() %>%
    mutate(K=K) %>%
    mutate(Tss= qnorm(1 - Eb1)) %>%
    mutate(Tos= qnorm(1 - Eb2)) %>% 
    mutate(Tincome= qnorm(1 - Ebincome)) %>% 
    mutate(Taqi= qnorm(1 - Ebaqi)) %>% 
    mutate(Ttemp= qnorm(1 - Ebtemp)) %>% 
    mutate(rliabss= ((T - Tss)*sqrt(1 - (T^2 - Tss^2)*(1 - T/i))) / (i + Tss^2*(i - T))) %>% 
    mutate(rliabos= ((T - Tos)*sqrt(1 - (T^2 - Tos^2)*(1 - T/i))) / (i + Tos^2*(i - T))) %>%
    mutate(rliabincome= ((T - Tincome)*sqrt(1 - (T^2 - Tincome^2)*(1 - T/i))) / (i + Tincome^2*(i - T))) %>%
    mutate(rliabaqi= ((T - Taqi)*sqrt(1 - (T^2 - Taqi^2)*(1 - T/i))) / (i + Taqi^2*(i - T))) %>%
    mutate(rliabtemp= ((T - Ttemp)*sqrt(1 - (T^2 - Ttemp^2)*(1 - T/i))) / (i + Ttemp^2*(i - T))) %>%
    rowwise() %>%
    mutate(h2_liab=heritability_from_tetrachor(rliabss,rliabos,p)) %>%
    rowwise() %>%
    mutate(c2_liab=shared_environment(rliabss,rliabos,p)) %>%
    mutate(rres01=varres/vartot) %>%
    mutate(Ebres=K + (varres/K)) %>%
    mutate(Tres= qnorm(1 - Ebres)) %>%
    mutate(rliabres= ((T - Tres)*sqrt(1 - (T^2 - Tres^2)*(1 - T/i))) / (i + Tres^2*(i - T))) 
    
    
  
  
  stdError <- tidyFrame[,c('term','std.error')]
  stdError <- stdError[complete.cases(stdError),]

  stdError$term <- paste0(stdError$term, '_SE')
  stdError <- stdError %>% spread(term,std.error)
  
  analysisDF <- analysisDF %>%  bind_cols(.,stdError)
  
  return(analysisDF)
  
  
}



FF <-	17919
MF <- 	20642
MM <- 	17835


p_mz <- 1-(2*MF/(MM+FF+MF))

p_ss <- (FF+MM)/(MM+FF+MF)

p <- p_mz/p_ss

p <- 0.4223689

phewas_code <- as.character(args[1])



df <- getPheWAS_data(phewas_code,cn, column_selection)



zipcodeData_corr <- readRDS('zipcode_twins_aqi_income_temp.rds')
 
zipcodeData_corr %>% summarise(income_temp=cor(family_income,avg_temp,method = "spearman"))
zipcodeData_corr %>% summarise(income_temp=cor(family_income,avg_aqi,method = "spearman"))
zipcodeData_corr %>% summarise(income_temp=cor(avg_temp,avg_aqi,method = "spearman"))


m <- zipcodeData_corr %>% group_by(income_group,temp_group,aqi_group) %>% tally() %>% ungroup()

saveRDS(m,file='count_zipcode_income_temp_aqi.rds')


zipcodeData <- readRDS('twinPair_environmentalData.rds') %>%
  mutate(aqi_group = as.factor(aqi_50_group)) %>%
  mutate(temp_group = as.factor(avgtemp_50_group)) %>%
  mutate(income_group = as.factor(family_income_group))



df <- df %>% inner_join(.,zipcodeData, by=c('individual_idT1','individual_idT2'))



df <- mixData(df)





df2 <- prepareDatasets_RF(df, T1_column, T2_column, names_RF, regression_column)

df3 <- mixDataCount(df, phewas_code)
h2_c2 <- h2_c2_calc(df)





fit <- randomEffects(df2, lmer_formula)
tidyDF <- tidy(fit)
glanceDF <- glance(fit)



K <- df3$prevalence
i <- df3$i

T <- df3$T



results <- boostrap_SE(df, 500, K, p, i, lmer_formula, T1_column, T2_column, names_RF, regression_column)




SE_Estimates <- data.frame(rliabss_SE_se=sd(results$rliabss_SE), 
                           rliabos_SE_se=sd(results$rliabos_SE), 
                           h2_liab_SE_se=sd(results$h2_liab_SE), 
                           c2_liab_SE_se=sd(results$c2_liab_SE),
                           lambda_OS_SE_se=sd(results$lambda_OS_SE),
                           lambda_SS_SE_se=sd(results$lambda_SS_SE),
                           rliabincome_SE_se=sd(results$rliabincome_SE),
                           rliabaqi_SE_se=sd(results$rliabaqi_SE),
                           rliabtemp_SE_se=sd(results$rliabtemp_SE), 
                           rliabres_SE_se=sd(results$rliabres_SE))


lmmModelEstimates <- flattenTidyFile(tidyDF,K,i,T) %>% bind_cols(.,SE_Estimates)
allData <- h2_c2 %>% bind_cols(., lmmModelEstimates)
allData <- allData %>% inner_join(., df3,by='phewas_code')
allData %>% select(phewas_code, h2,c2,h2_liab, c2_liab)



allData %>% select(h2_SE, c2_SE,T2_0,T2_1,T1_1,T1_0)



s <- df %>% dplyr::filter(.,PheWAS_Indicator_T1 == 1) %>% select(individual_id=individual_idT1, phewas_code)
sss <- df %>% dplyr::filter(.,PheWAS_Indicator_T2 == 1) %>% select(individual_id=individual_idT2, phewas_code)
ssss <- s %>% bind_rows(.,sss)



dd <- readRDS(file='costPheWAS.rds')

ddssss <-  ssss %>% inner_join(.,dd, by=c('individual_id','phewas_code') ) %>% 
  select(individual_id, costPheWAS_month,ratioTotal_PheWAS) %>% unique() %>% 
  summarise(avgCostPheWAS_month = mean(costPheWAS_month, na.rm='TRUE'),
            avgCostPheWAS_ratio = mean(ratioTotal_PheWAS,na.rm='TRUE'))


allData <- allData %>% bind_cols(.,ddssss)






filename <- paste0('phewas_zipcode_',phewas_code,'.csv')

write.table(allData,filename, row.names = FALSE, col.names = FALSE)

saveRDS(names(allData),file='names_binary_zipcode.rds')




