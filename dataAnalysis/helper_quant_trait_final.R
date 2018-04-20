#!/usr/bin/env Rscript

#args = commandArgs(trailingOnly=TRUE)



rm(list=ls())


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

library(lme4)

WD <- getwd()
setwd(WD)

library(SqlServerJtds)

cn <- connect.sql.server()


options(dplyr.width = Inf)


library(broom)
library(tidyr)
library(dplyr)
library(magrittr)
library(boot)



library(psych)

column_selection <- c('phewas_code','genderPairType','individual_idT1','individual_idT2','GenderT1','GenderT2','ageAverageT1','ageAverageT2', 
                      'PheWAS_Indicator_T1','PheWAS_Indicator_T2', 'ZipCodeT1', 
                      'ZipCodeT2','ageTest_T1','ageTest_T2')

T1_column <- c('phewas_code','individual_idT1','GenderT1','SubscriberIdT1','ageAverageT1',
               'PheWAS_Indicator_T1', 'genderConfig','extraSSIDT1','ageTest_T1','monthEnrollmentT1','monthEnrollmentT2' )



T2_column <- c('phewas_code','individual_idT2','GenderT2','SubscriberIdT2','ageAverageT2', 
               'PheWAS_Indicator_T2', 'genderConfig','extraSSIDT2','ageTest_T2','monthEnrollmentT1','monthEnrollmentT2')

names_RF <- c('phewas_code','MemberId','Gender', 'SubscriberId', 'MemberBirthYear',  
              'Indicator', 'genderConfig', 'extraSSID',  'ageTest','monthEnrollment')


regression_column <- c('phewas_code','MemberId', 'Indicator', 'MemberBirthYear', 'Gender',  
                       'SubscriberId','genderConfig', 'extraSSID', 'ageTest','monthEnrollment')


lmer_formula <- as.formula(Indicator ~  (1|SubscriberId) + (1|extraSSID)  + as.factor(Gender) + ageTest  )


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
  phewas_code <- paste0("select * from cl313.dbo.TwinPair_lab_results_all a where a.phewas_code = '", phewascode, "'" )
  sqlResults <- dbGetQuery(cn, phewas_code)
   sqlResults <- sqlResults[,column_selection]
    sqlResults <- sqlResults[complete.cases(sqlResults),]
    
    sqlResults$genderConfig <- ifelse(sqlResults$genderPairType == 'MF', 'OS', 'SS')
    
    
    
    return(sqlResults)
}



rankTransPheno <-  function(pheno,para_c){
    pheno<-qnorm((rank(pheno)-para_c)/(length(pheno)-2*para_c+1))
    return(pheno)
  }



getPheWAS_data_revised <- function(phewascode,cn,column_selection){

  month_enrollment <- paste0("select a.individual_idT1, a.individual_idT2, a.monthEnrollmentT1, a.monthEnrollmentT2 FROM cl313.[dbo].[twinPair_Full_info_single_pair] a" )
  sqlResults_enrollment <- dbGetQuery(cn, month_enrollment)
  
  
    phewas_code <- paste0("select a.* from cl313.dbo.TwinPair_lab_results_all a where a.phewas_code = '", phewascode, "'" )
  sqlResults <- dbGetQuery(cn, phewas_code)
  #zipCodeACS <- readRDS('/home/cl313/aetna_geo/zipCodeACS.rds')
  #sqlResults <- merge(sqlResults,zipCodeACS, by.x= "ZipCodeT1", by.y="zipcode", all.x=TRUE )
  #sqlResults <- sqlResults[,column_selection]
  
  sqlResults <- sqlResults %>% group_by(individual_idT1, individual_idT2) %>% sample_n(1) %>% ungroup()
  
  
  sqlResults <- sqlResults[complete.cases(sqlResults),]
  
  sqlResults_T1 <- sqlResults %>% select(individual_id = individual_idT1, PheWAS_Indicator = PheWAS_Indicator_T1, DateServiceStarted=DateServiceStarted_T1, ageTest=ageTest_T1, Gender=GenderT1 )
  sqlResults_T2 <- sqlResults %>% select(individual_id=individual_idT2, PheWAS_Indicator=PheWAS_Indicator_T2 ,DateServiceStarted = DateServiceStarted_T2,ageTest=ageTest_T2, Gender=GenderT2)
  
  sqlResults_both <- sqlResults_T1 %>% rbind(.,sqlResults_T2) 

  
  sqlResults_both <- sqlResults_both %>% cbind(.,data.frame(PheWAS_Indicator_normal=rankTransPheno(sqlResults_both$PheWAS_Indicator,para_c=3/8)))
  
  sqlResults_both <- sqlResults_both %>% mutate(PheWAS_Indicator=PheWAS_Indicator_normal) %>% select(-PheWAS_Indicator_normal)
  

  sqlResults_both <- sqlResults_both %>% select(-ageTest, -Gender)
  
  names(sqlResults_both) <- c('individual_id','PheWAS_Indicator','DateServiceStarted')
  
  
  sqlResults <- sqlResults %>% inner_join(.,sqlResults_both, by=c('individual_idT1' = 'individual_id', 'DateServiceStarted_T1'='DateServiceStarted') )
    
  sqlResults <- sqlResults %>% mutate(PheWAS_Indicator_T1=PheWAS_Indicator) %>% select(-PheWAS_Indicator)
  
  
  
  sqlResults <- sqlResults %>% inner_join(.,sqlResults_both, by=c('individual_idT2' = 'individual_id', 'DateServiceStarted_T2'='DateServiceStarted') )
  
  sqlResults <- sqlResults %>% mutate(PheWAS_Indicator_T2=PheWAS_Indicator) %>% select(-PheWAS_Indicator) %>% unique()
  
  
  sqlResults <- sqlResults %>% group_by(individual_idT1, individual_idT2) %>% sample_n(1) %>% ungroup()
  
  sqlResults <- sqlResults[, column_selection]
  
  sqlResults$genderConfig <- ifelse(sqlResults$genderPairType == 'MF', 'OS', 'SS')
  
  
  
  sqlResults <- sqlResults %>% inner_join(.,sqlResults_enrollment, by=c('individual_idT1','individual_idT2'))
  
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
  OSData$extraSSIDT1 <- OSData$extraSSIDT1 + 80000
  OSData$extraSSIDT2 <- OSData$extraSSIDT2 + 160000
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

mixDataCount <- function(sqlData){
  
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
  mixData_count$prevalence <- (mixData_count$T2_1 + mixData_count$T1_1) / mixData_count$total
  mixData_count$Evary01 <- mixData_count$prevalence * (1-mixData_count$prevalence)
  mixData_count$T <- qnorm(1-mixData_count$prevalence)
  mixData_count$z_height <- dnorm(mixData_count$T)
  mixData_count$i <- mixData_count$z_height/ifelse(mixData_count$prevalence == 0, .00000001,mixData_count$prevalence)
  mixData_count$total <- nrow(mixData)
  
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

h2_c2_calc <- function(mixData, p){
  
  phewas_code <- mixData %>% select(phewas_code) %>% unique()
  
  mixData$genderConfig <- ifelse(mixData$genderPairType == 'MF', 'OS', 'SS')
  
  data <- mixData
  
  #data[is.na(data)] <- 0
  
  
  dataSS <- data[data$genderConfig == 'SS',] %>% 
    do(glance(cor.test(.$PheWAS_Indicator_T1,.$PheWAS_Indicator_T2))) %>% 
    mutate(r_SS_SE=(sqrt((1 - estimate^2)/parameter))) %>% 
    select(r_SS=estimate, r_SS_SE)
  
  
  
  dataOS <- data[data$genderConfig == 'OS',] %>%
    do(glance(cor.test(.$PheWAS_Indicator_T1,.$PheWAS_Indicator_T2))) %>% 
    mutate(r_OS_SE=(sqrt((1 - estimate^2)/parameter))) %>% 
    select(r_OS=estimate, r_OS_SE)
  
  dataBoth <- dataSS %>% bind_cols(.,dataOS) %>%
    mutate(h2=heritability_from_tetrachor(r_SS, r_OS,p), c2=shared_environment(r_SS, r_OS,p)) %>%
    mutate(h2_SE=twin_heritability_SE(r_OS_SE, r_SS_SE,p), c2_SE=twin_environment_SE(r_OS_SE, r_SS_SE,p))
  
  dataBoth$phewas_code <- phewas_code$phewas_code
  
  return(dataBoth)  
  
}

randomEffects <- function(dataAnalysis, lmer_formula){
      fit <- lmer(lmer_formula, data=dataAnalysis)
      return(fit)
}

fit_h2_c2 <- function(d,p, formula_lmer){
  

  fit <- lmer(formula_lmer, data = d)
  
  
  tidyDF_boot <- tidy(fit)
  
  analysisDF <- tidyDF_boot[,c('term','estimate')] %>% 
    spread(term, estimate) %>% 
    mutate(varpair=`sd_(Intercept).SubscriberId` * `sd_(Intercept).SubscriberId`) %>%
    mutate(varextraSS=`sd_(Intercept).extraSSID` * `sd_(Intercept).extraSSID`) %>%
    mutate(varres=`sd_Observation.Residual` * `sd_Observation.Residual`) %>%
    mutate(vartot=varpair + varextraSS + varres) %>%
    mutate(varos=varpair) %>%
    mutate(varss=varpair + varextraSS) %>%                      
    mutate(rss01=varss/vartot) %>%
    mutate(ros01=varos/vartot) %>%
    mutate(h2_liab=heritability_from_tetrachor(rss01,ros01,p)) %>%
    mutate(c2_liab=shared_environment(rss01,ros01,p)) %>%
    mutate(rres01=varres/vartot) %>%
    mutate(e2_liab=1-h2_liab-c2_liab)
    
    
  
  cor_h2_c2_df <- data.frame(analysisDF$rss01,analysisDF$ros01,analysisDF$h2_liab,analysisDF$c2_liab, analysisDF$rres01, analysisDF$e2_liab)
  names(cor_h2_c2_df) <-c('rliabss_SE','rliabos_SE','h2_liab_SE','c2_liab_SE', 'rres01_SE', 'e2_liab_SE')
  return(cor_h2_c2_df)
  
}

boostrap_SE <- function(dataAnalysis,iterations,p,formula_lmer, T1_column, T2_column, names_RF, regression_column ){
  bootIter <- 0

  while(bootIter < iterations)
  {
    #print(bootIter)
    dataAnalysis_1 <- dataAnalysis %>% sample_n(., nrow(.),replace=TRUE)
    dataAnalysis2 <- prepareDatasets_RF(dataAnalysis_1,T1_column, T2_column, names_RF, regression_column)
    d <- dataAnalysis2
    
    
    
      data_h2_c2_row <- fit_h2_c2(d,p, formula_lmer)
      
      if(bootIter ==0){
        data_h2_c2_all <- data_h2_c2_row
        bootIter <- bootIter+1
        
      }

      data_h2_c2_all <- data_h2_c2_all %>% bind_rows(.,data_h2_c2_row)
      bootIter <- bootIter+1
      
        }
  return(data_h2_c2_all)
  
  
  
}

flattenTidyFile <- function(tidyFrame,p){
  #tidyFrame
  
  analysisDF <- tidyFrame[,c('term','estimate')] %>% 
    spread(term, estimate) %>% 
    mutate(varpair=`sd_(Intercept).SubscriberId` * `sd_(Intercept).SubscriberId`) %>% 
    mutate(varextraSS=`sd_(Intercept).extraSSID` * `sd_(Intercept).extraSSID`) %>%
    mutate(varres=`sd_Observation.Residual` * `sd_Observation.Residual`) %>%
    mutate(vartot=varpair + varextraSS + varres) %>%
    mutate(varos=varpair) %>%
    mutate(varss=varpair + varextraSS) %>%                      
    mutate(rss01=varss/vartot) %>%
    mutate(ros01=varos/vartot) %>%
    mutate(h2_liab=heritability_from_tetrachor(rss01,ros01,p)) %>%
    mutate(c2_liab=shared_environment(rss01,ros01,p)) %>%
    mutate(rres01=varres/vartot) %>%
    mutate(e2_liab=1-h2_liab-c2_liab)
  
  
  stdError <- tidyFrame[,c('term','std.error')]
  stdError <- stdError[complete.cases(stdError),]
  
  stdError$term <- paste0(stdError$term, '_SE')
  stdError <- stdError %>% spread(term,std.error)
  
  analysisDF <- analysisDF %>%  bind_cols(.,stdError)
  
  return(analysisDF)
  
  
}

getNumTests <- function(phewas_code,cn){
  
  phewas_code <- paste0("select * from cl313.dbo.TwinIndividual_quant_traits_numTest a where a.phewas_code = '", phewas_code, "'" )
  sqlResults <- dbGetQuery(cn, phewas_code)
  return(sqlResults)
  
}





phewas_code <- paste0("select distinct phewas_code FROM cl313.dbo.twinPair_PheWAS_quantitative_single_highPairCount a ")
phewas_list <- dbGetQuery(cn, phewas_code)


for (i in 1:nrow(phewas_list) ){

  print(phewas_code)
  
  phewas_code <- phewas_list[i,c("phewas_code")]
  

df <- getPheWAS_data_revised(phewas_code,cn, column_selection)


#numTests_df <- getNumTests(phewas_code,cn) %>% select(-SubscriberId,-phewas_code,-Gender)

FF <-	df %>% dplyr::filter(genderPairType == 'FF') %>% dplyr::tally()
MF <- 	df %>% dplyr::filter(genderPairType == 'MF') %>% dplyr::tally()
MM <- 	df %>% dplyr::filter(genderPairType == 'MM') %>% dplyr::tally()


p_mz <- 1-(2*MF/(MM+FF+MF))

p_ss <- (FF+MM)/(MM+FF+MF)

p <- p_mz/p_ss

p <- p[1,]





#df <- mixData(df)


df2 <- prepareDatasets_RF(df, T1_column, T2_column, names_RF, regression_column) 



df3 <- mixDataCount(df)
df4 <- createMixData(df)




h2_c2 <- h2_c2_calc(df, p)




lmer_formula <- as.formula(Indicator ~  (1|SubscriberId) + (1|extraSSID) + as.factor(Gender) + ageTest + monthEnrollment + monthEnrollment^2 )



fit <- lmer(lmer_formula, data=df2)


tidyDF <- tidy(fit)
glanceDF <- glance(fit)

results <- boostrap_SE(df, 500, p, lmer_formula, T1_column, T2_column, names_RF, regression_column)



SE_Estimates <- results %>% summarise_each(., funs(se=sd(., na.rm=TRUE)))



lmmModelEstimates <- flattenTidyFile(tidyDF,p) %>% bind_cols(.,SE_Estimates)

allData <- h2_c2 %>% bind_cols(., lmmModelEstimates)

allData <- allData %>% inner_join(., df3,by='phewas_code')

allData <- allData %>% select(phewas_code, everything())

allData$p <- p





rawData <- df 

filename <- paste0('',phewas_code,'.csv')

write.table(allData,filename, row.names = FALSE, col.names = FALSE)

saveRDS(names(allData),file='/home/cl313/aetna_manuscript_analysis/names_quant.rds')




filename_raw <- paste0('',phewas_code,'.rds')

saveRDS(rawData,file=filename_raw)


}



