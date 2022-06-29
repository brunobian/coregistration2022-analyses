setwd('/media/brunobian/DATABRUNO/Co-registro_2018/Data')
library('lme4')
library(ggplot2)
library(dplyr)
library(emmeans)
library(svglite)
source('/media/brunobian/DATABRUNO/corregistro/Analysis_bb/2020/my_functions/logit.R')

std   <- function(x) sd(x, na.rm = TRUE)/(sqrt(sum(!is.na(x))))
mean0 <- function(x,rm0) mean(x[x>50], na.rm=rm0)

df <- read.csv('DATA_EM.csv',header=TRUE,sep=",")
which(is.na(df$EM_subject))
df <- subset(df, EM_FPRT < 1000)
df <- subset(df, EM_FFD  < 1000)

ind <- df$catgram == "Q"
which(is.na(df$EM_subject))
for (cat in c("A","N","V","R")){ind <- ind + (cat==df$catgram)}
df <- df[ind==1,]

df$logFPRT <- log(df$EM_FPRT+1)
df$logFreq <- log(df$freq)

df <- df[df$pos > 1 & df$pos < 12,] 

df$logFPRT <- scale(df$logFPRT,center=TRUE,scale=FALSE)
df$logFreq <- scale(df$logFreq,center=TRUE,scale=FALSE)
df$pred    <- scale(df$pred,center=TRUE,scale=FALSE)
df$length  <- scale(1/df$length,center=TRUE,scale=FALSE)
df$SntType <- as.factor(df$SntType)

# Analisis por pos abs ####

# % FFD	(First Fixation Duration) Duration of the first fixation on a position if (and only if) the fixation was progressive. Zero otherwise.
# % SFD	(Single Fixation Duration) Duration of the fixation on a position if it was the *only* fixation on this region, i.e. if no subsequent fixation on this position followed. Zero if there were several fixations on this region.
# % FPRT  (First Pass Reading Time, Gaze Duration) Sum of all first-pass fixation durations on a region before *any* other region is fixated. (What exactly constitutes a first pass is determined by the parameter 'regressiveFirstPass'.)
# % RBRT  (Right Bounded Reading Time) Sum of all first-pass fixation durations on a position before another position to the *right* is fixated. (What exactly constitutes a first pass is determined by the parameter 'regressiveFirstPass'.)
# % RPD	(Regression Path Duration, Go-Past Duration) Sum of all first-pass fixation durations on a position n and all preceding positions in the time period between the first fixation on n and the first fixation on anything to the right of n.
# % CRPD  (Cumulative Regression Path Duration) The CRPD of position n is the total amount of time a participant spent reading the sentence until reaching region n+1. It is the sum of the RPDs of all preceding regions and the RPD of the current region.
# % TFT	(Total Fixation Time) Sum of all fixation durations on a region.
# % RRT	(Re-reading Time) Sum of all second-pass fixation durations. (RRT = TFT - FPRT)
# % RRTR  (ReReading Time Regressive) Sum of all second-pass fixation durations on a position that occured *after* a fixation on a region further to the right. (RRTR = TFT - RBRT)
# % RRTP  (ReReading Time Progressive) Sum of all second-pass fixation durations on a position that took place *before* a fixation on a region further to the right. (RRTP = RBRT - FPRT)
# % FFP	(First Fixation Progressive) 0 if material downstream was viewed before the first fixation on this position, 1 otherwise.
# % TRC	(Total Regression Count) Total number of regressions from this position.
# % RBRC  (Right-Bounded Regression Count) Number of regressions from this position given *before* any region further to the right has been fixated.                    
# % TRI   (Total Regressions In) Number of regressions TO this position.     
# % skipped 0 if word was fixated, 1 if skipped
# % refixed 0 if word was not re-fixated, 1 if refixed

campos <- c("EM_FPRT","EM_FFD","EM_CRPD","EM_TRC","EM_TRI","EM_NFT","EM_NFFP",
            "EM_RRT","EM_TFT","EM_skipped","EM_refixed")
campos <- c("EM_skipped")
for (f in campos){
  print(f)
  if (f=="EM_skipped"){
    df1 <- df %>%
      group_by(pos,SntType,EM_subject) %>%
      summarise(EM = mean(.data[[f]], na.rm=TRUE)) #interp(~mean(var, na.rm = T), var = as.name(f))) %>%
  }else{
    df1 <- df %>%
            group_by(pos,SntType,EM_subject) %>%
            summarise(EM = mean0(.data[[f]], rm0=TRUE)) #interp(~mean(var, na.rm = T), var = as.name(f))) %>%
  }
  
  df2 <- df1 %>% 
            group_by(pos,SntType) %>%
            dplyr::summarize(m = mean(EM, na.rm=TRUE),
                             e = std(EM))
  
  ggplot(df2, aes(x=pos, y=m, color=SntType)) + 
    geom_line() +
    geom_point()+
    geom_errorbar(aes(ymin=m-e, ymax=m+e), width=.2) +
    ylab(f)
  
  # basePath = "/media/brunobian/DATABRUNO/Co-registro_2018/Data/figs/EM_Analysis/"
  # filename = paste0(basePath, f, 'vsPos.svg')
  # ggsave(filename)
}

# pred no tiene sentido hacerlo por suj, lo hago afuera del for, eliminando duplicados
dfPred <- subset(df,select=c(pred,EM_roi,pos,SntType,EM_trial))
dfPred <- dfPred[!duplicated(dfPred),]

df2 <- dfPred %>% 
  group_by(pos,SntType) %>%
  dplyr::summarize(m = mean(pred, na.rm=TRUE),
                   e = std(pred))
ggplot(df2, aes(x=pos, y=m, color=SntType)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=m-e, ymax=m+e), width=.2) +
  ylab(f)
# basePath = "/media/brunobian/DATABRUNO/Co-registro_2018/Data/figs/EM_Analysis/"
# filename = paste0(basePath, 'PredvsPos.svg')
# ggsave(filename)

# Stats complex ####
# Pred
dfPred$posFac <- as.factor(dfPred$pos)
contrasts(dfPred$posFac)  <-  contr.treatment(unique(dfPred$posFac), contrasts = FALSE)
dfPred$SntType <- as.factor(dfPred$SntType)
contrasts(dfPred$SntType)  <- contr.treatment(unique(dfPred$SntType), base = 2)
  m4 <- lmer(pred ~  posFac*SntType + (1|EM_trial), data=dfPred)
  emm4 <- emmeans(m4, pairwise ~ SntType | posFac , type = "response")
  plot(emm4$contrasts)

df$posFac <- as.factor(df$pos)
contrasts(df$posFac)  <-  contr.treatment(unique(df$posFac), contrasts = FALSE)
df$SntType <- as.factor(df$SntType)
contrasts(df$SntType)  <- contr.treatment(unique(df$SntType), base = 2)

  # Gaze Duration
  dfFPRT <- subset(df, EM_FPRT < 1000 & EM_FPRT > 0)
  m1 <- lmer(EM_FPRT ~  posFac*SntType + (1|EM_subject)+ (1|EM_trial), data=dfFPRT)
  emm1.1 <- emmeans(m1, pairwise ~ SntType | posFac , type = "response") #entre tipos
  plot(emm1.1$contrasts)
  emm1.2 <- emmeans(m1, eff ~ posFac | SntType , type = "response") # dentro de tipos
  plot(emm1.2$contrasts)
  
  # Skips
  m2 <- lmer(EM_skipped ~  posFac*SntType + (1|EM_subject)+ (1|EM_trial), data=df)
  emm2 <- emmeans(m2, pairwise ~ SntType | posFac , type = "response")
  plot(emm2$contrasts)
  
  # Refix
  m3 <- lmer(EM_refixed ~  posFac*SntType + (1|EM_subject)+ (1|EM_trial), data=df)
  emm3 <- emmeans(m3, pairwise ~ SntType | posFac , type = "response")
  plot(emm3$contrasts)
  


# Analisis por pos rel ####
df_prov <- df[df$SntType == 0,]
df_prov <- df_prov[df_prov$posRelRP > -4 & df_prov$posRelRP < 6,]

df2 <- df_prov %>%
  group_by(posRelRP,EM_subject,SntType) %>%
  dplyr::summarize(EM = mean(EM_FPRT, na.rm=TRUE))

df3 <- df2 %>%
  group_by(posRelRP,SntType) %>%
  dplyr::summarize(m = mean(EM, na.rm=TRUE),
                   e = std(EM))

ggplot(df3, aes(x=posRelRP, y=m, color=SntType)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=m-e, ymax=m+e), width=.2) 


df_prov$posRelRP <- as.factor(df_prov$posRelRP)
contrasts(df_prov$posRelRP) <- contr.helmert(unique(df_prov$posRelRP))

m3 <- lmer(EM_FPRT ~  posRelRP + (1|EM_subject) + (1|EM_trial), 
           data=df_prov)#, contrasts = list(posRelRP = contrasts(df_prov$posRelRP)))
summary(m3, corr=F)

# Analisis LMM full ####
m4 <- lmer(EM_FPRT ~  freq + pos + SntType + pred + predNext + (1|EM_subject) + (1|EM_trial), data=df)
summary(m4, corr=F)

# Analisis skips 2
f = 'EM_skipped'
df <- df[df$pos < 11,] # & df$gramFilt==0
df$pos[df$pos==2|df$pos==3] = 2.5
df$pos[df$pos==7|df$pos==8] = 7.5
df$pos[df$pos==9|df$pos==10] = 9.5

df2 <- df %>%
  group_by(pos,SntType,EM_subject) %>%
  summarise(EM = mean(.data[[f]], na.rm = T)) %>% #interp(~mean(var, na.rm = T), var = as.name(f))) %>%
  group_by(pos,SntType) %>%
  dplyr::summarize(m = mean(EM, na.rm=TRUE),
                   e = std(EM))

ggplot(df2, aes(x=pos, y=m, color=SntType)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=m-e, ymax=m+e), width=.2) +
  ylab(f)

basePath = "/media/brunobian/DATABRUNO/Co-registro_2018/Data/figs/EM_Analysis/" 
filename = paste0(basePath,'SkipvsPosAgrupado.svg')
ggsave(filename)
