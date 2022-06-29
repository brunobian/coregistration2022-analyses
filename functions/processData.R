processData <- function(dataSet) {
  
  resamplear <- 1
  
  # dataSet <- dataSet[nchar(as.character(dataSet$pal))>2,]
  dataSet$pal <- gsub(' +','',dataSet$pal)
  dataSet$pal <- gsub('[!¡.,:;?¿"]','',dataSet$pal)
  dataSet$pal <- tolower(dataSet$pal)
  dataSet$pal <- as.factor(dataSet$pal)

  dataSet$suj_id <- as.factor(dataSet$suj_id)

  # Delete bad epochs 
  # dataSet <- dataSet[dataSet$filter == 1,] # Elimino los malos trials
  # dataSet <- dataSet[!is.na(dataSet$filter),] # Elimino los malos trials
  
  # Genero variable dummy para oraicones comunes y proverbios
  # dataSet$tipo_Dummy_common  = dataSet$tipo
  # dataSet$tipo_Dummy_common[dataSet$tipo_Dummy_common == 0] = NA
  # dataSet$tipo_Dummy_proverb = abs(1 - dataSet$tipo)
  # dataSet$tipo_Dummy_proverb[dataSet$tipo_Dummy_proverb == 0] = NA
  # 
  dataSet$sntType <- as.factor(dataSet$sntType-.5)
  dataSet$sntTypePrePost <- as.factor(dataSet$sntTypePrePost)
  
  # Transformo y centro variables
  # dataSet$freq   <- log(dataSet$freq+1)
  # dataSet$fixDur <- log(dataSet$fixDur)
  
  dataSet$freq    <- scale(dataSet$freq,   center = TRUE, scale = F)
  dataSet$pred    <- scale(dataSet$pred,   center = TRUE, scale = F)
  dataSet$pos     <- scale(dataSet$pos,    center = TRUE, scale = F)
  dataSet$lngth   <- scale(dataSet$lngth,  center = TRUE, scale = F)
  dataSet$fixDur  <- scale(dataSet$fixDur, center = TRUE, scale = F)
  dataSet$saccDur <- scale(dataSet$saccDur,center = TRUE, scale = F)
  
  #summary(dataSet)
  # Uso las dummy de arriba para generar las pred de prov y oraciones comunes por separado
  # De esta forma me ahorro el contraste y puedo ver N400 de ambos efectos po sep.
  # dataSet$pred_common  <- dataSet$pred * dataSet$tipo_Dummy_common
  # dataSet$pred_proverb <- dataSet$pred * dataSet$tipo_Dummy_proverb
  # 

  if (resamplear){
    inds <- resample(dataSet)
    dataSet <- dataSet[inds,]
  }
  
  return(dataSet)
}