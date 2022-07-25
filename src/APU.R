#!/usr/bin/env Rscript
#library(tictoc)
library(philentropy)
Diffusion<- function(nDat, vP, vRN, mW_norm, dAlpha){
  vG0 <- rep(0, nDat)
  vG0[vP] <- 1
  vG0[vRN] <- -length(vP)/length(vRN)
  dCheck <- 1
  vG <- vG0
  vG_tot <- c()
  while(dCheck>10^(-6)){
    vG1 <- (1-dAlpha)*t(mW_norm)%*%vG +dAlpha*vG0
    vG_tot <- cbind(vG_tot,vG1)
    dCheck <- sum(abs(vG1-vG))
    vG <- vG1
    print(dCheck)
  }
  DiffOut<- list()
  DiffOut$vG<- vG
  DiffOut$vG_tot<- vG_tot
  return(DiffOut)
}

dir_out<- '../NeDBIT_features/'
#dir_save<- '../data/APU_scores/'
setwd(dir_out)
disease_folder <- list.files(dir_out, full.names = TRUE)
disease_name <-list.files(dir_out, full.names = FALSE)
nDisease <- length(disease_folder)
vNeigh <- seq(2, 100, 2)
for(dis in 1:nDisease){
  #dir_dis <- disease_folder[dis]
  #setwd(dir_dis)
  #out_folder <- list.files(dir_dis, full.names = FALSE)
  #nOut <- length(out_folder)
  #for (sFold in 1:nOut){
    #print(paste(disease_name[dis], 'Fold', toString(sFold)))
    sNameData <- disease_name[dis]
    dat_ <- read.csv(sNameData, header = TRUE)
    vDim<- dim(dat_)
    NetMeasure <- dat_[,5:vDim[2]]
    vDim_meas <- dim(NetMeasure)
    vP <- which(dat_$class==1)
    nPos <- length(vP)
    mDat <- NetMeasure
    mDat <- as.matrix(mDat)
    vNames <- colnames(mDat)
    vDim <- dim(mDat)
    mData_norm <- matrix(0, vDim[1], vDim[2])
    for (i in 1:vDim[2]){
      mData_norm[, i] <- (mDat[,i]-min(mDat[,i]))/(max(mDat[,i])-min(mDat[,i]))
    }
    mEuclDist <-distance(mData_norm, method = "euclidean")
    dMin <- min(mEuclDist)
    dMax <- max(mEuclDist)
    nDat <- vDim_meas[1]
    mW_paola <- matrix(1, nDat, nDat)- (mEuclDist - matrix(dMin, nDat, nDat))/(dMax-dMin)
    rm(mEuclDist)
    vMeanFeat <- rep(0, vDim[2])
    vMeanFeat_P <- rep(0, vDim[2])
    for (i in 1:vDim[2]){
      vMeanFeat_P[i]<- mean(mData_norm[vP, i])
      vMeanFeat[i]<- mean(mData_norm[, i])
    }
    vDist <- rep(0, nDat)
    for (i in 1: nDat){
      vDist[i] <- sqrt(sum((mData_norm[i,]-vMeanFeat_P)^(2)))
    }
    hist(vDist)
    dq_neg <- quantile(vDist, 0.8)
    vRN <- which(vDist>= dq_neg)
    
    mW_sign <- mW_paola
    rm(mW_paola)
    vD <- rep(0, nDat)
    mW_norm <- matrix(0, nDat, nDat)
    for (i in 1: nDat){
      vD[i] <- sum(mW_sign[i, ])
      mW_norm[i,] <- mW_sign[i,]/vD[i]
    }
    #######--Start diffusion process--------------
    dAlpha <- 0.8
    DiffOut <- Diffusion(nDat, vP, vRN, mW_norm, dAlpha)
    #To save ordering through the score
    dat_prova <- data.frame(dat_$name, DiffOut$vG)
    colnames(dat_prova) <- c('name', 'out')
    sName <- paste(sNameData, '_Score.csv', sep='')
    #sName <- paste(dir_save, sNameData, '_Score.csv', sep='')
    write.csv(dat_prova, sName, row.names = FALSE)
  #}
}



