## collected plotting functions for project field analysis
library("HDInterval")







plot_PEf <- function(fit, sim, byGROUP, IO, c, nocontrol, plottofile) {
  # function to plot estimates/predictions of PEf, Protective efficacy in the field
  # 'fit' -- stanfit object
  # 'sim' -- simulation object, i.e. treatment and grouping information for all simulated PEf samples
  # 'byGROUP' -- grouping characteristic to be shown in plot
  # 'IO' -- 'IN'  or 'OUT' indicating indoor or outdoor data
  # 'c' -- real number that was used in the stan program to avoid the denominator in FEf being 0
  # 'nocontrol' -- logical: if TRUE: do not show contorl, if FALSE show control (as a check)
  # 'plottofile' -- character string of path to save plot to
  
  library(plyr)
  library(dplyr)
  library(reshape2)
  library(ggplot2)
  
  # prepare simulated data
  PEf_IO <- paste0("PEf_", IO) #value
  samplePEf_IO <- rstan::extract(fit, pars=PEf_IO) #extract
  samplePEf_IO<- as.data.frame(samplePEf_IO[[PEf_IO]]) #convert
  simdata_PEf <- melt(cbind(sim, t(samplePEf_IO)), id.vars = names(sim)) #add simulation information
  #name treatments, using factors
    simdata_PEf$TREAT <- as.factor(simdata_PEf$TREAT)
    simdata_PEf$TREAT <- plyr::revalue(simdata_PEf$TREAT, c("1"="control", "2"="pull", '3'='push', '4'='pushpull'))
  names(simdata_PEf)[names(simdata_PEf) == 'value'] <- 'PEf'  # rename variable
  names(simdata_PEf)[names(simdata_PEf) == 'variable'] <- 'iter'  # rename variable
  # remove control if nocontrol==TRUE
  if (nocontrol){
    simdata_PEf <- simdata_PEf[!simdata_PEf$TREAT=='control',]
  }
  
  
  # remove draws with PEf==NA (i.e. controls with zero count) 
  n_NA <- sum(is.na(simdata_PEf[,'PEf']))
  simdata_PEf_noNA <- simdata_PEf[complete.cases(simdata_PEf[,'PEf']),] 
  if (nrow(simdata_PEf) != n_NA + nrow(simdata_PEf_noNA)){stop('problem with complete.cases()')}
  
  # remove -Inf  HEURISTIC! NO IDEA WHY THERE ARE -Inf's!
  n_Inf <- sum(is.infinite(simdata_PEf_noNA[,'PEf']))
  simdata_PEf_noNA <- simdata_PEf_noNA[!is.infinite(simdata_PEf_noNA[,'PEf']),]
  
  if (!is.na(byGROUP)){
    
    #name byGROUP, using factors
    simdata_PEf_noNA[,byGROUP] <- as.factor(simdata_PEf_noNA[,byGROUP])
    simdata_PEf_noNA[,byGROUP] <- plyr::mapvalues(simdata_PEf_noNA[,byGROUP], from=levels(simdata_PEf_noNA[,byGROUP]), to=sapply(levels(simdata_PEf_noNA[,byGROUP]), FUN = function(s){paste0(byGROUP,' ',s)}))
    
    # rename variable byGROUP to generic 'GROUP' so that this stupid facet_grid can handle it!
    names(simdata_PEf_noNA)[names(simdata_PEf_noNA)==byGROUP] <- 'GROUP'
    
    # compute summary statistics
    stats_PEf_noNA <- ddply(simdata_PEf_noNA, c('GROUP', 'TREAT'), summarise, mean=mean(PEf), median=quantile(PEf, probs=0.5))
    
    # plot
    plot = ggplot(simdata_PEf_noNA, aes(x=PEf)) +
      geom_histogram(aes(y=stat(ncount)), position = 'identity', alpha=0.5, breaks=seq(-5,1, by=0.05)) +
      geom_vline(data = stats_PEf_noNA, aes(xintercept = mean), linetype="longdash", alpha=1, size=0.2) +
      geom_vline(data = stats_PEf_noNA, aes(xintercept = median), linetype="dotdash", alpha=1, size=0.2) +
      facet_grid(GROUP ~ TREAT) +
      coord_cartesian(xlim = c(-5,1)) + xlab(PEf_IO) + ylab('normalised counts')
    ggsave(plottofile, plot, height=12, width=18, dpi=300)
    
  }else{
    # prepare data
      groups <- names(sim)[!names(sim)=='TREAT'] # get grouping characteristic i.e. groups
      groups_new <- sapply(sim[,groups], max)
      simdata_PEf_noNA_pred <- simdata_PEf_noNA[simdata_PEf_noNA[,groups[1]]==groups_new[1]&simdata_PEf_noNA[,groups[2]]==groups_new[2],]
    
    # compute summary statistics
    stats_PEf_noNA_pred <- ddply(simdata_PEf_noNA_pred, 'TREAT', summarise, mean=mean(PEf), median=quantile(PEf, probs=0.5))
      
    # plot
    plot = ggplot(simdata_PEf_noNA_pred, aes(x=PEf)) +
      geom_histogram(aes(y=stat(ncount)), position = 'identity', alpha=0.5, breaks=seq(-5,1, by=0.05)) +
      geom_vline(data = stats_PEf_noNA_pred, aes(xintercept = mean), linetype="longdash", alpha=1, size=0.2) +
      geom_vline(data = stats_PEf_noNA_pred, aes(xintercept = median), linetype="dotdash", alpha=1, size=0.2) +
      facet_grid(~TREAT) +
      coord_cartesian(xlim = c(-5,1)) + xlab(PEf_IO) + ylab('normalised counts')
    ggsave(plottofile, plot, height=3, width=18, dpi=300)
  }
}




hdi_PEf <- function(fit, sim, IO, c, prob) {
  # function to compute highest density credible intervals of PEf, Protective efficacy in the field
  # 'fit' -- stanfit object
  # 'sim' -- simulation object, i.e. treatment and grouping information for all simulated PEf samples
  # 'IO' -- 'IN'  or 'OUT' indicating indoor or outdoor data
  # 'c' -- real number that was used in the stan program to avoid the denominator in FEf being 0
  # 'prob' -- is the probability for the credibility level
  
  library(plyr)
  library(dplyr)
  library(reshape2)
  
  # prepare simulated data
  PEf_IO <- paste0("PEf_", IO) #value
  samplePEf_IO <- rstan::extract(fit, pars=PEf_IO) #extract
  samplePEf_IO<- as.data.frame(samplePEf_IO[[PEf_IO]]) #convert
  simdata_PEf <- melt(cbind(sim, t(samplePEf_IO)), id.vars = names(sim)) #add simulation information
  #name treatments, using factors
    simdata_PEf$TREAT <- as.factor(simdata_PEf$TREAT)
    simdata_PEf$TREAT <- plyr::revalue(simdata_PEf$TREAT, c("1"="control", "2"="pull", '3'='push', '4'='pushpull'))
  names(simdata_PEf)[names(simdata_PEf) == 'value'] <- 'PEf'  # rename variable
  names(simdata_PEf)[names(simdata_PEf) == 'variable'] <- 'iter'  # rename variable
  simdata_PEf <- simdata_PEf[!simdata_PEf$TREAT=='control',]

  
  # remove draws with PEf==NA (i.e. controls with zero count) 
  n_NA <- sum(is.na(simdata_PEf[,'PEf']))
  simdata_PEf_noNA <- simdata_PEf[complete.cases(simdata_PEf[,'PEf']),] 
  if (nrow(simdata_PEf) != n_NA + nrow(simdata_PEf_noNA)){stop('problem with complete.cases()')}
  
  # remove -Inf  HEURISTIC! NO IDEA WHY THERE ARE -Inf's!
  n_Inf <- sum(is.infinite(simdata_PEf_noNA[,'PEf']))
  simdata_PEf_noNA <- simdata_PEf_noNA[!is.infinite(simdata_PEf_noNA[,'PEf']),]
  
  # prepare data
  groups <- names(sim)[!names(sim)=='TREAT'] # get grouping characteristic i.e. groups
  groups_new <- sapply(sim[,groups], max)
  simdata_PEf_noNA_pred <- simdata_PEf_noNA[simdata_PEf_noNA[,groups[1]]==groups_new[1] & simdata_PEf_noNA[,groups[2]]==groups_new[2],]
  simdata_PEf_noNA_pred <- as.data.frame(simdata_PEf_noNA_pred)
  simdata_PEf_noNA_pred <- simdata_PEf_noNA_pred[,c('TREAT','PEf')]
  
  # wrapper around hdi() as a workaround for this stupid ddply to accept my variable prob (which has nothing to do with the data supplied)
  hdi_prob <<- function(...){
    hdint <- hdi(..., credMass = get('prob', envir=parent.env(environment())))
    return(as.numeric(hdint))
  }
  
  hdi_prob_left <<- function(...){
    return(hdi_prob(...)[1])
  }
  
  hdi_prob_right <<- function(...){
    return(hdi_prob(...)[2])
  }


  # compute
  hdi <- ddply(simdata_PEf_noNA_pred, 'TREAT', summarise, hdi_left=hdi_prob(PEf)[1], hdi_right=hdi_prob(PEf)[2]) # ask Lars how to pass 'prob' directly to credMass in hdi(), solution above with defining functions in global environment is probably very risky

  return(hdi)
}




median_PEf <- function(fit, sim, IO, c) {
  # function to compute highest density credible intervals of PEf, Protective efficacy in the field
  # 'fit' -- stanfit object
  # 'sim' -- simulation object, i.e. treatment and grouping information for all simulated PEf samples
  # 'IO' -- 'IN'  or 'OUT' indicating indoor or outdoor data
  # 'c' -- real number that was used in the stan program to avoid the denominator in FEf being 0
  # 'prob' -- is the probability for the credibility level
  
  library(plyr)
  library(dplyr)
  library(reshape2)
  
  # prepare simulated data
  PEf_IO <- paste0("PEf_", IO) #value
  samplePEf_IO <- rstan::extract(fit, pars=PEf_IO) #extract
  samplePEf_IO<- as.data.frame(samplePEf_IO[[PEf_IO]]) #convert
  simdata_PEf <- melt(cbind(sim, t(samplePEf_IO)), id.vars = names(sim)) #add simulation information
  #name treatments, using factors
  simdata_PEf$TREAT <- as.factor(simdata_PEf$TREAT)
  simdata_PEf$TREAT <- plyr::revalue(simdata_PEf$TREAT, c("1"="control", "2"="pull", '3'='push', '4'='pushpull'))
  names(simdata_PEf)[names(simdata_PEf) == 'value'] <- 'PEf'  # rename variable
  names(simdata_PEf)[names(simdata_PEf) == 'variable'] <- 'iter'  # rename variable
  simdata_PEf <- simdata_PEf[!simdata_PEf$TREAT=='control',]
  
  
  
  # remove draws with PEf==NA (i.e. controls with zero count) 
  n_NA <- sum(is.na(simdata_PEf[,'PEf']))
  simdata_PEf_noNA <- simdata_PEf[complete.cases(simdata_PEf[,'PEf']),] 
  if (nrow(simdata_PEf) != n_NA + nrow(simdata_PEf_noNA)){stop('problem with complete.cases()')}
  
  # remove -Inf  HEURISTIC! NO IDEA WHY THERE ARE -Inf's!
  n_Inf <- sum(is.infinite(simdata_PEf_noNA[,'PEf']))
  simdata_PEf_noNA <- simdata_PEf_noNA[!is.infinite(simdata_PEf_noNA[,'PEf']),]
  
  # prepare data
  groups <- names(sim)[!names(sim)=='TREAT'] # get grouping characteristic i.e. groups
  groups_new <- sapply(sim[,groups], max)
  simdata_PEf_noNA_pred <- simdata_PEf_noNA[simdata_PEf_noNA[,groups[1]]==groups_new[1] & simdata_PEf_noNA[,groups[2]]==groups_new[2],]
  simdata_PEf_noNA_pred <- as.data.frame(simdata_PEf_noNA_pred)
  simdata_PEf_noNA_pred <- simdata_PEf_noNA_pred[,c('TREAT','PEf')]
  
  # compute 
  median <- ddply(simdata_PEf_noNA_pred, 'TREAT', summarise, median = median(PEf))
  
  return(median)
}




hdi_mu <- function(fit, sim, IO, prob) {
  # function to compute highest density credible intervals of PEf, Protective efficacy in the field
  # 'fit' -- stanfit object
  # 'sim' -- simulation object, i.e. treatment and grouping information for all simulated PEf samples
  # 'IO' -- 'IN'  or 'OUT' indicating indoor or outdoor data
  # 'prob' -- is the probability for the credibility level
  
  library(plyr)
  library(dplyr)
  library(reshape2)
  
  # prepare simulated data
  zmeanIO <- paste0('zmean', IO)
  sample_zmean <- rstan::extract(fit, pars=zmeanIO) #extract
  sample_zmean<- as.data.frame(sample_zmean[[zmeanIO]]) #convert
  simdata_zmean <- melt(cbind(sim, t(sample_zmean)), id.vars = names(sim)) #add simulation information
  #name treatments, using factors
  simdata_zmean$TREAT <- as.factor(simdata_zmean$TREAT)
  simdata_zmean$TREAT <- plyr::revalue(simdata_zmean$TREAT, c("1"="control", "2"="pull", '3'='push', '4'='pushpull'))
  names(simdata_zmean)[names(simdata_zmean) == 'value'] <- 'zmean'  # rename variable
  names(simdata_zmean)[names(simdata_zmean) == 'variable'] <- 'iter'  # rename variable
  
  # select predictive simulated data for further, previously unseen groups with respect to all available grouping characteristics
  groups <- names(sim)[!names(sim)=='TREAT'] # get grouping characteristic i.e. groups
  groups_new <- sapply(sim[,groups], max)
  simdata_zmean_pred <- simdata_zmean[simdata_zmean[,groups[1]]==groups_new[1] & simdata_zmean[,groups[2]]==groups_new[2],]
  simdata_zmean_pred <- as.data.frame(simdata_zmean_pred)
  simdata_zmean_pred <- simdata_zmean_pred[,c('TREAT','zmean')]
  
  # wrapper around hdi() as a workaround for this stupid ddply to accept my variable prob (which has nothing to do with the data supplied)
  hdi_prob <<- function(...){
    hdint <- hdi(..., credMass = get('prob', envir=parent.env(environment())))
    return(as.numeric(hdint))
  }
  
  hdi_prob_left <<- function(...){
    return(hdi_prob(...)[1])
  }
  
  hdi_prob_right <<- function(...){
    return(hdi_prob(...)[2])
  }
  
  
  # compute 
  hdi <- ddply(simdata_zmean_pred, 'TREAT', summarise, hdi_left=hdi_prob(zmean)[1], hdi_right=hdi_prob(zmean)[2])
  
  return(hdi)
}






hdi_psi <- function(fit, IO, prob) {
  # function to compute highest density credible intervals of PEf, Protective efficacy in the field
  # 'fit' -- stanfit object
  # 'sim' -- simulation object, i.e. treatment and grouping information for all simulated PEf samples
  # 'IO' -- 'IN'  or 'OUT' indicating indoor or outdoor data
  # 'prob' -- is the probability for the credibility level
  
  library(plyr)
  library(dplyr)
  library(reshape2)
  
  # prepare simulated data
  logpsiIO <- paste0('logpsi', IO)
  sample_logpsi <- rstan::extract(fit, pars=logpsiIO) #extract
  sample_logpsi<- as.data.frame(sample_logpsi[[logpsiIO]]) #convert
  names(sample_logpsi) <- c("control", "pull", 'push', 'pushpull')

  sample_logpsi_melt <- melt(sample_logpsi)
  names(sample_logpsi_melt)[names(sample_logpsi_melt)=='variable'] <- 'TREAT'
  
  # wrapper around hdi() as a workaround for this stupid ddply to accept my variable prob (which has nothing to do with the data supplied)
  hdi_prob <<- function(...){
    hdint <- hdi(..., credMass = get('prob', envir=parent.env(environment())))
    return(as.numeric(hdint))
  }
  
  hdi_prob_left <<- function(...){
    return(hdi_prob(...)[1])
  }
  
  hdi_prob_right <<- function(...){
    return(hdi_prob(...)[2])
  }
  
  
  # compute 
  hdi <- ddply(sample_logpsi_melt, 'TREAT', summarise, hdi_left=hdi_prob(value)[1], hdi_right=hdi_prob(value)[2])
  
  return(hdi)
}






plot_ppc <- function(file, species, suna, fit, sim, byGROUP, IO, pool, statistics, mean_dist, plottofile) {
  # function to plot estimates/predictions of PEf, Protective efficacy in the field
  # 'file' -- csv file with observed data
  # 'species' -- mosquito species
  # 'suna' -- switch for plotting suna trap catches instead of HLC or CDC
  # 'fit' -- stanfit object
  # 'sim' -- simulation object, i.e. treatment and grouping information for all simulated PEf samples
  # 'byGROUP' -- grouping characteristic to be shown in plot
  # 'IO' -- 'IN'  or 'OUT' indicating indoor or outdoor data
  # 'pool' -- switch for plotting additional last row with pooling all group-specific predictions and poolign all true data, in order to compare to predictin for further, previously unseen group
  # 'statistics' -- switch for plotting means and medians as vertical lines
  # 'mean_dist' -- switch for plotting predictive distribution of mean; The mean of this distribution is the mean plotted if statistics=TRUE (functionallity below to check that commented out below)
  # 'plottofile' -- character string of path to save plot to
  
  
  library(plyr)
  library(dplyr)
  library(reshape2)
  library(ggplot2)
  
  # get species specific measures of observed data
  if (species=='gamb'){
      if (IO=='IN'){
        spec <- 'An..gambiae.s.l..Unfed'
      } else if (IO=='OUT'){
        spec <- 'HLC_gambiae'
      }
  } else if (species=='fune'){
      if (IO=='IN'){
        spec <- 'An..funestus.s.l.Unfed'
      } else if (IO=='OUT'){
        spec <- 'HLC_funestus'
      }    
  } else if (species =='cule'){
      if (IO=='IN'){
        spec <- 'Culex.spp.Unfed'
      } else if (IO=='OUT'){
        spec <- 'HLC_culex'
    }
  } else if (species =='mans'){
    if (IO=='IN'){
      spec <- 'Mansonia.Unfed'
    } else if (IO=='OUT'){
      spec <- 'HLC_mansonia'
    }
  }
    
  
  # prepare observed data
    # extract 
    if (IO=='IN'){
      data_trap <- read.csv(file=file, header=TRUE, sep=",")
      names(data_trap)[names(data_trap) == 'house_ID'] <- 'houseID' # rename variable to match with HLC data
      names(data_trap)[names(data_trap) == 'week_no'] <- 'week' # rename variable to match with HLC data
      if (suna) {data_trap <- data_trap[data_trap$trap=='suna',]} 
      else {data_trap <- data_trap[data_trap$trap=='cdc',]}
      data <- data_trap
    }else if (IO=='OUT'){
      data_HLC <- read.csv(file=file, header=TRUE, sep=",")
      data_HLC <- subset(data_HLC, select = -c(treatID,treat) )
      names(data_HLC)[names(data_HLC) == 'treat_label'] <- 'treatment' # rename variable to match with HLC data
      data <- data_HLC
    }
  
    #remove NAs
    data <- data[complete.cases(data[[spec]]),] 
    #rename variables
    names(data)[names(data) == "treatment"] <- "TREAT"
    names(data)[names(data) == "week"] <- "WEEK"
    names(data)[names(data) == "houseID"] <- "HOUSE"
    #create arrays encoding features
    xHOUSE <- as.numeric(as.factor(data$HOUSE))
    xTREAT <- dplyr::recode(data$TREAT, control=1, pull=2, push=3, push_pull=4)
    xWEEK <- as.numeric(data$WEEK)
  
    #process observed data
    realdata <- data[,c('HOUSE','TREAT', 'WEEK', spec)]
    names(realdata)[names(realdata) == spec] <- 'value'  # rename variable
    realdata$data <- "true" # label observed data as such
    if (!is.numeric(realdata[,byGROUP])){ 
      realdata[,byGROUP] <- as.numeric(as.factor(realdata[,byGROUP])) #STRANGE! NEEDED?
    }
    realdata$TREAT <- dplyr::recode(data$TREAT, control=1, pull=2, push=3, push_pull=4)
    # add pooled data to compare with prediction for further, previously unseen group
    realdata_pooled <- realdata
    realdata_pooled[,byGROUP] <- max(sim[,byGROUP])
    realdata <- rbind(realdata, realdata_pooled)
    # add pooled data to compare with pooling all group-specific predictions
    if (pool){
      realdata_pool <- realdata
      realdata_pool[,byGROUP] <- max(sim[,byGROUP])+1
      realdata <- rbind(realdata, realdata_pool)      
    }
    #
    realdata$iter <- NA
    realdata[,'mean'] <- NA
    maxcount <- max(realdata$value)
    
    
    # prepare simdata
    zIO <- paste0('z', IO)
    zmeanIO <- paste0('zmean', IO)
    simfit <- rstan::extract(fit, pars=c(zIO, zmeanIO))
    #
    simfit_IO<- as.data.frame(simfit[[zIO]])
    simdata <- melt(cbind(sim, t(simfit_IO)), id.vars = names(sim))
    names(simdata)[names(simdata) == 'variable'] <- 'iter'  # rename variable
    #
    simfit_zmeanIO<- as.data.frame(simfit[[zmeanIO]])
    simdata_zmeanIO<- melt(cbind(sim, t(simfit_zmeanIO)), id.vars = names(sim))
    names(simdata_zmeanIO)[names(simdata_zmeanIO) == 'value'] <- 'mean' # rename variable
    names(simdata_zmeanIO)[names(simdata_zmeanIO) == 'variable'] <- 'iter'  # rename variable
    # 
    simdata <- merge(simdata, simdata_zmeanIO, by=c(names(sim),"iter"))
    simdata$data <- 'simulated' # label sim data as 'simulated'
    # create '14th house' by pooling of simulated data for houses 1-12, to be compared to house 13
    if (pool){
      simdata_pool <- simdata[simdata[,byGROUP]!=max(sim[,byGROUP]),]
      simdata_pool[,byGROUP] <- max(sim[,byGROUP])+1
      simdata <- rbind(simdata, simdata_pool)
    }
    
    # combine realdata and simdata
    combdata <- rbind(simdata, realdata)
    
    #name treatments, using factors
    combdata$TREAT <- as.factor(combdata$TREAT)
    combdata$TREAT <- plyr::revalue(combdata$TREAT, c("1"="control", "2"="pull", '3'='push', '4'='pushpull'))
    
    #name byGROUP, using factors
    combdata[,byGROUP] <- as.factor(combdata[,byGROUP])
    combdata[,byGROUP] <- plyr::mapvalues(combdata[,byGROUP], from=levels(combdata[,byGROUP]), to=sapply(levels(combdata[,byGROUP]), FUN = function(s){paste0(byGROUP,' ',s)}))
    
    # rename variable byGROUP to generic 'GROUP' so that this stupid facet_grid can handle it!
    names(combdata)[names(combdata)==byGROUP] <- 'GROUP'
    
    # compute statistics to compare prediction to data with
    stats <- ddply(combdata, c('data', 'GROUP', 'TREAT'), summarise, mean=mean(value), median=quantile(value, probs=0.5))#,  quants=quantile(value, probs=c(0.2, 0.8)))
    stats_melt <-  melt(stats, id.vars = c('data', 'GROUP', 'TREAT'))
    names(stats_melt)[names(stats_melt)=='variable'] <- 'statistics'
    
    stats3 <- ddply(combdata, c('data', 'GROUP', 'TREAT'), summarise, mean=mean(mean)) # only for double checking that the mean of the means given by the neg binomial distribution for each iter is the same as the mean of the simulated data
    
    ppc = ggplot(combdata, aes(x=value, fill=data)) +
      geom_histogram(aes(y=stat(ncount)), position = 'identity', alpha=0.5, breaks = seq(0,250, by=1))
      
    # add predictive distribution of mean of simulated data
    if (mean_dist){
      ppc = ppc + geom_histogram(data = combdata, aes(x=mean, y=stat(ncount)), position = 'identity', color='yellow', alpha=0.1, bins = 2000)
    }
    # geom_vline(data = stats3, aes(xintercept = mean, color = data), linetype="dotted", alpha=1, size=0.2) # only for double checking that the mean of the means given by the neg binomial distribution for each iter is the same as the mean of the simulated data
    
    # add vertical lines for means and medians  
    if (statistics){
      ppc = ppc + geom_vline(data = stats_melt, aes(xintercept = value, color = data, linetype=statistics), alpha=1, size=0.2) +
        scale_linetype_manual(values=c("longdash", "dotdash"))
    }
    
    # plot layout and saving
    ppc = ppc + facet_grid(GROUP ~ TREAT) +
      coord_cartesian(xlim = c(0,maxcount)) + xlab("mosquito count") + ylab('normalised counts')
    ggsave(plottofile, ppc, height=12, width=18, dpi=300)

}
