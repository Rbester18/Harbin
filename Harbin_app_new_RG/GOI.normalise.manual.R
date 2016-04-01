GOI.normalise.manual<-function(GOIfiles, refgenefiles, GOImin, GOImax, Refmin, Refmax, refdatabase=NULL, answer, num, htest){
  # Determine number of reference genes
  nrefgenes<-length(refgenefiles)
  
  # Read in gene of interest files
  tempout<-GOIfiles
  GOIdata<-tempout[order(tempout[,"Name"]),]
  n<-nrow(GOIdata)
  
  # Determine which of GOI Rep.Calc.Conc values are invalid
  GOI.min<-GOImin
  GOI.max<-GOImax
  GOI.valid.ind<-rep(1,times=n)
  GOI.valid.ind[GOIdata[,"Rep.Calc.Conc"]<GOI.min]<- -999
  GOI.valid.ind[GOIdata[,"Rep.Calc.Conc"]>GOI.max]<- 999
  GOI.valid.ind[is.na(GOIdata[,"Rep.Calc.Conc"])]<- 0
  
  # Read in reference gene files and determine validity (i.e. whether the rows correspond with gene of interest set)
  refgenedata<-vector("list",nrefgenes)
  refgenevals<-vector("list",nrefgenes)
  refgene.valid.ind<-rep(1,times=n)
  for(j in 1:nrefgenes){
    tempout<-refgenefiles[[j]]
    refgenedata[[j]]<-tempout[order(tempout[,"Name"]),]
    refgene.min<-Refmin[j]
    refgene.max<-Refmax[j]
    refgene.valid.ind[refgenedata[[j]][,"Rep.Calc.Conc"]<refgene.min]<- -999
    refgene.valid.ind[refgenedata[[j]][,"Rep.Calc.Conc"]>refgene.max]<- 999
    refgene.valid.ind[is.na(refgenedata[[j]][,"Rep.Calc.Conc"])]<- NA
    
    refgenevals[[j]]<-as.numeric(refgenedata[[j]][,"Rep.Calc.Conc"])
  }
  
  # Normalisation of gene of interest Rep.Calc.Conc values
  refindex<-(Reduce("*",refgenevals))^(1/nrefgenes) # geometric mean of reference gene values
  GOI.normdata<-as.numeric(GOIdata[,"Rep.Calc.Conc"])/refindex
  
  # Select valid gene of interest data only
  GOI.validdata<-GOI.normdata[((GOI.valid.ind==1) & (refgene.valid.ind==1))]
  #print(GOI.validdata)
  
  if(answer == 1){
    
    # Comparison to reference data base
    refdatabase.data<-refdatabase
    for(i in 1:n){
      if(GOIdata[i,"Name"] %in% refdatabase.data[,"Name"]){
        duplicated_names <- "Warning: GOI names already found in reference data set!"
      } else {
        duplicated_names <- "No duplicated names in reference data set"
      }
    }
    refdatabase.valid.ind<-rep(0,times=nrow(refdatabase.data))
    refdatabase.valid.ind[refdatabase.data$Interval!=0]<-1
    #boundvals<-convergence.check(c(refdatabase.data[refdatabase.valid.ind==1,"GOI.normalised"],GOI.validdata),graph=graph)
    boundvals <-convergence.check(c(refdatabase.data[refdatabase.valid.ind==1,"GOI.normalised"],GOI.validdata))
    refdatabase.index<-rep(0,times=nrow(refdatabase.data))
    refdatabase.index[((refdatabase.data[,"GOI.normalised"]>=0) & (refdatabase.data[,"GOI.normalised"]<=boundvals[1]))]<-1
    refdatabase.index[((refdatabase.data[,"GOI.normalised"]>boundvals[1]) & (refdatabase.data[,"GOI.normalised"]<=boundvals[2]))]<-2
    refdatabase.index[((refdatabase.data[,"GOI.normalised"]>boundvals[2]) & (refdatabase.data[,"GOI.normalised"]<=boundvals[3]))]<-3
    refdatabase.index[((refdatabase.data[,"GOI.normalised"]>boundvals[3]) & (refdatabase.data[,"GOI.normalised"]<=boundvals[4]))]<-4
    refdatabase.index[refdatabase.data[,"GOI.normalised"]>boundvals[4]]<-5
    refdatabase.index[refdatabase.valid.ind!=1]<-0
    refdatabase.data[,"Interval"]<-refdatabase.index
    
    # Calculate proportion of labels changing in reference data base
    
    if (htest==1){
      refdatabase_length <- nrow(refdatabase.data)
      old <- refdatabase$Interval
      new <- refdatabase.data$Interval
      number_change <- length(which(old != new))
      proportion <- number_change/refdatabase_length
      count.labels <- as.numeric(proportion)
      
      # Kolmogorov-Smirnov test: Reference data base vs. new data
      line1 <- "Kolmogorov-Smirnov Test: New dataset vs. Reference dataset"
      ks.out<-ks.test(x=refdatabase.data[,"GOI.normalised"],y=GOI.validdata,alternative='two.sided')
      numchange <- paste("Proportion of labels changing in reference data base: ",round(count.labels*100,2),"%",sep="")
      line2 <- "H0: New data originated from same distribution as reference data"
      line3 <- "H1: New data and reference data come from different distributions"
      line4 <- paste("Kolmogorov-Smirnov test p-value = ",round(ks.out$p.value,4),sep="")
      if(ks.out$p.value<=0.05){
        line5 <- "WARNING: Reference data and new data may not be compatible!"
      } else {
        line5 <- "Reference data and new data may be compatible!"
      }
      } else {
        # Harbin test: Reference data base vs. new data
        ## (NOTE: HARBIN TEST NOT YET WORKING AS INTENDED!!)
        line1 <- "Harbin Test: New dataset vs. Reference dataset"
        harbin.out<-harbin.test(x=refdatabase.data[,"GOI.normalised"],y=GOI.validdata,reps=1000)
        numchange <- paste("Proportion of labels changing in reference data base: ",round(harbin.out$statistic*100,2),"%",sep="")
        line2 <- "H0: New data originated from same distribution as reference data"
        line3 <- "H1: New data and reference data come from different distributions"
        line4 <- paste("Harbin test p-value = ",harbin.out$p.value,"\n",sep="")
        if(harbin.out$p.value<=0.05){
          line5 <- "WARNING: Reference data and new data may not be compatible!"
        } else {
          line5 <- "Reference data and new data may be compatible!"
        }
      }
    } else {
      duplicated_names <- "No reference data set to compare"
      line1 <- "No reference data set to compare"
      numchange <- "No reference data set to compare"
      line2 <- "No reference data set to compare"
      line3 <- "No reference data set to compare"
      line4 <- "No reference data set to compare"
      line5 <- "No reference data set to compare"
      
      boundvals <- convergence.check(GOI.validdata)
    }
  
  # Calculate categories for current data set
  interval.index<-rep(0,times=n)
  interval.index[((GOI.normdata>=0) & (GOI.normdata<=boundvals[1]))]<-1
  interval.index[((GOI.normdata>boundvals[1]) & (GOI.normdata<=boundvals[2]))]<-2
  interval.index[((GOI.normdata>boundvals[2]) & (GOI.normdata<=boundvals[3]))]<-3
  interval.index[((GOI.normdata>boundvals[3]) & (GOI.normdata<=boundvals[4]))]<-4
  interval.index[GOI.normdata>boundvals[4]]<-5
  interval.index[!((GOI.valid.ind==1) & (refgene.valid.ind==1))]<-0
  notescol<-rep("Pass",times=n)
  notescol[is.na(as.numeric(GOIdata[,"Rep.Calc.Conc"])/refindex)]<-"Fail"
  notescol[is.na(GOIdata[,"Rep.Calc.Conc"])]<-"No data"
  notescol[(GOI.valid.ind==-999)]<- "GOI below standard curve range"
  notescol[(GOI.valid.ind==999)]<- "GOI above standard curve range"
  notescol[(refgene.valid.ind==-999)]<- "Refgene below standard curve range"
  notescol[(refgene.valid.ind==999)]<- "Refgene above standard curve range"
  resultsmat<-data.frame(Name=GOIdata[,"Name"],
                         GOI.data=GOIdata[,"Rep.Calc.Conc"],
                         Reference.index=refindex,
                         GOI.normalised=as.numeric(GOIdata[,"Rep.Calc.Conc"])/refindex,
                         Interval=interval.index,
                         Date=GOIdata[,"Date"],
                         #Time=GOIdata[,"Time"],
                         Note=notescol)
  
  # Add new data to reference data base (if applicable)
  if(num == 1){
    new_db <- rbind(refdatabase.data,resultsmat)
  }  else{
    new_db <- NULL
  } 
  
  # Output current data set results
  GOI.results <- list("duplicated_names" = duplicated_names, "newdb" = new_db, "resultsmat"=resultsmat, "line1" = line1,  "numchange"=numchange, "line2"=line2, "line3"=line3, "line4"=line4, "line5"=line5)
  
  return(GOI.results)
  
}