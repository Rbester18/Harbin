
file.check <-function(GOIfiles, refgenefiles) {

nrefgenes<-length(refgenefiles)
tempout<-GOIfiles
GOIdata<-tempout[order(tempout[,"Name"]),]
n<-nrow(GOIdata)

refgenedata<-vector("list",nrefgenes)
NO_GOI_1<-vector("list",nrefgenes)
NO_GOI_2<-vector("list",nrefgenes)
NO_GOI_3<-vector("list",nrefgenes)

for(j in 1:nrefgenes){
  tempout<-refgenefiles[[j]]
  refgenedata[[j]]<-tempout[order(tempout[,"Name"]),]
refrows <- nrow(refgenedata[[j]])
if(n > refrows){
  NO_GOI_1[[j]] <- "Warning: More gene of interest rows. Look at number of replicates per sample"
  names(NO_GOI_1)[j] <- paste("Row number check for reference gene ",j," file(s)",sep="")
} else if(n < refrows){
  NO_GOI_1[[j]] <- "Warning: More reference gene rows. Look at number of replicates per sample"
  names(NO_GOI_1)[j] <- paste("Row number check for reference gene ",j," file(s)",sep="")
} else {
  NO_GOI_1[[j]] <- "Number of rows match gene of interest number of rows"
  names(NO_GOI_1)[j] <- paste("Row number check for reference gene ",j," file(s)",sep="")
}
for(i in 1:n){
  
  if(!(GOIdata[i,"Name"] %in% refgenedata[[j]][,"Name"])){
    NO_GOI_2[[j]] <- paste("Error: Missing sample name in gene of interest file(s)",sep="")
    names(NO_GOI_2)[j] <- paste("Name check for reference gene ",j," file(s)",sep="")
  } else if (!(refgenedata[[j]][i,"Name"] %in% GOIdata[,"Name"])){
    NO_GOI_2[[j]] <- paste("Error: Missing sample name in Reference gene ",j," file(s)",sep="")
    names(NO_GOI_2)[j] <- paste("Name check for reference gene ",j," file(s)",sep="")
  } else {
    NO_GOI_2[[j]] <- "Sample names match gene of interest sample names"
    names(NO_GOI_2)[j] <- paste("Name check for reference gene ",j," file(s)",sep="")
  }
}

GOI_complete <- length(complete.cases(GOIdata[,"Rep.Calc.Conc"]))
RG_complete <- length(complete.cases(refgenedata[[j]][,"Rep.Calc.Conc"]))
if (GOI_complete > RG_complete){
  NO_GOI_3[[j]] <- paste("Error: Missing concentration value in reference gene ",j," file(s)",sep="")
  names(NO_GOI_3)[j] <- paste("Sample value check for reference gene ",j," file(s)",sep="")
}

else if(RG_complete > GOI_complete){
  NO_GOI_3[[j]] <- "Error: Missing concentration value in gene of interest file(s)"
  names(NO_GOI_3)[j] <- paste("Sample value check for reference gene ",j," file(s)",sep="")
} else {
  NO_GOI_3[[j]] <- "All genes have equal number of concentration values"
  names(NO_GOI_3)[j] <- paste("Sample value check for reference gene ",j," file(s)",sep="")
}
}
file.results <- list("NO_GOI_1" = NO_GOI_1, "NO_GOI_2" = NO_GOI_2, "NO_GOI_3" = NO_GOI_3)
return(file.results)

}