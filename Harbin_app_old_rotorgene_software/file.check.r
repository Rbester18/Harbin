
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
} else if(n < refrows){
  NO_GOI_1[[j]] <- "Warning: More reference gene rows. Look at number of replicates per sample"
} else {
  NO_GOI_1[[j]] <- "Number of rows per gene match"
}
for(i in 1:n){
  
  if(!(GOIdata[i,"Name"] %in% refgenedata[[j]][,"Name"])){
    NO_GOI_2[[j]] <- paste("Error: Missing sample name in gene of interest file(s)!",sep="")
  } else if (!(refgenedata[[j]][i,"Name"] %in% GOIdata[,"Name"])){
    NO_GOI_2[[j]] <- paste("Error: Missing sample name in Reference gene ",j," file(s)!",sep="")
  } else {
    NO_GOI_2[[j]] <- "Sample names match"
  }
}

GOI_complete <- length(complete.cases(GOIdata[,"Rep.Calc.Conc"]))
RG_complete <- length(complete.cases(refgenedata[[j]][,"Rep.Calc.Conc"]))
if (GOI_complete > RG_complete){
  NO_GOI_3[[j]] <- paste("Error: Missing concentration value in reference gene ",j," file(s)!",sep="")
}

else if(RG_complete > GOI_complete){
  NO_GOI_3[[j]] <- "Error: Missing concentration value in gene of interest file(s)!"
} else {
  NO_GOI_3[[j]] <- "All genes have equal number of concentration values"
}
}
file.results <- list("NO_GOI_1" = NO_GOI_1, "NO_GOI_2" = NO_GOI_2, "NO_GOI_3" = NO_GOI_3)
return(file.results)

}