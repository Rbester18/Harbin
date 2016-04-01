library(shiny)
library(psych)
library(car)
library(beeswarm)

source("convergence.check.R")
source("convergence.plot.R")
source("GOI.normalise.manual.R")
source("GOI.normalise.RG.R")
source("harbin.test.old.R")
source("file.check.r")

shinyServer (function(input, output) {
  
  #####Example instructions download
  
  output$manual <- downloadHandler(
    filename="Harbin_instruction_manual.pdf",
    content=function(file) {
      file.copy("Harbin_manual.pdf", file)
    }
  )
  
  ####Number of reference genes file inputs
  
  output$input_ref <- renderUI({
    num <- as.integer(input$refgenenum.answer)
    if(is.null(num))
      return(NULL)
    lapply(1:num, function(i) {
      fileInput(paste0("Reference gene ", i), label = paste0("Reference gene ", i), multiple =TRUE)
    })
  })
  
  ####If applicable, reference database file input
  
  output$refdb <- renderUI({
    answer <- as.integer(input$Database)
    if(is.null(answer))
      return(NULL)
    if(answer == 1){
      fileInput("refdatabase", label="Select reference data set")
    }
  })
  
  #####Example file download
  
  output$Examplefile <- downloadHandler(
    filename="RotorGeneQ_example.csv",
    content=function(con) {
      file.copy("RotorGeneQ_example_input.csv", con)
    }
  )
 
   output$Examplefilem <- downloadHandler(
       filename="Manual_Cq_input_example.csv",
       content=function(con) {
         file.copy("Manual_Cq_example_input.csv", con)
       }
     ) 
     
     
     #####Example refdatebase file download
     
      output$RD_Examplefile <- downloadHandler(
       filename="Reference_dataset_example.csv",
       content=function(con) {
         file.copy("Reference_dataset_example_input.csv", con)
       }
     )
      
      output$RD_Examplefile_m <- downloadHandler(
          filename="Reference_dataset_example.csv",
          content=function(con) {
            file.copy("Reference_dataset_example_input.csv", con)
          }
      )
  
  
  output$Harbindescribe <- renderUI({
    answer <- as.integer(input$Database)
    if(is.null(answer))
      return(NULL)
    if(answer == 1){
      h5(strong("Choose between the Kolmogorov-Smirnov test or the Harbin test to compare normalised data to a reference data set."), em("Kolmogorov-Smirnov test is the default. The harbin test is more conservative and has a higher false negative rate. The harbin test is applicable if your reference data set has more samples than the test data."))
    }
  })
  
  
  output$Harbintesttype <- renderUI({
    answer <- as.integer(input$Database)
    if(is.null(answer))
      return(NULL)
    if(answer == 1){
      radioButtons("Harbintesttype", "Choose test for comparison", 
                   choices = list("Kolmogorov-Smirnov Test" = 1, "Harbin test (experimental)" = 2), selected=1)
    }
  })
 
  ####Function to load gene of interest (GOI) files
  
  GOI_o <- function(){
    new.unknown.frame<-NULL
    new.standard.frame<-NULL
    files_df <-input$GOI
    if(is.null(files_df))
      return(NULL)
    nfiles<-nrow(files_df)
    for(j in 1:nfiles){
      file_info <- files_df[j,]
      file <- file_info[,4]
      headerdata1 <- scan(file, what="character", sep="\n", quiet=TRUE)
      headerdata2 <- grep("No.",headerdata1, fixed = TRUE)
      headerdata3 <- headerdata2+4
      datestamp<-scan(file,what="character",nlines=1,skip=3,sep=",",quiet=TRUE)[2] # scan date stamp on 4th line
      #timestamp<-scan(file,what="character",nlines=1,skip=4,sep=",",quiet=TRUE)[2] # scan time stamp on 5th line
      varnames<-scan(file,what="character",nlines=1,skip=headerdata3,sep=",",quiet=TRUE) # scan variable names on 9th line
      tempdata<-scan(file,what="character",skip=headerdata3+1,sep=",",quiet=TRUE) # scan all data and store in vector
      nlines<-floor(length(tempdata)/14) # calculate number of lines in the file
      datamat<-NULL
      for(i in 1:nlines){
        temprow<-tempdata[((i-1)*14+1):(i*14)] # convert the data vector to lines again
        datamat<-rbind(datamat,temprow)
      }
      colnames(datamat)<-varnames
      allnames<-unique(datamat[datamat[,"Type"]=="Unknown","Name"])
      datamat<-datamat[datamat[,"Rep. Calc. Conc."]!="",,drop=FALSE] # select rows with non-empty cells in Rep. Calc. Conc. column
      non.empty.names<-unique(datamat[datamat[,"Type"]=="Unknown","Name"])
      names.NA<-allnames[!(allnames %in% non.empty.names)]
      new.unknown<-datamat[datamat[,"Type"]=="Unknown",c("Name","Type","Rep. Calc. Conc."),drop=FALSE] # select only rows of 'Unknown' type
      colnames(new.unknown)<-c("Name","Type","Rep.Calc.Conc.")
      if(length(names.NA)>0){
        new.unknown<-rbind(new.unknown,data.frame(Name=names.NA,Type="Unknown",Rep.Calc.Conc.=NA))
      }
      new.standard<-datamat[datamat[,"Type"]=="Standard",c("Name","Type","Rep. Calc. Conc."),drop=FALSE] # select only rows of 'Standard' type
      unknownvals<-as.numeric(as.character(new.unknown[,"Rep.Calc.Conc."]))
      standardvals<-as.numeric(as.character(new.standard[,"Rep. Calc. Conc."]))
      if(nrow(new.unknown)>0){
        #new.unknown.frame<-rbind(new.unknown.frame,data.frame(Name=new.unknown[,"Name"],Rep.Calc.Conc=unknownvals,Date=datestamp, Time=timestamp))
        new.unknown.frame<-rbind(new.unknown.frame,data.frame(Name=new.unknown[,"Name"],Rep.Calc.Conc=unknownvals,Date=datestamp))
      }
      if(nrow(new.standard)>0){
        new.standard.frame<-rbind(new.standard.frame,data.frame(Name=new.standard[,"Name"],Rep.Calc.Conc=standardvals))
      }
    }
    new.unknown.frame[,"Name"]<-as.character(new.unknown.frame[,"Name"])
    GOIfiles <- list(Unknown=new.unknown.frame,Min=min(new.standard.frame[,"Rep.Calc.Conc"]),Max=max(new.standard.frame[,"Rep.Calc.Conc"]))
  
    return(GOIfiles)
  }
  
  ####Function to load reference files
  
  RG_o <- function(){
    num <- as.integer(input$refgenenum.answer)
    refgenefiles <- vector("list",num)
    refmin <- vector("list",num)
    refmax <- vector("list",num)
    for(k in 1:num){
      new.unknown.frame <- NULL
      new.standard.frame <- NULL
      files_df <-input[[paste0("Reference gene ", k)]]
      if(is.null(files_df))
        return(NULL)
      nfiles<-nrow(files_df)
      for(j in 1:nfiles){
        file_info <- files_df[j,]
        file <- file_info[,4]
        headerdata1 <- scan(file, what="character", sep="\n", quiet=TRUE)
        headerdata2 <- grep("No.",headerdata1, fixed = TRUE)
        headerdata3 <- headerdata2+4
        datestamp<-scan(file,what="character",nlines=1,skip=3,sep=",",quiet=TRUE)[2] # scan date stamp on 4th line
        #timestamp<-scan(file,what="character",nlines=1,skip=4,sep=",",quiet=TRUE)[2] # scan time stamp on 5th line
        varnames<-scan(file,what="character",nlines=1,skip=headerdata3,sep=",",quiet=TRUE) # scan variable names on 9th line
        tempdata<-scan(file,what="character",skip=headerdata3+1,sep=",",quiet=TRUE) # scan all data and store in vector
        nlines<-floor(length(tempdata)/14) # calculate number of lines in the file
        datamat<-NULL
        for(i in 1:nlines){
          temprow<-tempdata[((i-1)*14+1):(i*14)] # convert the data vector to lines again
          datamat<-rbind(datamat,temprow)
        }
        colnames(datamat)<-varnames
        allnames<-unique(datamat[datamat[,"Type"]=="Unknown","Name"])
        datamat<-datamat[datamat[,"Rep. Calc. Conc."]!="",,drop=FALSE] # select rows with non-empty cells in Rep. Calc. Conc. column
        non.empty.names<-unique(datamat[datamat[,"Type"]=="Unknown","Name"])
        names.NA<-allnames[!(allnames %in% non.empty.names)]
        new.unknown<-datamat[datamat[,"Type"]=="Unknown",c("Name","Type","Rep. Calc. Conc."),drop=FALSE] # select only rows of 'Unknown' type
        colnames(new.unknown)<-c("Name","Type","Rep.Calc.Conc.")
        if(length(names.NA)>0){
          new.unknown<-rbind(new.unknown,data.frame(Name=names.NA,Type="Unknown",Rep.Calc.Conc.=NA))
        }
        new.standard<-datamat[datamat[,"Type"]=="Standard",c("Name","Type","Rep. Calc. Conc."),drop=FALSE] # select only rows of 'Standard' type
        unknownvals<-as.numeric(as.character(new.unknown[,"Rep.Calc.Conc."]))
        standardvals<-as.numeric(as.character(new.standard[,"Rep. Calc. Conc."]))
        if(nrow(new.unknown)>0){
          #new.unknown.frame<-rbind(new.unknown.frame,data.frame(Name=new.unknown[,"Name"],Rep.Calc.Conc=unknownvals,Date=datestamp, Time=timestamp))
          new.unknown.frame<-rbind(new.unknown.frame,data.frame(Name=new.unknown[,"Name"],Rep.Calc.Conc=unknownvals,Date=datestamp))
        }
        if(nrow(new.standard)>0){
          new.standard.frame<-rbind(new.standard.frame,data.frame(Name=new.standard[,"Name"],Rep.Calc.Conc=standardvals))
        }
      }
      new.unknown.frame[,"Name"]<-as.character(new.unknown.frame[,"Name"])

      refgenefiles[[k]] <- new.unknown.frame
      refmin[[k]] <- min(new.standard.frame[,"Rep.Calc.Conc"])
      refmax[[k]] <- max(new.standard.frame[,"Rep.Calc.Conc"])
      
    }
    refgene_con <- list(RG=refgenefiles, RGmin=refmin, RGmax=refmax)
    return(refgene_con)
  }
  
  ####Output of GOI file table
  
  output$GOI_o <- renderTable({
    GOI_df <- GOI_o()
    validate(
      need(GOI_df != "", label = "GOI data")
    )
    as.data.frame(GOI_df$Unknown)
  })
  
  ####Output of reference gene files table (concatenated)
  
  output$RG_o <- renderTable({
    RG_df <- RG_o()
    validate(
      need(RG_df != "", label = "Reference gene data")
    )
    RG_out <- cbind(as.data.frame(RG_df$RG))
    num <- as.integer(input$refgenenum.answer)
    if (num ==1){
     colnames(RG_out) <- c("Refgene1_Name", "Rep.Calc.Conc", "Date")
    }else if(num ==2) {
     colnames(RG_out) <- c("Refgene1_Name", "Rep.Calc.Conc", "Date", "Refgene2_Name", "Rep.Calc.Conc", "Date")
    }else if(num ==3) {
      colnames(RG_out) <- c("Refgene1_Name", "Rep.Calc.Conc", "Date", "Refgene2_Name", "Rep.Calc.Conc", "Date", "Refgene3_Name", "Rep.Calc.Conc", "Date")
    }else if(num ==4) {
      colnames(RG_out) <- c("Refgene1_Name", "Rep.Calc.Conc", "Date", "Refgene2_Name", "Rep.Calc.Conc", "Date", "Refgene3_Name", "Rep.Calc.Conc", "Date", "Refgene4_Name", "Rep.Calc.Conc", "Date")
    }else if(num ==5) {
      colnames(RG_out) <- c("Refgene1_Name", "Rep.Calc.Conc", "Date", "Refgene2_Name", "Rep.Calc.Conc", "Date", "Refgene3_Name", "Rep.Calc.Conc", "Date", "Refgene4_Name", "Rep.Calc.Conc", "Date", "Refgene5_Name", "Rep.Calc.Conc", "Date")
    }else if(num ==6) {
      colnames(RG_out) <- c("Refgene1_Name", "Rep.Calc.Conc", "Date", "Refgene2_Name", "Rep.Calc.Conc", "Date", "Refgene3_Name", "Rep.Calc.Conc", "Date", "Refgene4_Name", "Rep.Calc.Conc", "Date", "Refgene5_Name", "Rep.Calc.Conc", "Date", "Refgene6_Name", "Rep.Calc.Conc", "Date")
    }else if(num ==7) {
      colnames(RG_out) <- c("Refgene1_Name", "Rep.Calc.Conc", "Date", "Refgene2_Name", "Rep.Calc.Conc", "Date", "Refgene3_Name", "Rep.Calc.Conc", "Date", "Refgene4_Name", "Rep.Calc.Conc", "Date", "Refgene5_Name", "Rep.Calc.Conc", "Date", "Refgene6_Name", "Rep.Calc.Conc", "Date", "Refgene7_Name", "Rep.Calc.Conc", "Date")
    }else if(num ==8) {
      colnames(RG_out) <- c("Refgene1_Name", "Rep.Calc.Conc", "Date", "Refgene2_Name", "Rep.Calc.Conc", "Date", "Refgene3_Name", "Rep.Calc.Conc", "Date", "Refgene4_Name", "Rep.Calc.Conc", "Date", "Refgene5_Name", "Rep.Calc.Conc", "Date", "Refgene6_Name", "Rep.Calc.Conc", "Date", "Refgene7_Name", "Rep.Calc.Conc", "Date", "Refgene8_Name", "Rep.Calc.Conc", "Date")
    }else if(num ==9) {
      colnames(RG_out) <- c("Refgene1_Name", "Rep.Calc.Conc", "Date", "Refgene2_Name", "Rep.Calc.Conc", "Date", "Refgene3_Name", "Rep.Calc.Conc", "Date", "Refgene4_Name", "Rep.Calc.Conc", "Date", "Refgene5_Name", "Rep.Calc.Conc", "Date", "Refgene6_Name", "Rep.Calc.Conc", "Date", "Refgene7_Name", "Rep.Calc.Conc", "Date", "Refgene8_Name", "Rep.Calc.Conc", "Date", "Refgene9_Name", "Rep.Calc.Conc", "Date")
    }else{
      colnames(RG_out) <- c("Refgene1_Name", "Rep.Calc.Conc", "Date", "Refgene2_Name", "Rep.Calc.Conc", "Date", "Refgene3_Name", "Rep.Calc.Conc", "Date", "Refgene4_Name", "Rep.Calc.Conc", "Date", "Refgene5_Name", "Rep.Calc.Conc", "Date", "Refgene6_Name", "Rep.Calc.Conc", "Date", "Refgene7_Name", "Rep.Calc.Conc", "Date", "Refgene8_Name", "Rep.Calc.Conc", "Date", "Refgene9_Name", "Rep.Calc.Conc", "Date", "Refgene10_Name", "Rep.Calc.Conc", "Date")
    }
    RG_out
  })
  
  ####Function to load reference data set if applicable
  
  RD_o <- reactive({
    files_df <-input$refdatabase
    validate(
      need(files_df != "", label = "Reference data set")
    )
    refdatabase <-files_df$datapath
    RD_f <- read.csv(refdatabase)
    return(RD_f)
  })
  
  output$RD_o <- renderTable({
    RD_o()
  })
  
  ####Running function to normalise RG files
  
  GOI_norm <- reactive({
    answer <- as.integer(input$Database)
    num <- as.integer(input$Database_add)
    htest <- as.integer(input$Harbintesttype)
    GOI_df <- GOI_o()
    validate(
      need(GOI_df != "", label = "GOI data")
    )
    RG_df <- RG_o()
    validate(
      need(RG_df != "", label = "Data for number of reference genes selected")
    )
    GOIfiles <- GOI_df$Unknown
    RGfiles <- RG_df$RG
    GOImin <- GOI_df$Min
    GOImax <- GOI_df$Max
    RGmin <- RG_df$RGmin
    RGmax <- RG_df$RGmax
    answer <- as.integer(input$Database)
    num <- as.integer(input$Database_add)
    htest <- as.integer(input$Harbintesttype)
    if(answer == 2) {
      GOI.normalise.RG(GOIfiles,RGfiles, GOImin, GOImax, RGmin, RGmax, refdatabase=NULL, answer, num, htest)
    } else {
      GOI.normalise.RG(GOIfiles,RGfiles, GOImin, GOImax, RGmin, RGmax, refdatabase=RD_o(), answer, num, htest)
    }
  })
  
  output$GOI_norm_table <- renderTable({
    GOI_norm()$resultsmat
  })
  
  
  output$NO_GOI_1 <- renderPrint({
    GOI_df <- GOI_o()
    GOI_df <- GOI_o()
    validate(
      need(GOI_df != "", label = "GOI data")
    )
    RG_df <- RG_o()
    RG_df <- RG_o()
    validate(
      need(RG_df != "", label = "Data for number of reference genes selected")
    )
    GOIfiles <- GOI_df$Unknown
    RGfiles <- RG_df$RG
    file_check <- file.check(GOIfiles,RGfiles)
    file_check$NO_GOI_1
  })
  
  output$NO_GOI_2 <- renderPrint({
    GOI_df <- GOI_o()
    GOI_df <- GOI_o()
    validate(
      need(GOI_df != "", label = "GOI data")
    )
    RG_df <- RG_o()
    RG_df <- RG_o()
    validate(
      need(RG_df != "", label = "Data for number of reference genes selected")
    )
    GOIfiles <- GOI_df$Unknown
    RGfiles <- RG_df$RG
    file_check <- file.check(GOIfiles,RGfiles)
    file_check$NO_GOI_2
  })
  
  output$NO_GOI_3 <- renderPrint({
    GOI_df <- GOI_o()
    GOI_df <- GOI_o()
    validate(
      need(GOI_df != "", label = "GOI data")
    )
    RG_df <- RG_o()
    RG_df <- RG_o()
    validate(
      need(RG_df != "", label = "Data for number of reference genes selected")
    )
    GOIfiles <- GOI_df$Unknown
    RGfiles <- RG_df$RG
    file_check <- file.check(GOIfiles,RGfiles)
    file_check$NO_GOI_3
  })
  
  #####Duplicated names in data set
  
  output$duplicated_names <- renderText({
    GOI_norm()$duplicated_names
  }) 
  
  ####Kolmogorov-Smirnov test
  
  output$line1 <- renderText({
    GOI_norm()$line1
  })
  
  output$numchange <- renderText({
    GOI_norm()$numchange
  })
  
  output$line2 <- renderText({
    GOI_norm()$line2
  }) 
  
  output$line3 <- renderText({
    GOI_norm()$line3
  }) 
  
  output$line4 <- renderText({
    GOI_norm()$line4
  }) 
  
  output$line5 <- renderText({
    GOI_norm()$line5
  }) 
  
  ####Normalised output download
  
  output$downloadGOI <- downloadHandler(
    filename = function() {
      paste("Normalised_GOI.csv")
    },
    content = function(file) {
      write.table(GOI_norm()$resultsmat, file, sep=",",
                  row.names = FALSE)
    }
  )
  
  ####New database (Normalised output added to database) download
  
  
  output$downloadDB <- downloadHandler(
    filename = function() {
      paste("New_reference_dataset.csv")
    },
    content = function(file) {
      write.table(GOI_norm()$newdb, file, sep=",",
                  row.names = FALSE)
    }
  )
  
  ####Manual importing of Cq values
  
  
  output$Ref_slopes <- renderUI({
    numrefs <- as.integer(input$num_refgenes)
    if(is.null(numrefs))
      return(NULL)
    lapply(1:numrefs, function(i) {
      numericInput(paste0("SlopeRefgene", i), 
                   label = paste0("Slope of reference gene ", i),
                   value = -3.2
      )
    })
  })
  
  output$Ref_intercept <- renderUI({
    numrefs <- as.integer(input$num_refgenes)
    if(is.null(numrefs))
      return(NULL)
    lapply(1:numrefs, function(i) {
      numericInput(paste0("interceptRefgene", i), 
                   label = paste0("Y-intercept of reference gene ", i), 
                   value = 20
      )
    })
  })
  
  output$SD_min <- renderUI({
    numrefs <- as.integer(input$num_refgenes)
    if(is.null(numrefs))
      return(NULL)
    lapply(1:numrefs, function(i) {
      numericInput(paste0("SD_min_", i), 
                   label = paste0("Minimum Cq of standard curve for reference gene ", i),
                   value = 15
      )
    })
  })
  
  output$SD_max <- renderUI({
    numrefs <- as.integer(input$num_refgenes)
    if(is.null(numrefs))
      return(NULL)
    lapply(1:numrefs, function(i) {
      numericInput(paste0("SD_max_", i), 
                   label = paste0("Maximum Cq of standard curve for reference gene ", i),
                   value = 30
      )
    })
  })
  
  manualdata <- reactive({
    numrefs <- as.integer(input$num_refgenes)
    GOIy <- as.numeric(input$GOIintercept)
    GOIs <- as.numeric(input$GOIslope)
    GOI_SDmin <- as.numeric(input$GOI_SD_min)
    GOI_SDmax <- as.numeric(input$GOI_SD_max)
    GOI_SDmax2 <- 10^((GOI_SDmin-GOIy)/GOIs)
    GOI_SDmin2 <- 10^((GOI_SDmax-GOIy)/GOIs)
    Cqs <-input$Cqvals
    if(is.null(Cqs))
      return(NULL)
    Cqs_table <-Cqs$datapath
    CqsData <- read.csv(Cqs_table)
    GOI_conc <- sapply(CqsData[,2], function(x) 10^((x-GOIy)/GOIs))
    GOI_df <- cbind(CqsData[,1], GOI_conc)
    GOI_df2 <- aggregate( GOI_df[,2] ~ GOI_df[,1], GOI_df, geometric.mean)
    num_samples <- nrow(GOI_df2)
    expdate <- rep(as.character(input$expdate), num_samples)
    GOI_df3 <- cbind(GOI_df2, expdate)
    colnames(GOI_df3) <- c("Name","Rep.Calc.Conc", "Date")
    ref_result <- vector("list",numrefs)
    ref_result2 <- vector("list",numrefs)
    ref_result3 <- vector("list",numrefs)
    ref_result4 <- vector("list",numrefs)
    Refy <- vector("list",numrefs)
    Refs <- vector("list",numrefs)
    Refmin <- vector("list",numrefs)
    Refmax <- vector("list",numrefs)
    Refmin2 <- vector("list",numrefs)
    Refmax2 <- vector("list",numrefs)
    for (i in 1:numrefs){
      Refy[[i]] <- as.numeric(input[[paste0("interceptRefgene", i)]])
      Refs[[i]] <- as.numeric(input[[paste0("SlopeRefgene", i)]])
      Refmin[[i]] <- as.numeric(input[[paste0("SD_min_", i)]])
      Refmax[[i]] <- as.numeric(input[[paste0("SD_max_", i)]])
      Refmax2[[i]] <- 10^((Refmin[[i]]-Refy[[i]])/Refs[[i]])
      Refmin2[[i]] <- 10^((Refmax[[i]]-Refy[[i]])/Refs[[i]])
      ref_result[[i]] <- sapply(CqsData[,i+2], function(x) 10^((x-Refy[[i]])/Refs[[i]]))
      ref_result2[[i]] <- cbind(CqsData[,1], ref_result[[i]])
      ref_result3[[i]] <- aggregate( ref_result2[[i]][,2] ~ ref_result2[[i]][,1], ref_result2[[i]], geometric.mean)
      expdate <- rep(as.character(input$expdate), num_samples)
      ref_result4[[i]] <- cbind(ref_result3[[i]], expdate)
      colnames(ref_result4[[i]]) <- c("Name","Rep.Calc.Conc", "Date")
    }
    
    conc_files <- list(GOI_m=GOI_df3, GOImin_m=GOI_SDmin2, GOImax_m=GOI_SDmax2, RG_m=ref_result4, RGmin_m=Refmin2, RGmax_m=Refmax2)
    return(conc_files)
  })
  
  ####Output of manual GOI file table
  
  output$GOI_om <- renderTable({
    results_m <- manualdata()
    validate(
      need(results_m != "", label = "Cq values")
    )
    results_m$GOI_m
  })
  
  
  ####Output of manual reference gene files table (concatenated)
  
  output$RG_om <- renderTable({
    results_m <- manualdata()
    validate(
      need(results_m != "", label = "Cq values")
    )
    RG_out_m <- cbind(as.data.frame(results_m$RG_m))
    numrefs <- as.integer(input$num_refgenes)
    if (numrefs ==1){
      colnames(RG_out_m) <- c("Refgene1_Name", "Rep.Calc.Conc", "Date")
    }else if(numrefs ==2) {
      colnames(RG_out_m) <- c("Refgene1_Name", "Rep.Calc.Conc", "Date", "Refgene2_Name", "Rep.Calc.Conc", "Date")
    }else if(numrefs ==3) {
      colnames(RG_out_m) <- c("Refgene1_Name", "Rep.Calc.Conc", "Date", "Refgene2_Name", "Rep.Calc.Conc", "Date", "Refgene3_Name", "Rep.Calc.Conc", "Date")
    }else if(numrefs ==4) {
      colnames(RG_out_m) <- c("Refgene1_Name", "Rep.Calc.Conc", "Date", "Refgene2_Name", "Rep.Calc.Conc", "Date", "Refgene3_Name", "Rep.Calc.Conc", "Date", "Refgene4_Name", "Rep.Calc.Conc", "Date")
    }else if(numrefs ==5) {
      colnames(RG_out_m) <- c("Refgene1_Name", "Rep.Calc.Conc", "Date", "Refgene2_Name", "Rep.Calc.Conc", "Date", "Refgene3_Name", "Rep.Calc.Conc", "Date", "Refgene4_Name", "Rep.Calc.Conc", "Date", "Refgene5_Name", "Rep.Calc.Conc", "Date")
    }else if(numrefs ==6) {
      colnames(RG_out_m) <- c("Refgene1_Name", "Rep.Calc.Conc", "Date", "Refgene2_Name", "Rep.Calc.Conc", "Date", "Refgene3_Name", "Rep.Calc.Conc", "Date", "Refgene4_Name", "Rep.Calc.Conc", "Date", "Refgene5_Name", "Rep.Calc.Conc", "Date", "Refgene6_Name", "Rep.Calc.Conc", "Date")
    }else if(numrefs ==7) {
      colnames(RG_out_m) <- c("Refgene1_Name", "Rep.Calc.Conc", "Date", "Refgene2_Name", "Rep.Calc.Conc", "Date", "Refgene3_Name", "Rep.Calc.Conc", "Date", "Refgene4_Name", "Rep.Calc.Conc", "Date", "Refgene5_Name", "Rep.Calc.Conc", "Date", "Refgene6_Name", "Rep.Calc.Conc", "Date", "Refgene7_Name", "Rep.Calc.Conc", "Date")
    }else if(numrefs ==8) {
      colnames(RG_out_m) <- c("Refgene1_Name", "Rep.Calc.Conc", "Date", "Refgene2_Name", "Rep.Calc.Conc", "Date", "Refgene3_Name", "Rep.Calc.Conc", "Date", "Refgene4_Name", "Rep.Calc.Conc", "Date", "Refgene5_Name", "Rep.Calc.Conc", "Date", "Refgene6_Name", "Rep.Calc.Conc", "Date", "Refgene7_Name", "Rep.Calc.Conc", "Date", "Refgene8_Name", "Rep.Calc.Conc", "Date")
    }else if(numrefs ==9) {
      colnames(RG_out_m) <- c("Refgene1_Name", "Rep.Calc.Conc", "Date", "Refgene2_Name", "Rep.Calc.Conc", "Date", "Refgene3_Name", "Rep.Calc.Conc", "Date", "Refgene4_Name", "Rep.Calc.Conc", "Date", "Refgene5_Name", "Rep.Calc.Conc", "Date", "Refgene6_Name", "Rep.Calc.Conc", "Date", "Refgene7_Name", "Rep.Calc.Conc", "Date", "Refgene8_Name", "Rep.Calc.Conc", "Date", "Refgene9_Name", "Rep.Calc.Conc", "Date")
    }else{
      colnames(RG_out_m) <- c("Refgene1_Name", "Rep.Calc.Conc", "Date", "Refgene2_Name", "Rep.Calc.Conc", "Date", "Refgene3_Name", "Rep.Calc.Conc", "Date", "Refgene4_Name", "Rep.Calc.Conc", "Date", "Refgene5_Name", "Rep.Calc.Conc", "Date", "Refgene6_Name", "Rep.Calc.Conc", "Date", "Refgene7_Name", "Rep.Calc.Conc", "Date", "Refgene8_Name", "Rep.Calc.Conc", "Date", "Refgene9_Name", "Rep.Calc.Conc", "Date", "Refgene10_Name", "Rep.Calc.Conc", "Date")
    }
    RG_out_m
  })
  
  #####Reference database for manual input option
  
  output$refdbm <- renderUI({
    answer <- as.integer(input$Database_m)
    if (is.null(answer)) 
      return(NULL)
    if(answer == 1){
      fileInput("refdatabase_m", label="Select reference data set")
    }
  })
  
  ####Function to load reference database if applicable
  
  RD_om <- reactive({
    files_dfm <-input$refdatabase_m
    validate(
      need(files_dfm != "", label = "Reference data set")
    )
    refdatabase_m <-files_dfm$datapath
    RD_f <- read.csv(refdatabase_m)
    return(RD_f)
  })
  
  output$RD_om <- renderTable({
    RD_om()
  })
  
  ####Running function to normalise manual input
  
  GOI_norm_m <- reactive({
    answer <- as.integer(input$Database_m)
    num <- as.integer(input$Database_add_m)
    htest <- as.integer(input$Harbintesttype_m)
    results_m <- manualdata()
    if (is.null(results_m)) 
      return(NULL)
    GOIfiles <- results_m$GOI_m
    RGfiles <- results_m$RG_m
    GOImin <- results_m$GOImin_m
    GOImax <- results_m$GOImax_m
    RGmin <- results_m$RGmin_m
    RGmax <- results_m$RGmax_m
    if(answer == 2) {
      GOI.normalise.manual(GOIfiles,RGfiles, GOImin, GOImax, RGmin, RGmax,refdatabase=NULL, answer, num, htest)
    } else {
      GOI.normalise.manual(GOIfiles,RGfiles, GOImin, GOImax, RGmin, RGmax, refdatabase=RD_om(), answer, num, htest)
    }
  })
  
  output$GOI_norm_table_m <- renderTable({
    results_m <- manualdata()
    validate(
      need(results_m != "", label = "Cq values")
    )
    GOI_norm_m()$resultsmat
  })
  
  #####Duplicated names in database
  
  output$duplicated_names_m <- renderText({
    GOI_norm_m()$duplicated_names
  }) 
  
  output$Harbindescribe_m <- renderUI({
    answer <- as.integer(input$Database_m)
    if(is.null(answer))
      return(NULL)
    if(answer == 1){
      h5(strong("Choose between the Kolmogorov-Smirnov test or the Harbin test to compare normalised data to a reference data set."), em("Kolmogorov-Smirnov test is the default. The harbin test is more conservative and has a higher false negative rate. The harbin test is applicable if your reference data set has more samples than the test data."))
    }
  })
  
  
  output$Harbintesttype_m <- renderUI({
    answer <- as.integer(input$Database_m)
    if(is.null(answer))
      return(NULL)
    if(answer == 1){
      radioButtons("Harbintesttype_m", "Choose test for comparison", 
                   choices = list("Kolmogorov-Smirnov Test" = 1, "Harbin test (experimental)" = 2), selected=1)
    }
  })
  
  
  
  ####Kolmogorov-Smirnov test
  
  output$line1_2 <- renderText({
    GOI_norm_m()$line1
  }) 
  
  output$numchange_2 <- renderText({
    GOI_norm_m()$numchange
  })
  
  output$line2_2 <- renderText({
    GOI_norm_m()$line2
  }) 
  
  output$line3_2 <- renderText({
    GOI_norm_m()$line3
  }) 
  
  output$line4_2 <- renderText({
    GOI_norm_m()$line4
  }) 
  
  output$line5_2 <- renderText({
    GOI_norm_m()$line5
  }) 
  
  ####Normalised output download
  
  output$downloadGOIm <- downloadHandler(
    filename = function() {
      paste("Normalised_GOI.csv")
    },
    content = function(file) {
      write.table(GOI_norm_m()$resultsmat, file, sep=",",
                  row.names = FALSE)
    }
  )
  
  ####New database (Normalised output added to database) download
  
  
  output$downloadDBm <- downloadHandler(
    filename = function() {
      paste("New_reference_dataset.csv")
    },
    content = function(file) {
      write.table(GOI_norm_m()$newdb, file, sep=",",
                  row.names = FALSE)
    }
  )
  
  ######Plot to compare intervals of database and new data:
  output$intervalbarplot <- renderPlot({
    input_answer <- as.integer(input$Input_type)
    if (input_answer==1){
      newdb_answer <- as.integer(input$Database_add)
      olddb_answer <- as.integer(input$Database)
      if (newdb_answer==1 && olddb_answer==1){
        new_data <- GOI_norm()$resultsmat
        new_db <- GOI_norm()$newdb
        old_db <- RD_o()
        plotmain_old <- "Unchanged reference data set"
        plotmain_new <- "Normalised gene of interest data (new)" 
        plotmain_oldnew <- "New reference data set"
        par(mfrow=c(3,1))
        convergence.plot(old_db[,4], plotmain_old)
        convergence.plot(new_data[,4], plotmain_new)
        convergence.plot(new_db[,4], plotmain_oldnew)

        par(mfrow=c(1,1))
        
      } else if (newdb_answer==2 && olddb_answer==1) {
        new_data <- GOI_norm()$resultsmat
        old_db <- RD_o()
        plotmain_old <- "Unchanged reference data set"
        plotmain_new <- "Normalised gene of interest data (new)" 
        par(mfrow=c(2,1))
        convergence.plot(old_db[,4], plotmain_old)
        convergence.plot(new_data[,4], plotmain_new)
       
        par(mfrow=c(1,1))
        
      } else if (newdb_answer==2 && olddb_answer==2) {
        new_data <- GOI_norm()$resultsmat
        plotmain_new <- "Normalised gene of interest data (new)" 
        par(mfrow=c(1,1))
        convergence.plot(new_data[,4], plotmain_new)
        par(mfrow=c(1,1))
        
      } else {
        print("Can not plot new data and new reference data set if old reference data set was not loaded")
      }
    } else {
      newdb_answer <- as.integer(input$Database_add_m)
      olddb_answer <- as.integer(input$Database_m)
      if (newdb_answer==1 && olddb_answer==1){
        new_data <- GOI_norm_m()$resultsmat
        if (is.null(new_data)) 
          return(NULL)
        new_db <- GOI_norm_m()$newdb
        old_db <- RD_om()
        plotmain_old <- "Unchanged reference data set"
        plotmain_new <- "Normalised gene of interest data (new)" 
        plotmain_oldnew <- "New reference data set"
        par(mfrow=c(3,1))
        convergence.plot(old_db[,4], plotmain_old)
        convergence.plot(new_data[,4], plotmain_new)
        convergence.plot(new_db[,4], plotmain_oldnew)
        
        par(mfrow=c(1,1))
        
      } else if (newdb_answer==2 && olddb_answer==1) {
        new_data <- GOI_norm_m()$resultsmat
        if (is.null(new_data)) 
          return(NULL)
        old_db <- RD_om()
        plotmain_old <- "Unchanged reference data set"
        plotmain_new <- "Normalised gene of interest data (new)" 
        par(mfrow=c(2,1))
        convergence.plot(old_db[,4], plotmain_old)
        convergence.plot(new_data[,4], plotmain_new)
        
        par(mfrow=c(1,1))
        
        
      } else if (newdb_answer==2 && olddb_answer==2) {
        new_data <- GOI_norm_m()$resultsmat
        if (is.null(new_data)) 
          return(NULL)
        plotmain_new <- "Normalised gene of interest data (new)" 
        par(mfrow=c(1,1))
        convergence.plot(new_data[,4], plotmain_new)
        par(mfrow=c(1,1))
        
      } else {
        print("Can not plot new data and new reference data set if old reference data set was not loaded")
      }
    }
  })
 
  ####File input of samples (GOI) to compare and basic stats calculation
  
  output$NormNew <- renderUI({
    answer <- as.integer(input$file_choice)
    if(answer == 3){
      fileInput("Normdata", label = "Load data for testing")}
  })
  
  ####File input of samples (GOI) to compare table
  
  Compare <- reactive({
    input_answer2 <- as.integer(input$Input_type2)
    if (input_answer2==1){
      loganswer <- as.integer(input$logtransform)
      if(loganswer==3){
        answer <- as.integer(input$file_choice)
        if(answer == 1){
          data <- GOI_norm()$resultsmat
        } else if (answer == 2){
          data <- GOI_norm()$newdb
        } else {
          fileName <-input$Normdata
          validate(
            need(fileName != "", label = "'Other file' option selected. File")
          )
          data <- read.csv(fileName$datapath)
        }
        nolog <- data[,4]
        datalog <- cbind.data.frame(data[,1], nolog, data[,7])
        colnames(datalog) <- c("Name","Normalised GOI", "Note")
        return(datalog)
        
      } else if (loganswer==2){
        answer <- as.integer(input$file_choice)
        if(answer == 1){
          data <- GOI_norm()$resultsmat
        } else if (answer == 2){
          data <- GOI_norm()$newdb
        } else {
          fileName <-input$Normdata
          validate(
            need(fileName != "", label = "'Other file' option selected. File")
          )
          data <- read.csv(fileName$datapath)
        }
        logten <- log10(data[,4])
        datalog <- cbind.data.frame(data[,1], logten, data[,7])
        colnames(datalog) <- c("Name","Normalised GOI (log base 10)", "Note")
        return(datalog)
        
      } else {
        answer <- as.integer(input$file_choice)
        if(answer == 1){
          data <- GOI_norm()$resultsmat
        } else if (answer == 2){
          data <- GOI_norm()$newdb
        } else {
          fileName <-input$Normdata
          validate(
            need(fileName != "", label = "'Other file' option selected. File")
          )
          data <- read.csv(fileName$datapath)
        }
        lognat <- log(data[,4])
        datalog <- cbind.data.frame(data[,1], lognat, data[,7])
        colnames(datalog) <- c("Name","Normalised GOI (Natural log)", "Note")
        return(datalog)
        
      }
    } else {
      loganswer <- as.integer(input$logtransform)
      if(loganswer==3){
        answer <- as.integer(input$file_choice)
        if(answer == 1){
          data <- GOI_norm_m()$resultsmat
        } else if (answer == 2){
          data <- GOI_norm_m()$newdb
        } else {
          fileName <-input$Normdata
          validate(
            need(fileName != "", label = "'Other file' option selected. File")
          )
          data <- read.csv(fileName$datapath)
        }
        nolog <- data[,4]
        datalog <- cbind.data.frame(data[,1], nolog, data[,7])
        colnames(datalog) <- c("Name","Normalised GOI", "Note")
        return(datalog)
        
      } else if (loganswer==2){
        answer <- as.integer(input$file_choice)
        if(answer == 1){
          data <- GOI_norm_m()$resultsmat
        } else if (answer == 2){
          data <- GOI_norm_m()$newdb
        } else {
          fileName <-input$Normdata
          validate(
            need(fileName != "", label = "'Other file' option selected. File")
          )
          data <- read.csv(fileName$datapath)
        }
        logten <- log10(data[,4])
        datalog <- cbind.data.frame(data[,1], logten, data[,7])
        colnames(datalog) <- c("Name","Normalised GOI (log base 10)", "Note")
        return(datalog)
        
      } else {
        answer <- as.integer(input$file_choice)
        if(answer == 1){
          data <- GOI_norm_m()$resultsmat
        } else if (answer == 2){
          data <- GOI_norm_m()$newdb
        } else {
          fileName <-input$Normdata
          validate(
            need(fileName != "", label = "'Other file' option selected. File")
          )
          data <- read.csv(fileName$datapath)
        }
        lognat <- log(data[,4])
        datalog <- cbind.data.frame(data[,1], lognat, data[,7])
        colnames(datalog) <- c("Name","Normalised GOI (Natural log)", "Note")
        return(datalog)
        
      }
    }
    
  })
  
  output$table1  <-  DT::renderDataTable({
    Compare()
  })

  
  ####Output of rows selected per group
  
  output$table1_2 <-  renderPrint({
    group1_rows  <-  input$table1_rows_selected
    if (length(group1_rows)) {
      cat('These rows were selected for group 1:\n\n')
      cat(group1_rows, sep = ', ')
    }
  })
  
  ####File input of samples (GOI) to compare table
  
  output$table2 <-  DT::renderDataTable({
    Compare()
  })
  
  ####Output of rows selected per group
  
  output$table2_1 <-  renderPrint({
    group2_rows <-  input$table2_rows_selected
    if (length(group2_rows)) {
      cat('These rows were selected for group 2:\n\n')
      cat(group2_rows, sep = ', ')
    }
  })
  
  output$table3  <-  DT::renderDataTable({
    groupnum <- as.integer(input$num_groups)
    if (is.null(groupnum)) 
      return(NULL)
    if(groupnum >= 3){
    Compare()
    }
  })
  
  
  ####Output of rows selected per group
  
  output$table3_1 <-  renderPrint({
    groupnum <- as.integer(input$num_groups)
    if (is.null(groupnum)) 
      return(NULL)
    if(groupnum >= 3 ){
      group3_rows  <-  input$table3_rows_selected
      if (length(group3_rows)) {
        cat('These rows were selected for group 3:\n\n')
        cat(group3_rows, sep = ', ')
      }
    } else {
      print("Group not selected")
    }
  })
  
  output$table4  <-  DT::renderDataTable({
    groupnum <- as.integer(input$num_groups)
    if (is.null(groupnum)) 
      return(NULL)
    if(groupnum >= 4){
      Compare()
    }
  })
  
  
  ####Output of rows selected per group
  
  output$table4_1 <-  renderPrint({
    groupnum <- as.integer(input$num_groups)
    if (is.null(groupnum)) 
      return(NULL)
    if(groupnum >= 4 ){
      group4_rows  <-  input$table4_rows_selected
      if (length(group4_rows)) {
        cat('These rows were selected for group 4:\n\n')
        cat(group4_rows, sep = ', ')
      }
    } else {
      print("Group not selected")
    }
  })
  
  output$table5  <-  DT::renderDataTable({
    groupnum <- as.integer(input$num_groups)
    if (is.null(groupnum)) 
      return(NULL)
    if(groupnum == 5){
      Compare()
    }
  })
  
  
  ####Output of rows selected per group
  
  output$table5_1 <-  renderPrint({
    groupnum <- as.integer(input$num_groups)
    if (is.null(groupnum)) 
      return(NULL)
    if(groupnum == 5 ){
      group5_rows  <-  input$table5_rows_selected
      if (length(group5_rows)) {
        cat('These rows were selected for group 5:\n\n')
        cat(group5_rows, sep = ', ')
      }
    } else {
      print("Group not selected")
    }
  })
  
  ####Subsetting of values to work with
  
  bs <- reactive({
    data2 <- Compare()
    groupnum <- as.integer(input$num_groups)
    if (groupnum == 2 ){
    group1_rows  <-  input$table1_rows_selected
    if(is.null(group1_rows))
      return(NULL)
    g1 <- data2[group1_rows,2]
    group2_rows <-  input$table2_rows_selected
    g2 <- data2[group2_rows,2]
    result1 <- describe(g1)[c(2,3,4,5,8,9,10,13)]
    result2 <- describe(g2)[c(2,3,4,5,8,9,10,13)]
    row.names(result1) <- "Group 1  "
    row.names(result2) <- "Group 2  "
    return(list(result1, result2))
    } else if (groupnum == 3){
      group1_rows  <-  input$table1_rows_selected
      if(is.null(group1_rows))
        return(NULL)
      g1 <- data2[group1_rows,2]
      group2_rows <-  input$table2_rows_selected
      g2 <- data2[group2_rows,2]
      group3_rows <-  input$table3_rows_selected
      g3 <- data2[group3_rows,2]
      result1 <- describe(g1)[c(2,3,4,5,8,9,10,13)]
      result2 <- describe(g2)[c(2,3,4,5,8,9,10,13)]
      result3 <- describe(g3)[c(2,3,4,5,8,9,10,13)]
      row.names(result1) <- "Group 1  "
      row.names(result2) <- "Group 2  "
      row.names(result3) <- "Group 3  "
      return(list(result1, result2, result3))
    } else if (groupnum == 4){
      group1_rows  <-  input$table1_rows_selected
      if(is.null(group1_rows))
        return(NULL)
      g1 <- data2[group1_rows,2]
      group2_rows <-  input$table2_rows_selected
      g2 <- data2[group2_rows,2]
      group3_rows <-  input$table3_rows_selected
      g3 <- data2[group3_rows,2]
      group4_rows <-  input$table4_rows_selected
      g4 <- data2[group4_rows,2]
      result1 <- describe(g1)[c(2,3,4,5,8,9,10,13)]
      result2 <- describe(g2)[c(2,3,4,5,8,9,10,13)]
      result3 <- describe(g3)[c(2,3,4,5,8,9,10,13)]
      result4 <- describe(g4)[c(2,3,4,5,8,9,10,13)]
      row.names(result1) <- "Group 1  "
      row.names(result2) <- "Group 2  "
      row.names(result3) <- "Group 3  "
      row.names(result4) <- "Group 4  "
      return(list(result1, result2, result3, result4))
    } else if (groupnum == 5) {
      group1_rows  <-  input$table1_rows_selected
      if(is.null(group1_rows))
        return(NULL)
      g1 <- data2[group1_rows,2]
      group2_rows <-  input$table2_rows_selected
      g2 <- data2[group2_rows,2]
      group3_rows <-  input$table3_rows_selected
      g3 <- data2[group3_rows,2]
      group4_rows <-  input$table4_rows_selected
      g4 <- data2[group4_rows,2]
      group5_rows <-  input$table5_rows_selected
      g5 <- data2[group5_rows,2]
      result1 <- describe(g1)[c(2,3,4,5,8,9,10,13)]
      result2 <- describe(g2)[c(2,3,4,5,8,9,10,13)]
      result3 <- describe(g3)[c(2,3,4,5,8,9,10,13)]
      result4 <- describe(g4)[c(2,3,4,5,8,9,10,13)]
      result5 <- describe(g5)[c(2,3,4,5,8,9,10,13)]
      row.names(result1) <- "Group 1  "
      row.names(result2) <- "Group 2  "
      row.names(result3) <- "Group 3  "
      row.names(result4) <- "Group 4  "
      row.names(result5) <- "Group 5  "
      return(list(result1, result2, result3, result4, result5))
    }
  })
  
  groups_selected <- reactive({
    data2 <- Compare()
    groupnum <- as.integer(input$num_groups)
    if(groupnum == 2 ){
      group1_rows  <-  input$table1_rows_selected
      if(is.null(group1_rows))
        return(NULL)
      g1 <- data2[group1_rows,2]
      group2_rows <-  input$table2_rows_selected
      g2 <- data2[group2_rows,2]
      group_results <- list(g1=g1, g2=g2)
    } else if (groupnum == 3){
      group1_rows  <-  input$table1_rows_selected
      if(is.null(group1_rows))
        return(NULL)
      g1 <- data2[group1_rows,2]
      group2_rows <-  input$table2_rows_selected
      g2 <- data2[group2_rows,2]
      group3_rows <-  input$table3_rows_selected
      g3 <- data2[group3_rows,2]
      group_results <- list(g1=g1, g2=g2, g3=g3)
    } else if (groupnum == 4){
      group1_rows  <-  input$table1_rows_selected
      if(is.null(group1_rows))
        return(NULL)
      g1 <- data2[group1_rows,2]
      group2_rows <-  input$table2_rows_selected
      g2 <- data2[group2_rows,2]
      group3_rows <-  input$table3_rows_selected
      g3 <- data2[group3_rows,2]
      group4_rows <-  input$table4_rows_selected
      g4 <- data2[group4_rows,2]
      group_results <- list(g1=g1, g2=g2, g3=g3, g4=g4)
    } else if (groupnum == 5) {
      group1_rows  <-  input$table1_rows_selected
      if(is.null(group1_rows))
        return(NULL)
      g1 <- data2[group1_rows,2]
      group2_rows <-  input$table2_rows_selected
      g2 <- data2[group2_rows,2]
      group3_rows <-  input$table3_rows_selected
      g3 <- data2[group3_rows,2]
      group4_rows <-  input$table4_rows_selected
      g4 <- data2[group4_rows,2]
      group5_rows <-  input$table5_rows_selected
      g5 <- data2[group5_rows,2]
      group_results <- list(g1=g1, g2=g2, g3=g3, g4=g4, g5=g5)
    }
  })
  
  
  ####Function to plot distribution per group
  
  makedistPlot <- function(){
    
    groupnum <- as.integer(input$num_groups)
    if(groupnum == 2 ){
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      
      
      p1 <- density(g1)
      p2 <- density(g2)
      
      max_y <- max(c(p1$y,p2$y))
      min_y <- min(c(p1$y,p2$y))
      
      plot(p1, las=1, xlab = "Group 1 is expressed in blue and Group 2 in red. Vertical lines show the mean.",
           main = "", col = "blue", ylim=c(min_y,max_y))
      lines(p2, las=1, xlab = "", main = "", col = "red")
      
      abline(v = mean(g1), col = "blue", lwd = 2)
      abline(v = mean(g2), col = "red", lwd = 2)
      
    } else if (groupnum == 3){
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      g3 <- groups_selected()$g3
      
      p1 <- density(g1)
      p2 <- density(g2)
      p3 <- density(g3)
      
      max_y <- max(c(p1$y,p2$y,p3$y))
      min_y <- min(c(p1$y,p2$y,p3$y))
      
      plot(p1, las=1, xlab = "Group 1 is expressed in blue; Group 2 in red and Group 3 in green. Vertical lines show the mean.",
           main = "", col = "blue", ylim=c(min_y,max_y))
      lines(p2, las=1, xlab = "", main = "", col = "red")
      lines(p3, las=1, xlab = "", main = "", col = "green")
      
      abline(v = mean(g1), col = "blue", lwd = 2)
      abline(v = mean(g2), col = "red", lwd = 2)
      abline(v = mean(g3), col = "green", lwd = 2)
      
    } else if (groupnum == 4){
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      g3 <- groups_selected()$g3
      g4 <- groups_selected()$g4
    
      
      p1 <- density(g1)
      p2 <- density(g2)
      p3 <- density(g3)
      p4 <- density(g4)
      
      max_y <- max(c(p1$y,p2$y,p3$y,p4$y))
      min_y <- min(c(p1$y,p2$y,p3$y,p4$y))
      
      plot(p1, las=1, xlab = "Group 1 is expressed in blue; Group 2 in red; Group 3 in green and Group 4 in orange. Vertical lines show the mean.",
           main = "", col = "blue", ylim=c(min_y,max_y))
      lines(p2, las=1, xlab = "", main = "", col = "red")
      lines(p3, las=1, xlab = "", main = "", col = "green")
      lines(p4, las=1, xlab = "", main = "", col = "Orange")
      
      abline(v = mean(g1), col = "blue", lwd = 2)
      abline(v = mean(g2), col = "red", lwd = 2)
      abline(v = mean(g3), col = "green", lwd = 2)
      abline(v = mean(g4), col = "Orange", lwd = 2)
      
    } else if (groupnum == 5) {
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      g3 <- groups_selected()$g3
      g4 <- groups_selected()$g4
      g5 <- groups_selected()$g5
      
      
      p1 <- density(g1)
      p2 <- density(g2)
      p3 <- density(g3)
      p4 <- density(g4)
      p5 <- density(g5)
      
      max_y <- max(c(p1$y,p2$y, p3$y, p4$y, p5$y))
      min_y <- min(c(p1$y,p2$y, p3$y, p4$y, p5$y))
      
      plot(p1, las=1, xlab = "Group 1 is expressed in blue; Group 2 in red; Group 3 in green; Group 4 in orange and Group 5 in purple. Vertical lines show the mean.",
           main = "", col = "blue", ylim=c(min_y,max_y))
      lines(p2, las=1, xlab = "", main = "", col = "red")
      lines(p3, las=1, xlab = "", main = "", col = "green")
      lines(p4, las=1, xlab = "", main = "", col = "Orange")
      lines(p5, las=1, xlab = "", main = "", col = "Purple")
      
      abline(v = mean(g1), col = "blue", lwd = 2)
      abline(v = mean(g2), col = "red", lwd = 2)
      abline(v = mean(g3), col = "green", lwd = 2)
      abline(v = mean(g4), col = "Orange", lwd = 2)
      abline(v = mean(g5), col = "Purple", lwd = 2)
    }
  }
  
  

  
  ####Plotting distribution per group
  
  output$distPlot <- renderPlot({
    print(makedistPlot())
  })
  
  ####Function for boxplot   
  
  makeboxPlot <- function(){
    
    groupnum <- as.integer(input$num_groups)
    if(groupnum == 2 ){
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      
      score <- c(g1, g2)
      group <- factor(c(rep("Group 1", length(g1)), rep("Group 2", length(g2))))
      
      boxplot(score ~ group, las=1, ylab= "Means and +/-1 SDs are displayed in red.")
      
      beeswarm(score ~ group, col = 4, pch = 16, add = TRUE)
      
      points(1.2, mean(g1), pch = 18, col = "red", cex = 2)
      arrows(1.2, mean(g1), 1.2, mean(g1) + sd(g1), length = 0.1, angle = 45, col = "red")
      arrows(1.2, mean(g1), 1.2, mean(g1) - sd(g1), length = 0.1, angle = 45, col = "red")
      
      points(2.2, mean(g2), pch = 18, col = "red", cex = 2)
      arrows(2.2, mean(g2), 2.2, mean(g2) + sd(g2), length = 0.1, angle = 45, col = "red")
      arrows(2.2, mean(g2), 2.2, mean(g2) - sd(g2), length = 0.1, angle = 45, col = "red")
      
      
    } else if (groupnum == 3){
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      g3 <- groups_selected()$g3
      
      score <- c(g1, g2, g3)
      group <- factor(c(rep("Group 1", length(g1)), rep("Group 2", length(g2)), rep("Group 3", length(g3))))
      
      boxplot(score ~ group, las=1, ylab= "Means and +/-1 SDs are displayed in red.")
      
      beeswarm(score ~ group, col = 4, pch = 16, add = TRUE)
      
      points(1.2, mean(g1), pch = 18, col = "red", cex = 2)
      arrows(1.2, mean(g1), 1.2, mean(g1) + sd(g1), length = 0.1, angle = 45, col = "red")
      arrows(1.2, mean(g1), 1.2, mean(g1) - sd(g1), length = 0.1, angle = 45, col = "red")
      
      points(2.2, mean(g2), pch = 18, col = "red", cex = 2)
      arrows(2.2, mean(g2), 2.2, mean(g2) + sd(g2), length = 0.1, angle = 45, col = "red")
      arrows(2.2, mean(g2), 2.2, mean(g2) - sd(g2), length = 0.1, angle = 45, col = "red")
      
      points(3.2, mean(g3), pch = 18, col = "red", cex = 2)
      arrows(3.2, mean(g3), 3.2, mean(g3) + sd(g3), length = 0.1, angle = 45, col = "red")
      arrows(3.2, mean(g3), 3.2, mean(g3) - sd(g3), length = 0.1, angle = 45, col = "red")
      
    } else if (groupnum == 4){
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      g3 <- groups_selected()$g3
      g4 <- groups_selected()$g4
      
      score <- c(g1, g2, g3, g4)
      group <- factor(c(rep("Group 1", length(g1)), rep("Group 2", length(g2)), rep("Group 3", length(g3)), rep("Group 4", length(g4))))
      
      boxplot(score ~ group, las=1, ylab= "Means and +/-1 SDs are displayed in red.")
      
      beeswarm(score ~ group, col = 4, pch = 16, add = TRUE)
      
      points(1.2, mean(g1), pch = 18, col = "red", cex = 2)
      arrows(1.2, mean(g1), 1.2, mean(g1) + sd(g1), length = 0.1, angle = 45, col = "red")
      arrows(1.2, mean(g1), 1.2, mean(g1) - sd(g1), length = 0.1, angle = 45, col = "red")
      
      points(2.2, mean(g2), pch = 18, col = "red", cex = 2)
      arrows(2.2, mean(g2), 2.2, mean(g2) + sd(g2), length = 0.1, angle = 45, col = "red")
      arrows(2.2, mean(g2), 2.2, mean(g2) - sd(g2), length = 0.1, angle = 45, col = "red")
      
      points(3.2, mean(g3), pch = 18, col = "red", cex = 2)
      arrows(3.2, mean(g3), 3.2, mean(g3) + sd(g3), length = 0.1, angle = 45, col = "red")
      arrows(3.2, mean(g3), 3.2, mean(g3) - sd(g3), length = 0.1, angle = 45, col = "red")
      
      points(4.2, mean(g4), pch = 18, col = "red", cex = 2)
      arrows(4.2, mean(g4), 4.2, mean(g4) + sd(g4), length = 0.1, angle = 45, col = "red")
      arrows(4.2, mean(g4), 4.2, mean(g4) - sd(g4), length = 0.1, angle = 45, col = "red")
      
    } else if (groupnum == 5) {
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      g3 <- groups_selected()$g3
      g4 <- groups_selected()$g4
      g5 <- groups_selected()$g5
      
      score <- c(g1, g2, g3, g4, g5)
      group <- factor(c(rep("Group 1", length(g1)), rep("Group 2", length(g2)), rep("Group 3", length(g3)), rep("Group 4", length(g4)), rep("Group 5", length(g5))))
      
      boxplot(score ~ group, las=1, ylab= "Means and +/-1 SDs are displayed in red.")
      
      beeswarm(score ~ group, col = 4, pch = 16, add = TRUE)
      
      points(1.2, mean(g1), pch = 18, col = "red", cex = 2)
      arrows(1.2, mean(g1), 1.2, mean(g1) + sd(g1), length = 0.1, angle = 45, col = "red")
      arrows(1.2, mean(g1), 1.2, mean(g1) - sd(g1), length = 0.1, angle = 45, col = "red")
      
      points(2.2, mean(g2), pch = 18, col = "red", cex = 2)
      arrows(2.2, mean(g2), 2.2, mean(g2) + sd(g2), length = 0.1, angle = 45, col = "red")
      arrows(2.2, mean(g2), 2.2, mean(g2) - sd(g2), length = 0.1, angle = 45, col = "red")
      
      points(3.2, mean(g3), pch = 18, col = "red", cex = 2)
      arrows(3.2, mean(g3), 3.2, mean(g3) + sd(g3), length = 0.1, angle = 45, col = "red")
      arrows(3.2, mean(g3), 3.2, mean(g3) - sd(g3), length = 0.1, angle = 45, col = "red")
      
      points(4.2, mean(g4), pch = 18, col = "red", cex = 2)
      arrows(4.2, mean(g4), 4.2, mean(g4) + sd(g4), length = 0.1, angle = 45, col = "red")
      arrows(4.2, mean(g4), 4.2, mean(g4) - sd(g4), length = 0.1, angle = 45, col = "red")
      
      points(5.2, mean(g5), pch = 18, col = "red", cex = 2)
      arrows(5.2, mean(g5), 5.2, mean(g5) + sd(g5), length = 0.1, angle = 45, col = "red")
      arrows(5.2, mean(g5), 5.2, mean(g5) - sd(g5), length = 0.1, angle = 45, col = "red")
    }
    
  }
  
  ####Plotting boxplot   
  
  output$boxPlot <- renderPlot({
    print(makeboxPlot())
  })
  
  ####Function for normality testing
  
  testnorm <- reactive({
    groupnum <- as.integer(input$num_groups)
    if(groupnum == 2 ){
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      
      #group.1ks <- ks.test(scale(g1), "pnorm")
      group.1sh <- shapiro.test(g1)
      
      #group.2ks <- ks.test(scale(g2), "pnorm")
      group.2sh <- shapiro.test(g2)
      
      #return(list(Group.1 = group.1ks, Group.1 = group.1sh, Group.2 = group.2ks, Group.2 = group.2sh))
      return(list(Group.1 = group.1sh, Group.2 = group.2sh))
      
    } else if (groupnum == 3){
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      g3 <- groups_selected()$g3
      
      #group.1ks <- ks.test(scale(g1), "pnorm")
      group.1sh <- shapiro.test(g1)
      
      #group.2ks <- ks.test(scale(g2), "pnorm")
      group.2sh <- shapiro.test(g2)
      
      #group.3ks <- ks.test(scale(g3), "pnorm")
      group.3sh <- shapiro.test(g3)
      
      #return(list(Group.1 = group.1ks, Group.1 = group.1sh, Group.2 = group.2ks, Group.2 = group.2sh, Group.3 = group.3ks, Group.3 = group.3sh))
      return(list(Group.1 = group.1sh, Group.2 = group.2sh, Group.3 = group.3sh))
      
    } else if (groupnum == 4){
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      g3 <- groups_selected()$g3
      g4 <- groups_selected()$g4
      
      #group.1ks <- ks.test(scale(g1), "pnorm")
      group.1sh <- shapiro.test(g1)
      
      #group.2ks <- ks.test(scale(g2), "pnorm")
      group.2sh <- shapiro.test(g2)
      
      #group.3ks <- ks.test(scale(g3), "pnorm")
      group.3sh <- shapiro.test(g3)
      
      #group.4ks <- ks.test(scale(g4), "pnorm")
      group.4sh <- shapiro.test(g4)
      
      #return(list(Group.1 = group.1ks, Group.1 = group.1sh, Group.2 = group.2ks, Group.2 = group.2sh, Group.3 = group.3ks, Group.3 = group.3sh, Group.4 = group.4ks, Group.4 = group.4sh))
      return(list(Group.1 = group.1sh, Group.2 = group.2sh, Group.3 = group.3sh, Group.4 = group.4sh))
      
    } else if (groupnum == 5) {
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      g3 <- groups_selected()$g3
      g4 <- groups_selected()$g4
      g5 <- groups_selected()$g5
      
      #group.1ks <- ks.test(scale(g1), "pnorm")
      group.1sh <- shapiro.test(g1)
      
      #group.2ks <- ks.test(scale(g2), "pnorm")
      group.2sh <- shapiro.test(g2)
      
      #group.3ks <- ks.test(scale(g3), "pnorm")
      group.3sh <- shapiro.test(g3)
      
      #group.4ks <- ks.test(scale(g4), "pnorm")
      group.4sh <- shapiro.test(g4)
      
      #group.5ks <- ks.test(scale(g5), "pnorm")
      group.5sh <- shapiro.test(g5)
      
      #return(list(Group.1 = group.1ks, Group.1 = group.1sh, Group.2 = group.2ks, Group.2 = group.2sh, Group.3 = group.3ks, Group.3 = group.3sh, Group.4 = group.4ks, Group.4 = group.4sh, Group.5 = group.5ks, Group.5 = group.5sh))
      return(list(Group.1 = group.1sh, Group.2 = group.2sh, Group.3 = group.3sh, Group.4 = group.4sh, Group.5 = group.5sh))
      
      }
    
  })
  
  ####function for Levene test for homoscedasticity
  
  levene <- reactive({
    groupnum <- as.integer(input$num_groups)
    if(groupnum == 2 ){
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      
      score <- c(g1, g2)
      group <- factor(c(rep("Group 1", length(g1)), rep("Group 2", length(g2))))
      
      leveneTest(score, group, center=mean)
      
      
    } else if (groupnum == 3){
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      g3 <- groups_selected()$g3
      
      score <- c(g1, g2, g3)
      group <- factor(c(rep("Group 1", length(g1)), rep("Group 2", length(g2)), rep("Group 3", length(g3))))
      
      leveneTest(score, group, center=mean)
      
      
    } else if (groupnum == 4){
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      g3 <- groups_selected()$g3
      g4 <- groups_selected()$g4
      
      score <- c(g1, g2, g3, g4)
      group <- factor(c(rep("Group 1", length(g1)), rep("Group 2", length(g2)), rep("Group 3", length(g3)), rep("Group 4", length(g4))))
      
      leveneTest(score, group, center=mean)
      
    } else if (groupnum == 5) {
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      g3 <- groups_selected()$g3
      g4 <- groups_selected()$g4
      g5 <- groups_selected()$g5
      
      score <- c(g1, g2, g3, g4, g5)
      group <- factor(c(rep("Group 1", length(g1)), rep("Group 2", length(g2)), rep("Group 3", length(g3)), rep("Group 4", length(g4)), rep("Group 5", length(g5))))
      
      leveneTest(score, group, center=mean)
      
      
    } 
  })
  
  ####function for Bartlett's test for homoscedasticity
  
  bartlett <- reactive({
    groupnum <- as.integer(input$num_groups)
    if(groupnum == 2 ){
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      
      score <- c(g1, g2)
      group <- factor(c(rep("Group 1", length(g1)), rep("Group 2", length(g2))))
      
      bartlett.test(score~group)
      
    } else if (groupnum == 3){
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      g3 <- groups_selected()$g3
      
      score <- c(g1, g2, g3)
      group <- factor(c(rep("Group 1", length(g1)), rep("Group 2", length(g2)), rep("Group 3", length(g3))))
      
      bartlett.test(score~group)
      
    } else if (groupnum == 4){
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      g3 <- groups_selected()$g3
      g4 <- groups_selected()$g4
      
      score <- c(g1, g2, g3, g4)
      group <- factor(c(rep("Group 1", length(g1)), rep("Group 2", length(g2)), rep("Group 3", length(g3)), rep("Group 4", length(g4))))
      
      bartlett.test(score~group)
      
    } else if (groupnum == 5) {
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      g3 <- groups_selected()$g3
      g4 <- groups_selected()$g4
      g5 <- groups_selected()$g5
      
      score <- c(g1, g2, g3, g4, g5)
      group <- factor(c(rep("Group 1", length(g1)), rep("Group 2", length(g2)), rep("Group 3", length(g3)), rep("Group 4", length(g4)), rep("Group 5", length(g5))))
      
      bartlett.test(score~group)
      
    } 
  })
  
  ####Function for Parametric tests
  
  t <- reactive({
    groupnum <- as.integer(input$num_groups)
    if(groupnum == 2 ){
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      
      score <- c(g1, g2)
      group <- factor(c(rep("Group 1", length(g1)), rep("Group 2", length(g2))))
      
      normal.t <- t.test(score ~ group, var.equal=TRUE)
      Welch.t <- t.test(score ~ group, var.equal=FALSE)
      
      return(list(normal.t, Welch.t))
      
    } else {
      print("More than two groups selected for testing.")
    }
  })
  
  
  ####Output for basic statistics
  
  output$textarea.out <- renderPrint({
    bs()
    
  })
  
  
  ####Output for normality test
  
  output$testnorm.out <- renderPrint({
    testnorm()
    
  })
  
  ####Output for levene test
  
  output$levene.out <- renderPrint({
    levene()

  })
  
  output$Bartlett.out <- renderPrint({
    bartlett()
    
  })
  
  
  ####Output for parametric test
  output$t.out <- renderPrint({
    t()
    
  })
  
  output$KrWs <- renderPrint({
    groupnum <- as.integer(input$num_groups)
    if(groupnum == 2 ){
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      
      score <- c(g1, g2)
      group <- factor(c(rep("Group 1", length(g1)), rep("Group 2", length(g2))))
      
      kruskalW = kruskal.test(score ~ group)
      return(kruskalW)
      
      
    } else if (groupnum == 3){
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      g3 <- groups_selected()$g3
      
      score <- c(g1, g2, g3)
      group <- factor(c(rep("Group 1", length(g1)), rep("Group 2", length(g2)), rep("Group 3", length(g3))))
      
      kruskalW = kruskal.test(score ~ group)
      return(kruskalW)
      
      
    } else if (groupnum == 4){
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      g3 <- groups_selected()$g3
      g4 <- groups_selected()$g4
      
      score <- c(g1, g2, g3, g4)
      group <- factor(c(rep("Group 1", length(g1)), rep("Group 2", length(g2)), rep("Group 3", length(g3)), rep("Group 4", length(g4))))
      
      kruskalW = kruskal.test(score ~ group)
      return(kruskalW)
      
    } else if (groupnum == 5) {
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      g3 <- groups_selected()$g3
      g4 <- groups_selected()$g4
      g5 <- groups_selected()$g5
      
      score <- c(g1, g2, g3, g4, g5)
      group <- factor(c(rep("Group 1", length(g1)), rep("Group 2", length(g2)), rep("Group 3", length(g3)), rep("Group 4", length(g4)), rep("Group 5", length(g5))))
      
      kruskalW = kruskal.test(score ~ group)
      return(kruskalW)
    }

  })
  
  
  ####Output for non-parametric test
  
  output$wilcox <- renderPrint({
    groupnum <- as.integer(input$num_groups)
    if(groupnum == 2 ){
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      
      score <- c(g1, g2)
      group <- factor(c(rep("Group 1", length(g1)), rep("Group 2", length(g2))))
      
      wilcox <- wilcox.test(score ~ group, correct=FALSE)
      
      return(wilcox)
      
    } else {
      print("More than two groups selected for testing.")
    }
  })
  
  output$anova <- renderPrint({
    groupnum <- as.integer(input$num_groups)
    if(groupnum == 2 ){
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      
      score <- c(g1, g2)
      group <- factor(c(rep("Group 1", length(g1)), rep("Group 2", length(g2))))
      
      linmodel = lm(score ~ group)
      anova_results <- anova(linmodel)
      return(anova_results)
      
    }else if (groupnum == 3){
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      g3 <- groups_selected()$g3
      
      score <- c(g1, g2, g3)
      group <- factor(c(rep("Group 1", length(g1)), rep("Group 2", length(g2)), rep("Group 3", length(g3))))
      
      linmodel = lm(score ~ group)
      anova_results <- anova(linmodel)
      return(anova_results)
      
      
    } else if (groupnum == 4){
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      g3 <- groups_selected()$g3
      g4 <- groups_selected()$g4
      
      score <- c(g1, g2, g3, g4)
      group <- factor(c(rep("Group 1", length(g1)), rep("Group 2", length(g2)), rep("Group 3", length(g3)), rep("Group 4", length(g4))))
      
      linmodel = lm(score ~ group)
      anova_results <- anova(linmodel)
      return(anova_results)
      
    } else if (groupnum == 5) {
      g1 <- groups_selected()$g1
      if(is.null(g1))
        return(NULL)
      g2 <- groups_selected()$g2
      g3 <- groups_selected()$g3
      g4 <- groups_selected()$g4
      g5 <- groups_selected()$g5
      
      score <- c(g1, g2, g3, g4, g5)
      group <- factor(c(rep("Group 1", length(g1)), rep("Group 2", length(g2)), rep("Group 3", length(g3)), rep("Group 4", length(g4)), rep("Group 5", length(g5))))
      
      linmodel = lm(score ~ group)
      anova_results <- anova(linmodel)
      return(anova_results)
    }
  })
  
})

