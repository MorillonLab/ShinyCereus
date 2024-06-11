#!/usr/bin/env Rscript

library(shiny)
library(shinydashboard)
library(shinyFiles)
library(shinyBS)
library(shinycssloaders)
library(shinyDirectoryInput)
library(DT)
library(ggplot2)
library(tidyverse)
library(circlize)
library(cluster)
library(data.table)
library(dplyr)
library(reshape2)
library(plyr)
library(reactlog)
library(fst)
library(JBrowseR)
library(Biostrings)
library(seqinr)
library(Rsamtools)
library(markdown)


{
  ui <- dashboardPage(
    dashboardHeader(),
    
    dashboardSidebar(
      
      conditionalPanel(condition="output.loadDir!=0",style = "display: none;",
                       
                       h3(" Custom Filters :"),
                       uiOutput("filters_panel")
                       
      )
    ),
    
    dashboardBody(
      tabsetPanel(
        id = "tabs",
        tabPanel(
          title = "Main",
          #upload directory (server side)
          shinyDirButton("upload_directory",label="Upload your directory",title="Upload your directory", multiple=FALSE),br(),
          
          #if chosen_directory is not empty, show content
          conditionalPanel(condition="output.chosen_directory!=' '",
                           style = "display: none;",
                           div(style="width:100%",verbatimTextOutput("chosen_directory"))),
          br(),
          
          #if chosen_directory is not empty, and files are not yet validated, display the button that will allow the validation of the input directory and its content
          conditionalPanel(condition="output.loadDir==0 && output.chosen_directory!=' '",
                           style = "display: none;",
                           actionButton("validate_directory","validate",icon("thumbs-up"))),
          
          conditionalPanel(condition="(output.loadDir!=0 && output.chosen_directory!=' ') || (output.loadDir==0 && output.chosen_directory!=' ' && output.validation_message!='-> all files are validated ! waiting for application to launch...')",
                           style = "display: none;",
                           div(style="width:100%",verbatimTextOutput("validation_message"))),
          
          #if loadDir different from 0, load all the parts of the dashboard
          conditionalPanel(condition="output.loadDir!=0",style = "display: none;",
                           #tabPanel("Summary2",fluidRow(column(1,uiOutput("X_axis",inline = T),uiOutput("Y_axis"),uiOutput("colored_parameter",inline = T),uiOutput("showscatterplot"),inline = T)))
                           #column(1,uiOutput("scatterplot",inline = T))))
                           
                           uiOutput("scatterplot_panel"),
                           uiOutput("download_plots_panel"),
                           uiOutput("contig_table_panel"),
                           uiOutput("browser_panel")
          )
        ),
        tabPanel(
          title = "Readme",
            includeMarkdown("README.md")
        )
      )
    )
  )
}

server <- function(input, output, session) {
  
  plotsReady <- reactiveVal(FALSE)
  

  session$onSessionEnded(function() {
    data_server$stop_server()
    selected_directory <- isolate(directory_uploaded$value)
    subbam_files <- list.files(selected_directory, pattern = "^sub_bam.*$", full.names = TRUE)

    if (length(subbam_files) > 0) {
      file.remove(subbam_files)
      cat("Done: subbam files removed.\n")
    }
    else {
      cat("No subbam files found.\n")
    }

    rm(list = ls())
    gc()
    stopApp()
    .rs.restartR()
  })
  
  

  
  
  # onStop(function() {
  # 
  #   data_server$stop_server()
  #   selected_directory <- isolate(directory_uploaded$value)
  #   subbam_files <- list.files(selected_directory, pattern = "^sub_bam.*$", full.names = TRUE)
  # 
  #   if (length(subbam_files) > 0) {
  #     file.remove(subbam_files)
  #     cat("Done: subbam files removed.\n")
  #     }
  #   else {
  #     cat("No subbam files found.\n")
  #     }
  # 
  #   # rm(list = ls())
  #   # gc()
  # 
  # })
  
  
  #function to get peptides from nucleotide sequences
  get3framesPeptidesFromNuc<-function(myseq=""){
    
    #check if thesupplied sequence is not empty
    if(myseq!=""){
      
      #split it in characters in order to give it to seqinr::translate()
      #seqinr::translate() function needs lowercase characters...
      nuc<-tolower(unlist(strsplit(myseq,split="")))
      
      #translate the 1st frame, take the peptide until the stop ("*") after collapsing the characters, and convert it back to characters
      peptide_frame1<-unlist(strsplit(unlist(strsplit(paste(seqinr::translate(seq = nuc,frame=0),collapse=""),"\\*"))[1],""))
      #if we have a true stop, add it in the list of character
      if(grepl("\\*",paste(seqinr::translate(seq = tolower(nuc),frame=0),collapse=""))){
        peptide_frame1<-c(peptide_frame1,"*")
        
      }
      
      #translate the 2nd frame, take the peptide until the stop ("*") after collapsing the characters, and convert it back to characters
      peptide_frame2<-unlist(strsplit(unlist(strsplit(paste(seqinr::translate(seq = nuc,frame=1),collapse=""),"\\*"))[1],""))
      #if we have a true stop, add it in the list of character
      if(grepl("\\*",paste(seqinr::translate(seq = tolower(nuc),frame=1),collapse=""))){
        peptide_frame2<-c(peptide_frame2,"*")
        
      }
      
      #translate the 3rd frame, take the peptide until the stop ("*") after collapsing the characters, and convert it back to characters
      peptide_frame3<-unlist(strsplit(unlist(strsplit(paste(seqinr::translate(seq = nuc,frame=2),collapse=""),"\\*"))[1],""))
      #if we have a true stop, add it in the list of character
      if(grepl("\\*",paste(seqinr::translate(seq = tolower(nuc),frame=2),collapse=""))){
        peptide_frame3<-c(peptide_frame3,"*")
        
      }
      
      #dataframe that will have in column, the linear position of each nucleotide, the nucleotides, the peptides in the 3 frames, and additional positions values in order to plot
      matching_nuc_peptide<-data.frame(pos=1:length(nuc),nuc=toupper(nuc),frame1="",frame2="",frame3="")
      
      
      #for aa in frame1
      #initialize the position to display the aa
      j<-0
      for(i in 1:length(nuc)){
        
        #at each 3 nuc, insert the aa
        if(i%%3==0){
          j<-j+1
          
          if(!is.na(peptide_frame1[j])){
            
            matching_nuc_peptide$frame1[i]<-peptide_frame1[j]
            
          }
          
        }
      }
      
      #for aa in frame1
      #initialize the position to display the aa
      j<-0
      for(i in 2:length(nuc)){
        
        #at each 3 nuc, insert the aa
        if(i%%3==0){
          j<-j+1
          
          if(!is.na(peptide_frame2[j])){
            
            matching_nuc_peptide$frame2[i+1]<-peptide_frame2[j]
            
          }
          
        }
      }
      
      #for aa in frame3
      #initialize the position to display the aa
      j<-0
      for(i in 3:length(nuc)){
        
        #at each 3 nuc, insert the aa
        if(i%%3==0){
          j<-j+1
          
          if(!is.na(peptide_frame3[j])){
            
            matching_nuc_peptide$frame3[i+2]<-peptide_frame3[j]
            
          }
          
        }
      }
      
      
      
      #this column will contain values that will allow to break the sequence in "floors", as very long sequences can't be displayed in one line
      #initialization of the floor
      matching_nuc_peptide$floor<-0
      j<-0
      #loop across the number of nucleotides
      for(i in 1:nrow(matching_nuc_peptide)){
        
        matching_nuc_peptide[i,"floor"]<-j
        
        #each time we reach 35 nuleotides, we create a new floor (= new line), and we give it a value (lower value than the one before)
        if(i%%35==0){
          
          j=j-1
        }
        
      }
      
      #final dataframe that will be used by ggplot
      #here, we will display the nucleotide position (display_pos), at each 5 nucleotides, for each "floor"
      final_frame<-data.frame()
      #loop across the number of floors
      for(i in unique(matching_nuc_peptide$floor)){
        
        #take the sub dataframe corresponding to one floor
        tmp_frame<-matching_nuc_peptide[which(matching_nuc_peptide$floor==i),]
        
        #compute the position of each nucleotide for the floor (rel_pos)
        tmp_frame$rel_pos<-1:nrow(tmp_frame)
        
        #initialize the position to display
        tmp_frame$display_pos<-""
        
        #at each 5 nucleotides, store the position to display (real position on the total nucleotide sequence)
        for(j in 1:nrow(tmp_frame)){
          
          if(tmp_frame$pos[j]%%5==0){
            
            tmp_frame$display_pos[j]<-tmp_frame$pos[j]
            
          }
          
        }
        
        #store the result in the final dataframe
        final_frame<-rbind(final_frame,tmp_frame)
        
      }
      
      #function to obtain a white background
      white_background<-ggplot2::theme(axis.line = element_blank(),
                                       panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(),
                                       panel.background = element_blank(),
                                       axis.ticks = element_blank(),
                                       axis.text = element_blank(),
                                       axis.title = element_blank())
      
      ggplot(final_frame) +
        
        #plot the nucleotide sequence
        geom_text(data=final_frame,aes(x=rel_pos,y=0+floor,label=nuc),size=4,color="black")+
        
        #plot the peptide in frame1
        geom_text(data=final_frame,aes(x=rel_pos,y=-0.2+floor,label=frame1),size=3,color="blue")+
        #plot the peptide in frame2
        geom_text(data=final_frame,aes(x=rel_pos,y=-0.3+floor,label=frame2),size=3,color="red")+
        #plot the peptide in frame3
        geom_text(data=final_frame,aes(x=rel_pos,y=-0.4+floor,label=frame3),size=3,color="darkgreen")+
        #plot the nucleotide position
        geom_text(data=final_frame,aes(x=rel_pos,y=-0+floor+0.2,label=display_pos),size=3,color="purple")+
        
        #extend the limits according to the number of floors we have
        scale_y_continuous(limits = c(min(final_frame$floor)-2,0.2))+
        
        #add white background
        white_background
      
      #if there's no supplied sequence, say it
    }else{
      
      cat("no sequence supplied !\n")
      
    }
  }
  
  
  
  #reactive value when contigs data are loaded
  loadDir<-reactiveValues(value=0)
  directory_uploaded<-reactiveValues(value=" ")
  
  #return to UI the result of loadDir
  output$loadDir<-renderText({
    
    loadDir$value
    
    
  })
  
  outputOptions(output, 'loadDir', suspendWhenHidden=FALSE)
  
  #when upload_directory is triggered, put loadDir to 0
  observeEvent(input$upload_directory,{
    
    loadDir$value<-0
    
    
  })
  
  shinyDirChoose(input,"upload_directory",roots=getVolumes()(),session=session,restrictions=system.file(package='base'))
  
  #when upload is triggered, process the input
  observeEvent(input$upload_directory,{
    
    #store the selected directory
    directory_uploaded$value<-paste(as.character(parseDirPath(roots=getVolumes()(),input$upload_directory)),"/",sep="")
    
    if(directory_uploaded$value=="/" | directory_uploaded$value=="NA/"){
      
      directory_uploaded$value<-" "
    }
    
    
  })
  
  #show chosen directory
  output$chosen_directory<-renderText({
    
    if(directory_uploaded$value==" "){
      
      directory_uploaded$value
      
    }else{
      
      paste("directory : ",directory_uploaded$value)
      
    }
    
  })
  
  outputOptions(output, 'chosen_directory', suspendWhenHidden=FALSE)
  
  #When validate_directory is triggered, check the files in the input directory
  observeEvent(input$validate_directory,{
    
    
    
    if(length(directory_uploaded$value)<1){
      
      directory_uploaded$value<-" "
    }
    
    selected_directory<-directory_uploaded$value
    
    #if it exists, process the content
    if(dir.exists(selected_directory)){
      
      
      #initialize all the vectors
      bam_file=table_annotated=vcf_file=gff_file=reads_counts_cancer_Gtex_file=reads_counts_patients_file=cancer_normal_quantiles_file=genome_file=sample_order_file=cancer_cols_file=""
      
      #bam file of contigs
      bam_file<-list.files(selected_directory,pattern="Aligned-fixed.out.bam$")
      
      #annotated contigs
      table_annotated<-list.files(selected_directory,pattern="table_annotated.fst$")
      
      #vcf annotation (indel, insertions, snp)
      vcf_file<-list.files(selected_directory,pattern="clinvar_20210828_gencode_chromosomes.vcf.gz$")
      
      #annotation of known genes and transcripts
      gff_file<-list.files(selected_directory,pattern="gencode.v32.annotation_sorted_official_chromosomes.gff3.gz$")
      
      #counts for cancer and Gtex (normal tissue samples)
      reads_counts_cancer_Gtex_file<-list.files(selected_directory,pattern="reads_counts_cancer_Gtex.fst$")
      
      #counts for cancer patients
      reads_counts_patients_file<-list.files(selected_directory,pattern="reads_counts_patients.fst$")
      
      #quantiles comparing cancer cell lines samples vs normal samples
      cancer_normal_quantiles_file<-list.files(selected_directory,pattern="cancer_normal_quantiles.fst$")
      
      #quantiles comparing cancer patients samples vs normal samples
      patients_quantiles_file<-list.files(selected_directory,pattern="patients_quantiles.fst$")
      
      #human genome
      genome_file<-list.files(selected_directory,pattern="GRCh38.primary_assembly_genome_official_chromosomes.fa$")
      
      #cancer samples names (sorted according to the counts)
      sample_order_file<-list.files(selected_directory,pattern="sample_order.txt$")
      
      #prefix of the cancer samples
      cancer_cols_file<-list.files(selected_directory,pattern="cancer_cols.txt$")
      
      #json file for jbrowse
      jbrowse_gff_json_file<-list.files(selected_directory,pattern="_meta.json$")
      
      #vector that will store missing inputs
      missing_inputs<-c()
      
      if(length(bam_file)>0 ){
        
        #if files (and their index for some), are missing, store the vector name in missing_inputs
        if(length(bam_file)==""| !file.exists(paste(selected_directory,bam_file,".bai",sep=""))){
          
          missing_inputs<-c(missing_inputs,deparse(substitute(bam_file)))
          
          
        }
        
      }else{
        
        missing_inputs<-c(missing_inputs,"bam_file")
        
        
      }
      
      
      if(length(table_annotated)>0 ){
        
        if(table_annotated==""){
          
          missing_inputs<-c(missing_inputs,deparse(substitute(table_annotated)))
          
        }
        
      }else{
        
        missing_inputs<-c(missing_inputs,"table_annotated")
        
        
      }
      
      if(length(vcf_file)>0 ){
        
        if(vcf_file=="" | !file.exists(paste(selected_directory,vcf_file,".tbi",sep=""))){
          
          missing_inputs<-c(missing_inputs,deparse(substitute(vcf_file)))
          
        }
        
      }else{
        
        missing_inputs<-c(missing_inputs,"bam_file")
        
        
      }
      
      if(length(gff_file)>0 ){
        
        if(gff_file=="" | !file.exists(paste(selected_directory,gff_file,".tbi",sep=""))){
          
          missing_inputs<-c(missing_inputs,deparse(substitute(gff_file)))
          
        }
        
      }else{
        
        missing_inputs<-c(missing_inputs,"gff_file")
        
        
      }
      
      if(length(reads_counts_cancer_Gtex_file)>0 ){
        
        if(reads_counts_cancer_Gtex_file==""){
          
          missing_inputs<-c(missing_inputs,deparse(substitute(reads_counts_cancer_Gtex_file)))
          
        }
        
      }else{
        
        missing_inputs<-c(missing_inputs,"reads_counts_cancer_Gtex_file")
        
      }
      
      
      # for the patients
      
      if(length(reads_counts_patients_file)>0 ){
        
        if(reads_counts_patients_file==""){
          
          #missing_inputs<-c(missing_inputs,deparse(substitute(reads_counts_patients_file)))
          reads_counts_patients <- data.frame()
          
        }
        
      }else{
        
        #missing_inputs<-c(missing_inputs,"reads_counts_patients_file")
        reads_counts_patients <- data.frame()
        
      }
      
      
      
      if(length(cancer_normal_quantiles_file)>0 ){
        
        if(cancer_normal_quantiles_file==""){
          
          missing_inputs<-c(missing_inputs,deparse(substitute(cancer_normal_quantiles_file)))
          
        }
        
      }else{
        
        missing_inputs<-c(missing_inputs,"cancer_normal_quantiles_file")
        
      }
      
      # for the patients
      
      if(length(patients_quantiles_file)>0 ){
        
        if(patients_quantiles_file==""){
          
          #missing_inputs<-c(missing_inputs,deparse(substitute(patients_quantiles_file)))
          patients_quantiles <- data.frame()
          
        }
        
      }else{
        
        #missing_inputs<-c(missing_inputs,"patients_quantiles_file")
        patients_quantiles <- data.frame()
        
      }
      
      
      
      if(length(genome_file)>0 ){
        
        if(genome_file=="" | !file.exists(paste(selected_directory,genome_file,".fai",sep=""))){
          
          missing_inputs<-c(missing_inputs,deparse(substitute(genome_file)))
          
        }
        
      }else{
        
        missing_inputs<-c(missing_inputs,"genome_file")
        
        
      }
      
      if(length(sample_order_file)>0 ){
        
        if(sample_order_file==""){
          
          missing_inputs<-c(missing_inputs,deparse(substitute(sample_order_file)))
          
        }
        
      }else{
        
        missing_inputs<-c(missing_inputs,"sample_order_file")
        
        
      }
      
      if(length(cancer_cols_file)>0 ){
        
        if(cancer_cols_file==""){
          
          missing_inputs<-c(missing_inputs,deparse(substitute(cancer_cols_file)))
          
        }
        
      }else{
        
        missing_inputs<-c(missing_inputs,"cancer_cols_file")
        
        
      }
      
      
      #if missing_inputs is not empty, don't validate the input directory
      #otherwise, show a validation message, then initialize reactive values
      if(length(missing_inputs)>0){
        
        output$validation_message<-renderText({
          
          paste("-> some files are missing (or their index) :\n\n",paste(missing_inputs,collapse="\n"),sep="")
          
        })
        
      }else{
        
        output$validation_message<-renderText({
          
          "-> all files are validated ! waiting for application to launch..."
          
        })
        
        outputOptions(output, 'validation_message', suspendWhenHidden=FALSE)

        observeEvent(plotsReady(), {
          if (plotsReady()) {
            output$validation_message<-renderText({""})
          }
        })
        
        
        #show progress (initialization)
        progress1<-Progress$new(session,min=1,max=5)
        
        progress1$set(message=paste("step 1/5 : load annotated contigs..."),value=1)
        
        #readRDS replaced by read_fst, it's faster
        #https://github.com/fstpackage/fst
        my_data<-read_fst(paste(selected_directory,table_annotated,sep=""))
        
        progress1$set(message=paste("step 2/5 : load quantiles..."),value=2)
        
        cancer_normal_quantiles<-read_fst(paste(selected_directory,cancer_normal_quantiles_file,sep=""))
        #patients_quantiles<-read_fst(paste(selected_directory,patients_quantiles_file,sep=""))
        if(length(patients_quantiles_file)>0 ){
          patients_quantiles<-read_fst(paste(selected_directory,patients_quantiles_file,sep=""))
        }
        
        progress1$set(message=paste("step 3/5 : load counts..."),value=3)
        
        reads_counts_cancer_Gtex<-read_fst(paste(selected_directory,reads_counts_cancer_Gtex_file,sep=""))
        #reads_counts_patients<-read_fst(paste(selected_directory,reads_counts_patients_file,sep=""))
        if(length(reads_counts_patients_file)>0 ){
          reads_counts_patients<-read_fst(paste(selected_directory,reads_counts_patients_file,sep=""))
        }
        
        
        cancer_cols<-as.character(unlist(read.delim(paste(selected_directory,cancer_cols_file,sep=""),header=F)))
        
        sample_order<-unlist(strsplit(as.character(unlist(read.delim(paste(selected_directory,sample_order_file,sep=""),header=F))),","))
        sample_order<-data.frame(sample_name=sample_order,sample_order=1:length(sample_order))
        
        progress1$set(message=paste("step 4/5 : open port..."),value=4)
        
        
        #port for jbrowser
        port_number<-5000
        #initialize error vector
        error<-1
        
        
        #catch the result of server building
        while(error==1){
          
          #catch message when trying to create the server
          data_server <<- tryCatch(serve_data(selected_directory,port=port_number),error=function(e) e,warning=function(w) w)
          
          
          if(any(class(data_server)%in%c("error","simpleWarning"))){
            
            error<-1
            
            port_number<-port_number+1
            
            cat("port ", port_number,"failed\n")
            
          }else{
            
            
            error<-0
            
            cat("port ", port_number,"succeeded\n")
            
          }
          
        }
        
        progress1$set(message=paste("step 5/5 : build graphics..."),value=5)
        
        vals <- reactiveValues(scatterplot="",refined_table=my_data,boxplots="",peptide_plot="",selected_row_in_table=NULL,triggerscatterplot=FALSE,chromosome_start_end="chr12:6534512-6538374")
        
        
        #give to UI panel for the sidebar
        output$filters_panel<-renderUI({
          if (!is.null(my_data$recurrence_patients)) {
          #if (length(grep("recurrence", grep("patients", names(my_data), value = TRUE), value = TRUE)) > 0) {
            
            recurrence_patients_values <- my_data[, grep("recurrence", grep("patients", names(my_data), value = TRUE), value = TRUE)]
            recurrence_patients_numeric <- as.numeric(recurrence_patients_values[!is.na(recurrence_patients_values)])
            
            if (length(recurrence_patients_numeric) > 0) {
              max_recurrence_patients <- max(recurrence_patients_numeric)
              choices_recurrence_patients <- c("all", seq(1, max_recurrence_patients, 1))
              selected_recurrence_patients <- "all"
            } else {
              choices_recurrence_patients <- "no patients data"
              selected_recurrence_patients <- "no patients data"
            }
          } else {
            choices_recurrence_patients <- "no patients data"
            selected_recurrence_patients <- "no patients data"
          }
          
          #test1
          head(my_data)
          
          fluidRow(column(12,
                          selectInput("length_filter", "Select min contig length",choices = sort(unique(my_data$contig_length)),multiple =F,selected=31),
                          selectInput("gene_name_filter", "Select genes",choices = c("all",unique(my_data$gene_symbol)),selected="all",multiple =T),
                          selectInput("cell_line_filter","expression in one of those :",choices=c("any",sample_order$sample_name),selected="any",multiple=T),
                          selectInput("splice_filter", "spliced or not",choices = c("only_spliced","not_spliced","all"),selected="all",multiple =F),
                          selectInput("multihit_filter", "Select multihit or not",choices = c("single_hit","multihit_only","all"),selected="all",multiple =F),
                          selectInput("recurrence_CCLE_filter", "Min authorized recurrence in cancer cell lines",choices = c("all",seq(1,max(my_data[,grep("recurrence",grep("CCLE",names(my_data),value=T),value=T)]),1)),selected=max(my_data[,grep("recurrence",grep("CCLE",names(my_data),value=T),value=T)]),multiple =F),
                          selectInput("recurrence_patients_filter", "Min authorized recurrence in patients", choices = choices_recurrence_patients, selected = selected_recurrence_patients, multiple = F),
                          selectInput("recurrence_Gtex_filter", "Max authorized recurrence in Gtex",choices = c("all",seq(0,max(my_data[,grep("recurrence",grep("Gtex",names(my_data),value=T),value=T)]),1)),selected=max(my_data[,grep("recurrence",grep("Gtex",names(my_data),value=T),value=T)]),multiple =F),
                          selectInput("perct_filter_CCLE", "Select at least X % perct cancer cell lines samples to be higher than Gtex",choices = names(cancer_normal_quantiles)[grepl("CCLE",names(cancer_normal_quantiles))],selected="all_95perct_CCLE",multiple =F),
                          selectInput("perct_filter_Gtex", "Select at least X % perct Gtex samples to be lower than cancer",choices = names(cancer_normal_quantiles)[grepl("Gtex",names(cancer_normal_quantiles))],selected="all_99perct_Gtex",multiple =F),
                          selectInput("perct_filter_patients", "Select at least X % perct patients samples to be higher than Gtex",choices = names(patients_quantiles)[grepl("patients",names(patients_quantiles))],selected="all_95perct_patients",multiple =F),
                          actionButton("applyCustomFilters", HTML("Click to apply custom filters"),icon("paper-plane"),style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                          bsButton("applyNeoFilters", "Click to add Neo filters",icon("paper-plane"))


          ))
          
        })

          
        #give to UI panel for the body
        output$scatterplot_panel <- renderUI({
          div(
            style = "border: 4px double black; margin-left:3px; margin-right:3px; padding:4px;",
            tagList(
              fluidRow(
                style = "margin-left:3px; margin-right:3px;",
                h3("Visualization parameters:"),
                column(width=12,
                       div(style = "flex: 1 1 50%;",
                           box(
                             title = "Scatterplot",
                             status = "primary",
                             solidHeader = TRUE,
                             collapsible = TRUE,
                             selectInput("X_axis", "Select X axis", width = "100%", choices = c(colnames(my_data)), multiple = FALSE, selected = grep("^median", grep("CCLE", names(my_data), value = TRUE), value = TRUE)),
                             selectInput("Y_axis", "Select Y axis", width = "100%", choices = c(colnames(my_data)), multiple = FALSE, selected = grep("^median", grep("Gtex", names(my_data), value = TRUE), value = TRUE)),
                             selectInput("colored_parameter", "Parameter to color", width = "100%", choices = c(colnames(my_data)), multiple = FALSE, selected = "gene_biotype"),
                             actionButton("showscatterplot", "Show scatterplot", icon("paper-plane"), style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                             plotOutput("scatterplot", height = "394px", width = "100%")
                           )
                       ),
                       div(style = "flex: 1 1 50%;",
                           box(
                             title = "Summary Table",
                             status = "warning",
                             solidHeader = TRUE,
                             collapsible = TRUE,
                             DTOutput("summary_table", height = "651px", width = "100%")
                           )
                       ))
              ),
              fluidRow(
                style = "margin-left:3px; margin-right:3px;",
                column(width=12,
                       div(style = "flex: 1 1 50%;",
                           box(
                             title = "Boxplot",
                             status = "info",
                             solidHeader = TRUE,
                             collapsible = TRUE,
                             plotOutput("boxplots", height = "300px", width = "100%")
                           )),
                       div(style = "flex: 1 1 50%;",
                           box(
                             title = "Peptides",
                             status = "danger",
                             solidHeader = TRUE,
                             collapsible = TRUE,
                             plotOutput("peptide_plot", height = "300px", width = "100%")
                           ))
                       
                       
                )
              )
            )
          )
        })
        
        
        
        #give to UI panel for downlaod
        output$download_plots_panel<-renderUI({
          
          fluidRow(style = "margin-left:3px; margin-right:3px;",
                   column(3,downloadButton("downloadscatterplot", "download scatterplot")),
                   column(3,downloadButton("downloadsummary", "download summary")),
                   column(3,downloadButton("downloadboxplot", "download boxplot")),
                   column(3,downloadButton("downloadpeptides", "download peptides")))
          
        })
        
        #give to UI panel for contig table
        output$contig_table_panel<-renderUI({
          
          fluidRow(style = "border: 4px double black; margin-left:3px; margin-right:3px; padding:4px;",
                   htmlOutput("selection_title"),
                   downloadButton("downloadcontigstable", "download selected contigs"),
                   DTOutput("refined_table"))
          
        })
        
        
        # create the necessary JB2 assembly configuration
        assembly <- assembly(
          paste(data_server$url,"/",basename(genome_file),sep="")
        )
        
        # create configuration for a JB2 GFF FeatureTrack
        annotations_track <- track_feature(
          paste(data_server$url,"/",basename(gff_file),sep="") ,
          assembly
        )
        
        # create configuration for a JB2 GFF alignmenttrack
        # alignments_track <- track_alignments(
        #   paste(data_server$url,"/",basename(bam_file),sep="") ,
        #   assembly
        # )
        
        # create configuration for a JB2 GFF alignmenttrack
        variants_track <- track_variant(
          paste(data_server$url,"/",basename(vcf_file),sep="") ,
          assembly
        )
        

        
        ########################################################################
        
        ################## VERSION AVEC LA RECHERCHE PAR GENE (LENTE) ##################

        gff_json <- paste(
          "{\"dateCreated\":\"2024-04-25T09:33:18.394Z\",",
          "\"tracks\":[",
          "{\"trackId\":\"", gff_file, "\",",
          "\"attributesIndexed\":[\"gene_name\", \"gene_id\", \"transcript_name\"]",
          "\"excludedTypes\":[\"CDS\", \"exon\"],",
          "\"adapterConf\":{",
          "\"type\":\"Gff3TabixAdapter\",",
          "\"gffGzLocation\":{",
          "\"localPath\":\"", data_server$url, "/", gff_file, "\",",
          "\"locationType\":\"LocalPathLocation\"",
          "}",
          "}",
          "}],",
          "\"assemblyNames\":[]}",
          sep = ""
        )



        cat("dir is : ",selected_directory," ; url is : ",data_server$url,"; json is : ",gff_json,"\n")

        cat(gff_json,file=paste(selected_directory,"/gff_json.json",sep=""))

        gff_json<-paste(data_server$url,"/gff_json.json",sep="")
        ix_file<-paste(data_server$url,"/",gff_file,".ix",sep="")
        ixx_file<-paste(data_server$url,"/",gff_file,".ixx",sep="")


        # Messages de debug
        cat("gff_json path: ", gff_json, "\n")
        cat("ix_file path: ", ix_file, "\n")
        cat("ixx_file path: ", ixx_file, "\n")

        output$browserOutput <- renderJBrowseR(

          if(!is.null(vals$selected_row_in_table)){

            row_num= as.numeric(vals$selected_row_in_table)

            sub_table<-vals$refined_table

            contig_name<-as.character(sub_table[row_num,1])

            sub_bam<-paste(selected_directory,"/","sub_bam_",contig_name,".bam",sep="")

            filter <- FilterRules(list(iWantOneContig=function(x) startsWith(x$qname, contig_name)))

            filterBam(paste(selected_directory,"/",bam_file,sep=""), sub_bam, filter=filter)

            alignments_track_selected_contig <- track_alignments(
              paste(data_server$url,"/","sub_bam_",contig_name,".bam",sep="") ,
              assembly
            )

            # create the tracks array to pass to browser
            tracks <- tracks(

              annotations_track,
              variants_track,
              alignments_track_selected_contig
              # alignments_track


            )

            default_session <- default_session(
              assembly,
              c(annotations_track,variants_track,alignments_track_selected_contig)

            )

            JBrowseR::JBrowseR(
              view="View",
              #view="ViewHg38",
              assembly = assembly,
              tracks = tracks,
              #chr12:6534512-6538374
              location = vals$chromosome_start_end,
              defaultSession = default_session,

              #https://github.com/GMOD/JBrowseR/blob/main/example_apps/basic_usage_with_text_index/app.R
              #assembly name should be the one you see when you click on "about" for the annotation track ! and it is inherited from the function assembly()
              text_index = text_index(ix_file,ixx_file,gff_json,unlist(strsplit(basename(genome_file),"\\."))[1])

            )

          }else{
            
            # create the tracks array to pass to browser
            tracks <- tracks(
              
              annotations_track,
              #alignments_track,
              variants_track
              
            )
            
            # set up the default session for the browser
            default_session <- default_session(
              assembly,
              c(annotations_track,variants_track)
              
            )
            

            JBrowseR::JBrowseR(
              view="View",
              #view="ViewHg38",
              assembly = assembly,
              tracks = tracks,
              #vals$chromosome_start_end,
              location = vals$chromosome_start_end,
              defaultSession = default_session,

              #https://github.com/GMOD/JBrowseR/blob/main/example_apps/basic_usage_with_text_index/app.R
              #assembly name should be the one you see when you click on "about" for the annotation track ! and it is inherited from the function assembly()
              text_index = text_index(ix_file,ixx_file,gff_json,unlist(strsplit(basename(genome_file),"\\."))[1])

            )


          }

        )
        
        #give to UI panel panel for browser
        output$browser_panel<-renderUI({
          
          JBrowseROutput("browserOutput")
          
        })
        ########################################################################
        
        ################## ANCIENNE VERSION ##################
        
        # output$browserOutput <- renderJBrowseR(
        #   JBrowseR::JBrowseR(
        #     #view="View",
        #     view="ViewHg38",
        #     #assembly = assembly,
        #     tracks = tracks,
        #     location = vals$chromosome_start_end,
        #     defaultSession = default_session
        #     #theme = my_theme
        #   )
        # )
        # 
        # #give to UI panel panel for browser
        # output$browser_panel<-renderUI({
        #   
        #   JBrowseROutput("browserOutput")
        #   
        # })
        
        
        
        
        
        

        
        
        
        progress1$close()
        
        
        loadDir$value<-1
        
        
        contigs_stats <- function (all){
          stats<-c(contigs=nrow(all),
                   mapped=sum(!is.na(all$mapped_to)),
                   intergenic=sum(!is.na(all$mapped_to) & is.na(all$gene_id)),
                   genic=sum(!is.na(all$gene_id)),
                   exonic=sum(all$is_exonic,na.rm=T),
                   intronic=sum(all$is_intronic,na.rm=T),
                   spliced=sum(all$nb_splice>0,na.rm=T),
                   indel=sum(all$nb_del>0 | all$nb_ins>0,na.rm=T),
                   snv=sum(all$nb_snv>0,na.rm=T),
                   clipped=sum(all$clipped_5p>0 | all$clipped_3p >0, na.rm=T),
                   clippedsup5=sum(all$clipped_5p>5 | all$clipped_3p >5, na.rm=T),
                   circular=sum(all$is_circ,na.rm=T),
                   chimeric=sum(all$is_chimeric,na.rm=T),
                   repeats=sum(!is.na(all$HumanRepeats))
          )
          statf <- data.frame(feature=names(stats),counts=stats)
          statf$percent= round(100*statf[,"counts"]/statf[1,"counts"],2)
          if("meanA" %in% colnames(all)) {
            stats<-c(Nups=sum(all$meanB>all$meanA,na.rm=T),
                     Ndowns=sum(all$meanA>all$meanB,na.rm=T),
                     mean_meanA=mean(all$meanA),
                     mean_meanB=mean(all$meanB),
                     median_meanA=median(all$meanA),
                     median_meanB=median(all$meanB))
            statf2 <- data.frame(feature=names(stats),counts=stats)  
            statf2$percent=NA
            statf=rbind(statf,statf2)
          }
          statf$counts=round(statf$counts)
          return(statf)
        }
        
        load_contigs2 <- function (contigfile) {
          
          all<-as_tibble(contigfile)
          all$length=nchar(all$contig)
          all$neo =
            (all$is_exonic==0 | is.na(all$is_exonic))  | # not exonic
            (all$is_exonic==1 & (all$nb_del>0 | all$nb_ins>0)) | # exonic & indel
            (all$is_circ==1 | all$is_chimeric==1) |
            (all$clipped_5p>5 | all$clipped_3p >5)
          # replace chr names with chr by just number
          all$chromosome=str_replace(all$chromosome, "chr", "")
          # create output object
          CTG<-list()
          CTG$all<-all
          CTG$stats<-contigs_stats(all)
          return(CTG)
        }
        
        
        
        observeEvent(c(input$applyCustomFilters),{
          
          #if original dataframe exists, apply filters
          #if(is.data.frame(refined_table)){
          
          updateButton(session,"applyNeoFilters",label="Click to add Neo filters")
          
          tmp_table1<-my_data[which(my_data$contig_length>=as.numeric(input$length_filter)),]
          
          vals$refined_table<-tmp_table1
          
          
          if(input$splice_filter=="only_spliced"){
            
            tmp_table1<-tmp_table1[which(tmp_table1$nb_splice>=1),]
            
          }else if(input$splice_filter=="not_spliced"){
            
            tmp_table1<-tmp_table1[which(tmp_table1$nb_splice==0),]
            
          }
          vals$refined_table<-tmp_table1
          
          
          if(input$multihit_filter=="multihit_only"){
            
            tmp_table1<-tmp_table1[which(tmp_table1$nb_hit>1),]
            
          }else if(input$multihit_filter=="single_hit"){
            
            tmp_table1<-tmp_table1[which(tmp_table1$nb_hit==1),]
            
          }
          
          vals$refined_table<-tmp_table1
          
          #patients
          if (!is.null(my_data$recurrence_patients)) {
          #if (length(grep("recurrence", grep("patients", names(my_data), value = TRUE), value = TRUE)) > 0) {
          #if ("recurrence_patients" %in% colnames(my_data)) {
            if(input$recurrence_patients_filter!="all"){
              
              tmp_table1<-tmp_table1[which(tmp_table1[,grep("recurrence",grep("patients",names(tmp_table1),value=T),value=T)]>=as.numeric(input$recurrence_patients_filter)),]
              
            }else{
              
              tmp_table1<-tmp_table1[which(tmp_table1[,grep("recurrence",grep("patients",names(tmp_table1),value=T),value=T)]>=max(tmp_table1[,grep("recurrence",grep("patients",names(tmp_table1),value=T),value=T)])),]
            }
            vals$refined_table<-tmp_table1
          }
          
          
          
          # if(input$recurrence_patients_filter!="all"){
          #   
          #   
          #   tmp_table1<-tmp_table1[which(tmp_table1[,grep("recurrence",grep("patients",names(tmp_table1),value=T),value=T)]>=as.numeric(input$recurrence_patients_filter)),]
          #   
          # }else{
          #   
          #   tmp_table1<-tmp_table1[which(tmp_table1[,grep("recurrence",grep("patients",names(tmp_table1),value=T),value=T)]>=max(tmp_table1[,grep("recurrence",grep("patients",names(tmp_table1),value=T),value=T)])),]
          # }
          
          
          
          if(input$recurrence_CCLE_filter!="all"){


            tmp_table1<-tmp_table1[which(tmp_table1[,grep("recurrence",grep("CCLE",names(tmp_table1),value=T),value=T)]>=as.numeric(input$recurrence_CCLE_filter)),]

          }else{

            tmp_table1<-tmp_table1[which(tmp_table1[,grep("recurrence",grep("CCLE",names(tmp_table1),value=T),value=T)]>=max(tmp_table1[,grep("recurrence",grep("CCLE",names(tmp_table1),value=T),value=T)])),]
          }

          vals$refined_table<-tmp_table1
          
          
          
          if(input$recurrence_Gtex_filter!="all"){
            
            tmp_table1<-tmp_table1[which(tmp_table1[,grep("recurrence",grep("Gtex",names(tmp_table1),value=T),value=T)]<=as.numeric(input$recurrence_Gtex_filter)),]
            
          }else{
            
            tmp_table1<-tmp_table1[which(tmp_table1[,grep("recurrence",grep("CCLE",names(tmp_table1),value=T),value=T)]<1e100),]
          }
          
          vals$refined_table<-tmp_table1
          
          #Filter for (X % CCLE > X % Gtex) && (X % patients > X % Gtex)
          selected_quantiles_cancer<-as.character(cancer_normal_quantiles[which(cancer_normal_quantiles[,input$perct_filter_CCLE]>cancer_normal_quantiles[,input$perct_filter_Gtex] ),]$tag)
          
          if (nrow(patients_quantiles) == 0) {
            selected_quantiles <- selected_quantiles_cancer
          } else {
            selected_quantiles_patients<-as.character(patients_quantiles[which(patients_quantiles[,input$perct_filter_patients]>cancer_normal_quantiles[,input$perct_filter_Gtex]),]$tag)
            selected_quantiles <- intersect (selected_quantiles_cancer, selected_quantiles_patients)
          }
          # selected_quantiles_patients<-as.character(patients_quantiles[which(patients_quantiles[,input$perct_filter_patients]>cancer_normal_quantiles[,input$perct_filter_Gtex]),]$tag)
          # selected_quantiles <- intersect (selected_quantiles_cancer, selected_quantiles_patients)
          
          tmp_table1<-tmp_table1[which(tmp_table1$tag%in%selected_quantiles),]
          
          vals$refined_table<-tmp_table1
          
          if(any(as.character(unlist(input$gene_name_filter))%in%c("all"))==F){
            
            
            tmp_table1<-tmp_table1[which(tmp_table1$gene_symbol%in%c(as.character(unlist(input$gene_name_filter)))),]
            
          }
          
          vals$refined_table<-tmp_table1

          
          
          
          #filter on wanted min samples in cancer
          
          if(any(as.character(unlist(input$cell_line_filter))%in%c("any"))==F){
            
            sample_order_selected<-sample_order[which(sample_order$sample_name%in%as.character(unlist(input$cell_line_filter))),]
            
            reads_counts_cancer_Gtex_subset<-reads_counts_cancer_Gtex[,c(1,sample_order_selected$sample_order+1)]
            
            selected_expression_in_samples<-as.character(reads_counts_cancer_Gtex_subset[apply(reads_counts_cancer_Gtex_subset[,2:ncol(reads_counts_cancer_Gtex_subset)], 1, function(row) any(row >0 )),]$tag)
            
            tmp_table1<-tmp_table1[which(tmp_table1$tag%in%selected_expression_in_samples),]
            
            
            
            vals$refined_table<-tmp_table1
            
            
          }
          
          #}
          
        })
        
        
        observeEvent(input$applyNeoFilters,{
          
          #if original dataframe exists, apply filters
          #if(is.data.frame(refined_table)){
          
          updateButton(session,"applyNeoFilters",style = "warning",label="Neo filters are applied !")
          
          tmp_table1<-vals$refined_table
          
          vals$refined_table<-tmp_table1[which((tmp_table1$is_exonic==0 | is.na(tmp_table1$is_exonic))| # not exonic
                                                 tmp_table1$is_exonic==1 & (tmp_table1$nb_del>0 | tmp_table1$nb_ins>0| tmp_table1$nb_snv>0 )|# exonic & (indel|insertion|snv)
                                                 (tmp_table1$is_circ==1 | tmp_table1$is_chimeric==1)|
                                                 (tmp_table1$clipped_5p>5 | tmp_table1$clipped_3p >5)
          ),]
          
          #}
          
        })
        
        
        
        #when you trigger the scatter, put the value to TRUE to plot it
        observeEvent(input$showscatterplot,{
          
          #if(is.data.frame(refined_table)){
          
          vals$triggerscatterplot<-TRUE
          
          #}
          
        })
        
        #each time an other reactive value change, put trigger scatterplot to FALSE in order to have the scatter at null
        observeEvent(c(input$applyCustomFilters,input$X_axis,input$Y_axis,input$colored_parameter),{
          
          vals$triggerscatterplot<-FALSE
          
        })
        
        
        
        #https://community.rstudio.com/t/problem-in-downloading-the-plot-with-downloadhandler-for-shiny-app/44876
        output$scatterplot <- renderPlot({
          
          if(vals$triggerscatterplot==T){
            
            sub_table<-vals$refined_table
            
            # patients_data <- reads_counts_patients[which(reads_counts_patients$feature %in% sub_table$tag),]
            # 
            # combined_table <- merge(sub_table, patients_data, by.x = "tag", by.y = "feature", all = TRUE)
            # 
            cat("classes : ",class(sub_table[,input$Y_axis]),class(sub_table[,input$X_axis]),class(sub_table[,input$colored_parameter]))
            
            
            #if x discrete, take n rows, only
            if(class(sub_table[,input$X_axis])%in%c("numeric","integer")){
              
              x_scale<-scale_x_continuous()
              
            }else{
              
              x_scale<-scale_x_discrete()
              
              
              #if it's greater than the threshold, do a subset
              if(length(unique(sub_table[,input$X_axis]))>10){
                
                
                #if it's chromosome parameter, take all (~20)
                if(grepl("chrom",input$X_axis,ignore.case = T)){
                  
                  sub_x<-head(unique(sub_table[,input$X_axis]),25)
                  
                }else{
                  
                  sub_x<-head(unique(sub_table[,input$X_axis]),10)
                }
                
                
                sub_table<-sub_table[which(sub_table[,input$X_axis]%in%sub_x),]
                
              }
              
            }
            
            #if y discrete, take n rows, only
            if(class(sub_table[,input$Y_axis])%in%c("numeric","integer")){
              
              y_scale<-scale_y_continuous()
              
              
            }else{
              
              y_scale<-scale_y_discrete()
              
              if(length(unique(sub_table[,input$Y_axis]))>10){
                
                
                if(grepl("chrom",input$Y_axis,ignore.case = T)){
                  
                  sub_y<-head(unique(sub_table[,input$Y_axis]),25)
                  
                }else{
                  
                  sub_y<-head(unique(sub_table[,input$Y_axis]),10)
                }
                
                
                
                sub_table<-sub_table[which(sub_table[,input$Y_axis]%in%sub_y),]
                
                
              }
              
            }
            
            
            #if colored parameter discrete, color only n features, the rest is seen as "others"
            if(class(sub_table[,input$colored_parameter])%in%c("numeric","integer")){
              
              color_scale<-scale_colour_gradient(low =
                                                   "yellow", high = "red")
              
              
              
            }else{
              
              sub_table[,input$colored_parameter]<-as.character(sub_table[,input$colored_parameter])
              
              color_scale<-scale_color_discrete()
              
              if(length(unique(sub_table[,input$colored_parameter]))>15){
                
                if(grepl("chrom",input$colored_parameter,ignore.case = T)){
                  
                  sub_colored<-head(unique(sub_table[,input$colored_parameter]),25)
                  
                }else{
                  
                  sub_colored<-head(unique(sub_table[,input$colored_parameter]),15)
                }
                
                
                
                sub_table[,input$colored_parameter][!sub_table[,input$colored_parameter] %in% sub_colored] <- "others"
                
              }
              
            }
            
            if(input$Y_axis%in%c(grep("median|mean",names(sub_table),value=T)) & input$X_axis%in%c(grep("median|mean",names(sub_table),value=T,perl=T))){
              
              diag_line<-geom_abline(intercept=0,size=2,color="black")
              
              all_max<-max(c(max(sub_table[,input$X_axis]),max(sub_table[,input$Y_axis])))
              
              all_min<-min(c(min(sub_table[,input$X_axis]),min(sub_table[,input$Y_axis])))
              
              x_scale<-scale_x_continuous(limits=c(all_min,all_max))
              
              y_scale<-scale_y_continuous(limits=c(all_min,all_max))
              
            }else{
              
              diag_line<-geom_abline(intercept=0,size=2,color="transparent")
              
              
            }
            
            #vals$refined_table<-sub_table
            
            my_title<-paste(nrow(sub_table)," contigs",sep="")
            
            
            white_background<-ggplot2::theme(axis.line = element_line(colour = "black",size=2),
                                             panel.grid.major = element_blank(),
                                             panel.grid.minor = element_blank(),
                                             panel.background = element_blank(),
                                             axis.ticks.length.x = unit(0.5, "cm"),
                                             axis.ticks.length.y = unit(0.5, "cm"),
                                             axis.ticks.x=element_line(size=2,color="black"),
                                             axis.ticks.y=element_line(size=2,color="black"),
                                             legend.key = element_rect(fill = "white"))
            
            
            vals$scatterplot<-ggplot(sub_table,aes_string(x=input$X_axis,y=input$Y_axis,color=input$colored_parameter))+
              ggtitle(my_title)+
              diag_line+
              geom_point(size=2,pch=19)+
              color_scale+
              #scale_x_continuous(expand=c(0,0),limits=c(min(sub_table[input$X_axis])-1,max(sub_table[input$X_axis])+1),breaks=seq(min(sub_table[input$X_axis]),max(sub_table[input$X_axis]),round(max(sub_table[input$X_axis])/8)))+
              x_scale+
              #scale_y_continuous(expand=c(0,0),limits=c(min(sub_table[input$Y_axis])-1,max(sub_table[input$Y_axis])+1),breaks=seq(min(sub_table[input$Y_axis]),max(sub_table[input$Y_axis]),round(max(sub_table[input$Y_axis])/8)))+
              y_scale+
              white_background+
              ggplot2::theme(plot.margin = margin(t=0,r=0,b=0,l=0,unit="cm"),axis.text.x= element_text(color="black",size=12,face="bold",hjust =1,vjust=0.5,angle=90),axis.text.y= element_text(color="black",size=12,face="bold",hjust = 0.5),axis.title.x =element_text(color="black",size=14,face="bold"),axis.title.y =element_text(color="black",size=14,face="bold"),legend.text=element_text(size=12),legend.title=element_text(size=14,face="bold"),plot.title=element_text(color="black",size=16,face="bold",hjust = 0.5),aspect.ratio=1)
            
            print(vals$scatterplot)
            
          }else{
            
            ggplot() +
              ggplot2::theme_void() +
              ggplot2::geom_text(aes(0,0,label='Click on "Show the scatterplot"'), size = 4) +
              xlab(NULL)
            
            
          }
          
          
        })
        
        output$summary_table <- renderDT({
          
          
          contigs_stats<-load_contigs2(vals$refined_table)$stats
          
          
          return(contigs_stats)
          
          #https://stackoverflow.com/questions/55152819/r-shiny-datatable-pagination-and-show-all-rows-as-options
        },options = list(scrollX = T,autoWidth = TRUE,searchHighlight = TRUE,columnDefs = list(
          list(width = 'auto', targets ="_all"),
          list(className = "dt-center", targets = "_all")
        ),dom ="ft",
        pageLength = 10000),selection="single",rownames=FALSE,filter = "top")
        
        
        
        #download the scatterplot
        output$downloadscatterplot <- downloadHandler(
          filename = function() {
            
            #give the used parameters as file name
            paste(input$X_axis,"_vs_",input$Y_axis,"_vs_",input$colored_parameter,".png",sep="")
          },
          content = function(file) {
            
            png(file,width = 1000, height = 800)
            
            print(vals$scatterplot)
            
            dev.off()
          }
        )
        
        #download the boxplot
        output$downloadboxplot <- downloadHandler(
          filename = function() {
            
            contig_name<-as.character(vals$refined_table[as.numeric(vals$selected_row_in_table),]$tag)
            
            
            #give the used parameters as file name
            paste("boxplot_",contig_name,".png",sep="")
          },
          content = function(file) {
            
            png(file,width = 1000, height = 800)
            
            print(vals$boxplots)
            
            dev.off()
          }
        )
        
        
        #download the peptides
        output$downloadpeptides <- downloadHandler(
          filename = function() {
            contig_name<-as.character(vals$refined_table[as.numeric(vals$selected_row_in_table),]$tag)
            
            
            #give the used parameters as file name
            paste("peptides_",contig_name,".pdf",sep="")
          },
          content = function(file) {
            
            pdf(file)
            print(vals$peptide_plot)
            
            dev.off()
          }
        )
        
        #download the contigs table
        output$downloadcontigstable <- downloadHandler(
          filename = function() {
            
            table_name<-paste("table_of_",nrow(vals$refined_table),"_selected_contigs", ".tsv", sep = "")
            
            #give the used parameters as file name
            paste(table_name,sep="")
          },
          #filename = paste("table_of_",nrow(vals$refined_table),"_selected_contigs", ".tsv", sep = ""),
          content = function(file) {
            
            write.table(my_data[which(my_data$tag%in%vals$refined_table$tag),],file,row.names=F,col.names=T,quote=F,sep="\t")
          }
        )
        
        
        
        output$selection_title<-renderText({
          
          paste("<h3>Select a row to visualize the corresponding result :<h3>",sep="")
          
        })
        
        #show the table
        # output$refined_table <-renderDT({

          # datatable(
          #   vals$refined_table[,c("tag","contig","contig_length","chromosome","start","end","strand","nb_hit","gene_symbol","gene_biotype","is_chimeric", "is_circ","is_exonic","is_intronic","HumanRepeats","nb_snv","clipped_5p","clipped_3p",grep("^median",grep("CCLE",names(vals$refined_table),value=T),value=T,perl=T),grep("^median",grep("Gtex",names(vals$refined_table),value=T),value=T,perl=T),grep("^median",grep("patients",names(vals$refined_table),value=T),value=T,perl=T),grep("median",grep("log2FC",names(vals$refined_table),value=T),value=T,perl=T),grep("recurrence",grep("CCLE",names(vals$refined_table),value=T),value=T,perl=T),grep("recurrence",grep("Gtex",names(vals$refined_table),value=T),value=T,perl=T),grep("recurrence",grep("patients",names(vals$refined_table),value=T),value=T,perl=T))],

          # table <- summary_info(state)
        #     options = list(
        #       scrollX = T,
        #       autoWidth = TRUE,
        #       searchHighlight = TRUE,
        #       columnDefs = list(
        #         list(width = '300px', targets = 2),
        #         list(targets = 2, createdCell = JS(
        #           "function(cell) {",
        #           "  $(cell).css('max-width', '300px');",
        #           "  $(cell).css('overflow', 'hidden');",
        #           "  $(cell).css('text-overflow', 'ellipsis');",
        #           "  $(cell).css('white-space', 'nowrap');",
        #           "}"
        #         ))
        #       ),
        #       paging = TRUE,
        #       pageLength = 5
        #     ),
        #     selection="single",rownames=TRUE,filter = "top")
        # })
      
        
        output$refined_table <- renderDT({

          fixed_columns <- c("tag", "contig", "contig_length", "chromosome", "start", "end", "strand",
                             "nb_hit", "gene_symbol", "gene_biotype", "is_chimeric", "is_circ",
                             "is_exonic", "is_intronic", "HumanRepeats", "nb_snv", "clipped_5p",
                             "clipped_3p")

          # optional columns
          dynamic_columns <- c(
            grep("^median", grep("CCLE", names(vals$refined_table), value = TRUE), value = TRUE, perl = TRUE),
            grep("^median", grep("Gtex", names(vals$refined_table), value = TRUE), value = TRUE, perl = TRUE),
            grep("^median", grep("patients", names(vals$refined_table), value = TRUE), value = TRUE, perl = TRUE),
            grep("median", grep("log2FC", names(vals$refined_table), value = TRUE), value = TRUE, perl = TRUE),
            grep("recurrence", grep("CCLE", names(vals$refined_table), value = TRUE), value = TRUE, perl = TRUE),
            grep("recurrence", grep("Gtex", names(vals$refined_table), value = TRUE), value = TRUE, perl = TRUE),
            grep("recurrence", grep("patients", names(vals$refined_table), value = TRUE), value = TRUE, perl = TRUE)
          )

          # only include the optional columns that exist
          dynamic_columns <- dynamic_columns[dynamic_columns %in% names(vals$refined_table)]

          selected_columns <- c(fixed_columns, dynamic_columns)

          datatable(
            vals$refined_table[, selected_columns, drop = FALSE],
            options = list(
              scrollX = TRUE,
              autoWidth = TRUE,
              searchHighlight = TRUE,
              columnDefs = list(
                list(width = '300px', targets = 2),
                list(targets = 2, createdCell = JS(
                  "function(cell) {",
                  "  $(cell).css('max-width', '300px');",
                  "  $(cell).css('overflow', 'hidden');",
                  "  $(cell).css('text-overflow', 'ellipsis');",
                  "  $(cell).css('white-space', 'nowrap');",
                  "}"
                ))
              ),
              paging = TRUE,
              pageLength = 5
            ),
            selection = "single",
            rownames = TRUE,
            filter = "top"
          )
        })
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        #########################
        observeEvent(input$refined_table_rows_selected,{
          
          vals$selected_row_in_table<-input$refined_table_rows_selected
          
          vals$chromosome_start_end<-paste(as.character(vals$refined_table[input$refined_table_rows_selected,]$chromosome),":",as.integer(vals$refined_table[input$refined_table_rows_selected,]$start)-50,"-",as.integer(vals$refined_table[input$refined_table_rows_selected,]$end)+50,sep="")
        
          
        },ignoreNULL=F)
        #######################
        
        
        #PEPTIDES
        
        #plot 3-frames peptides when you click on one line, otherwise, print empty plot
        
        output$peptide_plot <- renderPlot(
          
          if (!is.null(vals$selected_row_in_table)) {
            selected_contig_sequence <- vals$refined_table[vals$selected_row_in_table, "contig"]
            
            a <- get3framesPeptidesFromNuc(selected_contig_sequence)
            
            vals$peptide_plot<-a
            
            print(a)
          }else{
              a <- ggplot() +
                ggplot2::theme_void() +
                ggplot2::geom_text(aes(0,0,label='Click on a line in the table below to show\nthe peptides corresponding to the contig'), size = 4) +
                xlab(NULL)
              
              vals$peptide_plot<-a
              
              print(a)
            })
 
        
        
        #BOXPLOT
        
        #plot boxplot when you click on one line, otherwise, print empty plot
        output$boxplots <- renderPlot(
          
         
          if(!is.null(vals$selected_row_in_table)){
            
            row_num= as.numeric(vals$selected_row_in_table)
            
            sub_table<-vals$refined_table
            
            one_kmer<-as.character(sub_table[row_num,1])
            
            
            reads_counts_c_G<-reads_counts_cancer_Gtex[which(reads_counts_cancer_Gtex$tag==one_kmer),]
            reads_counts_p<-reads_counts_patients[which(reads_counts_patients$feature==one_kmer),]
            
            
            reads_counts_c_G <- reshape2::melt(reads_counts_c_G,measure.vars=c(names(reads_counts_c_G)[2:ncol(reads_counts_c_G)]))
            names(reads_counts_c_G)[3]<-"counts"
            reads_counts_c_G$counts<-log10(reads_counts_c_G$counts+1)
            reads_counts_c_G$condition<-gsub("_[0-9]+$","",reads_counts_c_G$variable)
            
            if (nrow(reads_counts_p) > 0) {
              reads_counts_p <- reshape2::melt(reads_counts_p,measure.vars=c(names(reads_counts_p)[2:ncol(reads_counts_p)]))
              names(reads_counts_p)[3]<-"counts"
              reads_counts_p$counts<-log10(reads_counts_p$counts+1)
              reads_counts_p$condition<- "Patients"
              
              names(reads_counts_p) <- names(reads_counts_c_G)
              combined_counts <- rbind(reads_counts_c_G, reads_counts_p)
              
            } else {
              combined_counts <- reads_counts_c_G
            }
            
            # reads_counts_p <- reshape2::melt(reads_counts_p,measure.vars=c(names(reads_counts_p)[2:ncol(reads_counts_p)]))
            # names(reads_counts_p)[3]<-"counts"
            # reads_counts_p$counts<-log10(reads_counts_p$counts+1)
            # reads_counts_p$condition<- "Patients"


            # names(reads_counts_p) <- names(reads_counts_c_G)
            # combined_counts <- rbind(reads_counts_c_G, reads_counts_p)
            
            combined_counts$condition <- factor(combined_counts$condition, levels = unique(combined_counts$condition))
            
            
            white_background<-ggplot2::theme(axis.line = element_line(colour = "black",size=1),
                                             panel.grid.major = element_blank(),
                                             panel.grid.minor = element_blank(),
                                             panel.background = element_blank(),
                                             panel.border = element_blank())
            
            
            p<-ggplot(combined_counts, aes(x = condition, y =counts)) +
              geom_boxplot(color="black",outlier.color = "transparent",outlier.size = 1,lwd=1)+
              
              #xlab(one_gene)+
              ylab(paste("counts (log10 CPM)",sep=""))+
              scale_y_continuous(expand = c(0,0))+
              #ggtitle(global_title)+
              ggplot2::theme_bw()+
              
              facet_wrap(tag~.,scales = c("free_y"), dir = "v",ncol =2)+
              ggplot2::theme(axis.text.x= element_text(color="black",size=14,face="bold",angle = 90,hjust = 1,vjust=0.5),axis.text.y = element_text(color="black",size=14,face="bold"),plot.title = element_text(hjust = 0.5,size=10),axis.title=element_text(size=18,face="bold"),legend.text=element_text(size=18),legend.title=element_text(size=18,face="bold"),
                             axis.ticks.length.x = unit(0.25, "cm"),axis.ticks.length.y = unit(0.25, "cm"),axis.ticks.x=element_line(size=2,color="black"),axis.ticks.y=element_line(size=2,color="black"),legend.position="none",
                             strip.text.x = element_text(size = 10))+
              white_background
            #facet_grid(feature~.,scales = c("free"))
            
            
            dat <- ggplot_build(p)$data[[1]]
            
            meds <- ddply(combined_counts, .(tag,condition), summarize, med = median(counts))
            
            dat<-cbind(dat,meds)
            
            p<-p +geom_point(data =combined_counts,aes(colour=condition),size=1.5,pch = 19, position = position_jitterdodge(5))
            
            p<-p+geom_segment(data=dat, aes(x=xmin, xend=xmax,y=med, yend=med), colour="black", size=1)+ggplot2::theme(strip.text = element_text(size=4),strip.background = element_rect(colour="black", fill="grey80",size=1))
            
            vals$boxplots<-p
            
            print(p)
            
            plotsReady(TRUE)
            
          }else{
            
            p<-ggplot() +
              ggplot2::theme_void() +
              geom_text(aes(0,0,label='Click on a line in the table below to show\nthe contig expression in CCLE vs GTEX vs Patients'), size = 4) +
              xlab(NULL)
            
            vals$boxplots<-p
            
            print(p)
            
            plotsReady(TRUE)
          }
        )
        
      }
      
      
      
    }
    
  })
  
  
  
  
  ########### browser ##############
  
  
  
  # {
  #   
  #   
  #   # create the necessary JB2 assembly configuration
  #   assembly <- assembly(
  # 	paste(data_server$url,"/",basename(genome),sep="")
  #   )
  #   
  #   # create configuration for a JB2 GFF FeatureTrack
  #   annotations_track <- track_feature(
  # 	paste(data_server$url,"/",basename(ref_gff),sep="") ,
  # 	assembly
  #   )
  #   
  #   # create configuration for a JB2 GFF alignmenttrack
  #   alignments_track <- track_alignments(
  # 	paste(data_server_alignment$url,"/",basename(alignment),sep="") ,
  # 	assembly
  #   )
  #   
  #   # create configuration for a JB2 GFF alignmenttrack
  #   variants_track <- track_variant(
  # 	paste(data_server_variant$url,"/",basename(variant),sep="") ,
  # 	assembly
  #   )
  #   
  #   # create the tracks array to pass to browser
  #   tracks <- tracks(
  #	 
  # 	annotations_track,
  # 	alignments_track,
  # 	variants_track
  #	 
  #   )
  #   
  #   # set up the default session for the browser
  #   default_session <- default_session(
  # 	assembly,
  # 	c(annotations_track,variants_track,alignments_track)
  #	 
  #   )
  #   
  #   my_theme<-json_config(my_config_theme)
  #   
  #   #link the UI with the browser widget
  #   output$browserOutput <- renderJBrowseR(
  # 	JBrowseR::JBrowseR(
  #   	view="View",
  #   	#view="ViewHg38",
  #   	assembly = assembly,
  #   	tracks = tracks,
  #   	location = vals$chromosome_start_end,
  #   	defaultSession = default_session,
  #   	theme = my_theme
  # 	)
  #   )
  #   
  #   # options <- parseAndValidateGenomeSpec(genomeName="hg38",  initialLocus="GAPDH")
  #   # output$browserOutput <- renderIgvShiny({
  #   #   igvShiny(options)
  #   # })
  #   
  # }
  #
  ###############################
  
  
  
  

  
  
  
  
}

#shinyApp(ui, server, options = list(launch.browser = F) )
shinyApp(ui, server, options = list(shiny.maxRequestSize=Inf))
