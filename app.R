# Scott Schumacker
# 3 Prime Analytics Application

# Load Packages
library(shinydashboard)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(sleuth)
library(biomaRt)
library(shinybusy)
library(plotly)
library(dplyr)
library(fastqcr)

ui <- dashboardPage(skin = "purple",
                    
  Header <- dashboardHeader(title = "3 Prime Analytics",
                            dropdownMenu(type = "messages",
                                         messageItem(
                                           from = "3 Prime",
                                           message = "Welcome to 3 Prime Analytics!"
                                         )
                            )),
  
  #################################################################################
  Sidebar <- dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "Home", icon = icon("home")),
      menuItem("Test", 
               menuSubItem("TestSub", tabName = "TestSub", icon = icon("fire"))),
      menuItem("FastQC", tabName = "FastQC", icon = icon("fire")),
      menuItem("Alignment Data Input", tabName = "RawDataInput", icon = icon("fire")),
      menuItem("QC", tabName = "QC", icon = icon("fire")),
      menuItem("Sleuth Object", tabName = "Sleuth", icon = icon("fire")),
      menuItem("PCA", tabName = "PCA", icon = icon("fire")),
      menuItem("Statistical Analysis", tabName = "Computation", icon = icon("fire")),
      menuItem("Volcano Plot", tabName = "Volcano", icon = icon("fire")),
      menuItem("Transcript Abundance", tabName = "Transcript", icon = icon("fire"))
    )
  ),
  
  #################################################################################
  Body <- dashboardBody(
    tabItems(
      
      # Home Tab
      tabItem(tabName = "Home",
              
              tags$h1("Welcome to 3 Prime Analytics!"),
              tags$br(),
              tags$h3("Page under construction. Coming soon...")
              
      ),
      
      # FastQC Tab
      tabItem(tabName = "FastQC",
              add_busy_spinner(spin = "fading-circle"),
              textInput("FastQCInput", "FastQ File Path Input", value = "", 
                        width = NULL, 
                        placeholder = "Please paste file path here"),
              tags$br(),
              textInput("ResultsInput", "Results Desired File Path Input", value = "", 
                        width = NULL, 
                        placeholder = "Please paste file path here"),
              tags$br(),
              fluidRow(
                actionButton("FastQCButton", "Generate FastQC Report")
              ),
              tags$br(),
              tags$h3("Page under construction. Coming soon...")
              
      ),
      
      # Raw Data Input Tab
      tabItem(tabName = "RawDataInput",
              
              textInput("DataInput", "Kallisto Data Input", value = "", 
                        width = NULL, 
                        placeholder = "Please paste file path here"),
              
              textInput("ConditionInput", "Condition Table.txt", value = "", 
                        width = NULL, 
                        placeholder = "Please paste file path here")
      ),
      
      # QC Tab
      tabItem(tabName = "QC",
              
              add_busy_spinner(spin = "fading-circle"),
              
              fluidRow(
              actionButton("s2cButton", "Show Samples"),
              align="center"
              ),
              dataTableOutput("Table"),
              
              fluidRow(actionButton("t2gButton", 
                                    "Generate Human Transcript Data from biomaRt"), 
                       align="center"),
              dataTableOutput("Table2")
      ),
      
      # Sleuth Object Tab
      tabItem(tabName = "Sleuth",
              
              add_busy_spinner(spin = "fading-circle"),
              actionButton("goButton3", "Generate Sleuth Object"),
              dataTableOutput("Maps")
      ),
      
      # PCA Tab
      tabItem(tabName = "PCA",
              
              add_busy_spinner(spin = "fading-circle"),
              actionButton("PCAButton", "Generate PCA Plot"),
              plotlyOutput("plot")
      ),
      
      # Statistical Analysis Tab
      tabItem(tabName = "Computation",
              
              add_busy_spinner(spin = "fading-circle"),
              actionButton("CompButton", "Generate models and run Wald's Test"),
              dataTableOutput("CompTable")
      ),
      
      # Volcano Plot Tab
      tabItem(tabName = "Volcano",
              
              add_busy_spinner(spin = "fading-circle"),
              actionButton("VolcanoButton", "Generate Volcano Plot"),
              plotOutput("Volc")
      ),
      
      # Transcript Abundance Plot Tab
      tabItem(tabName = "Transcript",
              
              add_busy_spinner(spin = "fading-circle"),
              textInput("TranscriptInput", "Transcript Input", value = "", 
                        width = NULL, 
                        placeholder = "Enter transcript number"),
              actionButton("TranscriptButton", "Generate Transcript Abundance Plot"),
              plotOutput("Abundance")
      )
    )
  )
)

#################################################################################
server <- function(input, output) {
    
  # Generate FastQC Report
  #observeEvent(input$FastQCButton, {
    
    #id <- showNotification(paste("Loading (this may take a few minutes)..."), duration = NULL)
    #fastqc(fq.dir = "~/Documents/FASTQ", # FASTQ files directory
           #qc.dir = "~/Documents/FASTQC", # Results direcory
           #threads = 4                    # Number of threads
    #)
    #qc.dir <- system.file("fastqc_results", package = "fastqcr")
    #qc.dir
    #qc_report(qc.dir, result.file = "~/Desktop/multi-qc-result",
              #experiment = "Exome sequencing of colon cancer cell lines")
    #qc.file <- system.file("fastqc_results", "S1_fastqc.zip",  package = "fastqcr")
    #qc.file
    #qc_report(qc.dir, result.file = "~/Desktop/multi-qc-result",
              #experiment = "Exome sequencing of colon cancer cell lines")
    #qc <- qc_read(qc.file)
    #qc_plot(qc, "Per sequence GC content")
    #qc_plot(qc, "Per base sequence quality")
    #qc_plot(qc, "Per sequence quality scores")
    #qc_plot(qc, "Per base sequence content")
    #qc_plot(qc, "Sequence duplication levels")
    
    #})
    
     # Data Acquisition
      sample_id <- dir(file.path("/Users/ScottSchumacker/Desktop/Research/SciDataPaper/MacRetData"))
      
      kal_dirs <- file.path("/Users/ScottSchumacker/Desktop/Research/SciDataPaper/MacRetData",sample_id)
      
      s2c <- read.table(file.path("/Users/ScottSchumacker/Desktop/Research/SciDataPaper/MacRetTable.txt"),
                        header = TRUE,
                        stringsAsFactors = FALSE,
                        sep = "\t")
      s2c <- dplyr::mutate(s2c, path = kal_dirs)
    
      # Generate transcript table from BiomaRt
      observeEvent(input$t2gButton, {
        
        id <- showNotification(paste("Loading (this may take a few minutes)..."), duration = NULL)
        ensembl = useMart("ENSEMBL_MART_ENSEMBL")
        mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                                 dataset = "hsapiens_gene_ensembl",
                                 host = "www.ensembl.org")
        t2G <<- biomaRt::getBM(
          attributes = c("ensembl_transcript_id", "transcript_version",
                         "ensembl_gene_id", "external_gene_name", "description",
                         "transcript_biotype"), mart = mart)
        t2G <<- dplyr::mutate(t2G, target_id = paste(ensembl_transcript_id, transcript_version, sep = "."))
        t2G <<- dplyr::select(t2G, target_id, ens_gene = ensembl_gene_id, Gene_Name = external_gene_name, description, transcript_biotype)
        output$Table2 <- renderDataTable({t2G})
        id <- showNotification(paste("Completed"), duration = NULL)
        
      })

      # Generate Sleuth Object
      observeEvent(input$goButton3, {
        
        id <- showNotification(paste("Loading (this may take 10-20 minutes)..."), duration = NULL)
        so <<- sleuth_prep(s2c,
                          target_mapping = t2G,
                          extra_bootstrap_summary = TRUE,
                          read_bootstrap_tpm = TRUE)
        Mapping <- so$target_mapping
        View(Mapping)
        output$Maps <- renderDataTable({Mapping})
        id <- showNotification(paste("Model Creation Complete"), duration = NULL)
        
      })
  
      # Generate PCA
      observeEvent(input$PCAButton,{
        
      boldTextPlot <- element_text(face = "bold", color = "black", size = 19)
      boldTextAxis <- element_text(face = "bold", color = "black", size = 19)
      PCA <- plot_pca(so, color_by = 'tissue', point_size = 7)
      PCA <- ggplotly(PCA)
      output$plot <- renderPlotly({PCA})
      id <- showNotification(paste("PCA Completed"), duration = NULL)
      
      })

      # Generate Sample Table
      observeEvent(input$s2cButton, {
        
        output$Table <- renderDataTable({s2c})
        id <- showNotification(paste("Completed"), duration = NULL)
        
      })
  
      # Generate Models and Run Statistical Tests
      observeEvent(input$CompButton, {
        
        id <- showNotification(paste("Loading (this may take 5-10 minutes)..."), duration = NULL)
        so <- sleuth_fit(so, ~tissue, 'full')
        so <- sleuth_fit(so, ~1, 'reduced')
        so <- sleuth_wt(so, which_beta = 'tissueretina', which_model = 'full')
        retMacSheet <<- sleuth_results(so, 'tissueretina', 'wt')
        output$CompTable <- renderDataTable({retMacSheet})
        id <- showNotification(paste("Model Creation Complete"), duration = NULL)
        
      })
  
      # Generate Volcano Plot
      observeEvent(input$VolcanoButton,{
        
        inputData <- retMacSheet
    
        inputData <- mutate(inputData, sig=ifelse(inputData$qval<0.05, "FDR < 0.05", "Not Significant"))
    
    
        customPlot <<- ggplot(inputData, aes(b, -log10(qval))) + geom_point(aes(x=b, y=-log10(qval), 
                                                                           color = ifelse(-log10(qval)>-log10(0.05), "Statistically Significant", "no")), size=1.5)
    
        customPlot <<- customPlot + geom_point(aes(x=b, y = -log10(qval), color = ifelse(b < -1.5 & -log10(qval) > -log10(0.05) | b > 1.5 & -log10(qval) > -log10(0.05), "Biologically and Statistically Significant" , NA)))
        customPlot <<- customPlot + geom_point(aes(x=b, y = -log10(qval), color = ifelse(-log10(qval) < -log10(0.05) & (b) < -1.5 | (-log10(qval) < -log10(0.05) & (b) > 1.5), "Biologically Significant", NA)))
        customPlot <<- customPlot + scale_color_manual(values = c("#44AA99", "orange", "black","#332288"), label = c("no" = "Not Signficant"))
    
        customPlot <<- customPlot + geom_vline(xintercept=c(-1.5, 1.5), linetype="dashed")
        customPlot <<- customPlot + geom_hline(yintercept=c(1.30), linetype="dashed")
        output$Volc <<- renderPlot({customPlot})
    
        id <- showNotification(paste("Volcano Plot Completed"), duration = NULL)
        
      })
  
      # Generate Transcript Abundance Plot
      observeEvent(input$TranscriptButton,{
       
        AbundancePlot <- plot_bootstrap(so, input$TranscriptInput, color_by = 'tissue')
        output$Abundance <- renderPlot({AbundancePlot})
    
        id <- showNotification(paste("Abundance Plot Completed"), duration = NULL)
      })
  
}

shinyApp(ui, server)