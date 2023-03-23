library(singscore)
library(Biobase)
library(openxlsx)
library(ggbeeswarm)
library(ggpubr)
library(gridExtra)
library(shiny)
library(shinythemes)

#######################################

#Must change the eset within the folder since the eset file size is too large for RShiny

#######################################
create_IQR_filter_RNASeq_eset = function(eset, samplenames){
  
  expM = exprs(eset)
  # Select required samples
  expM = expM[,samplenames]
  # Apply IQR filtration
  IQR_affy = as.data.frame(apply(expM,1,function(x) IQR(as.vector(x))))
  colnames(IQR_affy) = "IQR_value"
  
  # Read fData
  fdat = fData(eset)
  colnames(fdat)[which(colnames(fdat) == "ensembl")] = "probe"
  df = cbind(IQR_affy,fdat)
  
  dup_geneSym = df[!is.na(df$geneSymbol),]
  dup_geneSym = dup_geneSym[dup_geneSym$geneSymbol != "",]
  dup_geneSym = dup_geneSym[!duplicated(dup_geneSym$geneSymbol),]
  
  # Select the genes with IQR > 0
  IQR_zero = dup_geneSym[which(dup_geneSym$IQR_value > 0),] 
  IQR_zero = IQR_zero[order(IQR_zero$probe),]
  rownames(IQR_zero) = IQR_zero$probe
  
  # Use this as your fData now
  fdata = IQR_zero
  rm(IQR_zero,fdat)
  
  pdat = pData(eset) 
  pdat_sub = pdat[rownames(pdat) %in% colnames(expM),]
  pdat_sub = pdat_sub[order(rownames(pdat_sub)),]
  
  # Use this as your pData now
  pdata = pdat_sub
  rm(pdat_sub,pdat)
  
  expM = data.frame(expM) # Convert matrix to data frame
  expM$probe = rownames(expM)
  expM_new = expM[expM$probe %in% fdata$probe,]
  expM_new = expM_new[order(expM_new$probe),]
  expM_new$probe = NULL
  expM_new = expM_new[,order(colnames(expM_new))]
  
  # Use this as your assay data
  expM = expM_new
  rm(expM_new)
  
  # Create new eset with IQR filter applied
  eset_new = ExpressionSet(assayData = as.matrix(expM),
                           phenoData = AnnotatedDataFrame(data=pdata),
                           featureData = AnnotatedDataFrame(data = fdata))
  return(eset_new)
}

total_rank=function(ranked,num){
  genes=c()
  for(k in 1:ncol(ranked)){
    r1<-ranked[,k]
    r2<-r1[order(r1,decreasing=TRUE)]
    r3<-names(r2)
    r4<-r3[1:num]
    genes<-cbind(genes,r4)
  }
  return(genes)
}

gene_ranks <- function(ranked,Genes){
  rank <- c()
  names <- c()
  for(i in 1:length(Genes)){
    for(j in 1:nrow(ranked)){
      if (rownames(ranked)[j]==Genes[i]){
        index <- which(rownames(ranked)==Genes[i])
        FTL <- ranked[index,]
        m <- mean(FTL)
        rank <- append(rank,m)
        names <- append(names,Genes[i])
      }
    }
  }
  ranks <- cbind(names,rank)
  return(ranks)
  }
# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  theme=shinytheme("slate"),
  # App title ----
  titlePanel("Singscore Script"),
  
  p("Given an eset this script returns a list of the top ranked genes for each patient
     of an amount that the user specifies. It also generates the average rank of a specified list of genes that the user provides. The list of genes will be of every 
     all the genes in the spreadsheet, so if the user wants to specify a gene module they should just input those genes."),
  
  downloadButton("download_button", label = h5("Download the total ranked genes"), "Download Data"),
  
  downloadButton("download_button2", label = h5("Download the gene ranks"), "Download Data"),
  
  fileInput("eset", label = h3("Eset File Input"), accept = ".RData"),
  
  fileInput("file", label = h3("Gene Set File Input"), accept = ".xlsx"),
  
  
  fluidRow(
      column(4,
             
             # Copy the line below to make a slider bar 
             sliderInput("slider1", label = h3("Pick the top _ genes that you want to see"), min = 0, 
                         max = 1000, value = 50)
      )
    ),
    
    hr(),
    
    fluidRow(
      column(4, verbatimTextOutput("value"))
    ),
  
  mainPanel(
    tableOutput("table")
  )
  
)
  
server <- function(input, output) {
  
  options(shiny.maxRequestSize=100*1024^2)
    
  reactiveData = reactiveValues()
  
  tab1 <- observeEvent(input$eset,{
    if(is.null(input$eset)) return(NULL)
    
    inFile = isolate({input$eset})
    
    eset_file = inFile$datapath
    eset1 = get(load(eset_file, envir = .GlobalEnv))
    
    s = sampleNames(eset1)
    expS = create_IQR_filter_RNASeq_eset(eset = eset1, samplenames = s)
    
    pdat = pData(expS)
    fdat = fData(expS)
    exp = data.frame(exprs(expS))
    
    rownames(exp) <- fdat$geneSymbol
    
    ranked = rankGenes(exp)
    
    ranked_small = ranked[1:50,]
    
    n <- input$slider1
    
    ranky <- total_rank(ranked,n)
    
    colnames(ranky) <- colnames(exp)
    
    reactiveData$pd = pdat
    reactiveData$fd = fdat
    reactiveData$ex = exp
    reactiveData$ra = ranky
    reactiveData$r  = ranked_small
  })
  
    excel <- reactive({
      file1 <- input$file
      if(is.null(file1)){
      return()
      }
      Celgene_GSVA <- read.xlsx(file1$datapath)
      Celgene_GSVA_stacked <- stack(Celgene_GSVA)
      if(ncol(Celgene_GSVA)==1){
        return(Celgene_GSVA)
      }
      index<-which(is.na(Celgene_GSVA_stacked))
      Celgene_GSVA_stacked <- Celgene_GSVA_stacked[-index,]
      return(unique(Celgene_GSVA_stacked[,1]))
    })
  
  ran <- reactive({reactiveData$ra}) 
  
  p <- reactive({
    if(is.null(input$file) | is.null(input$eset)){
      return()
    }
    return(gene_ranks(ran(),excel()))
  })
    
  output$table <- renderTable({
    if(is.null(excel())){
      return()
    }
    p()
  }) 
  
  output$value <- renderPrint({
    input$slider1
  })
  
  output$download_button <- downloadHandler(
    filename = "top_genes.csv",
    content = function(file) {
      write.csv(reactiveData$ra, file)
    }
  )
  
  output$download_button2 <- downloadHandler(
    filename = "gene_ranks.csv",
    content = function(file) {
      write.csv(p(),file)
    }
  )
}

shinyApp(ui = ui, server = server)
