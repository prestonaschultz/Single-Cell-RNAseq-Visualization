source("package_install.R")

library(shiny)
library(ggplot2)
library(Seurat)
library(dplyr)
library(labeling)


options(shiny.maxRequestSize=10000*1024^2)
ui=fluidPage(
  
  titlePanel("Single Cell RNA-Sequencing Visualization Tools"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Upload Seurat R File", multiple=F, accept=".Robj, .RData"),
      
      selectInput("Module", "Chose a graph type:", choices=c("UMAP", "tSNE", "FeaturePlot", "Violin Plot", "Dot Plot", "Heatmap", "Ridge Plot")),
      
      conditionalPanel(
        condition="input.Module=='FeaturePlot', 'Violin Plot', 'Dot Plot', 'Heatmap', 'Ridge Plot'",
        selectInput("scGene",
                    label="Select Genes",
                    choices=NULL,
                    selected=NULL,
                    multiple=T))
    ),
    mainPanel(
      tabsetPanel(type="tabs",
                  tabPanel("About",
                           h3("How to upload data:"),
                           p("The file type used on this site is one generated on the", a("Computational Suite for Bioinformaticians and Biologists (CSBB)", href="https://github.com/praneet1988/CSBB-Shiny", target="_blank"), ". Users will create a single cell file based on their own paramters using CSBB. The file type created on CSBB is an R object, or .Robj, file. That is the only file that can be uploaded on this site. The purpose of this site is to create an easy way for plot analysis with previously generated single cell files. "),
                           h3("Common parameters used by our lab on CSBB:"),
                           p("First, change the Module to SingleCell RNASeq Analysis. The Raw Counts Matrix will be your file type if you are using bulk RNA-seq data or a publicly available dataset. Your Expression File is a .txt file. Next, select your species."),
                           p("General parameters are set by default on the CSBB site. These parameters are meant to be adjusted and worked with. There is no right answer on what numbers you should use since each dataset is unique. However, below are some common parameters used by our lab:"),
                           p("1. Input minimum number of cells to express all genes: Try 3 at first but it is acceptable to go to 5."),
                           p("2. Input minimum number of features all cells should express: Try 100 at first but it is acceptable to go up to 200."),
                           p("3. Use LogNormalizaiton"),
                           p("4. Dimension reduction: UMAP is generally more preferred nowadays but it will be dependent upon your analysis."),
                           p("5. Regress cell cycle effect: Yes = removes cell cycle effect; No = keeps cell cycle effect. Our lab is typically not interested in this effect as we are primarily interested in how cell fates change due to loss of cilium."),
                           p("6. Input resolution for clustering: Try 0.8 or 1 but it is acceptable to use 0.5."),
                           p("Your .Robj file type can now be downloaded and used on this site. Again, these parameters are just suggestions for creating stringent analysis. These parameters are meanted to be played with until you believe a conclusion can be drawn from the data. If you are unsure how to use these numbers, there a link on the CSBB site (which is linked above) to a Youtube video for how to use the site."),
                           h3("Which plot should I use?"),
                           p(strong("UMAP:"), "This plot preserves both local AND global structure in the data. This means you can make comparisons on the distance and positions between points within a single cluster AND make comparisons between the distance and positions of two different clusters. This is a newer method for analysis and is typically the most accepted plot type nowadays between UMAP or tSNE."),
                           p(strong("tSNE:"), "This plot preserves the local structure of the data. This means you can compare the distance between points within a single cluster. You cannot make comparisons between the distance of 2 different clusters."),
                           p(strong("Feature Plot:"), "This plot will show you the expression of your gene of interest by coloring single cells. With this plot, you can see much your gene of interest is expressed and if it localizes to any specific clusters. To see if it localizes to any clusters, compare the expression levels to the cluster number in a tSNE plot."),
                           p(strong("Violin Plot:"), "Similar to a box plot but it shows the full distribution of the data for your gene of interest. The black dots represent individual cells for the gene and the colored violin shows the distribution of the data. The cells do not correlate to the actual violin."),
                           p(strong("Dot Plot:"), "Similar to a heatmap but shows the average expression of genes in the cells for your gene of interest. Darker expression color represents higher average gene expression from the cells. Larger dots indicate that gene was detected in a greater proportion of cells."),
                           p(strong("Top 10 Markers Heatmap:"), "Creates a heatmap of the top 10 marker expressions in each cluster. Marker genes are shown on the left and cluster numbers are at the top."),
                           p(strong("Ridge Plot:"), "Similar to a violin plot but does not show individual cells. A ridge plot shows the average expression density of your gene of interest in each cluster. Clusters are numbered on the left and expression levels are at the bottom."),
                           h3("Questions to ask yourself while examining plots:"),
                           p("Are the parameters you used on CSBB the same between wild type and mutant?"),
                           p("Are your parameters stringent enough to draw conclusions from these plots?"),
                           p("Are the genes you'd expect to be present in this tissue and at this time point expressed?"),
                           p("In your UMAP/tSNE plots, are there differences in cluster spacing and orientation? Are there differences between how clusters look in wild type and mutant?"),
                           p("Violin plots, Dot plots, Heatmaps, and Ridge plots essentially show the same thing: how is your gene of interest expressed in each cluster? Which one of these plots is the best way to represent your data and the conclusions you're trying to draw?"),
                           br(),
                           p("Large file sizes are normal for single cell analysis. If the application is not responding (especially when looking at Top 10 Markers Heatmap and the Results table), give it a few minutes."),
                           br(),
                           p("Developed by Preston Schultz as a senior capstone project by the University of Cincinnati and Dr. Samantha Brugmann's lab at Cincinnati Children's Hospital.")),
                          p("For questions, please contact Preston using the following email: prestonaschultz@gmail.com"),
                  tabPanel("Plots", downloadButton('downloadPlot', 'Download Plot'), plotOutput("contents")),
                  tabPanel("Results/Data Table", DT::dataTableOutput("result"), downloadButton('downloadResult', 'Download Results'))
      )
    )
  )
)


server=function(input,output,session) {
  
  observe({
    inFile=input$file1
    if (is.null(inFile))
      return(NULL)
    load(inFile$datapath)
    seurat.object=UpdateSeuratObject(seurat.object)
    genes=rownames(seurat.object)
    genes=data.frame(genes)
    genes=genes$genes
    updateSelectInput(session, "scGene",
                      label="Select Genes",
                      choices=genes)
  })
  
  SingleCell=reactive({
    inFile=input$file1
    if (is.null(inFile))
      return(NULL)
    load(inFile$datapath)
    seurat.object=UpdateSeuratObject(seurat.object)
    DefaultAssay(seurat.object)="RNA"
    Idents(seurat.object)=seurat.object@meta.data$seurat_clusters
    return(seurat.object)
  })
  
  SingleCellMarkers=reactive({
    seurat.object=SingleCell()
    if (is.null(seurat.object))
      return(NULL)
    marker_clusters=FindAllMarkers(seurat.object)
    return(marker_clusters)
  })
  
  SingleCellPlot=reactive({
    seurat.object=SingleCell()
    if(input$Module=="UMAP"){
      UMAPPlot(seurat.object, label=T, pt.size=2)
    }
    else if(input$Module=="tSNE"){
      TSNEPlot(seurat.object, label=T, pt.size=2)
    }
    else if(input$Module=="FeaturePlot"){
      FeaturePlot(seurat.object, features=input$scGene, pt.size=2, order=T)
    }
    else if(input$Module=="Violin Plot"){
      VlnPlot(seurat.object, features=input$scGene) + RotatedAxis()
    }
    else if(input$Module=="Dot Plot"){
      DotPlot(seurat.object, features=input$scGene) + RotatedAxis()
    }
    else if(input$Module=="Heatmap"){
      top10=SingleCellMarkers() %>% group_by(cluster)	%>% top_n(n=10, wt=avg_logFC)
      DoHeatmap(seurat.object, features=top10$gene)
    }
    else if(input$Module=="Ridge Plot"){
      RidgePlot(seurat.object, features=input$scGene, ncol=2)
    }
  })
  
  output$contents=renderPlot({
    SingleCellPlot()
  })
  
  output$result=DT::renderDataTable({
    DT::datatable(SingleCellMarkers())
  })
  
  output$downloadResult=downloadHandler(
    filename=function(){
      paste0(input$Module, "_Analysis_Result", "-", Sys.Date(), ".txt")
    },
    content=function(file){
      write.table(SingleCellMarkers(), file, sep="\t", quote=F)
    }
  )
  
  output$downloadPlot=downloadHandler(
    filename=function(){
      paste0(input$Module, "_Plot", "-", Sys.Date(), ".png")
    },
    content=function(file){
      SingleCellPlot()
      ggsave(file, width=15, height=15)
    },
    contentType='image/png'
  )
}

shinyApp(ui = ui, server = server)
