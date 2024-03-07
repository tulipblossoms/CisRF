#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(DT)
library(shinythemes)
library(plotly)
library(pheatmap)
organs<-c("Breast","Breast","Endocrine","Endocrine","Endocrine","Endocrine","Endocrine","Endocrine","Endocrine","Endocrine","Endocrine","Endocrine","Endocrine","Endocrine","Endocrine","Endocrine","Endocrine","Endocrine","Endocrine","Endocrine","Endocrine","Endocrine","Endocrine","Endocrine","Endocrine","Endocrine","Endocrine","Endocrine","Fibroblast","Fibroblast","Fibroblast","Fibroblast","Fibroblast","Fibroblast","Fibroblast","Vessel","Vessel","Vessel","Vessel","Vessel","Vessel","Nervous","Nervous","Nervous","Nervous","Nervous","Nervous","Nervous","Nervous","Nervous","Nervous","Nervous","Lymphocyte","Lymphocyte","Lymphocyte","Lymphocyte","Lymphocyte","Lymphocyte","Lymphocyte","Lymphocyte","Lymphocyte","Lymphocyte","Lymphocyte","Lymphocyte","Lymphocyte","Lymphocyte","Lymphocyte","Skin","Skin","Skin","Skin","Skin","Skin","Skin","Skin","Skin","Cardiac Muscle","Cardiac Muscle","Myeloid","Myeloid","Myeloid","Colon","Colon","Colon","Colon","Colon","Colon","Colon","Colon","Colon","Colon","Colon","Colon","Colon","Colon","Colon","Embryo","Embryo","Embryo","Embryo","Embryo","Embryo","Embryo","Embryo","Embryo","Embryo","Embryo","Embryo","Embryo","Embryo","Embryo","Embryo","Embryo","Embryo","Embryo","Embryo","Embryo","Embryo","Embryo","Embryo","Embryo","Embryo","Embryo","Embryo","Lung","Lung","Lung","Lung","Lung","Lung","Lung","Lung","Lung","Lung","Lung","Lung","Lung","Lung","Lung","Lung","Lung","Lung","Lung","Lung","Lung","Lung","Lung","Upper Muscle","Upper Muscle","Upper Muscle","Upper Muscle","Upper Muscle","Upper Muscle","Upper Muscle","Upper Muscle","Upper Muscle","Upper Muscle","Upper Muscle","Skeletal Muscle","Skeletal Muscle","Skeletal Muscle","Skeletal Muscle","Skeletal Muscle","Heart","Heart","Ventricle","Ventricle","Ventricle","Ventricle","Ventricle","Ventricle","Ventricle","Ventricle","Immune","Immune","Immune","Immune","Immune","Immune","Immune","Immune","Immune","Immune","Immune","Immune","Immune","Liver","Liver","Liver","Liver","Liver","Liver","Liver","Muscle","Muscle","Muscle","Muscle","Muscle","Muscle","Muscle","Muscle","Muscle","Muscle","Muscle","Muscle","Muscle","Muscle","Muscle","Muscle","Muscle","Muscle","Kidney","Kidney","Kidney","Kidney","Kidney","Kidney","Kidney","Kidney","Kidney","Kidney","Kidney","Kidney","Kidney","Kidney","Kidney","Kidney","Kidney","Kidney","Kidney","Kidney","Kidney","Kidney","Kidney","Kidney","Kidney","Kidney","Atrium","Atrium","Atrium","Atrium","Atrium","ES","ES","Lower Muscle","Lower Muscle","Lower Muscle","Lower Muscle","Lower Muscle","Lower Muscle","Lower Muscle","Lower Muscle","Gut","Gut","Gut","Gut","Gut","Gut","Gut","Gut","Gut","Gut","Gut","Gut","Gut","Gut","Gut","Gut","Gut","Gut","Reproductive","Reproductive","Reproductive")
snp_map<-read.csv("C:/Users/kelly/R/GenePrediction/training_data_gwas_genes/training_data_gwas_genes/full_snp_map_1881.csv")
colnames(snp_map)<-snp_map[1,]
snp_map<-snp_map[,2:ncol(snp_map)]
qualified_genes<-read.csv("C:/Users/kelly/R/GenePrediction/matched_snp_ccres/qualified_genes2.csv")
qualified_genes<-qualified_genes[,5:ncol(qualified_genes)]
setwd("C:/Users/kelly/R/GenePrediction/training_data_gwas_genes/training_data_gwas_genes")

# Define UI for application that draws
ui <- dashboardPage(skin="green",
                    dashboardHeader(title="CisRF Shiny",
                                    titleWidth=200
                    ),
                    dashboardSidebar(
                      sidebarMenu(
                        menuItem("home", tabName = "Home", icon = icon("tree")),
                        menuItem("snp_explore", tabName = "snp_explore", icon = icon("car")),
                        menuItem("user_guide",tabName="user_guide",icon=icon("dashboard"))    
                      )
                    ),
                    dashboardBody(
                      tabItems(
                        tabItem("Home",
                                fluidPage(
                                fluidRow(
                                    column(width=12,
                                           h1("Welcome to CisRF Shiny!"),
                                    h1("
                                        Predictions of effects of regulatory-element SNPs on gene expression."),
                                    box(width = 12,solidHeader=TRUE,title="CisRF Overview",status="primary",
                                    imageOutput("cover",height="auto")
                                    ))
                                    ),
                                  fluidRow(
                                    column(width=12,
                                           box(width=12,solidHeader=TRUE,title="Heatmap: CisRF predictions",collapsible=TRUE,status="primary",
                                               imageOutput("heatmap",height="auto"),
                                               textOutput("heatmap_description")))
                                  )
                                  ,fluidRow(column(width=4,
                                                   valueBoxOutput("valueBox1",width=12)),
                                            column(width=4,valueBoxOutput("valueBox2",width=12)),
                                            column(width=4,valueBoxOutput("valueBox3",width=12))
                                            )
                                )

                                  ),
                          tabItem("snp_explore",
                            fluidPage(
                            h1("explore the effects of SNPs on context-dependent gene expression."),
                            box(width=12,solidHeader=TRUE,status="primary",title="select a snp.",collapsible = TRUE,
                            dataTableOutput("snp_list"),
                            downloadButton("download_data", "Download Table")
                            ),
                            fluidRow(
                              column(width=8,
                                     box(width=12,title="Data Visualization",
                                     plotOutput("violin_plot")),
                                     height = "500px"
                                     ),
                              column(width=4,
                                     box(width=12,title="Mean Expression Across Clusters",
                                     plotOutput("barChart")),
                                     verbatimTextOutput("stats"))
                            )
                            )
                        ),
                        tabItem("user_guide",
                                fluidPage(
                                  fluidRow(
                                    box(width=12,solidHeader=TRUE,status="primary",title="User Instructions",collapsible=TRUE,
                                        imageOutput("user_instructions",height="auto"))
                                  ),
                                  fluidRow(
                                    column(offset=1,width=12,
                                    textOutput("github_info")
                                  ))
                                )
                                )
                        )
                      
                    
))





# Define server logic required to draw a histogram
server <- function(input, output) {
  output$cover<-renderImage({
    img_source<-"C:/Users/kelly/R/GenePrediction/CisRFDashboard/cover_image.jpg"
    width<-"100%"
    list(src=img_source,
         contentType="image/jpg",
         width=width,
         height="auto")
  },deleteFile=FALSE)
  
  output$heatmap<-renderImage({
  img_source<-"C:/Users/kelly/R/GenePrediction/CisRFDashboard/1881_pheatmap.jpg"
    width<-"100%"
  list(src=img_source,
       contentType="image/jpg",
       width=width,
       height="auto")
  },deleteFile=FALSE)

  output$heatmap_description<-renderText({
    "CisRF predictions on the impact of gene expression across 25 clusters of contexts, for 1881 SNPs. Each row is a SNP, and each column is a cluster of contexts."
  })
  output$valueBox1<-renderValueBox({
    valueBox(
      value=20,
      subtitle="Traits",
      icon=icon("heart"),
      color="blue"
    )
  })
  output$valueBox2<-renderValueBox({
    valueBox(
      value=1881,
      subtitle="SNPs",
      icon=icon("dna"),
      color="maroon"
    )
  })
  output$valueBox3<-renderValueBox({
    valueBox(
      value=1103,
      subtitle="Genes",
      icon=icon("flask"),
      color="purple"
    )
  })
  output$snp_list<-DT::renderDataTable(qualified_genes,server=FALSE)
  output$barChart<-renderPlot({
    if(length(input$snp_list_rows_selected)<=1){
      values<-as.numeric(snp_map[input$snp_list_rows_selected,])
      dataFrame<-data.frame(
        category=organs,
        value =values
        )
      ggplot(dataFrame,aes(x=category,y=value))+
        geom_bar(stat="identity",fill="skyblue")+
        labs(title=paste0("Chromosome #: ",qualified_genes$Chr.Number[input$snp_list_rows_selected],"        SNP Position: ",qualified_genes$SNP.Position[input$snp_list_rows_selected]), x="Clusters",y="Mean Expression")+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
    }else{
      plot(NULL,xlim = c(0, 1), ylim = c(0, 1))
    }
  })
  
  output$violin_plot<-renderPlot({
    if(length(input$snp_list_rows_selected)<=1){
    row<-input$snp_list_rows_selected
    dataSet<-data.frame(
      category=organs,
      values = abs(as.numeric(snp_map[row+1,]))
    )

    dataSet$category <- reorder(dataSet$category, -dataSet$values,fun=mean)
    ggplot(dataSet,mapping=aes(x=category,y=values,fill=category))+geom_violin(scale="width",width=0.5)+geom_boxplot(width=0.1,fill="white",color="black")+labs(title=paste0("Chromosome #: ",qualified_genes$Chr.Number[input$snp_list_rows_selected],"        SNP Position: ",qualified_genes$SNP.Position[input$snp_list_rows_selected]),x="",y="Difference")+theme_minimal()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
      theme(axis.title=element_text(size=25))+theme(plot.title=element_text(size=35,hjust=0.5))+theme(axis.text.x=element_text(size=c(20,20)))+theme(axis.text.y=element_text(size=c(15,15))+theme(legend.position="none")+theme(axis.text.x = element_text(angle = 45, hjust = 1)))+
      theme(axis.text.x = element_text(angle = 60, hjust = 1))+theme(legend.position="none")+stat_summary(fun=mean,geom="point",shape=18,size=3,color="red")
    
     }else if(length(input$snp_list_rows_selected)>1){
      s = input$snp_list_rows_selected
      heatmapData<-c()
      dataLabels<-c()
      clusters<-unique(organs)
      for(i in s){
        meanNums<-c()
        for(j in clusters){
          meanNums<-c(meanNums,mean(as.numeric(snp_map[i+1,which(organs==j)])))
        }
        heatmapData<-rbind(heatmapData,meanNums)
      }
      pheat<-pheatmap(heatmapData,main="Heatmap",fontsize=20,labels_col = c(clusters),labels_row=c(qualified_genes$snp_ids[s]))
    }
    })
  output$stats<-renderPrint({
    if(length(input$snp_list_rows_selected>=1)){
      s = input$snp_list_rows_selected
      paste0("Gene Used:" ,qualified_genes$Genes[s])
    }else{
      paste0("Genes Used: none")
    }
  })
  output$download_data <- downloadHandler(
    filename=function(){
      paste0("snp_data",".csv")
    },
    content=function(file){
      write.csv(qualified_genes, file)
    }
  )
  output$user_instructions<-renderImage({
    img_source<-"C:/Users/kelly/R/GenePrediction/CisRFDashboard/CisRF_Shiny_Guide.jpg"
    width<-"100%"
    list(src=img_source,
         contentType="image/jpg",
         width=width,
         height="auto")
  },deleteFile=FALSE)
  output$github_info<-renderText({
    "All code used to train CisRF and create CisRF Shiny is available on GitHub at this link: https://github.com/tulipblossoms/CisRF.git. Thank you for your interest in CisRF Shiny! 
"
  })
}
# Run the application 
shinyApp(ui = ui, server = server)

