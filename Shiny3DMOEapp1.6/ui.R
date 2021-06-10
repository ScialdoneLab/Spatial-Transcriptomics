
ui<-secure_app(head_auth = tags$script(inactivity), fluidPage( #App structure, where you define the title, the buttons or searchers and declare the outputs (plots) that will appear
  #titlePanel("Mouse Main Olfactory Epithelium Atlas"),
  theme = shinytheme("lumen"),
  
  tabsetPanel(type = "tabs",
              
    #tabPanel("Main",
  
      #fluidRow(wellPanel("Mouse Olfactory Mucosa Atlas", style = "background-color:coral;font-size:40px")),
  
      #fluidRow(wellPanel("Manual", style = "background-color:gold;font-size:20px")),
      #fluidRow(wellPanel("", style = "background-color:white;font-size:20px", 
                        #tags$iframe(src="manual-1.0.pdf", width="900", height="300"))),
    #),
    
    tabPanel("3D Zones",
             
      fluidRow(wellPanel("Mouse Olfactory Mucosa Atlas", style = "background-color:coral;font-size:40px")),
  
      fluidRow(wellPanel("3D Zones", style = "background-color:gold;font-size:20px")),
  
      fluidRow(wellPanel("", style = "background-color:white", actionButton("topics3D", "see 3D zones"))),
  
      
    #mainPanel(plotlyOutput(outputId = "ThreeDStructure")),
      plotlyOutput(outputId = "ThreeDStructure"),
  
      radioButtons("projectionZones", "Slices Projection:",
                c("LML x DV" = "LMLxDV",
                  "AP x DV" = "APxDV",
                  "AP x LML" = "APxLML")),
      plotlyOutput(outputId = "ThreeDSlicesZones", height = "25%")
    ),
    
    tabPanel("Gene search",
    
      fluidRow(wellPanel("Gene search", style = "background-color:gold;font-size:20px")),
      #selectizeInput(inputId = "gene", label = "Introduce a gene name", choices=allAxes_detected, multiple=T),
      selectizeInput(inputId = "gene", label = "Introduce a gene name", choices=NULL, selected=1),
      plotOutput(outputId = "OneDPlots"),
      downloadButton('download1D_DV',"Download DV 1D data (all genes)"),
      downloadButton('download1D_AP',"Download AP 1D data (all genes)"),
      downloadButton('download1D_LML',"Download LML 1D data (all genes)"),
      fluidRow(wellPanel("3D Patterns", style = "background-color:gold;font-size:15px")),
      actionButton("Calculate3D", "Calculate 3D pattern"),
      plotlyOutput(outputId = "ThreeDPlot"),
      radioButtons("projection", "Slices Projection:",
                  c("LML x DV" = "LMLxDV",
                    "AP x DV" = "APxDV",
                    "AP x LML" = "APxLML")),
      plotlyOutput(outputId = "ThreeDSlices", height = "25%"),
      downloadButton('download3D',"Download 3D data"),
      plotOutput(outputId = "geneZoneInds"),
      downloadButton('downloadDoBs',"Download Degrees of belonging (all non-random genes)"),
      plotOutput(outputId = "cellTypes"),
      downloadButton('downloadCellTypes',"Download Single cell data"),
    ),
    
    tabPanel("Correlation between 2 genes",
      fluidRow(wellPanel("Correlation between 2 genes", style = "background-color:gold;font-size:15px")),
      selectizeInput(inputId = "gene1", label = "Gene 1", choices=NULL, selected=1),
      selectizeInput(inputId = "gene2", label = "Gene 2", choices=NULL, selected=1),
      plotOutput(outputId = "twoGenesCor")
    ),
    
    tabPanel("Manual",
             
             #fluidRow(wellPanel("Mouse Olfactory Mucosa Atlas", style = "background-color:coral;font-size:40px")),
             
             fluidRow(wellPanel("Manual", style = "background-color:gold;font-size:20px")),
             fluidRow(wellPanel("", style = "background-color:white;font-size:20px", 
                                tags$iframe(src="manual-1.0.pdf", width="900", height="300"))),
    )
  )
  #fluidRow(wellPanel("Odorant search", style = "background-color:gold;font-size:20px")),
  #selectizeInput(inputId = "odorant", label = "Introduce an odorant CAS number", choices=rownames(odorantsZones), multiple=T),
  #actionButton("searchOdorants", "load odorants patterns data (takes aprox. 1 min)"),
  #uiOutput("Calculate3Dod"),
  #plotlyOutput(outputId = "ThreeDPlotOdorants"),
  #radioButtons("projectionOdorants", "Slices Projection:",
  #             c("LML x DV" = "LMLxDV",
  #               "AP x DV" = "APxDV",
  #               "AP x LML" = "APxLML")),
  #plotlyOutput(outputId = "ThreeDSlicesOdorants", height = "25%"),
  #plotOutput(outputId = "odorantZoneInds")
))