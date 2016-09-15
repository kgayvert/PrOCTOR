# Define UI for application that draws a histogram
shinyUI(fluidPage(
  #  useShinyjs(),
  # Application title
  titlePanel("Interactive PrOCTOR Model Feature Testing Tool"),
  
  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      h5("Combined Model Prediction"),
      verbatimTextOutput("overallPredictionInfo"),
      h5("Structural Model Prediction"),
      verbatimTextOutput("structuralPredictionInfo"),
      h5("Target Model Prediction"),
      verbatimTextOutput("targetPredictionInfo"),
      br(),
      selectInput("var", 
                  label = "Select a variable to change",
                  choices = c("Original",unique(colnames(sample_drugs))),
                  selected = "Original"),
      textInput(inputId = "newval", label = h5("Enter a new value:"), value = ""),
      actionButton("goButton", "Go!"),
      checkboxInput("correlates", "Change Using Correlations", FALSE),
      br(),
      selectInput("drug",label = "Or select a known drug to use",choices = c("Median Values",sort(rownames(sample_drugs))),selected = "Median Values"),
      actionButton("preselectButton", "Choose Selection!"),
      br(),
      p("Click the button to update the value displayed in the main panel.")
    ),
    
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel('Feature Quantile Plot', plotOutput("featurePlot", height=600,click = "plot_click",dblclick = dblclickOpts(id = "plot_dblclick"))),
        tabPanel('Predictions', plotOutput("predictionPanel")),
        tabPanel('Structure Feature Values', dataTableOutput("structureFeatureValues")),
        tabPanel('Target Feature Values', dataTableOutput("targetFeatureValues")))#,
    )
  )
))


