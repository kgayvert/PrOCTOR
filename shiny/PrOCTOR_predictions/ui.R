# Define UI for application that draws a histogram
shinyUI(fluidPage(
  #  useShinyjs(),
  # Application title
  titlePanel("Interactive PrOCTOR Prediction Tool"),
  
  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      h5("Inputted Drug"),
      verbatimTextOutput("drugInfo"),           
      h5("Combined Model Prediction"),
      verbatimTextOutput("overallPredictionInfo"),
      h5("Structural Model Prediction"),
      verbatimTextOutput("structuralPredictionInfo"),
      h5("Target Model Prediction"),
      verbatimTextOutput("targetPredictionInfo"),
      br(),
      selectizeInput("target", 
                  label = "Select a drug target",
                  choices = sort(rownames(GTEx_targets)),
                  selected = NULL,
                  multiple=TRUE),
      textInput(inputId = "smiles", label = h5("Enter the SMILES:"), value = ""),
      actionButton("goButton", "Go!"),
      p("Click the button to update the value displayed in the main panel.")
    ),
    
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel('Predictions', plotOutput("predictionPanel")),
        tabPanel('Structure Feature Values', dataTableOutput("structureFeatureValues")),
        tabPanel('Target Feature Values', dataTableOutput("targetFeatureValues")))
    )
  )
))


