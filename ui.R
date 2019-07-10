library(shiny)

# Define UI for application of pk1c ----

navbarPage(
  title = "Rivastigmine: Dataset Simulation (1-Comp model)",
  
  ## Data ----
  tabPanel(
    title = "Data",
    sidebarLayout(
      sidebarPanel(
        sliderInput("nSubj", "Subjects", min = 1, max = 100, value = 10),
        sliderInput("CL", "Clearance", min = 10, max = 300, value = 100, step = 10),
        sliderInput("V", "Volume of Distribution", min = 1, max = 500, value = 240, step = 10),
        sliderInput("Ka", "Absorption", min = 1, max = 20, value = 3, step=0.5),
        textInput("Time", "Sampling Time", value = "0, 0.25, 0.5, 1, 2, 4, 5, 7, 9, 12, 24"),
        textInput("DH1", "Dose Hx(t, Dose, t, Dose, ...)", value = "0, 100000, 12, 100000"),
        textInput("FullCov", "Full Covariance (2x2 matrix)", value = "0,0,0,0,0.49,0,0,0,0"), # CL, V, Ka
        sliderInput("PropE", "Proportional Error (SD)", min = 0, max = 1, value = 0, step = 0.1),
        sliderInput("AddE", "Additive Error (SD)", min = 0, max = 1, value = 0, step = 0.1),
        sliderInput("Jitter", "Time Jitter", min = 0, max = 10, value = 0, step = 1)
      ),
      
      mainPanel(
        tags$h2("Concentration-time Curves"),
        plotOutput("concTimePlot"),
        checkboxInput(inputId = "concLog", label = "Log scale"),
        tags$h2("NCA"),
        tags$h4("Noncompartmental Analysis of Simulated Dataset using `NonCompart` package"),
        tableOutput("ncarTable")
      )
    )
  ),
  
  ## Plot ----
  tabPanel(
    title = "Details",
    tags$h2("Individual Plots"),
    plotOutput("concTimeFacet", height = "800px"),
    tags$h2("Raw Data"),
    tableOutput("concTable")
  ),
  
  tabPanel(
    title = "NCA",
    includeMarkdown("parameters.md")
  ),
  
  ### 90 ###
  tabPanel(
    title = "Help", 
    includeMarkdown("README.md")
    #htmlOutput('README')
  ),
  ### 99 ###
  tabPanel(
    title = "Contact", 
    includeMarkdown("CONTACT.md")
  )
)

