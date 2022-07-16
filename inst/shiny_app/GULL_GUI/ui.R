library(shiny)
library(shinythemes)
customDownloadbutton <- function(outputId, label = "Download"){
  tags$a(id = outputId, class = "btn btn-default shiny-download-link", href = "",
         target = "_blank", download = NA, icon("accessible-icon"), label)
}

## Intro Panel
intro_panel <- tabPanel("Introduction",
                       h4("Workflow"),
                       HTML('<center><img src="Congruence.svg" width = "60%"></center>'),
                       br(),
                       column(8, offset = 2,
                              tabPanel("Introduction",
                                 "BioModelSelect is an interactive tool to select the most appropritate biological models based on genome-wide and pathway-specific considerations.
                                 Considering the running time, the data pre-processing step is needed beforehand. An RData file should be prepared prior to use this shiny app for exploration.
                                 In the RData file, the aligned tumor data matrix with sample labels, biological models data matrix, an differential expression analysis object, and a tuned SDA model should be included.
                                 Details of the file preparation can be referred to our tutorial.")))


## Genome Panel
genome_input <- sidebarPanel(
  h4("Data input"),
  fileInput("dataset", "Choose .RData including an GULL object.",
            accept = ".RData"),
  uiOutput('select_GULL_object'),

  actionButton("genome_analysis_start", "Submit")
)
genome_figure <- mainPanel(
  h3("Genome-wide selection"),
  textOutput("genome_text"),
  br(),
  column(9, align="center",
         plotOutput("genome_visualize", width = "120%"))
)
genome_panel <- tabPanel("Genome-wide pre-selection",
                         sidebarLayout(genome_input, genome_figure))


## Pathway Panel
pathway_input <- sidebarPanel(
  conditionalPanel(condition = "input.pathway == 1",uiOutput('download')),
  conditionalPanel(condition = "input.pathway == 2",uiOutput('select_interested_pathway')),
  conditionalPanel(condition = "input.pathway == 3",uiOutput('select_interested_pathway2'),uiOutput('select_interested_cell')),
  conditionalPanel(condition = "input.pathway == 4",uiOutput('select_interested_pathway3'), uiOutput('select_interested_cell2'))
)
pathway_figure <- mainPanel(
  h3("Pathway specific analysis"),
  tabsetPanel(type = "tabs", id = "pathway",
              tabPanel("Pathway congruence heatmap", value = 1, br(), textOutput("pathway_congruence_heatmap_text"), br(),
                       column(10, align="center", imageOutput("pathway_congruence_heatmap", width = "120%"))),
              tabPanel("Pathway specific heatmap", value = 2, br(), textOutput("pathway_specific_heatmap_text"), br(),
                       column(8, align="center", imageOutput("pathway_specific_heatmap", width = "120%"))),
              tabPanel("Pathway specific distribution", value = 3, column(8, align="center",
                       imageOutput("pathway_specific_distribution", width = "120%"))),
              tabPanel("KEGG PathView", value = 4, column(8, align="center", br(),
                        imageOutput("pathway_specific_pathview", width = "100%")))
              )
)
pathway_panel <- tabPanel("Pathway specific analysis",
                         sidebarLayout(pathway_input, pathway_figure))

ui <- fluidPage(
  theme = shinytheme("sandstone"),

  navbarPage(
    "Biological Model Selection",
    intro_panel,
    genome_panel,
    pathway_panel)
)

