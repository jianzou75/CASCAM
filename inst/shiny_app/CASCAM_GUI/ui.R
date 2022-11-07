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
                                 "CASCAM is an interactive tool to select the most appropriate cancer models to mimick the target tumor cohorts based on genome-wide and pathway-specific considerations.
                                 Considering the running time, particularly the Celligner alignment and the cross-validation for SDA, the data pre-processing step is needed beforehand.
                                 An RData file including a CASCAM object should be prepared before using this shiny app for exploration.
                                 Details of the file preparation can be referred to our tutorial.")))


## Genome Panel
genome_input <- sidebarPanel(
  h4("Data input"),
  fileInput("dataset", "Choose .RData including an CASCAM object.",
            accept = ".RData"),
  uiOutput('select_CASCAM_object'),

  h4("Criteria for genomic preselection"),
  numericInput("assignment_prob_cutoff", "Cutoff for assignment probability", 0.8, min = 0.5, max = 1),
  numericInput("min_sda_ds_pval", "Minimum value for the p-value of SDA based devaince score", 0.05, min = 0, max = 0.5),

  actionButton("genome_analysis_start", "Submit")
)
genome_figure <- mainPanel(
  h3("Genome-wide selection"),
  textOutput("genome_text"),
  br(),
  column(9, align="center",
         plotOutput("genome_visualize", width = "120%"))
)
genome_panel <- tabPanel("Genome-wide pre-selection by machine learning",
                         sidebarLayout(genome_input, genome_figure))


## Pathway Panel
pathway_input <- sidebarPanel(
  conditionalPanel(condition = "input.pathway == 1",uiOutput('download')),
  conditionalPanel(condition = "input.pathway == 2",uiOutput('select_interested_pathway')),
  conditionalPanel(condition = "input.pathway == 3",uiOutput('select_interested_pathway2'),uiOutput('select_interested_camod')),
  conditionalPanel(condition = "input.pathway == 4",uiOutput('select_interested_pathway3'), uiOutput('select_interested_camod2'))
)
pathway_figure <- mainPanel(
  h3("Pathway and mechanistic analysis"),
  tabsetPanel(type = "tabs", id = "pathway",
              tabPanel("Pathway congruence heatmap", value = 1, br(), textOutput("pathway_congruence_heatmap_text"), br(),
                       column(10, align="center", imageOutput("pathway_congruence_heatmap", width = "120%"))),
              tabPanel("Pathway specific heatmap", value = 2, br(), textOutput("pathway_specific_heatmap_text"), br(),
                       column(8, align="center", imageOutput("pathway_specific_heatmap", width = "120%"))),
              tabPanel("Pathway specific violin", value = 3, column(8, align="center",
                       imageOutput("pathway_specific_violin", width = "120%"))),
              tabPanel("KEGG PathView", value = 4, column(8, align="center", br(),
                        imageOutput("pathway_specific_pathview", width = "100%")))
              )
)
pathway_panel <- tabPanel("Pathway and mechanistic analysis",
                         sidebarLayout(pathway_input, pathway_figure))

ui <- fluidPage(
  theme = shinytheme("sandstone"),

  navbarPage(
    "CASCAM",
    intro_panel,
    genome_panel,
    pathway_panel)
)

