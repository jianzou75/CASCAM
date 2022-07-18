requiredPackages = c('shiny','boot','tidyverse', 'ggrepel', 'patchwork',
                     'cowplot', 'reshape2', 'ggridges', 'sparseLDA', 'pheatmap')
for(p in requiredPackages){
  if(!require(p, character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

requiredBioconductor = c("BiocManager", 'pathview')
for(p in requiredBioconductor){
  if(!require(p, character.only = TRUE)) BiocManager::install(p)
  library(p, character.only = TRUE)
}

library(CASCAM)
set.seed(12345)
load("./www/KEGG_name_ID_match.RData")

options(shiny.maxRequestSize=3000*1024^2)
lsdata <- function (fnam = ".Rdata"){
  x <- load(fnam, envir = environment())
  return(x)
}


server <- function(input, output, session) {
  output$select_CASCAM_object <- renderUI({
    req(input$dataset$datapath)
    selectInput(inputId = 'CASCAM_object', "Select the CASCAM object", choices = c("", lsdata(input$dataset$datapath)))
  })

  observe({
    req(input$dataset$datapath, input$CASCAM_object)
    e <- new.env()
    name <- load(input$dataset$datapath, envir = e)
    CASCAM <- e[[input$CASCAM_object]]

    output$select_interested_pathway <- renderUI({
      req(pathway_result())
      selectInput(inputId = 'interested_pathway', "Select the interested pathway",
                  choices =  c("", pathway_result()@available_pathways))
    })

    output$select_interested_pathway2 <- renderUI({
      req(pathway_result(), input$interested_pathway)
      selectInput(inputId = 'interested_pathway2', "Select the interested pathway",
                  choices =  c("", pathway_result()@available_pathways), selected = input$interested_pathway)
    })

   output$select_interested_camod <- renderUI({
      req(pathway_result(), input$interested_pathway)
      selectizeInput(inputId = 'interested_camod', "Select the interested cancer models",
                     choices =  c("", pathway_result()@selected_camods), multiple = TRUE, options = list(maxItems = 5))
    })

   output$select_interested_pathway3 <- renderUI({
     req(pathway_result(), input$interested_pathway)
     selectInput(inputId = 'interested_pathway3', "Select the interested pathway",
                 choices =  c("", pathway_result()@available_pathways), selected = input$interested_pathway)
   })

   output$select_interested_camod2 <- renderUI({
     req(pathway_result(), input$interested_pathway)
     selectInput(inputId = 'interested_camod2', "Select the interested cancer model",
                 choices =  c("", pathway_result()@selected_camods))
   })

    ## Genome analysis
    genome_result <- eventReactive(input$genome_analysis_start, {
      input$genome_analysis_start

      genome_selection(CASCAM)
    })
    output$genome_visualize <- renderPlot({
      CASCAM <- genome_selection_visualize(genome_result())
      plot_grid(CASCAM@genome_figure$sda_project_position,
                CASCAM@genome_figure$sda_ds_ci_rank,
                align = "h", axis = "bt", rel_widths = c(1.1, 1))
    }, height = session$clientData$output_genome_visualize_width/2)
    output$genome_text <- renderText({
      req(genome_result())
      "The left figure is for SDA projected deviance score visualization.
      Red circles are the one classified as interested subtype by the combination of SDA classification (assignment probability) and deviance score (p-value).
      Color represents the classification results and point shape represents the deviance score results.
      The SDA projected tumor distribution is drawn on the right.
      The right figure is for SDA projected deviance score with confidence interval on the genome-wide pre-selected biological models.
      The bootstrap with 1,000 times on the tumor data is performed to obtain the confidence interval, and the 95% confidence intervals are presented."
    })


    ## pathway analysis
    pathway_result <- eventReactive(input$genome_analysis_start,{
      req(genome_result())
      pathway_analysis(genome_result())
    })

    datasetOutput <- reactive({
      req(pathway_result())
      interested_subtype <- pathway_result()@interested_subtype
      uninterested_subtype <- setdiff(pathway_result()@sda_model$classes, pathway_result()@interested_subtype)

      genome_frame <- data.frame(rownames(pathway_result()@sda_ds),
                                 pathway_result()@sda_ds[,interested_subtype],
                                 pathway_result()@sda_ds[,uninterested_subtype],
                                 pathway_result()@sda_predict_prob[,interested_subtype],
                                 pathway_result()@sda_predict_prob[,uninterested_subtype],
                                 pathway_result()@sda_ds_pval[,interested_subtype],
                                 pathway_result()@sda_ds_pval[,uninterested_subtype])
      colnames(genome_frame) <- c("camods", paste0("SDA_DS_",interested_subtype),
                                  paste0("SDA_DS_",uninterested_subtype), paste0("SDA_P_",interested_subtype),
                                  paste0("SDA_P_",uninterested_subtype), paste0("SDA_DS_pval_",interested_subtype),
                                  paste0("SDA_DS_pval_",uninterested_subtype))

      pathway_interested_frame <- data.frame(pathway = pathway_result()@available_pathways,
                                             pathway_result()@pathway_ds[[interested_subtype]])
      pathway_uninterested_frame <- data.frame(pathway = pathway_result()@available_pathways,
                                               pathway_result()@pathway_ds[[uninterested_subtype]])

      gene_interested_frame <- data.frame(gene = rownames(pathway_result()@gene_ds[[interested_subtype]]),
                                          pathway_result()@gene_ds[[interested_subtype]])
      gene_uninterested_frame <- data.frame(gene = rownames(pathway_result()@gene_ds[[interested_subtype]]),
                                            pathway_result()@gene_ds[[uninterested_subtype]])

      list_of_output_data <- list("Genome" = genome_frame,
                                  "Pathway_interested_subtype" = pathway_interested_frame,
                                  "Pathway_uninterested_subtype" = pathway_uninterested_frame,
                                  "Gene_interested_subtype" = gene_interested_frame,
                                  "Gene_uninterested_subtype" = gene_uninterested_frame)
    })

    output$downloadData <- downloadHandler(
      filename = "analysis_results.xlsx",
      content = function(file) {
        openxlsx::write.xlsx(datasetOutput(), file)
      }
    )

    output$download <- renderUI({
      req(pathway_result())
      downloadButton('downloadData', 'Download Analysis Results')
    })

    output$pathway_congruence_heatmap_text <- renderText({
      req(pathway_result())
      "The heatmap below shows a general overview of the congruence of the cancer models to the tumor centers in different avaiable pathways.
       The first row represents the genome-wide deviance score (the smaller, the better).
       Pathway_Size shows the number of genes included for each pathway.
       NES shows the normalized enrichment score obtained from GSEA by inputting the log fold change of the tumor data. Positive means interested subtype is upregulated comparing to uninterested one; negative means interested subtype is downregulated.
       In the main figure, the color represents the pathway specific deviance score, and smaller values (more red) indicates better congruence."
    })

    output$pathway_congruence_heatmap <- renderImage({
      req(pathway_result())

      width  <- session$clientData$output_pathway_congruence_heatmap_width
      height <- session$clientData$output_pathway_congruence_heatmap_height
      pixelratio <- session$clientData$pixelratio

      outfile <- tempfile(fileext = '.png')

      ### Generate the png
      png(outfile, width = width*pixelratio, height = width/50 * length(pathway_result()@available_pathways) * pixelratio,
          res = 72*pixelratio)
      path_heatmap  <- pathway_congruence_heatmap(pathway_result())
      dev.off()

      ### Return a list containing the filename
      list(src = outfile,
           contentType = 'image/png',
           width = width,
           height = width/50 * length(pathway_result()@available_pathways),
           alt = "This is alternate text")
    }, deleteFile = TRUE)

    output$pathway_specific_heatmap_text <- renderText({
      req(input$interested_pathway)
      paste0("The heatmap below shows the gene specific deviance score for the differential expression genes, which is the smaller (more red) the better, in ", input$interested_pathway,
             " among the genome-wide pre-filtered cancer models. The first row shows the pathway specific deviance score.")
    })

    output$pathway_specific_heatmap <- renderImage({
      req(input$interested_pathway)

      width  <- session$clientData$output_pathway_specific_heatmap_width
      height <- session$clientData$output_pathway_specific_heatmap_height
      pixelratio <- session$clientData$pixelratio

      pathways <- c(qusage::read.gmt(system.file(package = 'CASCAM', 'extdata/c2.cp.kegg.v7.4.symbols.gmt')),
                    qusage::read.gmt(system.file(package = 'CASCAM', 'extdata/h.all.v7.4.symbols.gmt')))
      gene_pathway <- intersect(rownames(pathway_result()@gene_ds[[pathway_result()@interested_subtype]]),
                                pathways[[input$interested_pathway]])

      outfile <- tempfile(fileext = '.png')

      ### Generate the png
      png(outfile, width = width*pixelratio, height = width/25 * length(gene_pathway) * pixelratio,
          res = 72*pixelratio)
      path_specific_heatmap  <- pathway_specific_heatmap(pathway_result(), input$interested_pathway)
      dev.off()

      ### Return a list containing the filename
      list(src = outfile,
           contentType = 'image/png',
           width = width,
           height = width/25 * length(gene_pathway),
           alt = "This is alternate text")
    }, deleteFile = TRUE)

   output$pathway_specific_distribution_text <- renderText({
     req(input$interested_camod, input$interested_pathway2)
      paste0("The ridgeline figure below shows the detailed distribution on the differential expression genes in ", input$interested_pathway,
             " among the genome-wide pre-filtered cancer models. The genes are sorted by the adjusted p-value in DEA from the smallest to the largest.")
    })

  output$pathway_specific_distribution <- renderImage({
    req(input$interested_camod, input$interested_pathway2)

    width  <- session$clientData$output_pathway_specific_distribution_width
    height <- session$clientData$output_pathway_specific_distribution_height
    pixelratio <- session$clientData$pixelratio

    pathways <- c(qusage::read.gmt(system.file(package = 'CASCAM', 'extdata/c2.cp.kegg.v7.4.symbols.gmt')),
                  qusage::read.gmt(system.file(package = 'CASCAM', 'extdata/h.all.v7.4.symbols.gmt')))
    gene_pathway <- intersect(rownames(pathway_result()@gene_ds[[pathway_result()@interested_subtype]]),
                              pathways[[input$interested_pathway]])

    outfile <- tempfile(fileext = '.png')

    png(outfile, width = width*pixelratio, height = width/20 * length(gene_pathway) * pixelratio,
        res = 72*pixelratio)
    plot(pathway_specific_ridgeline(pathway_result(), pathway_name = input$interested_pathway2,
                                    interested_camods = input$interested_camod)@pathway_gene_ridgeline)
    dev.off()

    ### Return a list containing the filename
    list(src = outfile,
         contentType = 'image/png',
         width = width,
         height = width/20 * length(gene_pathway),
         alt = "This is alternate text")
  }, deleteFile = TRUE)


  output$pathway_specific_pathview <- renderImage({
    req(input$interested_camod2, input$interested_pathway3, pathway_result())
    pathview_input <- pathway_result()@gene_ds[[pathway_result()@interested_subtype]][,input$interested_camod2]
    kegg_id = KEGG_name_ID_match$gs_exact_source[KEGG_name_ID_match$gs_name == input$interested_pathway3]


    #Save the file in the temp folder
    p = pathview(gene.data = pathview_input, pathway.id = kegg_id,
                 species = "hsa", out.suffix = "pathview", gene.idtype = "SYMBOL",
                 limit=list(gene=c(-3,3)), bins = list(gene = 16),
                 low = list(gene = "#B2E6F0"), mid = list(gene = "#FF0018"), high = list(gene = "#FFF485"),
                 node.sum = "mean", kegg.dir = tempdir())

    list(src = paste0("./hsa", kegg_id,".pathview.png"),
         contentType = 'image/png',
         alt = "This is alternate text",
         width = "120%")
  }, deleteFile = TRUE)

})


}
