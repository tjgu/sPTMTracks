library(shiny)
library(plotly)
library(dplyr)

# credentials_db <- shinymanager::create_db(credentials = data.frame(user = "amhawkridge",
#                                                                    password = "mds_tracks",
#                                                                    stringsAsFactors = FALSE),
#                                           sqlite_path = "/scratch/g/tgu/svenkatesan/Hoffmeister_lab/MDSProteomics/app_credentials.sqlite")


info_card <- function(text) {
  bslib::card(style = "margin-top: 10px; margin-bottom: 12px;",
              bslib::card_body(style = "padding: 10px 14px; font-size: 13px; color: #333;",
                               text))
}

app_ui <- function(dataset_name) {
  lollipop_inst <- glue::glue("Each circle along the protein in the top lane represents a glycosylation site.
  There may be multiple glycans which can be bound to the protein at the site.
  The height of the lollipop is the mean log fold change of all glycans which pass both the log fold change and p value threshold at the site as set by the sliders on the left.
  The topology present in the second lane displays chains, helices and other features.
  The annotation lane at the bottom displays features annotated in UniProt such as disulfide bond locations, already known glycosylation, phosphorylation sites, motifs, etc.") %>%
    stringr::str_squish()

  alphafold_inst <- glue::glue("The color of the protein represents the confidence scores of the protein as predicted by AlphaFold.
  The spheres are glycosylation sites found in data and the color (red for up, blue for down, gray for non-significant) represents mean log fold change (same as in lollipop plot).") %>%
    stringr::str_squish()

  ui <- shiny::fluidPage(theme = bslib::bs_theme(version = 5),
                         shiny::div(style = paste("background: #1e3a5f; color: white;",
                                                  "padding: 12px 24px; border-radius: 8px;",
                                                  "margin-bottom: 20px; display: flex;",
                                                  "justify-content: space-between; align-items: center;"),
                                    shiny::div(shiny::span("sPTMTracks",
                                                           style = "font-size: 18px; font-weight: 600;"),
                                               shiny::span(glue::glue(" · {dataset_name}"),
                                                           style = "font-size: 14px; opacity: 0.75;"))),
                         shiny::sidebarLayout(sidebarPanel = shiny::sidebarPanel(shinyWidgets::pickerInput(inputId = "comparison", label = "Comparison", choices = NULL,
                                                                                                           options = shinyWidgets::pickerOptions(liveSearch   = TRUE,
                                                                                                                                                 liveSearchPlaceholder = "Search...",
                                                                                                                                                 style = "btn-outline-secondary btn-sm")),
                                                                                 shinyWidgets::pickerInput(inputId = "accession", label = "Gene", choices = NULL,
                                                                                                           options = shinyWidgets::pickerOptions(liveSearch   = TRUE,
                                                                                                                                                 liveSearchPlaceholder = "Search...",
                                                                                                                                                 style = "btn-outline-secondary btn-sm")),
                                                                                 shiny::sliderInput(inputId = "p_thresh", label = "Adjusted P-value threshold", min = 0, max = 1, value = 0.2, step = 0.01),
                                                                                 shiny::sliderInput(inputId = "fc_thresh", label = "Absolute Log2 Fold change threshold", min = 0, max = 5, value = 1, step = 0.1)),
                                              mainPanel = shiny::mainPanel(shiny::tabsetPanel(shiny::tabPanel(title = "Lollipop plot",
                                                                                                              info_card(lollipop_inst),
                                                                                                              bslib::card(style = "margin-top: 16px;",
                                                                                                                          bslib::card_body(shinycssloaders::withSpinner(plotly::plotlyOutput("protein_plot", height = "850px"), color = "#1e3a5f"),
                                                                                                                                           padding = 0))),
                                                                                              shiny::tabPanel(title = "Volcano",
                                                                                                              bslib::card(style = "margin-top: 16px;",
                                                                                                                          bslib::card_body(shinycssloaders::withSpinner(plotly::plotlyOutput("volcano_plot", height = "850px"), color = "#1e3a5f"),
                                                                                                                                           padding = 0))),
                                                                                              shiny::tabPanel(title = "Alphafold Viewer",
                                                                                                              info_card(alphafold_inst),
                                                                                                              bslib::card(style = "margin-top: 16px;",
                                                                                                                          bslib::card_body(shiny::uiOutput(outputId = "alphafold_legend"),
                                                                                                                                           shinycssloaders::withSpinner(ui_element = NGLVieweR::NGLVieweROutput(outputId = "alphafold_viewer", height = "850px"), color = "#1e3a5f"),
                                                                                                                                           padding = 0))),
                                                                                              shiny::tabPanel(title = "Functional annotation",
                                                                                                              shiny::div(style = "padding: 20px 10px;",
                                                                                                                         shiny::uiOutput("func_annotation_card"))),
                                                                                              shiny::tabPanel(title = "Heatmap",
                                                                                                              bslib::card(style = "margin-top: 16px;",
                                                                                                                          bslib::card_body(shinycssloaders::withSpinner(plotly::plotlyOutput(outputId = "protein_heatmap", height = "auto"), color = "#1e3a5f"),
                                                                                                                                           padding = 0)))))))

    # ui <- shinymanager::secure_app(ui)
    return(ui)
}

app_server <- function(sample_meta, res_list, site_summaries, sig_genes_list) {
  server <- function(input, output, session) {
    # shinymanager::secure_server(check_credentials = shinymanager::check_credentials(credentials_db))
    shiny::observe({
      comparisons <- sort(intersect(names(sig_genes_list), names(site_summaries)))
      shinyWidgets::updatePickerInput(session = session, inputId = "comparison", choices = comparisons, selected = comparisons[[1]])
    })

    res_data <- shiny::reactive({
      shiny::req(input$comparison)

      res_list[[input$comparison]]
    })

    sig_genes <- shiny::reactive({
      req(res_data(), input$fc_thresh, input$p_thresh)

      temp_df <- res_data() %>%
        dplyr::filter(adj.P.Val <= input$p_thresh & abs(logFC) >= input$fc_thresh) %>%
        tibble::rownames_to_column(var = "glycoform_id")

      temp_df
    })

    site_summary <- shiny::reactive({
      req(res_data(), input$fc_thresh, input$p_thresh)

      temp_df <- res_data() %>%
        dplyr::select(adj.P.Val, logFC, protein, gene, protein_site)

      temp_df <- temp_df %>%
        dplyr::mutate(is_sig = adj.P.Val <= input$p_thresh & abs(logFC) > input$fc_thresh) %>%
        dplyr::group_by(protein, gene, protein_site) %>%
        dplyr::summarise(n_glycoforms = dplyr::n(),
                         n_sig_glycoforms = sum(is_sig, na.rm = TRUE),
                         pct_sig_glycoforms = (n_sig_glycoforms/n_glycoforms)*100,
                         is_hit = any(is_sig, na.rm = TRUE),
                         mean_logFC = ifelse(is_hit, mean(logFC[is_sig], na.rm = TRUE), 0),
                         .groups = "drop") %>%
        dplyr::mutate(direction = dplyr::case_when(is_hit & mean_logFC > 0 ~ "up",
                                                   is_hit & mean_logFC < 0 ~ "down",
                                                   TRUE ~ "ns")) %>%
        dplyr::mutate(direction = factor(direction, levels = c("up", "down", "ns"))) %>%
        dplyr::group_by(protein) %>%
        dplyr::mutate(n_sites = dplyr::n(),
                      n_hit_sites = sum(is_hit),
                      pct_hit_sites = (n_hit_sites/n_sites)*100) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(dplyr::desc(is_hit), protein) %>%
        as.data.frame()

      temp_df
    })

    prot_info <- shiny::reactive({
      shiny::req(input$accession)

      httr::GET(url = glue::glue("https://rest.uniprot.org/uniprotkb/{input$accession}.json")) %>%
        httr::content()
    })

    shiny::observe({
      shiny::req(res_data())

      max_fc <- max(abs(res_data()$logFC), na.rm = TRUE)
      shiny::updateSliderInput(session = session, inputId = "fc_thresh", max = ceiling(max_fc), value = min(1, max_fc))
    })

    shiny::observeEvent(list(input$comparison, site_summary()), {
      gene_to_accession <- site_summary() %>%
        dplyr::distinct(gene, protein) %>%
        dplyr::arrange(gene, protein)

      current_accession <- shiny::isolate(input$accession)
      selected_accession <- if (!is.null(current_accession) && current_accession %in% gene_to_accession$protein) {
        current_accession
        } else {
          gene_to_accession$protein[[1]]
        }

      shinyWidgets::updatePickerInput(session = session, inputId = "accession", label = "Gene",
                                      choices = stats::setNames(gene_to_accession$protein, gene_to_accession$gene),
                                      selected = selected_accession)
    })

    output$protein_plot <- plotly::renderPlotly({
      shiny::req(input$accession)
      plot_protein_sites_interactive(prot_info = prot_info(), site_summary = site_summary(), sig_genes = sig_genes())
    })

    output$volcano_plot <- plotly::renderPlotly({
      shiny::req(input$comparison, input$accession)
      plot_volcano(accession = input$accession, res_data = res_data(), fc_thresh = input$fc_thresh, p_thresh = input$p_thresh)
    })

    output$alphafold_viewer <- NGLVieweR::renderNGLVieweR({
      shiny::req(input$accession)
      view_alphafold(accession = input$accession, site_summary = site_summary())
    })

    output$alphafold_legend <- shiny::renderUI({
      shiny::req(input$accession)
      htmltools::tags$div(style = "font-family: sans-serif; font-size: 12px; margin: 8px 12px; display: flex; align-items: center; gap: 16px;",
                          htmltools::tags$span(style = "font-weight: bold;",
                                               "Model confidence:"),
                          htmltools::tags$span(htmltools::tags$span(style = "display:inline-block; width:14px; height:14px; background:#2152CE; border-radius:2px; margin-right:4px; vertical-align:middle;"),
                                               "Very high (pLDDT > 90)"),
                          htmltools::tags$span(htmltools::tags$span(style = "display:inline-block; width:14px; height:14px; background:#7FC9EF; border-radius:2px; margin-right:4px; vertical-align:middle;"),
                                               "High (90 > pLDDT > 70)"),
                          htmltools::tags$span(htmltools::tags$span(style = "display:inline-block; width:14px; height:14px; background:#F9DC4D; border-radius:2px; margin-right:4px; vertical-align:middle;"),
                                               "Low (70 > pLDDT > 50)"),
                          htmltools::tags$span(htmltools::tags$span(style = "display:inline-block; width:14px; height:14px; background:#EE8453; border-radius:2px; margin-right:4px; vertical-align:middle;"),
                                               "Very low (pLDDT < 50)"))
    })

    output$func_annotation_card <- shiny::renderUI({
      shiny::req(input$accession)
      ann <- parse_protein_summary(prot_info = prot_info())

      shiny::div(style = "max-width: 700px;",
                 shiny::div(style = paste("background: #f8f9fa; border: 1px solid #dee2e6;",
                                          "border-radius: 8px; padding: 20px 24px;"),
                            shiny::div(style = "display: flex; justify-content: space-between; align-items: flex-start; margin-bottom: 16px;",
                                       shiny::div(shiny::h4(style = "margin: 0 0 4px 0; font-size: 16px;",
                                                            ann$protein_name),
                                                  shiny::span(style = "font-size: 13px; color: #6c757d; font-weight: 600;",
                                                              ann$gene_name)),
                                       shiny::tags$a(href = glue::glue("https://www.uniprot.org/uniprot/{ann$uniprot_id}"),
                                                     target = "_blank",
                                                     style  = "font-size: 12px; color: #0d6efd; white-space: nowrap; margin-left: 16px;",
                                                     ann$uniprot_id,
                                                     "\u2197")),
                            shiny::tags$hr(style = "margin: 0 0 16px 0; border-color: #dee2e6;"),
                            shiny::div(shiny::span(style = "font-size: 11px; font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em; color: #6c757d;",
                                                   "Function"),
                                       shiny::p(style = "font-size: 14px; line-height: 1.6; margin: 6px 0 0 0; color: #212529;",
                                                if (is.null(ann$function_text)) "No functional annotation available." else ann$function_text))))
    })

    output$protein_heatmap <- plotly::renderPlotly({
      shiny::req(input$accession)
      plot_protein_heatmap(accession = input$accession, res_data = res_data(), sample_meta = sample_meta, fc_thresh = input$fc_thresh, p_thresh = input$p_thresh)
    })
  }

  return(server)
}

run_app <- function(sample_meta, res_list, site_summaries, sig_genes_list, dataset_name) {
  app <- shiny::shinyApp(ui = app_ui(dataset_name = dataset_name),
                         server = app_server(sample_meta = sample_meta, res_list = res_list, site_summaries = site_summaries, sig_genes_list = sig_genes_list))

  return(app)
}
