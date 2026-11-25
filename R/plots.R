plot_volcano <- function(accession, res_data, fc_thresh, p_thresh) {
  hover_spec <- tibble::tibble(column = c("protein", "gene", "protein_site", "ModificationsNameList", "GlycanList", "logFC", "adj.P.Val"),
                               label = c("Protein", "Gene", "Site", "Modifications", "Glycans", "logFC", "adj.P.Val"),
                               expr = c("protein", "gene", "protein_site", "ModificationsNameList", "GlycanList", "round(logFC, 3)", "signif(adj.P.Val, 3)")) %>%
    dplyr::filter(column %in% names(res_data)) %>%
    dplyr::mutate(piece = glue::glue("{label}: {{{expr}}}"))
  hover_template <- glue::glue_collapse(hover_spec$piece, sep = "<br>")

  plot_data <- res_data %>%
    dplyr::select(protein, logFC, adj.P.Val, dplyr::any_of(hover_spec$column))

  plot_data <- plot_data %>%
    dplyr::mutate(neglogP = -log10(adj.P.Val),
                  fc_flag = abs(logFC) >= fc_thresh,
                  p_flag = adj.P.Val <= p_thresh,
                  group = dplyr::case_when(fc_flag & p_flag ~ "both",
                                           fc_flag & !p_flag ~ "log2FC",
                                           p_flag & !fc_flag ~ "p-value",
                                           TRUE ~ "none"),
                  hover = glue::glue(hover_template))

  cols <- c("none" = "#4C4C4C",
            "log2FC" = "#9BC495",
            "p-value" = "#4967DA",
            "both" = "#E78582")

  main_plot <- plotly::plot_ly(data = plot_data,
                               x = ~logFC, y = ~neglogP,
                               type = "scatter", mode = "markers",
                               marker = list(size = 6, opacity = 1, color = cols[plot_data$group]),
                               text = ~hover,
                               name = "all", showlegend = FALSE,
                               hovertemplate = glue::glue("<b>%{{text}}</b>",
                                                          "<br>logFC: %{{x:.3f}}",
                                                          "<br>-log10(adj.P.Val): %{{y:.3g}}",
                                                          "<extra></extra>")) %>%
    plotly::add_trace(data = plot_data %>% dplyr::filter(protein == accession),
                      x = ~logFC, y = ~neglogP,
                      type = "scatter", mode = "markers",
                      marker = list(size = 10, color = "#E41A1C",
                                    line = list(width = 1, color = "black")),
                      showlegend = FALSE) %>%
    plotly::layout(xaxis = list(title = "logFC"),
                   yaxis = list(title = "-log10(adj. P-value)"),
                   shapes = list(list(type = "line", x0 = -fc_thresh, x1 = -fc_thresh, y0 = 0, y1 = 1,
                                      xref = "x", yref = "paper",
                                      line = list(dash = "dash", color = "black")),
                                 list(type = "line", x0 = fc_thresh, x1 = fc_thresh, y0 = 0, y1 = 1,
                                      xref = "x", yref = "paper",
                                      line = list(dash = "dash", color = "black")),
                                 list(type = "line", x0 = min(plot_data$logFC, na.rm = TRUE), x1 = max(plot_data$logFC, na.rm = TRUE),
                                      y0 = -log10(p_thresh), y1 = -log10(p_thresh),
                                      xref = "x", yref = "y", line = list(dash = "dash", color = "black"))),
                   annotations = list(list(x = 0, y = -log10(p_thresh),
                                           xref = "x", yref = "y", text = "", showarrow = FALSE)))

  return(main_plot)
}

plot_protein_heatmap <- function(accession, res_data, sample_meta, fc_thresh, p_thresh) {
  sample_meta <- sample_meta %>%
    dplyr::select(sample, group_colors)

  sample_cols <- sample_meta$sample

  plot_data <- res_data %>%
    dplyr::select(protein, logFC, adj.P.Val, id, dplyr::any_of(sample_cols))

  plot_data <- plot_data %>%
    dplyr::filter(protein == accession) %>%
    dplyr::arrange(adj.P.Val) %>%
    tibble::rownames_to_column(var="glycoid") %>%
    tibble::column_to_rownames(var="id")

  if (nrow(plot_data) == 0) return(NULL)

  mat_unscaled <- plot_data %>%
    dplyr::select(dplyr::any_of(sample_cols)) %>%
    as.matrix()

  mat <- t(scale(t(mat_unscaled), center = TRUE, scale = TRUE))
  # mat <- mat_unscaled
  limit <- max(abs(mat), na.rm = TRUE)

  sig_flags <- plot_data$adj.P.Val <= p_thresh & abs(plot_data$logFC) >=fc_thresh

  group_colors <- sample_meta %>%
    dplyr::slice(match(colnames(mat), sample)) %>%
    dplyr::pull(group_colors)

  x_annotations <- lapply(seq_along(colnames(mat)), function(i) {
    list(x = i, y = min(1.002, 1 + 6 / nrow(mat)),
         text = colnames(mat)[i],
         xref = "x", yref = "paper",
         showarrow = FALSE,
         textangle = -90,
         font = list(color = group_colors[i], size = 11),
         xanchor = "center", yanchor = "bottom")
  })

  mat_unscaled_rev <- mat_unscaled[nrow(mat_unscaled):1, , drop = FALSE]
  main_plot <- heatmaply::heatmaply(mat,
                                    colors = rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")),
                                    limits = c(-limit, limit),
                                    dendrogram = "none",
                                    main = "A",
                                    show_dendrogram = c(FALSE, FALSE),
                                    grid_color = "grey60",
                                    plot_method = "plotly",
                                    height = max(400, nrow(mat) * 30)) %>%
    plotly::style(customdata = mat_unscaled_rev,
                  hovertemplate = "Row: %{y}<br>Col: %{x}<br>Intensity: %{customdata:.3f}<extra></extra>") %>%
    plotly::layout(title = list(text = glue::glue("{accession}"),
                                x = 0.25, y = 0.98,
                                xanchor = "left"),
                   autosize=TRUE,
                   margin = list(t = 100),
                   xaxis = list(showticklabels = FALSE, side = "top", automargin=TRUE),
                   yaxis = list(ticktext = dplyr::if_else(sig_flags, as.character(glue::glue("<b>* {rownames(plot_data)}</b>")), rownames(plot_data)),
                                tickvals = rev(seq_len(nrow(mat))),
                                tickmode = "array"),
                   annotations = x_annotations,
                   coloraxis = list(colorbar = list(title = "Z-score")))

  return(main_plot)
}

parse_protein_summary <- function(prot_info) {
  return(list(protein_name = prot_info$proteinDescription$recommendedName$fullName$value,
              gene_name = prot_info$genes[[1]]$geneName$value,
              uniprot_id = prot_info$primaryAccession,
              function_text = prot_info$comments %>%
                purrr::keep(\(c) identical(c$commentType, "FUNCTION")) %>%
                purrr::pluck(1, "texts", 1, "value", .default = NA_character_)))
}

view_alphafold <- function(accession, site_summary) {
  site_summary <- site_summary %>%
    dplyr::select(protein, protein_site, is_hit, mean_logFC, direction)

  meta <- httr2::request(glue::glue("https://alphafold.ebi.ac.uk/api/prediction/{accession}")) %>%
    httr2::req_perform() %>%
    httr2::resp_body_json()

  pdb_url <- meta[[1]]$pdbUrl
  sequence   <- meta[[1]]$uniprotSequence

  pdb_string <- httr2::request(pdb_url) %>%
    httr2::req_user_agent("sPTMTracks/0.1.0") %>%
    httr2::req_perform() %>%
    httr2::resp_body_string()

  plot_data <- site_summary %>%
    dplyr::filter(protein == accession) %>%
    dplyr::mutate(protein_site = as.integer(protein_site),
                  aa = stringr::str_sub(sequence, protein_site, protein_site),
                  site_label = glue::glue("{aa}{protein_site}"),
                  mean_logFC = dplyr::if_else(is_hit, mean_logFC, 0)) %>%
    dplyr::mutate(resno = protein_site,
                  color = dplyr::case_when(direction == "up" ~ "#d73027",
                                           direction == "down" ~ "#4575b4",
                                           TRUE ~ "#b0b0b0"))

  main_plot <- NGLVieweR::NGLVieweR(data = pdb_string, format = "pdb") %>%
    NGLVieweR::stageParameters(backgroundColor = "white") %>%
    NGLVieweR::setQuality(quality = "high") %>%
    NGLVieweR::setSpin(spin = FALSE) %>%
    NGLVieweR::addRepresentation(type = "cartoon",
                                 param = list(name = "cartoon", colorScheme = "bfactor", opacity = 1))

  # add a ball+stick sphere at each glycosite Cα
  for (i in seq_len(nrow(plot_data))) {
    main_plot <- main_plot %>%
      NGLVieweR::addRepresentation(type = "spacefill",
                                   param = list(name = glue::glue("site_{plot_data$protein_site[i]}"),
                                                sele = glue::glue("{plot_data$protein_site[i]} and .CA"),
                                                colorValue = plot_data$color[i],
                                                colorScheme = "uniform",
                                                radius = 1.5))
  }

  return(main_plot)
}

plot_protein_sites_interactive <- function(prot_info, site_summary, sig_genes) {
  sig_genes <- sig_genes %>%
    dplyr::select(protein, protein_site, logFC, adj.P.Val, GlycanList)

  site_summary <- site_summary %>%
    dplyr::select(protein, protein_site, is_hit, mean_logFC, n_sig_glycoforms, n_glycoforms, direction)

  # fetch from UniProt
  accession <- prot_info$primaryAccession
  sequence   <- prot_info$sequence$value
  seq_length <- prot_info$sequence$length
  gene_name  <- prot_info$genes[[1]]$geneName$value
  prot_name  <- prot_info$proteinDescription$recommendedName$fullName$value

  # parse features
  features_df <- purrr::map_dfr(prot_info$features, ~ data.frame(type = .x$type,
                                                                 description = .x$description,
                                                                 start = .x$location$start$value,
                                                                 end = .x$location$end$value))

  track1_types <- c("Signal", "Transmembrane", "Chain")
  track2_types <- c("Glycosylation", "Disulfide bond", "Binding site", "Motif", "Modified residue")
  track1_colors <- c("Signal" = "#b2abd2", "Chain" = "#aec7e8",
                     "Transmembrane" = "#d6604d")
  track2_colors <- c("Glycosylation" = "#4575b4", "Disulfide bond" = "#d73027", "Modified residue" = "green",
                     "Binding site" = "#f46d43", "Motif" = "#74add1")

  track1 <- features_df %>% dplyr::filter(type %in% track1_types)
  track2_other <- features_df %>% dplyr::filter(type %in% track2_types & type != "Disulfide bond")
  track2_disulfide <- features_df %>% dplyr::filter(type == "Disulfide bond")

  top_glycoforms <- sig_genes %>%
    dplyr::filter(protein == accession) %>%
    dplyr::mutate(protein_site = as.integer(protein_site)) %>%
    dplyr::group_by(protein_site) %>%
    dplyr::slice_min(adj.P.Val, n = 2, with_ties = FALSE) %>%
    dplyr::mutate(glycoform_label = glue::glue("{GlycanList} (logFC={round(logFC, 2)}, adj.P={round(adj.P.Val, 3)})")) %>%
    dplyr::summarise(top_glycoforms_text = glue::glue_collapse(glycoform_label, sep = "<br>") %>% as.character(), .groups = "drop")

  # plot data
  plot_data <- site_summary %>%
    dplyr::filter(protein == accession) %>%
    dplyr::mutate(protein_site = as.integer(protein_site),
                  aa = stringr::str_sub(sequence, protein_site, protein_site),
                  site_label = glue::glue("{aa}{protein_site}"),
                  mean_logFC = dplyr::if_else(is_hit, mean_logFC, 0)) %>%
    dplyr::left_join(top_glycoforms, by = "protein_site") %>%
    dplyr::mutate(top_glycoforms_text = tidyr::replace_na(top_glycoforms_text, "None")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(overlapping_features = features_df %>%
                    dplyr::filter(start <= protein_site & end >= protein_site) %>%
                    dplyr::mutate(label = dplyr::if_else(description == "", type, as.character(glue::glue("{type}: {description}")))) %>%
                    dplyr::pull(label) %>%
                    glue::glue_collapse(sep = "<br>") %>%
                    as.character()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(hover_text = glue::glue("<b>{site_label}</b><br>",
                                          "Mean logFC: {round(mean_logFC, 3)}<br>",
                                          "Sig glycoforms: {n_sig_glycoforms} / {n_glycoforms}<br>",
                                          "Direction: {direction}<br>",
                                          "<br><b>Top glycoforms:</b><br>{top_glycoforms_text}<br>",
                                          "<br><b>Features:</b><br>{overlapping_features}"))

  direction_colors <- c("#d73027", "#4575b4", "grey70")

  # --- plot 1: lollipop ---
  p1 <- plotly::plot_ly() %>%
    plotly::add_segments(x = 1, xend = seq_length, y = 0, yend = 0,
                         line = list(color = "grey40", width = 3),
                         showlegend = FALSE, hoverinfo = "none") %>%
    plotly::add_segments(data = plot_data,
                         x = ~protein_site, xend = ~protein_site,
                         y = 0, yend = ~mean_logFC,
                         line = list(color = "grey60", width = 1, dash = "dash"),
                         showlegend = FALSE, hoverinfo = "none") %>%
    plotly::add_markers(data = plot_data,
                        x = ~protein_site, y = ~mean_logFC,
                        color = ~direction, colors = direction_colors,
                        marker = list(size = 12, line = list(color = "white", width = 1.5)),
                        text = ~hover_text, hoverinfo = "text", name = ~direction) %>%
    plotly::add_text(data = plot_data,
                     x = ~protein_site, y = ~mean_logFC,
                     text = ~site_label,
                     textposition = ~ifelse(mean_logFC >= 0, "top center", "bottom center"),
                     textfont = list(size = 10, color = "black"),
                     showlegend = FALSE, hoverinfo = "none") %>%
    plotly::layout(yaxis = list(title = "Mean logFC", zeroline = TRUE),
                   xaxis = list(title = "", showticklabels = FALSE),
                   legend = list(title = list(text = "Direction")),
                   annotations = list(list(x = seq_length,
                                           y = 0,
                                           text = glue::glue("{seq_length} aa"),
                                           xanchor = "right",
                                           yanchor = "bottom",
                                           showarrow = FALSE,
                                           font = list(size = 11, color = "grey40"))))

  # --- plot 2: topology track ---
  p2 <- plotly::plot_ly()

  # backbone
  p2 <- p2 %>%
    plotly::add_segments(x = 1, xend = seq_length, y = 0.5, yend = 0.5,
                         line = list(color = "grey70", width = 2),
                         showlegend = FALSE, hoverinfo = "none")

  # one trace per feature type for legend
  for (ft in unique(track1$type)) {
    ft_data <- track1 %>% dplyr::filter(type == ft)
    for (i in seq_len(nrow(ft_data))) {
      row <- ft_data[i, ]
      hover <- ifelse(row$description == "", row$type,
                      glue::glue("{row$type}: {row$description}"))
      p2 <- p2 %>%
        plotly::add_polygons(x = c(row$start, row$end, row$end, row$start, row$start),
                             y = c(0.2, 0.2, 0.8, 0.8, 0.2),
                             color = I(track1_colors[[ft]]),
                             line = list(color = "black", width = 1),
                             text = hover,
                             hoverinfo = "text",
                             hoveron = "fills",
                             name = ft, legendgroup = ft, showlegend = (i == 1))
    }
  }

  p2 <- p2 %>%
    plotly::layout(yaxis = list(title = "Topology", showticklabels = FALSE,
                                range = c(0, 1), zeroline = FALSE),
                   xaxis = list(title = "", showticklabels = FALSE),
                   legend = list(title = list(text = "Feature")))

  # --- plot 3: site annotations track ---
  p3 <- plotly::plot_ly() %>%
    plotly::add_segments(x = 1, xend = seq_length, y = 0.5, yend = 0.5,
                         line = list(color = "grey70", width = 2),
                         showlegend = FALSE, hoverinfo = "none")

  # other point features
  if (nrow(track2_other) > 0) {
    for (ft in unique(track2_other$type)) {
      ft_data <- track2_other %>% dplyr::filter(type == ft)
      p3 <- p3 %>%
        plotly::add_markers(data = ft_data,
                            x = ~((start + end) / 2), y = 0.5,
                            marker = list(color = track2_colors[[ft]], size = 10, symbol = "diamond"),
                            name = ft,
                            legendgroup = ft,
                            text = ~ifelse(description == "", type, glue::glue("{type}: {description}")),
                            hoverinfo = "text")
    }
  }

  # disulfide bond arcs as bezier via many segments
  if (nrow(track2_disulfide) > 0) {
    first_disulfide <- TRUE
    for (i in seq_len(nrow(track2_disulfide))) {
      row <- track2_disulfide[i, ]
      # generate arc points
      t_seq <- seq(0, 1, length.out = 50)
      arc_x <- row$start + (row$end - row$start) * t_seq
      arc_y <- 0.5 + ((row$end - row$start) / (seq_length * 0.2)) * sin(pi * t_seq)
      hover_txt <- glue::glue("Disulfide bond: {row$start}-{row$end}")
      p3 <- p3 %>%
        plotly::add_lines(x = arc_x, y = arc_y,
                          line = list(color = "#d73027", width = 2),
                          name = "Disulfide bond",
                          legendgroup = "Disulfide bond", showlegend = first_disulfide,
                          text = hover_txt, hoverinfo = "text")
      first_disulfide <- FALSE
    }
  }

  p3 <- p3 %>%
    plotly::layout(yaxis = list(title = "Annotations", showticklabels = FALSE,
                                range = c(0, 2), zeroline = FALSE),
                   xaxis = list(title = "Sequence position"),
                   legend = list(title = list(text = "Annotation")))

  # combine with subplot
  plotly::subplot(p1, p2, p3,
                  nrows = 3,
                  shareX = TRUE,
                  heights = c(0.6, 0.2, 0.2),
                  titleY = TRUE) %>%
    plotly::layout(title = glue::glue("{accession} — {gene_name} — {prot_name}"),
                   legend = list(tracegroupgap = 20))
}
