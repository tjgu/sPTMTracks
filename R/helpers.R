standardize_data <- function(sample_meta, res_list, pval_col = "adj.P.Val", logfc_col = "logFC", sample_col = "sample_id_unique",
                             protein_col = "protein", gene_col = "gene", site_col = "protein_site", glycan_col = "GlycanList") {

  res_list <- lapply(res_list, function(df) {
    df <- df %>%
      dplyr::rename(adj.P.Val = dplyr::all_of(pval_col),
                    logFC = dplyr::all_of(logfc_col),
                    protein = dplyr::all_of(protein_col),
                    gene = dplyr::all_of(gene_col),
                    protein_site = dplyr::all_of(site_col),
                    GlycanList = dplyr::all_of(glycan_col)) %>%
      tibble::rowid_to_column(var = "glycoform_id") %>%
      dplyr::mutate(glycoform_id = glue::glue("glycoform_{glycoform_id}")) %>%
      tibble::column_to_rownames(var = "glycoform_id")

    return(df)
  })

  sample_meta <- sample_meta %>%
    dplyr::rename(sample = dplyr::all_of(sample_col))

  return(list(sample_meta = sample_meta,
              res_list = res_list))
}

prepare_site_summaries <- function(res_list) {
  site_summaries <- lapply(names(res_list), function(comp) {
    temp_df <- res_list[[comp]] %>%
      dplyr::select(adj.P.Val, logFC, protein, gene, protein_site)

    site_summary <- temp_df %>%
      dplyr::mutate(is_sig = adj.P.Val <= 0.2 & abs(logFC) > 1) %>%
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
      dplyr::mutate(direction = forcats::fct(direction, levels = c("up", "down", "ns"))) %>%
      dplyr::group_by(protein) %>%
      dplyr::mutate(n_sites = dplyr::n(),
                    n_hit_sites = sum(is_hit),
                    pct_hit_sites = (n_hit_sites/n_sites)*100) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(dplyr::desc(is_hit), protein) %>%
      as.data.frame()

    return(site_summary)
  })
  names(site_summaries) <- names(res_list)

  return(site_summaries)
}


prepare_sig_genes_list <- function(res_list) {
  sig_genes_list <- lapply(names(res_list), function(comp) {
    temp_df <- res_list[[comp]]

    temp_df <- temp_df %>%
      dplyr::filter(adj.P.Val <= 0.2 & abs(logFC) >= 1) %>%
      tibble::rownames_to_column(var = "glycoform_id")

    return(temp_df)
  })
  names(sig_genes_list) <- names(res_list)

  return(sig_genes_list)
}


map_ensembl_to_uniprot <- function(ids, poll_interval = 3, timeout = 600) {

  ids <- unique(trimws(as.character(ids)))
  ids <- ids[!is.na(ids) & nzchar(ids)]

  if (length(ids) == 0) {
    stop("No Ensembl protein IDs supplied.")
  }

  if (length(ids) > 100000) {
    stop("UniProt accepts a maximum of 100,000 IDs per mapping job.")
  }

  submission <- httr2::request("https://rest.uniprot.org/idmapping/run") %>%
    httr2::req_method("POST") %>%
    httr2::req_body_multipart(from = "Ensembl_Protein",
                              to = "UniProtKB",
                              ids = glue::glue_collapse(ids, sep = ",")) %>%
    httr2::req_perform() %>%
    httr2::resp_body_json()

  job_id <- submission$jobId

  if (is.null(job_id)) {
    stop("UniProt did not return a mapping job ID.")
  }

  start_time <- Sys.time()

  repeat {
    status <- httr2::request(base_url = glue::glue("https://rest.uniprot.org/idmapping/status/{job_id}")) %>%
      httr2::req_perform() %>%
      httr2::resp_body_json()

    if (!is.null(status$jobStatus)) {
      if (status$jobStatus == "FINISHED") {
        break
      }

      if (status$jobStatus == "ERROR") {
        stop("UniProt mapping job failed.")
      }

      if (!status$jobStatus %in% c("NEW", "RUNNING")) {
        stop(glue::glue("Unexpected UniProt job status: {status$jobStatus}"))
      }
    } else if (!is.null(status$results) || !is.null(status$failedIds)) {
      break
    } else {
      stop("Unexpected response from UniProt mapping status endpoint.")
    }

    if (as.numeric(difftime(Sys.time(), start_time, units = "secs")) > timeout) {
      stop("Timed out waiting for UniProt mapping job.")
    }

    Sys.sleep(poll_interval)
  }

  results <- httr2::request(base_url = glue::glue("https://rest.uniprot.org/idmapping/stream/{job_id}")) %>%
    httr2::req_url_query(format = "json") %>%
    httr2::req_perform() %>%
    httr2::resp_body_json()

  mapped_ids <- lapply(results$results, function(result) {
    tibble::tibble(ensembl_protein_id = result$from,
                   uniprot_id = as.character(result$to))
  }) %>%
    dplyr::bind_rows()

  unmapped_ids <- setdiff(ids, mapped_ids$ensembl_protein_id)

  unmapped_ids <- tibble::tibble(ensembl_protein_id = unmapped_ids,
                                 uniprot_id = NA_character_)

  mapping <- dplyr::bind_rows(mapped_ids, unmapped_ids) %>%
    dplyr::arrange(match(ensembl_protein_id, ids), uniprot_id) %>%
    dplyr::distinct(ensembl_protein_id, .keep_all = TRUE)

  return(mapping)
}
