library(httr2)
library(jsonlite)
library(dplyr)
library(purrr)

# Helper for NULL handling
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# 1. Get target ChEMBL ID from gene symbol
get_target_id <- function(gene_symbol) {
  url <- paste0(
    "https://www.ebi.ac.uk/chembl/api/data/target/search?q=",
    gene_symbol, "&format=json"
  )
  
  
  resp <- request(url) |> req_perform() |> resp_body_json()
  
  human_targets <- Filter(function(x) {
    x$organism == "Homo sapiens" && x$target_type == "SINGLE PROTEIN"
  }, resp$targets)
  
  if (length(human_targets) == 0) return(NULL)
  
  list(
    target_chembl_id = human_targets[[1]]$target_chembl_id,
    pref_name = human_targets[[1]]$pref_name
  )
}

# 2. Get drugs targeting a ChEMBL target
get_drugs_for_target <- function(target_chembl_id) {
  url <- paste0(
    "https://www.ebi.ac.uk/chembl/api/data/mechanism?",
    "target_chembl_id=", target_chembl_id,
    "&format=json"
  )
  
  resp <- request(url) |> req_perform() |> resp_body_json()
  
  if (length(resp$mechanisms) == 0) return(NULL)
  
  map_df(resp$mechanisms, function(m) {
    data.frame(
      molecule_chembl_id = m$molecule_chembl_id %||% NA,
      drug_name = m$molecule_name %||% NA,
      mechanism = m$mechanism_of_action %||% NA,
      action_type = m$action_type %||% NA,
      max_phase = m$max_phase %||% NA,
      stringsAsFactors = FALSE
    )
    
  }) |> distinct()
}

# 3. Get indications for a molecule
get_indications <- function(molecule_chembl_id) {
  url <- paste0(
    "https://www.ebi.ac.uk/chembl/api/data/drug_indication?",
    "molecule_chembl_id=", molecule_chembl_id,
    "&format=json"
  )
  
  resp <- tryCatch(
    request(url) |> req_perform() |> resp_body_json(),
    error = function(e) list(drug_indications = list())
  )
  
  if (length(resp$drug_indications) == 0) {
    return(data.frame(
      molecule_chembl_id = molecule_chembl_id,
      indication = NA,
      mesh_id = NA,
      efo_id = NA,
      max_phase_for_ind = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  map_df(resp$drug_indications, function(ind) {
    data.frame(
      molecule_chembl_id = molecule_chembl_id,
      indication = ind$mesh_heading %||% ind$efo_term %||% NA,
      mesh_id = ind$mesh_id %||% NA,
      efo_id = ind$efo_id %||% NA,
      max_phase_for_ind = ind$max_phase_for_ind %||% NA,
      stringsAsFactors = FALSE
    )
  })
}

# 4. Main function: gene -> drugs -> indications
get_drug_target_indications <- function(gene_symbol, sleep_time = 0.3) {
  
  message("Processing: ", gene_symbol)
  
  # Get target
  target <- get_target_id(gene_symbol)
  if (is.null(target)) {
    message("  No target found")
    return(NULL)
  }
  
  Sys.sleep(sleep_time)
  
  # Get drugs
  
  drugs <- get_drugs_for_target(target$target_chembl_id)
  if (is.null(drugs) || nrow(drugs) == 0) {
    message("  No drugs found")
    return(NULL)
  }
  
  message("  Found ", nrow(drugs), " drugs")
  
  # Get indications for each drug
  indications <- map_df(unique(drugs$molecule_chembl_id), function(mol_id) {
    Sys.sleep(sleep_time)
    get_indications(mol_id)
  })
  
  # Combine everything
  result <- drugs |>
    left_join(indications, by = "molecule_chembl_id", relationship = "many-to-many") |>
    mutate(
      gene = gene_symbol,
      target_chembl_id = target$target_chembl_id,
      target_name = target$pref_name
    ) |>
    select(gene, target_name, target_chembl_id, everything())
  
  message("  Found ", sum(!is.na(result$indication)), " indications")
  result
}




# Null coalescing operator: returns 'b' if 'a' is NULL or empty, otherwise 'a'
# Useful for handling missing fields in API responses
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

get_opentargets_comprehensive <- function(gene_symbol) {
  
  message("Processing: ", gene_symbol)
  
  
  # ---- Step 1: Resolve gene symbol to internal ID ----
  # Open Targets requires Ensembl IDs internally, but we can
  # use their search endpoint to convert gene symbols automatically
  search_query <- sprintf('
  {
    search(queryString: "%s", entityNames: ["target"]) {
      hits {
        id
      }
    }
  }', gene_symbol)
  
  search_resp <- request("https://api.platform.opentargets.org/api/v4/graphql") |>
    req_body_json(list(query = search_query)) |>
    req_perform() |>
    resp_body_json()
  
  hits <- search_resp$data$search$hits
  
  # Exit early if gene not found
  if (length(hits) == 0) {
    message("  Gene not found")
    return(NULL)
  }
  
  # Take first hit (best match)
  target_id <- hits[[1]]$id
  
  # ---- Step 2: Query drugs and indications for this target ----
  # GraphQL query to get:
  #   - drug name, type, max clinical phase
  #   - mechanism of action
  #   - disease indications and their trial phases
  query <- sprintf('
  {
    target(ensemblId: "%s") {
      knownDrugs {
        uniqueDrugs
        rows {
          drug {
            id
            name
            drugType
            maximumClinicalTrialPhase
          }
          mechanismOfAction
          disease {
            id
            name
          }
          phase
        }
      }
    }
  }', target_id)
  
  resp <- request("https://api.platform.opentargets.org/api/v4/graphql") |>
    req_body_json(list(query = query)) |>
    req_perform() |>
    resp_body_json()
  
  rows <- resp$data$target$knownDrugs$rows
  
  # Exit if no drugs found for this target
  if (length(rows) == 0) {
    message("  No drugs found")
    result <- data.frame(
      gene = gene_symbol,
      drug_name = NA,
      drug_type = NA,
      max_phase = NA,
      mechanism = NA,
      indication = NA,
      phase_for_indication = NA,
      stringsAsFactors = FALSE
    )
    return(result)
  }
  
  # ---- Step 3: Parse response into a data frame ----
  # Each row is a drug-indication pair
  # %||% NA handles missing fields gracefully
  result <- map_df(rows, function(r) {
    data.frame(
      gene = gene_symbol,
      drug_name = r$drug$name %||% NA,
      drug_type = r$drug$drugType %||% NA,
      max_phase = r$drug$maximumClinicalTrialPhase %||% NA,
      mechanism = r$mechanismOfAction %||% NA,
      indication = r$disease$name %||% NA,
      phase_for_indication = r$phase %||% NA,
      stringsAsFactors = FALSE
    )
  })
  
  message("  Found ", n_distinct(result$drug_name), " drugs")
  result
}



#' Get genes with genetic associations from Open Targets
#'
#' @param efo_codes Character vector of EFO disease codes (e.g., "EFO_0000305")
#' @param score_threshold Minimum association score (0-1), default 0
#' @param size Maximum number of associations to return per disease, default 500
#' @param genetic_only If TRUE, filter for genetic_association datatype only
#' @param verbose Print progress messages, default TRUE
#'
#' @return A data.frame with gene-disease associations
#'
#' @examples
#' results <- get_opentargets_genetic_associations("EFO_0000384")

get_opentargets_genetic_associations <- function(efo_codes, 
                                                 score_threshold = 0, 
                                                 size = 500,
                                                 genetic_only = TRUE,
                                                 verbose = TRUE) {
  
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required. Install with: install.packages('httr')")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required. Install with: install.packages('jsonlite')")
  }
  
  api_url <- "https://api.platform.opentargets.org/api/v4/graphql"
  
  # Simple query without datasource filtering (works reliably)
  query <- '
  query DiseaseAssociations($efoId: String!, $size: Int!, $index: Int!) {
    disease(efoId: $efoId) {
      id
      name
      associatedTargets(page: {size: $size, index: $index}) {
        count
        rows {
          target {
            id
            approvedSymbol
            approvedName
            biotype
          }
          score
          datatypeScores {
            id
            score
          }
        }
      }
    }
  }
  '
  
  all_results <- list()
  
  for (i in seq_along(efo_codes)) {
    efo_id <- efo_codes[i]
    
    if (verbose) {
      message(sprintf("[%d/%d] Querying: %s", i, length(efo_codes), efo_id))
    }
    
    # Prepare request body
    body <- list(
      query = query,
      variables = list(
        efoId = efo_id,
        size = as.integer(size),
        index = 0L
      )
    )
    
    response <- tryCatch({
      httr::POST(
        url = api_url,
        httr::content_type_json(),
        httr::accept_json(),
        body = jsonlite::toJSON(body, auto_unbox = TRUE),
        encode = "raw"
      )
    }, error = function(e) {
      warning(sprintf("Request failed for %s: %s", efo_id, e$message))
      return(NULL)
    })
    
    if (is.null(response)) next
    
    if (httr::status_code(response) != 200) {
      content <- httr::content(response, as = "text", encoding = "UTF-8")
      warning(sprintf("API returned status %d for %s: %s", 
                      httr::status_code(response), efo_id, content))
      next
    }
    
    content <- httr::content(response, as = "text", encoding = "UTF-8")
    result <- jsonlite::fromJSON(content, flatten = TRUE)
    
    if (!is.null(result$errors)) {
      warning(sprintf("GraphQL error for %s: %s", 
                      efo_id, result$errors$message[1]))
      next
    }
    
    disease_data <- result$data$disease
    
    if (is.null(disease_data)) {
      warning(sprintf("No data found for %s", efo_id))
      next
    }
    
    associations <- disease_data$associatedTargets$rows
    
    if (is.null(associations) || length(associations) == 0 || nrow(associations) == 0) {
      if (verbose) message(sprintf("  No associations found for %s", efo_id))
      next
    }
    
    # Extract genetic association scores from datatypeScores
    genetic_scores <- sapply(associations$datatypeScores, function(x) {
      if (is.null(x) || nrow(x) == 0) return(NA_real_)
      genetic_idx <- which(x$id == "genetic_association")
      if (length(genetic_idx) > 0) x$score[genetic_idx[1]] else NA_real_
    })
    
    df <- data.frame(
      disease_id = disease_data$id,
      disease_name = disease_data$name,
      gene_id = associations$target.id,
      symbol = associations$target.approvedSymbol,
      gene_name = associations$target.approvedName,
      gene_biotype = associations$target.biotype,
      overall_score = associations$score,
      genetic_association_score = genetic_scores,
      stringsAsFactors = FALSE
    )
    
    # Filter for genetic associations only if requested
    if (genetic_only) {
      df <- df[!is.na(df$genetic_association_score) & df$genetic_association_score > 0, ]
    }
    
    # Apply score threshold
    df <- df[df$overall_score >= score_threshold, ]
    
    if (nrow(df) > 0) {
      all_results[[efo_id]] <- df
      if (verbose) {
        message(sprintf("  Found %d associations (%d with genetic evidence)", 
                        nrow(df), sum(!is.na(df$genetic_association_score))))
      }
    }
    
    Sys.sleep(0.3)
  }
  
  if (length(all_results) == 0) {
    warning("No associations found for any of the provided EFO codes")
    return(data.frame())
  }
  
  final_df <- do.call(rbind, all_results)
  rownames(final_df) <- NULL
  
  # Sort by genetic association score
  final_df <- final_df[order(-final_df$genetic_association_score, na.last = TRUE), ]
  
  if (verbose) {
    message(sprintf("\nTotal: %d gene-disease associations across %d diseases", 
                    nrow(final_df), length(unique(final_df$disease_id))))
  }
  
  return(final_df)
}


#' Get detailed genetic evidence from Open Targets
#' 
#' Returns individual evidence records with variant-level information
#'
#' @param efo_code Single EFO disease code
#' @param ensembl_ids Optional vector of Ensembl gene IDs to filter
#' @param size Maximum number of evidence records, default 100
#'
#' @return A data.frame with detailed genetic evidence

get_opentargets_genetic_evidence <- function(efo_code,
                                             ensembl_ids = NULL,
                                             size = 100,
                                             verbose = TRUE) {
  
  if (!requireNamespace("httr", quietly = TRUE)) stop("Package 'httr' is required.")
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Package 'jsonlite' is required.")
  
  api_url <- "https://api.platform.opentargets.org/api/v4/graphql"
  
  # If no specific genes provided, first get associated genes
  if (is.null(ensembl_ids)) {
    if (verbose) message("Getting associated genes first...")
    assoc <- get_opentargets_genetic_associations(efo_code, genetic_only = TRUE, verbose = FALSE)
    if (nrow(assoc) == 0) {
      warning("No genetic associations found")
      return(data.frame())
    }
    ensembl_ids <- head(assoc$gene_id, 50)  # Limit to top 50 genes
  }
  
  # Query for genetic evidence
  query <- '
  query GeneticEvidence($efoId: String!, $ensemblIds: [String!]!, $size: Int!) {
    disease(efoId: $efoId) {
      id
      name
      evidences(
        ensemblIds: $ensemblIds
        datasourceIds: ["ot_genetics_portal", "gene_burden", "eva", "eva_somatic", 
                        "gene2phenotype", "orphanet", "genomics_england", "clingen",
                        "uniprot_literature", "uniprot_variants"]
        size: $size
      ) {
        count
        rows {
          id
          score
          datasourceId
          datatypeId
          variantRsId
          studyId
          publicationYear
          resourceScore
          oddsRatio
          beta
          pValueMantissa
          pValueExponent
          clinicalSignificances
          alleleOrigins
          target {
            id
            approvedSymbol
          }
          disease {
            id
            name
          }
          variantFunctionalConsequence {
            id
            label
          }
        }
      }
    }
  }
  '
  
  body <- list(
    query = query,
    variables = list(
      efoId = efo_code,
      ensemblIds = as.list(ensembl_ids),
      size = as.integer(size)
    )
  )
  
  if (verbose) message(sprintf("Querying genetic evidence for %s...", efo_code))
  
  response <- httr::POST(
    url = api_url,
    httr::content_type_json(),
    httr::accept_json(),
    body = jsonlite::toJSON(body, auto_unbox = TRUE),
    encode = "raw"
  )
  
  if (httr::status_code(response) != 200) {
    content <- httr::content(response, as = "text", encoding = "UTF-8")
    stop(sprintf("API error %d: %s", httr::status_code(response), content))
  }
  
  content <- httr::content(response, as = "text", encoding = "UTF-8")
  result <- jsonlite::fromJSON(content, flatten = TRUE)
  
  if (!is.null(result$errors)) {
    stop(sprintf("GraphQL error: %s", result$errors$message[1]))
  }
  
  evidences <- result$data$disease$evidences$rows
  
  if (is.null(evidences) || nrow(evidences) == 0) {
    warning("No genetic evidence found")
    return(data.frame())
  }
  
  if (verbose) message(sprintf("Found %d evidence records", nrow(evidences)))
  
  return(evidences)
}




get_opentargets_overall_association <- function(efo_codes, 
                                                 score_threshold = 0, 
                                                 page_size = 500,
                                                 max_targets = Inf,
                                                 verbose = TRUE) {
  
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required. Install with: install.packages('httr')")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required. Install with: install.packages('jsonlite')")
  }
  
  api_url <- "https://api.platform.opentargets.org/api/v4/graphql"
  
  query <- '
  query DiseaseAssociations($efoId: String!, $size: Int!, $index: Int!) {
    disease(efoId: $efoId) {
      id
      name
      associatedTargets(page: {size: $size, index: $index}) {
        count
        rows {
          target {
            id
            approvedSymbol
            approvedName
            biotype
          }
          score
          datatypeScores {
            id
            score
          }
          datasourceScores {
            id
            score
          }
        }
      }
    }
  }
  '
  
  # Helper to extract score by ID
  extract_score <- function(score_list, score_id) {
    sapply(score_list, function(x) {
      if (is.null(x) || nrow(x) == 0) return(NA_real_)
      idx <- which(x$id == score_id)
      if (length(idx) > 0) x$score[idx[1]] else NA_real_
    })
  }
  
  # Helper to fetch a single page
  fetch_page <- function(efo_id, page_index, page_size) {
    body <- list(
      query = query,
      variables = list(
        efoId = efo_id,
        size = as.integer(page_size),
        index = as.integer(page_index)
      )
    )
    
    response <- httr::POST(
      url = api_url,
      httr::content_type_json(),
      httr::accept_json(),
      body = jsonlite::toJSON(body, auto_unbox = TRUE),
      encode = "raw"
    )
    
    if (httr::status_code(response) != 200) {
      return(NULL)
    }
    
    content <- httr::content(response, as = "text", encoding = "UTF-8")
    jsonlite::fromJSON(content, flatten = TRUE)
  }
  
  all_results <- list()
  
  for (i in seq_along(efo_codes)) {
    efo_id <- efo_codes[i]
    
    if (verbose) {
      message(sprintf("[%d/%d] Querying: %s", i, length(efo_codes), efo_id))
    }
    
    # Fetch first page to get total count
    result <- tryCatch({
      fetch_page(efo_id, page_index = 0, page_size = page_size)
    }, error = function(e) {
      warning(sprintf("Request failed for %s: %s", efo_id, e$message))
      return(NULL)
    })
    
    if (is.null(result) || !is.null(result$errors) || is.null(result$data$disease)) {
      warning(sprintf("No data found for %s", efo_id))
      next
    }
    
    disease_data <- result$data$disease
    total_count <- disease_data$associatedTargets$count
    
    if (total_count == 0) {
      if (verbose) message(sprintf("  No associations found"))
      next
    }
    
    # Calculate number of pages needed
    targets_to_fetch <- min(total_count, max_targets)
    n_pages <- ceiling(targets_to_fetch / page_size)
    
    if (verbose) {
      message(sprintf("  Total associations: %d (fetching %d across %d pages)", 
                      total_count, targets_to_fetch, n_pages))
    }
    
    # Collect all pages
    disease_results <- list()
    
    for (page_idx in 0:(n_pages - 1)) {
      if (page_idx > 0) {
        # Already have first page, fetch subsequent pages
        Sys.sleep(0.2)  # Rate limiting
        result <- tryCatch({
          fetch_page(efo_id, page_index = page_idx, page_size = page_size)
        }, error = function(e) {
          warning(sprintf("Failed to fetch page %d for %s", page_idx, efo_id))
          return(NULL)
        })
        
        if (is.null(result) || is.null(result$data$disease)) {
          warning(sprintf("  Page %d returned no data", page_idx))
          next
        }
      }
      
      associations <- result$data$disease$associatedTargets$rows
      
      if (is.null(associations) || nrow(associations) == 0) {
        next
      }
      
      df <- data.frame(
        disease_id = disease_data$id,
        disease_name = disease_data$name,
        gene_id = associations$target.id,
        gene_symbol = associations$target.approvedSymbol,
        gene_name = associations$target.approvedName,
        gene_biotype = associations$target.biotype,
        overall_score = associations$score,
        genetic_association_score = extract_score(associations$datatypeScores, "genetic_association"),
        ot_genetics_portal_score = extract_score(associations$datasourceScores, "ot_genetics_portal"),
        gene_burden_score = extract_score(associations$datasourceScores, "gene_burden"),
        eva_score = extract_score(associations$datasourceScores, "eva"),
        stringsAsFactors = FALSE
      )
      
      disease_results[[page_idx + 1]] <- df
      
      if (verbose && n_pages > 1) {
        message(sprintf("    Page %d/%d: %d genes", page_idx + 1, n_pages, nrow(df)))
      }
    }
    
    # Combine all pages for this disease
    if (length(disease_results) > 0) {
      disease_df <- do.call(rbind, disease_results)
      
      # Apply filtering: genetic_association > 0 OR overall_score >= threshold
      has_genetic <- !is.na(disease_df$genetic_association_score) & 
        disease_df$genetic_association_score > 0
      meets_threshold <- disease_df$overall_score >= score_threshold
      
      disease_df <- disease_df[has_genetic | meets_threshold, ]
      
      if (nrow(disease_df) > 0) {
        all_results[[efo_id]] <- disease_df
        if (verbose) {
          n_genetic <- sum(!is.na(disease_df$genetic_association_score) & 
                             disease_df$genetic_association_score > 0)
          message(sprintf("  Retained: %d genes (%d with genetic evidence)", 
                          nrow(disease_df), n_genetic))
        }
      }
    }
    
    # Rate limiting between diseases
    Sys.sleep(0.3)
  }
  
  if (length(all_results) == 0) {
    warning("No associations found for any of the provided EFO codes")
    return(data.frame())
  }
  
  final_df <- do.call(rbind, all_results)
  rownames(final_df) <- NULL
  
  # Sort by genetic score, then overall score
  final_df <- final_df[order(-final_df$genetic_association_score, 
                             -final_df$overall_score, 
                             na.last = TRUE), ]
  
  if (verbose) {
    message(sprintf("\n=== Summary ==="))
    message(sprintf("Total: %d gene-disease associations", nrow(final_df)))
    message(sprintf("Unique genes: %d", length(unique(final_df$gene_id))))
    message(sprintf("Diseases: %d", length(unique(final_df$disease_id))))
  }
  
  return(final_df)
}
