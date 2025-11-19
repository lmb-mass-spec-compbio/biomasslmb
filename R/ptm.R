#' Filter PTM data from MaxQuant to retain those with a given PTM
#'
#' @description This function takes a `Qfeatures` object with the evidence.txt
#' output from MaxQuant and filters it to retain peptides with a given PTM.
#'
#' @param obj `QFeatures`. Proteomics dataset
#' @param i `string`. Index for the SummarizedExperiment with data without imputation
#' @param ptm_prob_col `character` column name for PTMs.
#' @return Returns a `SummarizedExperiment` where the value for the PTM probability column is not empty
#' @export
filter_maxquant_ptm <- function(obj, i, ptm_prob_col='Phospho..STY..Probabilities'){
  return(obj[[i]][rowData(obj[[i]])[[ptm_prob_col]]!='',])
}

#' Add rowData columns with positions of PTMs with respect to peptide sequence
#'
#' @description This function takes a `SummarizedExperiment` object with the evidence.txt
#' output from MaxQuant and adds rowData columns to describe the PTMs present in the peptide and their positions within the peptide' '
#' @param obj `SummarizedExperiment`. Proteomics dataset
#' @param ptms_to_retain `character` vector with PTMs to retain. Inspect values in rowData(obj)$Modified.sequence to get correct PTM name
#' @param ptm_encoding_pos `character` name vector describing whether modifification comes before (-1) or after (1) amino acid in Modified.sequence column values.
#' @return Returns a `SummarizedExperiment` with an additional column in the RowData describing the position of the PTMs with respect to the peptide sequence
#' @export
add_ptm_pos_rowdata_mq <- function(obj,
                                   ptms_to_retain = c('p', '(ph)', '(ac)p'),
                                   ptm_encoding_pos = c('p'=-1, '(ph)'=1, '(ac)p'=-1, '(ox)'=1),
                                   verbose=TRUE){


  check_se(obj)

  ptms_without_defined_encoding_pos <- setdiff(ptms_to_retain, names(ptm_encoding_pos))

  if(length(ptms_without_defined_encoding_pos)>0){
    stop(sprintf('Need to define position of modification with respect to amino acid in
  sequence with argument ptm_encoding_pos. Not defined for PTM: %s', paste(ptms_without_defined_encoding_pos, collapse='; ')))
  }
  if(any(!ptm_encoding_pos%in% c(-1,1))){
    stop(sprintf('values for ptm_encoding_pos must be 1 or -1, values are: %s', paste(unname(ptm_encoding_pos), collapse=', ')))
  }

  seqs_split <- strsplit(rowData(obj)$Sequence, '')
  mod_seqs_split <- strsplit(gsub('_$|^_', '', rowData(obj)$Modified.sequence), '')


  ptms_detected <- list()
  ptm_pos <- lapply(1:length(seqs_split), function(seq_ix){
    seq_split <- seqs_split[[seq_ix]]
    mod_seq_split <- mod_seqs_split[[seq_ix]]

    ix2_nudge <- 0
    ptm_positions <- c()
    ptms <- c()
    ptm_amino_acids <- c()

    for(ix in 1:length(seq_split)){

      if(seq_split[[ix]]!=mod_seq_split[[ix+ix2_nudge]]){

        mod_seq <- ''
        while(seq_split[[ix]]!=mod_seq_split[[ix+ix2_nudge]]){
          mod_seq <- paste0(mod_seq, mod_seq_split[[ix+ix2_nudge]])

          ix2_nudge <- ix2_nudge + 1

        }
        if(mod_seq %in% names(ptms_detected)){
          ptms_detected[[mod_seq]] <<- ptms_detected[[mod_seq]] + 1
        } else{
          ptms_detected[[mod_seq]] <<- 1
        }

        if(mod_seq %in% ptms_to_retain){
          ptms <- c(ptms, mod_seq)

          if(ptm_encoding_pos[[mod_seq]] == 1){
            ptm_positions <- c(ptm_positions, (ix-1))
            ptm_amino_acids <- c(ptm_amino_acids, seq_split[[ix-1]])
          }
          else{
            ptm_positions <- c(ptm_positions, (ix))
            ptm_amino_acids <- c(ptm_amino_acids, seq_split[[ix]])
          }

        }

      }

    }
    ptms <- paste(ptms, collapse = '; ')
    ptm_positions <- paste(ptm_positions, collapse = '; ')
    ptm_amino_acids <- paste(ptm_amino_acids, collapse = '; ')

    return(c('ptms'=ptms, 'ptm_positions'=ptm_positions, 'ptm_amino_acids'=ptm_amino_acids))
  }) %>% bind_rows()


  if(verbose){
    ptms_detected <- ptms_detected[order(-as.numeric(ptms_detected))]
    ptms_detected_str <- paste(sapply(names(ptms_detected), function(x) sprintf('%s: %s', x, ptms_detected[[x]])), collapse='; ')

    message(sprintf('The following PTMs were detected: %s', ptms_detected_str))
  }

  rowData(obj) <- bind_cols(data.frame(rowData(obj)), ptm_pos)

  return(obj)
}


#' Add rowData columns with positions of PTMs with respect to peptide sequence
#'
#' @description This function takes a `SummarizedExperiment` object with the evidence.txt
#' output from MaxQuant and adds rowData columns to describe the PTMs present in the peptide and their positions within the peptide' '
#' @param obj `SummarizedExperiment`. Proteomics dataset
#' @param ptms_to_retain `character` vector with PTMs to retain. Inspect values in rowData(obj)$Modified.sequence to get correct PTM name
#' @param ptm_encoding_pos `character` name vector describing whether modifification comes before (-1) or after (1) amino acid in Modified.sequence column values
#' @param prob_col `character` name of column containing PTM probabilities
#' @param min_prob `numeric` Minimum acceptable probability for PTM localisation
#' @param filter_pep_by_prob `logical` Filter the output to only return cases where the number of sites passing the probability threshold
#' equals the number of PTMs in the peptide
#' @param verbose `logical` Describe the number of PTMs detected
#' @return Returns a `SummarizedExperiment` with an additional column in the RowData describing the position of the PTMs with respect to the peptide sequence
#' @export
#' @export
add_filter_ptm_pos_rowdata_mq <- function(obj,
         ptms_to_retain = c('Phospho (STY)'),
         prob_col = 'Phospho..STY..Probabilities',
         min_prob=0.501,
         filter_pep_by_prob=TRUE,
         verbose=TRUE){

  check_se(obj)

  mod_seq <- gsub("_$|^_", "", rowData(obj)$Modified.sequence)
  mod_matches <- gregexpr('\\(.*?\\)\\)', mod_seq)
  mods <- regmatches(mod_seq, mod_matches)
  pep_split <- regmatches(mod_seq, mod_matches, invert = TRUE)

  prob_matches <- gregexpr('\\((0.\\d+|1)\\)', rowData(obj)[[prob_col]])
  mod_probs <- regmatches(rowData(obj)[[prob_col]], prob_matches)
  pep_split_mod <- regmatches(rowData(obj)[[prob_col]], prob_matches, invert = TRUE)

  ptms_detected <- list()

  ptm_pos <- lapply(1:length(mod_seq), function(seq_ix){

    seq_mods <- mods[[seq_ix]]
    names(seq_mods) <- head(cumsum(nchar(pep_split[[seq_ix]])), -1)

    seq_probs <- as.numeric(gsub('(\\(|\\))', '', mod_probs[[seq_ix]]))
    names(seq_probs) <- head(cumsum(nchar(pep_split_mod[[seq_ix]])), -1)

    ptm_positions <- c()
    ptms <- c()
    ptm_amino_acids <- c()
    n_ptms_detected <- 0

    for(ptm_pos in names(seq_mods)){

      seq_mod <- gsub('(^\\(|\\)$)', '', seq_mods[[ptm_pos]])

      if(seq_mod %in% names(ptms_detected)){
        ptms_detected[[seq_mod]] <<- ptms_detected[[seq_mod]] + 1
      } else{
        ptms_detected[[seq_mod]] <<- 1
      }

      if(seq_mod %in% ptms_to_retain){
        n_ptms_detected <- n_ptms_detected + 1

        mod_prob <- seq_probs[ptm_pos]

        if(mod_prob>=min_prob){
          ptms <- c(ptms, seq_mod)
          ptm_positions <- c(ptm_positions, ptm_pos)
          ptm_amino_acids <- c(ptm_amino_acids, substr(rowData(obj)$Sequence[[seq_ix]], as.numeric(ptm_pos), as.numeric(ptm_pos)))
        }
      }
    }

    n_ptms = length(ptms)
    ptms <- paste(ptms, collapse = '; ')
    ptm_positions <- paste(ptm_positions, collapse = '; ')
    ptm_amino_acids <- paste(ptm_amino_acids, collapse = '; ')

    return(c('ptms'=ptms, 'ptm_positions'=ptm_positions, 'ptm_amino_acids'=ptm_amino_acids,
             'n_ptms'=n_ptms, 'n_ptms_detected'=n_ptms_detected))
  }) %>% bind_rows()


  if(verbose){
    ptms_detected <- ptms_detected[order(-as.numeric(ptms_detected))]
    ptms_detected_str <- paste(sapply(names(ptms_detected), function(x) sprintf('%s: %s', x, ptms_detected[[x]])), collapse='; ')

    message(sprintf('The following PTMs were detected: %s', ptms_detected_str))
  }

  rowData(obj) <- bind_cols(data.frame(rowData(obj)), ptm_pos)

  if(filter_pep_by_prob){
    obj <- obj[rowData(obj)$n_ptms==rowData(obj)$n_ptms_detected,]
  }
  return(obj)
}


#' Add peptide positions with respect to the protein
#'
#' @description This function takes a `SummarizedExperiment` object and adds rowData
#' columns to describe the positions of the peptides within the proteins', taking
#' into account the possible peptides according to the digestion enzyme used.
#'
#' @param obj `SummarizedExperiment`. Proteomics dataset
#' @param proteome_fasta `character` filepath to fasta with protein sequences
#' digest_enzyme = "trypsin-simple",
#' @param digest_enzyme `character`. Enzyme used. See `?cleaver::cleave`
#' @param missed_cleavages `numeric`. Vector of allowed number of missed cleavages
#' @param master_protein_col `character`. Name of column containing master proteins
#' @param sequence_col `character`. Name of column containing peptide sequences
#' @return Returns a `SummarizedExperiment` with an additional column in the RowData describing the position of the PTMs with respect to the protein
#' @export
add_peptide_positions_from_cleavage <- function(obj,
                                                proteome_fasta,
                                                digest_enzyme = "trypsin-simple",
                                                missed_cleavages = c(0,1,2),
                                                master_protein_col = "Leading.razor.protein",
                                                sequence_col = "Sequence") {
  check_se(obj)

  proteins <- Biostrings::readAAStringSet(proteome_fasta)

  # Update name to be just UniprotID
  names(proteins) <- gsub('(sp|tr)\\|(\\S*)\\|.*', '\\2', names(proteins))

  digest_pep <- cleaver::cleave(proteins, enzym=digest_enzyme, unique=FALSE, missedCleavages = missed_cleavages)

  pep_ranges <- cleaver::cleavageRanges(proteins, enzym=digest_enzyme, missedCleavages = missed_cleavages)

  pep_loc <- unique(names(pep_ranges)) %>% lapply(function(x){
    .df <- data.frame(pep_ranges[[x]]) %>%
      mutate(sequence=as.character(digest_pep[[x]]),
             width=width(digest_pep[[x]]))

    if(nrow(filter(.df, width!=(end+1)-start))!=0){

      print(filter(.df, width!=(end+1)-start))
      stop(sprintf('Widths do not agree with start and end positions for protein: %s', x))
    }

    .df_trim_leading_M <- .df %>%
      filter(start==1, grepl('^M', sequence)) %>%
      mutate(sequence=gsub('^M', '', sequence), start=2, width=width-1)

    .df <- bind_rows(.df, .df_trim_leading_M)

    .df <- .df %>% group_by(sequence) %>%
      summarise(start=paste(start, collapse=';'),
                end=paste(end, collapse=';')) %>%
      mutate(protein=x)
    return(.df)
  }) %>% bind_rows()

  rowData(obj) <- rowData(obj) %>%
    data.frame() %>%
    mutate(original_order=1:nrow((rowData(obj)))) %>%
    merge(pep_loc,
          by.x=c(master_protein_col, sequence_col),
          by.y=c('protein', 'sequence'),
          all.x=TRUE) %>%
    arrange(original_order) %>%
    select(-original_order)

  return(obj)
}

#' Add rowData columns with details of PTMs positions
#'
#' @description This function takes a `SummarizedExperiment` object ...
#'
#' @param obj `SummarizedExperiment`. Proteomics dataset
#' @param proteome_fasta `character` filepath to fasta with protein sequences
#' digest_enzyme = "trypsin-simple",
#' @param digest_enzyme `character`. Enzyme used. See `?cleaver::cleave`
#' @param missed_cleavages `numeric`. Vector of allowed number of missed cleavages
#' @param master_protein_col `character`. Name of column containing master proteins
#' @param sequence_col `character`. Name of column containing peptide sequences
#' @return Returns a `SummarizedExperiment` with an additional column in the RowData describing the position of the PTMs with respect to the protein
#' @export
add_ptm_positions <- function(obj,
                              proteome_fasta,
                              digest_enzyme = "trypsin-simple",
                              missed_cleavages = c(0,1,2),
                              master_protein_col = "Leading.razor.protein",
                              sequence_col = "Sequence") {

  obj <- add_peptide_positions_from_cleavage(obj,
                                             proteome_fasta,
                                             digest_enzyme = digest_enzyme,
                                             missed_cleavages = missed_cleavages,
                                             master_protein_col = master_protein_col,
                                             sequence_col = sequence_col)

  rowData(obj)$ptm_positions_prot <- mapply(FUN=function(x, y) paste(as.numeric(x) + as.numeric(y) - 1, collapse = '; '),
                                            strsplit(rowData(obj)$ptm_positions, '; '),
                                            rowData(obj)$start)

  rowData(obj)$ptm_name <- mapply(FUN=function(x, y) paste(paste0(x,y), collapse='; '),
                                  strsplit(rowData(obj)$ptm_amino_acids, '; '),
                                  strsplit(rowData(obj)$ptm_positions_prot, '; '))


  return(obj)
}


#' Parse the PTM probabilities from Proteome Discoverer and add new columns with PTM information
#'
#' @description Extract PTM information from the PTM probabilities column in
#' the output from PD and add new columns with this information. Also optionally
#' prints a summary of how many features (PSMs/peptides) pass a given probability
#' threshold
#'
#' @param obj `SummarizedExperiment`. Proteomics dataset
#' @param threshold `numeric` If any score is below a set threshold, disregard all putative PTM sites
#' @param ptm_col `character` Columm name for PTM probabilities
#' @param prob_split `character` regex to split PTM probabilities
#' @param collapse_delimiter `character` delimiter for multiple values in output columns
#' @param verbose Set TRUE to print log of events to console
#'
#' @return `SummarizedExperiment`
#' @export
parse_PTM_scores_pd <- function(obj,
                                threshold = 95,
                                ptm_col = "ptmRS.Best.Site.Probabilities",
                                prob_split = '; |: ',
                                collapse_delimiter = "; ",
                                verbose = TRUE) {

  check_se(obj)
  new_obj <- obj

  if(verbose){
    message(sprintf(
      "Removed %s Features where the ptm_col value == `Inconclusive data`",
      sum(rowData(new_obj)[[ptm_col]] == "Inconclusive data")
    ))
  }

  new_obj <- new_obj[rowData(new_obj)[[ptm_col]] != "Inconclusive data", ]

  # initiate vectors with empty string
  # where the PTM(s) pass the threshold, these will be updated
  filtered_ptm_desc <- filtered_ptm_res <- filtered_ptm_pos <- filtered_ptm <- filtered_ptm_score <- rep("", nrow(rowData(new_obj)))

  split_probabities <- strsplit(rowData(new_obj)[[ptm_col]], split = prob_split)

  log <- list(
    "Total Features" = 0,
    "Total detected PTMFeatures" = 0,
    "Features passing filter" = 0,
    "Features failing filter" = 0,
    "BiPTM/multiPTM Features where some sites fail filter" = 0,
    "Total detected sites" = 0,
    "Sites passing filter" = 0,
    "Sites failing filter" = 0,
    "monoPTM passing filter" = 0,
    "biPTM passing filter" = 0,
    "multiPTM passing filter" = 0,
    "Too many isoforms" = 0
  )

  for (i in seq_along(split_probabities)) {
    log[["Total Features"]] <- log[["Total Features"]] + 1

    peptide_ptm_scores <- split_probabities[[i]]
    if (length(peptide_ptm_scores) == 0) {
      next()
    } # no PTM detected
    if (is.na(peptide_ptm_scores[[1]])) {
      next()
    } # no PTM detected

    log[["Total detected PTMFeatures"]] <- log[["Total detected PTMFeatures"]] + 1

    if (peptide_ptm_scores[[1]] == "Too many isoforms") {
      log[["Too many isoforms"]] <- log[["Too many isoforms"]] + 1
      log[["Features failing filter"]] <- log[["Features failing filter"]] + 1
      next()
    }

    log[["Total detected sites"]] <- log[["Total detected sites"]] + length(peptide_ptm_scores) / 2

    scores <- peptide_ptm_scores[seq(2, length(peptide_ptm_scores), 2)]

    # if any score is below threshold, disregard all putative ptm sites
    if (any(as.numeric(scores) < threshold)) {
      log[["Sites failing filter"]] <- log[["Sites failing filter"]] + length(peptide_ptm_scores) / 2
      log[["Features failing filter"]] <- log[["Features failing filter"]] + 1
      if (any(as.numeric(scores) >= threshold)) {
        log[["BiPTM/multiPTM Features where some sites fail filter"]] <- log[[
          "BiPTM/multiPTM Features where some sites fail filter"]] + 1
      }
      # if we want to handle this differently, can implement an alternative approach here
      # and move the rest of the code below into an else clause
      next()
    }

    log[["Sites passing filter"]] <- log[["Sites passing filter"]] + length(peptide_ptm_scores) / 2
    log[["Features passing filter"]] <- log[["Features passing filter"]] + 1

    if (length(scores) == 1) {
      log[["monoPTM passing filter"]] <- log[["monoPTM passing filter"]] + 1
    }
    else if (length(scores) == 2) {
      log[["biPTM passing filter"]] <- log[["biPTM passing filter"]] + 1
    }
    else {
      log[["multiPTM passing filter"]] <- log[["multiPTM passing filter"]] + 1
    }

    ptms <- peptide_ptm_scores[seq(1, length(peptide_ptm_scores), 2)] # extract the PTMs info
    split_ptms <- unlist(strsplit(ptms, split = '\\(|\\)')) # split to remove parantheses
    modifications <- split_ptms[seq(2, length(split_ptms), 2)] # extract modifications, e.g "phospho"
    positions <- split_ptms[seq(1, length(split_ptms), 2)] # extract the positions, e.g "S6"
    residues <- substr(positions, 1, 1) # extract first element, e.g S
    positions <- sub('.', '', positions) # remove first element and leave position, e.g 6

    # paste together the value, separated by option(collapse_delimiter) and update vectors which will become columns
    #filtered_ptm_desc[[i]] <- paste(peptide_ptm_scores, collapse = collapse_delimiter)
    filtered_ptm_res[[i]] <- paste(residues, collapse = collapse_delimiter)
    filtered_ptm_pos[[i]] <- paste(positions, collapse = collapse_delimiter)
    filtered_ptm[[i]] <- paste(modifications, collapse = collapse_delimiter)
    filtered_ptm_score[[i]] <- paste(scores, collapse = collapse_delimiter)
  }

  # add columns
  #rowData(new_obj)['filtered_PTM_desc'] = filtered_ptm_desc
  rowData(new_obj)['ptm_amino_acids'] = filtered_ptm_res
  rowData(new_obj)['ptm_positions'] = filtered_ptm_pos
  rowData(new_obj)['ptms'] = filtered_ptm
  rowData(new_obj)['ptm_scores'] = filtered_ptm_score

  if (verbose) {
    for (event in names(log)) {
      message(sprintf("%s: %i", event, log[[event]]))
    }
  }

  return(new_obj)
}



#' Get the amino acid sequence around a PTM
#'
#' @description Get the amino acid sequence around a PTM. Will return NA if
#' peptide maps to multiple proteins or has multiple PTMs. If padding extends
#' outside the protein AA sequence, padding will be extended with '_'.
#'
#' @param proteome `XStringSet` as generated by `Biostrings::readAAStringSet`
#' @param protein `character` protein ID
#' @param ptm_position `numeric` position of the PTM in the protein
#' @param pad `numeric` Number of amino acids around PTM
#'
#' @return `character` PSM-centered amino acid sequence
#' @export
get_sequence <- function(proteome, protein, ptm_position, pad = 7) {
  ptm_position <- suppressWarnings(as.numeric(as.character(ptm_position)))

  if (is.na(ptm_position)) {
    return(NA)
  }

  # if (grepl("; ", ptm_position)) {
  #   return(NA)
  # }

  if (!protein %in% names(proteome)) {
    return(NA)
  }

  protein_length <- length(proteome[[protein]])

  if (ptm_position > protein_length) {
    warning(sprintf(
      "PTM positions is outside protein sequence! Returning NA. %s: [-%s], PTM: %s",
      protein, protein_length, ptm_position
    ))
    return(NA)
  }

  start_pad <- end_pad <- ""

  start <- ptm_position - pad
  if (start <= 0) {
    start_pad <- paste0(rep("_", (start * -1) + 1), collapse = "")
    start <- 0
  }

  end <- ptm_position + pad
  if (end > protein_length) {
    end_pad <- paste0(rep("_", (end - protein_length)), collapse = "")
    end <- protein_length
  }

  mod_position <- pad + 1

  sequence <- as.character(proteome[[protein]][start:end])
  sequence <- paste0(start_pad, sequence, end_pad)

  sequence <- paste(base::substr(sequence, 1, pad),
                    tolower(base::substr(sequence, pad + 1, pad + 1)),
                    base::substr(sequence, pad + 2, pad + pad + 1),
                    sep = ""
  )

  return(sequence)
}

#' Add a column with amino acid sequence around a PTM
#'
#' @description Add a rowData column with amino acid sequence around a PTM. Value will be
#' NA if peptide maps to multiple proteins or has multiple PTMs. If padding extends
#' outside the protein AA sequence, padding will be extended with '_'. The PTM-
#' centered AA sequence is useful to integrate external databases. Input
#' should have been been passed through `add_ptm_positions`
#' to add the 'ptm_positions_prot' column to the rowData.
#'
#' @param obj `SummarizedExperiment`. Proteomics dataset
#' @param proteome_fasta `character` Filepath for proteome fasta
#' @param master_protein_col `character` Column name for master protein
#' @param sequence_pad `numeric` Number of amino acids to include on either side of PTM
#' @return `SummarizedExperiment`
#' @export
add_site_sequence <- function(obj,
                              proteome_fasta,
                              master_protein_col = "Master.Protein.Accessions",
                              sequence_pad=7) {

  check_se(obj)

  proteome <- Biostrings::readAAStringSet(proteome_fasta)
  # Update name to be just UniprotID
  names(proteome) <- gsub('(sp|tr)\\|(\\S*)\\|.*', '\\2', names(proteome))

  new_obj <- obj

  rowData(new_obj) <- rowData(new_obj) %>%
    data.frame() %>%
    rowwise() %>%
    mutate(site_seq = get_sequence(
      proteome,
      !!sym(master_protein_col),
      .data$ptm_positions_prot,
      pad = sequence_pad
    )) %>% ungroup()

  return(new_obj)
}
