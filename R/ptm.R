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
#' @param ptm_encoding_pos `character` name vector describing whether modifification comes before (-1) or after (1) amino acid in Modified.sequence column values.
#' @return Returns a `SummarizedExperiment` with an additional column in the RowData describing the position of the PTMs with respect to the peptide sequence
#' @export
#' @export
add_filter_ptm_pos_rowdata_mq <- function(obj,
         ptms_to_retain = c('Phospho (STY)'),
         prob_col = 'Phospho..STY..Probabilities',
         min_prob=0.501,
         filter_pep_by_prob=TRUE,
         verbose=TRUE){

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

