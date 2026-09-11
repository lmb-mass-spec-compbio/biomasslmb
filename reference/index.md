# Package index

## Protein annotation

Retrieving UniProt annotations, GO terms and membrane topology for the
accessions in your data. Worked through in the [protein
annotation](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/articles/protein_annotation.md)
article.

- [`get_uniprot_details()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_uniprot_details.md)
  : Get UniProt annotation details, with UniParc fallback for retired
  accessions
- [`collapse_uniprot_details_multi_accession()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/collapse_uniprot_details_multi_accession.md)
  : Map per-accession UniProt details onto protein IDs with multiple
  accessions
- [`get_uniparc_fallback()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_uniparc_fallback.md)
  : Recover annotations for retired UniProtKB accessions via UniParc
- [`check_uniprot_release()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/check_uniprot_release.md)
  : Check the current UniProt release
- [`add_gene_long_protein_name_pd()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/add_gene_long_protein_name_pd.md)
  : Extract the gene name and long-form protein name from the master
  protein descriptions column
- [`get_go_terms()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_go_terms.md)
  : Obtain GO term annotations for proteins
- [`get_all_mappings()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_all_mappings.md)
  : Get all mappings for GO terms
- [`get_ancestor_go()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_ancestor_go.md)
  : Get all ancestor GO terms
- [`expand_terms()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/expand_terms.md)
  : Expand data.frame GO terms
- [`determine_ancestor_function()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/determine_ancestor_function.md)
  : Determine GO ancestor object
- [`determine_offspring_function()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/determine_offspring_function.md)
  : Determine GO offspring object
- [`get_protein_tm_topology()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_protein_tm_topology.md)
  : Obtain protein transmembrane-domains and topology information
- [`query_protein_tm_topology()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/query_protein_tm_topology.md)
  : Query UniProt to determine protein transmembrane-domains and
  topology
- [`add_tm_info()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/add_tm_info.md)
  : Add transmembrane domain details
- [`add_topology_info()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/add_topology_info.md)
  : Add topology details

## Reading and filtering search engine output

Getting features out of Proteome Discoverer, MaxQuant, DIA-NN or
Spectronaut output and removing the ones that should not be quantified.
One function per search engine, plus the filters common to all of them
and the helper that reports how many features and master proteins
survive a step.

- [`filter_features_pd_dda()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/filter_features_pd_dda.md)
  : Filter Proteome Discoverer DDA output
- [`filter_features_mq_dda()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/filter_features_mq_dda.md)
  : Filter MaxQuant DDA output
- [`filter_features_diann()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/filter_features_diann.md)
  : Filter DIA-NN output
- [`filter_features_sn()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/filter_features_sn.md)
  : Filter Spectronaut output
- [`readDIANNFilterQJoin()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/readDIANNFilterQJoin.md)
  : Read in data from DIA-NN
- [`filter_TMT_PSMs()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/filter_TMT_PSMs.md)
  : Filter a PSM-level summarizedExperiment to remove low quality PSMs
- [`filter_features_per_protein()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/filter_features_per_protein.md)
  : Remove features which are assigned to a protein with too few
  supporting features in total
- [`filter_complete_groups()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/filter_complete_groups.md)
  : Filter a QFeatures assay to samples with complete group annotations
- [`filter_by_protein_fdr()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/filter_by_protein_fdr.md)
  : Filter to remove peptides from proteins failing the protein FDR
  threshold
- [`remove_redundant_psm_quant()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/remove_redundant_psm_quant.md)
  : Remove redundant Peptide Spectrum Matches
- [`update_peptide_assignments()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/update_peptide_assignments.md)
  : Update PSM level protein/master protein assignments using
  peptide-level assignments
- [`update_average_sn()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/update_average_sn.md)
  : Adds new feature describing the average reporter Signal/Noise ratio.
- [`message_parse()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/message_parse.md)
  : Report how many features and master proteins remain

## Contaminants

Building the accession lists the filtering functions match against. See
[contaminants and protein
FDR](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/articles/gotcha_contaminants_and_FDR.md)
for what each defence actually catches.

- [`get_contaminant_fasta_accessions()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_contaminant_fasta_accessions.md)
  : Extract the contaminants protein accessions from a contaminants
  fasta file
- [`get_maxquant_cont_accessions()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_maxquant_cont_accessions.md)
  : Extract the contaminants protein accessions from a MaxQuant
  contaminants fasta file
- [`get_crap_fasta_accessions()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_crap_fasta_accessions.md)
  : Extract the contaminants protein accessions from a cRAP fasta file
- [`sub_crap()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/sub_crap.md)
  : Insert cRAP numbers into a character vector

## Summarisation and normalisation

Turning features into protein-level abundance, and putting samples or
plexes onto a common scale.

- [`maxlfq_wrapper()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/maxlfq_wrapper.md)
  : Wrapper for iq::maxLFQ to use with QFeatures::aggregateFeatures
- [`get_protein_no_quant_mask()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_protein_no_quant_mask.md)
  : Identify proteins which have too few features to quantify protein
  abundance in each sample
- [`mask_protein_level_quant()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/mask_protein_level_quant.md)
  : Replace protein-level quantifications with NA if they derive from
  too few lower feature level quantifications
- [`get_n_feature_per_prot()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_n_feature_per_prot.md)
  : Identify how many features (PSMs/Peptides) are quantified for each
  protein
- [`get_medians()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_medians.md)
  : Extract the assay column medians from an MSnSet
- [`center_normalise_to_ref()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/center_normalise_to_ref.md)
  : Center-median normalise the expression matrix in an MSnSet using
  medians from a reference dataset
- [`bridge_normalise()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/bridge_normalise.md)
  : Put TMT plexes onto a common scale using bridge (pooled reference)
  channels

## Working with QFeatures objects

Helpers for the object every analysis is built in. See [working with
QFeatures
objects](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/articles/qfeatures_objects.md)
for the object model and the design-table contract.

- [`sync_coldata()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/sync_coldata.md)
  : Copy the experimental design onto one or more assays

## Missing values

Diagnosing how much missingness there is and whether it is explained by
abundance or by experimental condition, and imputing where that is
defensible. See [handling missing
values](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/articles/handling_missing_values.md).

- [`condition_miss_score()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/condition_miss_score.md)
  : Score features by how well missingness can be predicted by
  experimental condition
- [`condition_miss_index()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/condition_miss_index.md)
  : Summarise per-feature condition missingness scores into a
  dataset-level index
- [`global_condition_miss_score()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/global_condition_miss_score.md)
  : Compute a dataset-level score of condition-predictable missingness
  using logistic regression
- [`restrict_imputation()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/restrict_imputation.md)
  : Restrict imputed values to specific conditions
- [`create_long_form_imputed_data()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/create_long_form_imputed_data.md)
  : Create long format data with column defining if quantification value
  is imputed
- [`plot_missing_upset()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/plot_missing_upset.md)
  : Plot the most common missing value patterns
- [`plot_missing_SN()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/plot_missing_SN.md)
  : Plot the missing values vs signal:noise
- [`plot_missing_SN_per_sample()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/plot_missing_SN_per_sample.md)
  : Plot the missing values vs signal:noise for each sample
- [`mnar_score()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/mnar_score.md)
  : Deprecated: use condition_miss_score()
- [`mnar_index()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/mnar_index.md)
  : Deprecated: use condition_miss_index()
- [`mnar_global_score()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/mnar_global_score.md)
  : Deprecated: use global_condition_miss_score()

## Post-translational modifications

Localising a modification to a residue, placing it within the protein
and summarising to sites. Worked through in [PTM site
quantification](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/articles/PTM_site_quantification.md).

- [`parse_PTM_scores_pd()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/parse_PTM_scores_pd.md)
  : Parse the PTM probabilities from Proteome Discoverer and add new
  columns with PTM information
- [`filter_maxquant_ptm()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/filter_maxquant_ptm.md)
  : Filter PTM data from MaxQuant to retain those with a given PTM
- [`add_ptm_positions()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/add_ptm_positions.md)
  : Add rowData columns with details of PTMs positions
- [`add_ptm_pos_rowdata_mq()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/add_ptm_pos_rowdata_mq.md)
  : Add rowData columns with positions of PTMs with respect to peptide
  sequence
- [`add_filter_ptm_pos_rowdata_mq()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/add_filter_ptm_pos_rowdata_mq.md)
  : Add rowData columns with positions of PTMs with respect to peptide
  sequence
- [`add_site_sequence()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/add_site_sequence.md)
  : Add a column with amino acid sequence around a PTM

## Peptides, sequences and FASTA files

In silico digestion, peptide positions within a protein, and retrieving
or assembling sequence databases.

- [`make_fasta()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/make_fasta.md)
  : Make a FASTA using UniProt accessions
- [`append_fasta()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/append_fasta.md)
  : Append sequences to end of a FASTA
- [`get_sequence()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_sequence.md)
  : Get the amino acid sequence around a PTM
- [`add_peptide_positions()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/add_peptide_positions.md)
  : Add a column describing the position of the peptide sequence with
  respect to the protein
- [`add_peptide_positions_from_cleavage()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/add_peptide_positions_from_cleavage.md)
  : Add peptide positions with respect to the protein
- [`extract_peptide_positions()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/extract_peptide_positions.md)
  : Extract Peptide Position Start and End Columns from a QFeatures
  Object
- [`add_tryptic_termini_flags()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/add_tryptic_termini_flags.md)
  : Evaluate peptide termini against theoretical cleavage positions

## Exploration and plotting

Checking that an experiment worked, and displaying a result. Every
function here returns a `ggplot` object, so it can be modified in the
usual way.

- [`plot_quant()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/plot_quant.md)
  : Plot distributions for feature intensities per sample.
- [`plot_pca()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/plot_pca.md)
  : Create a Principal Component plot from the feature quantification
- [`plot_cor_samples()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/plot_cor_samples.md)
  : Plot the correlation between sample
- [`plot_protein_assays()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/plot_protein_assays.md)
  : Plot quantification at multiple levels
- [`plot_volcano()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/plot_volcano.md)
  : Volcano plot for differential abundance testing
- [`plot_rt_dist()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/plot_rt_dist.md)
  : Plot the retention time distribution
- [`plot_rt_vs_delta()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/plot_rt_vs_delta.md)
  : Plot the retention time vs delta precursor mass
- [`get_samples_present()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_samples_present.md)
  : Extract the number of samples each feature was detected in for each
  experiment in a Qfeatures object
- [`plot_samples_present()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/plot_samples_present.md)
  : Plots the number of samples each feature was detected in for each
  experiment in a Qfeatures object
- [`theme_biomasslmb()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/theme_biomasslmb.md)
  : A ggplot2 theme for the package
- [`get_cat_palette()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_cat_palette.md)
  : Generate a colour-blind friendly palette for categorical colour
  encoding

## Functional enrichment

Over-representation testing against the proteins that were tested, with
a correction for abundance bias. See [functional
enrichment](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/articles/functional_enrichment.md)
and [ORA or
GSEA?](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/articles/ORA_vs_GSEA.md).

- [`get_enriched_go()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/get_enriched_go.md)
  : GO term enrichment using goseq
- [`estimate_overrep()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/estimate_overrep.md)
  : Estimate effect size of over-representation
- [`add_independent_filtering_padj()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/add_independent_filtering_padj.md)
  : Perform independent filtering to threshold the over-representation
  testing on the number of features in each category to limit the
  multiple testing burden
- [`remove_redundant_go()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/remove_redundant_go.md)
  : Remove redundant GO terms
- [`plot_go()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/plot_go.md)
  : Plot selected GO terms
- [`plot_go_terms_upset()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/plot_go_terms_upset.md)
  : Plot selected GO terms

## Statistics

Fitting and adjustment helpers used by the articles.

- [`tls()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/tls.md)
  : Total least squares regression (for use with geom_smooth)
- [`predict(`*`<tls_model>`*`)`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/predict.tls_model.md)
  : Predict method for tls_model
- [`bb2014()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/bb2014.md)
  : Benjamini & Bogomolov (2014) selective-inference adjustment using
  external screening p-values (e.g., F-test p-values)

## Example data

Search engine output and processed `QFeatures` objects used throughout
the articles, so that every example can be run without your own data.

- [`psm_tmt_clock`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/psm_tmt_clock.md)
  : PSM-level PD output for a TMT12plex Control vs Mutant experiment

- [`tmt_clock_design`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/tmt_clock_design.md)
  :

  Experimental design for `psm_tmt_clock`

- [`psm_tmt_2plex`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/psm_tmt_2plex.md)
  : PSM-level PD output for a two-plex TMTpro 18plex experiment

- [`tmt_2plex_design`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/tmt_2plex_design.md)
  :

  Experimental design for `psm_tmt_2plex`

- [`psm_tmt_per2_mq`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/psm_tmt_per2_mq.md)
  : PSM-level MaxQuant output for a TMT18plex IP vs control experiment

- [`tmt_per2_mq_design`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/tmt_per2_mq_design.md)
  :

  Experimental design for `psm_tmt_per2_mq`

- [`psm_tmt_factorial`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/psm_tmt_factorial.md)
  : PSM-level MaxQuant output for a TMT18plex whole-proteome factorial
  experiment

- [`tmt_factorial_design`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/tmt_factorial_design.md)
  :

  Experimental design for `psm_tmt_factorial`

- [`psm_tmt_phospho`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/psm_tmt_phospho.md)
  : PSM-level PD output for the phospho-enriched fraction of a TMTpro
  experiment

- [`psm_tmt_phospho_total`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/psm_tmt_phospho_total.md)
  :

  PSM-level PD output for the total fraction matching `psm_tmt_phospho`

- [`tmt_phospho_design`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/tmt_phospho_design.md)
  :

  Experimental design for `psm_tmt_phospho` and `psm_tmt_phospho_total`

- [`psm_tmt_total`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/psm_tmt_total.md)
  : PSM-level PD output for total proteome TMT10-plex data

- [`tmt_qf`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/tmt_qf.md)
  : TMT data

- [`tmt_qf_mq`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/tmt_qf_mq.md)
  : TMT data (MaxQuant input)

- [`tmt_qf_factorial`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/tmt_qf_factorial.md)
  : Protein-level abundances for the TMT factorial experiment

- [`lfq_qf`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/lfq_qf.md)
  : LFQ-DDA data

- [`lfq_qf_turboid`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/lfq_qf_turboid.md)
  : LFQ-DDA TurboID pulldown data

- [`dia_qf`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/dia_qf.md)
  : LFQ-DIA data

## Utilities

Small helpers for tidying column names and identifiers.

- [`make_unique_all()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/make_unique_all.md)
  : Make Elements of a Character Vector Unique (Numbering From 1)
- [`remove_dots()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/remove_dots.md)
  : Remove duplicated full stops
- [`remove_x()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/remove_x.md)
  : Remove leading X
- [`qfeatures_long()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/qfeatures_long.md)
  : Compatibility wrapper for QFeatures long-format extraction
