#' Encode DNA 5-mer sequences as one-hot factors
#'
#' Given a character vector of fixed-length DNA (or RNA) sequences
#' (for example 5-mers), this helper converts each nucleotide into
#' a separate column and returns a data frame whose columns are
#' factors with levels A, T, C and G.
#'
#' @param dna_strings A character vector containing equal-length
#'   nucleotide strings.
#'
#' @return A data frame with one row per input sequence and one
#'   factor column per nucleotide position.
#'
#' @examples
#' \dontrun{
#' dna_encoding(c("ATCGA", "GGGTT"))
#' }
dna_encoding <- function(dna_strings) {
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(
    unlist(strsplit(dna_strings, "")),
    ncol = nn,
    byrow = TRUE
  )
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m, stringsAsFactors = FALSE)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Predict m6A status for multiple sequences
#'
#' This function prepares feature columns and encoded 5-mer DNA
#' sequences, then uses a fitted random forest model to obtain
#' m6A probabilities and Positive/Negative calls for many rows
#' at once.
#'
#' @param ml_fit A fitted \code{randomForest} classification model
#'   for m6A prediction.
#' @param feature_df A data frame containing at least the columns
#'   \code{"gc_content"}, \code{"RNA_type"}, \code{"RNA_region"},
#'   \code{"exon_length"}, \code{"distance_to_junction"},
#'   \code{"evolutionary_conservation"} and \code{"DNA_5mer"}.
#' @param positive_threshold Numeric value between 0 and 1 giving
#'   the decision threshold on the Positive-class probability.
#'   Defaults to 0.5.
#'
#' @return A copy of \code{feature_df} with two extra columns:
#'   \code{predicted_m6A_prob} (numeric probabilities) and
#'   \code{predicted_m6A_status} (factor or character, "Positive"
#'   or "Negative").
#'
#' @examples
#' \dontrun{
#'   # Load model and example feature data shipped with the package
#'   rf_model <- readRDS(system.file("extdata", "rf_fit.rds",
#'                                   package = "m6APrediction"))
#'   example_df <- read.csv(system.file("extdata", "m6A_input_example.csv",
#'                                      package = "m6APrediction"))
#'
#'   pred_df <- prediction_multiple(rf_model, example_df,
#'                                  positive_threshold = 0.6)
#'   head(pred_df)
#' }
#'
#' @export
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5) {
  required_cols <- c(
    "gc_content",
    "RNA_type",
    "RNA_region",
    "exon_length",
    "distance_to_junction",
    "evolutionary_conservation",
    "DNA_5mer"
  )
  stopifnot(all(required_cols %in% colnames(feature_df)))

  df <- feature_df

  # ensure factor levels match training
  df$RNA_type <- factor(df$RNA_type,
                        levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  df$RNA_region <- factor(df$RNA_region,
                          levels = c("CDS", "intron", "3'UTR", "5'UTR"))

  enc <- dna_encoding(df$DNA_5mer)

  pred_input <- cbind(
    gc_content = df$gc_content,
    RNA_type = df$RNA_type,
    RNA_region = df$RNA_region,
    exon_length = df$exon_length,
    distance_to_junction = df$distance_to_junction,
    evolutionary_conservation = df$evolutionary_conservation,
    enc
  )

  prob_mat <- predict(ml_fit, newdata = pred_input, type = "prob")
  pos_prob <- prob_mat[, "Positive"]

  status <- ifelse(pos_prob >= positive_threshold, "Positive", "Negative")

  out <- feature_df
  out$predicted_m6A_prob   <- pos_prob
  out$predicted_m6A_status <- status

  return(out)
}

#' Predict m6A status for a single sequence
#'
#' This convenience wrapper creates a one-row feature data frame
#' from user-supplied covariates and calls
#' \code{\link{prediction_multiple}} to obtain a probability and a
#' Positive/Negative label for one candidate site.
#'
#' @param ml_fit A fitted \code{randomForest} classification model
#'   for m6A prediction.
#' @param gc_content Numeric GC content of the transcript region.
#' @param RNA_type Character or factor giving the transcript type
#'   (e.g. \code{"mRNA"}, \code{"lincRNA"}, \code{"lncRNA"},
#'   \code{"pseudogene"}).
#' @param RNA_region Character or factor describing the region
#'   (e.g. \code{"CDS"}, \code{"intron"}, \code{"3'UTR"},
#'   \code{"5'UTR"}).
#' @param exon_length Numeric exon length.
#' @param distance_to_junction Numeric distance to the nearest
#'   exon junction.
#' @param evolutionary_conservation Numeric conservation score.
#' @param DNA_5mer Character string containing the 5-mer centred
#'   on the candidate m6A site.
#' @param positive_threshold Numeric decision threshold on the
#'   Positive-class probability. Defaults to 0.5.
#'
#' @return A named vector of length two with elements
#'   \code{predicted_m6A_prob} and \code{predicted_m6A_status}.
#'
#' @examples
#' \dontrun{
#'   rf_model <- readRDS(system.file("extdata", "rf_fit.rds",
#'                                   package = "m6APrediction"))
#'   prediction_single(
#'     rf_model,
#'     gc_content = 0.55,
#'     RNA_type   = "mRNA",
#'     RNA_region = "CDS",
#'     exon_length          = 2200,
#'     distance_to_junction = 150,
#'     evolutionary_conservation = 0.8,
#'     DNA_5mer = "ATGCG",
#'     positive_threshold = 0.5
#'   )
#' }
#'
#' @import randomForest
#' @export
prediction_single <- function(ml_fit,
                              gc_content,
                              RNA_type,
                              RNA_region,
                              exon_length,
                              distance_to_junction,
                              evolutionary_conservation,
                              DNA_5mer,
                              positive_threshold = 0.5) {
  single_df <- data.frame(
    gc_content               = gc_content,
    RNA_type                 = RNA_type,
    RNA_region               = RNA_region,
    exon_length              = exon_length,
    distance_to_junction     = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer                 = DNA_5mer,
    stringsAsFactors         = FALSE
  )

  pred_df <- prediction_multiple(ml_fit, single_df, positive_threshold)
  returned_vector <- c(
    pred_df$predicted_m6A_prob[1],
    as.character(pred_df$predicted_m6A_status[1])
  )
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  return(returned_vector)
}
