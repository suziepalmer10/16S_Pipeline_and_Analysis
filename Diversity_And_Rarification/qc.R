#' Filter phyloseq object
#' @param ps a phyloseq object
#' @param taxa.min.abund retain taxa with minimial mean abundance >= taxa.min.abund
#' @param taxa.min.sample retain taxa with at least `taxa.min.sample` have reads 
#' @param sample.min.read retain sample having >= `sample.min.read`
#' @return a phyloseq object
qc <- function(ps,
               taxa.min.abund = 1e-5,
               taxa.min.sample = 1,
               sample.min.read = 5,
               sample.name.exclude = "Zymopositive") {

    ## remove strange/benchmark sample
    old = nrow(sample_data(ps))
    idx <- !grepl(pattern = sample.name.exclude, rownames(sample_data(ps)))
    sm <- rownames(sample_data(ps))[idx]
    ps <- phyloseq::prune_samples(sm, ps)
    # ps <- phyloseq::prune_samples(ps, sm)
    new = nrow(sample_data(ps))
    if (new != old) {
        cat(old - new, "samples filtered by sample names\n")
    }
    
    ## filter taxa
    ##   (mean fractional abundance > 1e-5)
    ##   taxa must have >= 1 samples
    old = nrow(tax_table(ps))
    ps = phyloseq::filter_taxa(ps, function(x) {mean(x, na.rm = TRUE) > 1e-5}, prune = TRUE)
    ps = phyloseq::filter_taxa(ps, function(x) {sum(!is.na(x)) > 1}, prune = TRUE)
    new = nrow(tax_table(ps))
    if (new != old) {
        cat(old - new, "taxa filtered by minimal abundance ", taxa.min.abund,
            " or number of samples ", taxa.min.abund, "\n")
    }
    
    ## filter samples
    ##   sample has >0 read
    ##   per-sample sum of read > 5
    old = nrow(sample_data(ps))    
    ps = prune_samples(!is.na(sample_sums(ps)), ps)
    ps = prune_samples(sample_sums(ps) > sample.min.read, ps)
    new = nrow(sample_data(ps))
    if (new != old) {
        cat(old -  new, "samples filtered by minimal read ", sample.min.read, "\n")
    }

    ps
}
