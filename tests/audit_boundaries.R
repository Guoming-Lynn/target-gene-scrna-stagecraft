args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
test_root <- normalizePath(file.path(dirname(file_arg), ".."), winslash = "/")
setwd(test_root)
source(file.path(test_root, 'scripts/part5_model_audit.R'))
cfg <- list(design=list(block='technical_block'), eligibility=list(n_formal=12, min_datasets_formal=3, min_rdf_formal=6), robustness=list(dominant_block='A'))
meta <- data.frame(dataset_donor_id=paste0('d',1:24), technical_block=paste0('t',1:24), source_block=rep(c('A','B','C'),8), dataset=rep(c('A','B','C'),8))
meta$technical_block <- rep(c('t1','t2'),12)
tab <- data.frame(subset=c('FULL','FULL','LA','LB','LC'), model=c('limma_primary','edgeR_QL',rep('limma_primary',3)), gene='G', status='SUCCESS', q_bh=.001, logFC=.5, rdf=20)
run <- function(t=tab, m=meta) audit_models(t,m,cfg,list(),c(A='LA',B='LB',C='LC'),identity)$audit
a <- run()
stopifnot(a$n_donors==24, a$n_blocks==2, a$formal_gates_pass)
t <- tab[tab$subset!='LB',]
stopifnot(!run(t)$formal_gates_pass)
t <- tab; t$status[t$subset=='LA'] <- 'NOT_ESTIMABLE'
stopifnot(run(t)$source_dependent, !run(t)$formal_gates_pass)
t <- tab; t$logFC[t$subset=='LA'] <- -.5
stopifnot(run(t)$source_dependent, !run(t)$formal_gates_pass)
m <- meta; m$dataset_donor_id <- NULL
stopifnot(inherits(try(run(m=m),silent=TRUE),'try-error'))
cat('AUDIT_BOUNDARIES_OK\n')
