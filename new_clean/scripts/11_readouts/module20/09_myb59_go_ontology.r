#!/usr/bin/env Rscript
# =============================================================================
# 09_myb59_go_ontology.r -- GO id -> term name AND ontology, dumped once.
#
# WHY THIS EXISTS. Two lookups already in the tree are each missing half of what a
# reader needs. `gene2go_<study>.tsv` carries GO IDs with no name and no ontology;
# `results/conservation/go_term_names.tsv` carries names with no ontology. Splitting
# a gene's terms into BP / MF / CC columns -- which is the form a bench biologist can
# actually read -- needs both.
#
# GO.db has both, but lives in topGO_env, and the steps that consume this are Python.
# So it is dumped ONCE to a small TSV rather than making every downstream step carry
# an R dependency and an env switch.
#
# RUN: ./run_all.sh 09
# =============================================================================

suppressMessages({ library(GO.db); library(AnnotationDbi) })

OUT <- "/dados04/jorge/comparative_saccharum/new_clean/docs/data/go_terms.tsv"
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)

k <- keys(GO.db)
d <- AnnotationDbi::select(GO.db, keys = k, columns = c("TERM", "ONTOLOGY"))
d <- d[!is.na(d$TERM) & !is.na(d$ONTOLOGY), c("GOID", "TERM", "ONTOLOGY")]
names(d) <- c("GO.ID", "Term", "Ontology")

# Obsolete/odd entries carry an ontology outside the three; drop them loudly rather
# than letting a fourth value appear in a column a reader expects to hold BP/MF/CC.
bad <- !d$Ontology %in% c("BP", "MF", "CC")
if (any(bad)) {
  cat(sprintf("dropping %d term(s) with ontology outside BP/MF/CC\n", sum(bad)))
  d <- d[!bad, ]
}

write.table(d, OUT, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("wrote %s\n  %d terms  (BP %d, MF %d, CC %d)\n", OUT, nrow(d),
            sum(d$Ontology == "BP"), sum(d$Ontology == "MF"),
            sum(d$Ontology == "CC")))
