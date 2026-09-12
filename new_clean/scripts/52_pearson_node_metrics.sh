#!/usr/bin/env bash
# =============================================================================
# 52_pearson_node_metrics.sh -- the gene universe for the Pearson-only networks.
#
# WHY THIS EXISTS AND WHY IT IS SO SMALL. Two module stages take a node-metrics
# table, and BOTH read only the `gene` column from it:
#
#   15_module_profile.r:82   universe <- unique(strip_version(fread(NODES, select = "gene")$gene))
#   18_module_go.r:118       node_ids <- read.table(NODES, ...)
#
# It is the test universe -- the set of genes that were in the network and could
# therefore have been called -- not a metrics table anything computes on. The
# merged network's table cannot serve: it contains genes that exist only through
# MI edges, and scoring the Pearson-only modules against a universe that includes
# genes this network never had would test a design that was never run.
#
# degree and strength are carried anyway, because they cost nothing: `mcx query`
# reports per-node degree and mean edge weight, so strength = degree * mean.
#
# transitivity is NA, deliberately. No consumer reads it (checked), and the only
# way to compute it is the 9-column edge table, which for these networks was
# deleted -- 70 GB whose sole consumer had been mcxload. Writing NA is honest;
# writing a number from a different network would not be.
#
# RUN: through run.sh  ->  ./run.sh nodemetrics sugarcane
# =============================================================================
set -euo pipefail

STUDY="${CLEAN_STUDY:?}"
WORK="${CLEAN_WORK_DIR:?}"
BIN="${CLEAN_MCL_BIN_DIR:?}"
OUT="${CLEAN_OUT_FILE:?}"
CORES="${CLEAN_CORES:-16}"

MCI="${WORK}/${STUDY}.mci"
TAB="${WORK}/${STUDY}.tab"
say() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }

for f in "$MCI" "$TAB"; do
  [ -s "$f" ] || { echo "FATAL: missing $f -- run ./run.sh pearsonmci $STUDY" >&2; exit 1; }
done
mkdir -p "$(dirname "$OUT")"

N_TAB=$(wc -l < "$TAB")
say "$STUDY: $(printf "%'d" "$N_TAB") genes in $(basename "$TAB")"

# mcx query columns: node degree mean min max median iqr
TMP="${OUT}.part"
"$BIN/mcx" query -imx "$MCI" -t "$CORES" 2>/dev/null \
  | awk -v tab="$TAB" '
      BEGIN { FS = OFS = "\t"
              while ((getline line < tab) > 0) { split(line, a, "\t"); g[a[1]] = a[2] } }
      NR == 1 { next }                       # mcx header
      { idx = $1
        if (!(idx in g)) { print "FATAL: index " idx " not in tab file" > "/dev/stderr"; exit 1 }
        printf "%s\t%d\t%.6g\tNA\n", g[idx], $2, $2 * $3 }
    ' > "$TMP"

N_OUT=$(wc -l < "$TMP")
if [ "$N_OUT" != "$N_TAB" ]; then
  echo "FATAL: wrote $N_OUT rows for $N_TAB genes -- the matrix and the tab disagree" >&2
  rm -f "$TMP"; exit 1
fi

{ printf 'gene\tdegree\tstrength\ttransitivity\n'; cat "$TMP"; } > "$OUT"
rm -f "$TMP"

# --- global metrics, for the topology figure's legend -------------------------
# The figure quotes nodes/edges/density/components. Those must describe THIS
# network: the table left over from the merged (Pearson+MI) build says purple is a
# single connected component, which is true of that graph and not of this one.
# Components come from `clm close --write-sizes`, which reads the matrix directly.
if [ -n "${CLEAN_GLOBAL_OUT:-}" ]; then
  say "computing connected components"
  CC=$(mktemp); "$BIN/clm" close -imx "$MCI" --write-sizes -o "$CC" 2>/dev/null
  NCC=$(tr ' ' '\n' < "$CC" | grep -cE '^[0-9]+$')
  GIANT=$(tr ' ' '\n' < "$CC" | grep -E '^[0-9]+$' | sort -rn | head -1)
  rm -f "$CC"
  NODES=$N_TAB
  ARCS=$(awk -F'\t' 'NR>1{s+=$2} END{print s+0}' "$OUT")
  EDGES=$(( ARCS / 2 ))
  MEANW=$(awk -F'\t' -v a="$ARCS" 'NR>1{s+=$3} END{printf "%.6g", s/a}' "$OUT")
  DENS=$(awk -v e="$EDGES" -v n="$NODES" 'BEGIN{printf "%.7g", 2*e/(n*(n-1))}')
  { printf 'Metric\tValue\n'
    printf 'Study\t%s\n' "$STUDY"
    printf 'Nodes\t%s\n' "$NODES"
    printf 'Edges\t%s\n' "$EDGES"
    printf 'Edge_density\t%s\n' "$DENS"
    printf 'Global_transitivity\t\n'
    printf 'Mean_normalised_weight\t%s\n' "$MEANW"
    printf 'N_connected_components\t%s\n' "$NCC"
    printf 'Giant_component_size\t%s\n' "$GIANT"
    printf 'Edge_file\t%s\n' "$(basename "$MCI") (Pearson-only, |r| >= 0.8)"
  } > "$CLEAN_GLOBAL_OUT"
  say "wrote $(basename "$CLEAN_GLOBAL_OUT"): $NCC components, giant $(printf "%'d" "$GIANT"), density $DENS"
fi

say "wrote $(basename "$OUT"): $(printf "%'d" "$N_OUT") genes"
say "  degree/strength from mcx query; transitivity NA (no consumer reads it)"
head -3 "$OUT" | column -t
