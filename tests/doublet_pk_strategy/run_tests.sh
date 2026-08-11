#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(
  cd "$(dirname "${BASH_SOURCE[0]}")/../.." &&
  pwd
)"
cd "$PROJECT_ROOT"

Rscript - <<'RS'
expressions <- parse(file = "scripts/01g_doublet_detection_scoring.R")
targets <- c(
  "resolve_sweep_size",
  "allocate_stratified_quota",
  "resolve_pk"
)

for (expression in expressions) {
  if (
    is.call(expression) &&
      identical(as.character(expression[[1]]), "<-") &&
      is.symbol(expression[[2]]) &&
      as.character(expression[[2]]) %in% targets
  ) {
    eval(expression, envir = .GlobalEnv)
  }
}

stopifnot(
  resolve_sweep_size(7785L, 2000L, 0.5, 5000L) == 3892L,
  resolve_sweep_size(20000L, 2000L, 0.5, 5000L) == 5000L,
  resolve_sweep_size(2500L, 2000L, 0.5, 5000L) == 2000L,
  resolve_sweep_size(1500L, 2000L, 0.5, 5000L) == 1500L
)

fixed <- resolve_pk(
  obj = NULL,
  pcs_use = 1:10,
  pN = 0.25,
  sweep_cores = 1L,
  pk_mode = "fixed",
  fixed_pk = 0.09,
  sweep_min_cells = 2000L,
  sweep_fraction = 0.5,
  sweep_max_cells = 5000L,
  sweep_seed_primary = 1234L,
  sweep_seed_secondary = 4321L,
  pk_tolerance = 0.01,
  sweep_scale_factor = 1.5
)
stopifnot(
  fixed$best_pk == 0.09,
  fixed$source == "configured_fixed_pk",
  fixed$attempts == 0L,
  fixed$sweep_cells == 0L
)

sizes <- c(alpha = 100L, beta = 50L, gamma = 25L)
quota <- allocate_stratified_quota(sizes, 80L)
stopifnot(
  identical(names(quota), names(sizes)),
  sum(quota) == 80L,
  all(quota >= 1L),
  all(quota <= sizes),
  quota[["alpha"]] > quota[["beta"]],
  quota[["beta"]] > quota[["gamma"]]
)

cat("Adaptive pK helper contracts: VALID\n")
RS

grep -q 'get_config_value(config, "doublet.pk_mode", "auto")' \
  scripts/01g_doublet_detection_scoring.R
grep -q '"auto", "sweep", "fixed"' \
  scripts/01g_doublet_detection_scoring.R

if grep -q 'pk_mode' configs/project_config_doublet_test.yaml; then
  echo "ERROR: legacy fixture unexpectedly requires pk_mode." >&2
  exit 1
fi

echo "Adaptive pK backwards-compatibility contract: VALID"
