# Persistence-image tests for the torus radius-shift experiment.
#
# Power comparison b uses:
#   null shift 0.00, block b, diagrams 1:50
#   alternative shift, block b, diagrams 1:50
#
# Validity uses 250 disjoint null comparisons:
#   blocks 1:250 versus blocks 251:500,
#   again using diagrams 1:50 in every block.

get_script_directory <- function() {
  command_line <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", command_line, value = TRUE)

  if (length(file_argument) > 0L) {
    script_path <- sub("^--file=", "", file_argument[[1L]])
    return(dirname(normalizePath(script_path, mustWork = FALSE)))
  }

  if (
    requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()
  ) {
    script_path <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(script_path)) {
      return(dirname(script_path))
    }
  }

  getwd()
}

setwd(get_script_directory())

script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)
source("functions.R")

# ---------------------------------------------------------------------------
# Configuration: these values agree with the Python notebook.
# ---------------------------------------------------------------------------
RADIUS_SHIFTS <- c(0.01, 0.016, 0.026, 0.036, 0.0)
ALTERNATIVE_SHIFTS <- RADIUS_SHIFTS[RADIUS_SHIFTS > 0]
NULL_SHIFT <- 0.00

NPC <- 50L
NSET <- 500L
NPAIR <- NSET %/% 2L
PI_RESOLUTION <- 40L
PI_RANGE <- c(0, 4)  # retained for the interface of ts_main_power()
ALPHA <- 0.05
C_VALUES <- c(0,0.2,0.4,0.6,0.8)
REPORTED_C <- 0.2

PI_INPUT_DIR <- file.path("PI_Torus_data", "torus_radius_shift")
RESULT_DIR <- file.path("..", "results", "simulation_results",
  "Torus_results"
)

dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

pi_file_for_shift <- function(radius_shift) {
  file.path(
    PI_INPUT_DIR,
    sprintf("torus_pi_linear_shift_%.3f.rds", radius_shift)
  )
}

load_shift_persistence_images <- function(radius_shift) {
  input_path <- pi_file_for_shift(radius_shift)
  if (!file.exists(input_path)) {
    stop(
      "Persistence-image file not found: ",
      normalizePath(input_path, mustWork = FALSE),
      ". Run Torus_data_PI.R first."
    )
  }

  record <- readRDS(input_path)
  required_fields <- c(
    "radius_shift", "persistence_images", "nset", "sample_size",
    "diagram_indices", "pi_resolution"
  )
  missing_fields <- setdiff(required_fields, names(record))
  if (length(missing_fields) > 0L) {
    stop(
      "PI file is missing metadata fields: ",
      paste(missing_fields, collapse = ", ")
    )
  }

  if (!isTRUE(all.equal(as.numeric(record$radius_shift), radius_shift))) {
    stop("The radius-shift metadata do not agree with the requested file.")
  }
  if (as.integer(record$nset) != NSET) {
    stop("The PI file does not contain the expected 500 blocks.")
  }
  if (as.integer(record$sample_size) != NPC) {
    stop("The PI file was not constructed with diagrams 1:50 per block.")
  }
  if (!identical(as.integer(record$diagram_indices), seq_len(NPC))) {
    stop("The PI file uses unexpected point-cloud indices.")
  }
  if (as.integer(record$pi_resolution) != PI_RESOLUTION) {
    stop("The PI resolution does not agree with PI_RESOLUTION.")
  }

  persistence_images <- record$persistence_images
  if (length(persistence_images) != NSET * NPC) {
    stop(
      "Expected ", NSET * NPC,
      " persistence images, but found ", length(persistence_images), "."
    )
  }

  valid_dimensions <- vapply(
    persistence_images,
    function(image) {
      is.matrix(image) &&
        identical(dim(image), c(PI_RESOLUTION, PI_RESOLUTION))
    },
    logical(1)
  )
  if (!all(valid_dimensions)) {
    stop("At least one stored persistence image has invalid dimensions.")
  }

  persistence_images
}

reported_c_index <- match(REPORTED_C, C_VALUES)
if (is.na(reported_c_index)) {
  stop("REPORTED_C must be one of the values in C_VALUES.")
}

# ---------------------------------------------------------------------------
# Power: the same replication block and the same first 50 diagrams used by
# the Python Bonferroni, PD, and PL tests.
# ---------------------------------------------------------------------------
message("Loading null persistence images.")
null_images <- load_shift_persistence_images(NULL_SHIFT)

power_rejections <- array(
  NA_integer_,
  dim = c(length(ALTERNATIVE_SHIFTS), NSET, length(C_VALUES)),
  dimnames = list(
    radius_shift = sprintf("%.3f", ALTERNATIVE_SHIFTS),
    replication = as.character(seq_len(NSET)),
    C = as.character(C_VALUES)
  )
)

for (shift_index in seq_along(ALTERNATIVE_SHIFTS)) {
  radius_shift <- ALTERNATIVE_SHIFTS[[shift_index]]
  message(sprintf("Power: radius_shift=%.3f", radius_shift))

  shifted_images <- load_shift_persistence_images(radius_shift)

  # ts_main_power() reads consecutive groups of NPC images.  Since both
  # lists were saved in replication-major order, replication b compares
  # null block b with shifted block b without any resampling.
  one_shift_result <- ts_main_power(
    onepi = list(null_images),
    twopi = list(shifted_images),
    sig = radius_shift,
    npc = NPC,
    nset = NSET,
    range = PI_RANGE,
    res = PI_RESOLUTION,
    alpha = ALPHA
  )

  power_rejections[shift_index, , ] <- one_shift_result[1L, , ]

  rm(shifted_images, one_shift_result)
  invisible(gc())
}

power_by_c <- apply(power_rejections, c(1L, 3L), mean)
power_at_reported_c <- power_by_c[, reported_c_index]

power_table <- data.frame(
  radius_shift = ALTERNATIVE_SHIFTS,
  empirical_power = as.numeric(power_at_reported_c),
  C = REPORTED_C,
  n = NPC,
  m = NPC,
  replications = NSET
)

print(power_table)

# This compact file preserves the format of the previous analysis: one power
# value per positive radius shift, evaluated at C = 0.2.
saveRDS(
  setNames(as.numeric(power_at_reported_c), sprintf("%.3f", ALTERNATIVE_SHIFTS)),
  file.path(RESULT_DIR, "Torus_PI_power_result.rds")
)

# ---------------------------------------------------------------------------
# Validity: form 250 disjoint pairs from the 500 null blocks.
#
# Pair 1 uses blocks 1 and 251, ..., pair 250 uses blocks 250 and 500.
# Calling ts_main_power() here is intentional: it applies exactly the same
# pixelwise test, low-variance handling, lower-triangle restriction, and BH
# correction as in the power calculation.  The older ts_main_fpr() call did
# not exist in the supplied functions file.
# ---------------------------------------------------------------------------
first_null_blocks_zero_based <- 0:(NPAIR - 1L)
second_null_blocks_zero_based <- NPAIR:(2L * NPAIR - 1L)

flat_indices_for_blocks <- function(blocks_zero_based) {
  as.vector(vapply(
    blocks_zero_based,
    function(block_index) {
      as.integer(block_index * NPC + seq_len(NPC))
    },
    integer(NPC)
  ))
}

validity_x <- null_images[
  flat_indices_for_blocks(first_null_blocks_zero_based)
]
validity_y <- null_images[
  flat_indices_for_blocks(second_null_blocks_zero_based)
]

message("Validity: 250 disjoint null-block pairs.")
validity_rejections <- ts_main_power(
  onepi = list(validity_x),
  twopi = list(validity_y),
  sig = NULL_SHIFT,
  npc = NPC,
  nset = NPAIR,
  range = PI_RANGE,
  res = PI_RESOLUTION,
  alpha = ALPHA
)

validity_matrix <- validity_rejections[1L, , , drop = TRUE]
validity_by_c <- colMeans(validity_matrix)
names(validity_by_c) <- as.character(C_VALUES)
type1_error_at_reported_c <- unname(validity_by_c[[reported_c_index]])

validity_table <- data.frame(
  C = C_VALUES,
  empirical_type1_error = as.numeric(validity_by_c),
  alpha = ALPHA,
  n = NPC,
  m = NPC,
  null_pairs = NPAIR
)

print(validity_table)

saveRDS(
  type1_error_at_reported_c,
  file.path(RESULT_DIR, "Torus_PI_validity_result.rds")
)


message("PI power and validity simulations are complete.")
message("Results saved under: ", normalizePath(RESULT_DIR))

#readRDS("../../results/simulation_results/Torus_results/torus_radius_shift/Torus_PI_power_result.rds")
