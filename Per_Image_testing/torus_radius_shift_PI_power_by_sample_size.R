# Persistence-image power for the torus radius-shift experiment.
# Each stored block has 50 images. Power uses the first n images from the
# matching null and shifted blocks, for n = 15, 25, and 40.
# The separate Type-I error experiment still uses n = m = 50.

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
      return(dirname(normalizePath(script_path, mustWork = FALSE)))
    }
  }

  getwd()
}

setwd(get_script_directory())
source("functions.R")

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
RADIUS_SHIFTS <- c(0.01, 0.016, 0.026, 0.036, 0.0)
ALTERNATIVE_SHIFTS <- RADIUS_SHIFTS[RADIUS_SHIFTS > 0]
NULL_SHIFT <- 0.0

BLOCK_SIZE <- 50L
POWER_SAMPLE_SIZES <- c(15L, 25L, 40L)
VALIDITY_SAMPLE_SIZE <- BLOCK_SIZE
NSET <- 500L
NPAIR <- NSET %/% 2L
CHECKPOINT_BATCH <- 25L
PI_RESOLUTION <- 40L
PI_RANGE <- c(0, 4)  # Retained for ts_main_power()'s interface.
ALPHA <- 0.05
C_VALUES <- c(0, 0.2, 0.4, 0.6, 0.8)
REPORTED_C <- 0.2

PI_INPUT_DIR <- file.path("PI_Torus_data", "torus_radius_shift")
RESULT_DIR <- file.path(
  "..", "results", "simulation_results", "Torus_results"
)
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

if (NSET %% 2L != 0L || any(POWER_SAMPLE_SIZES < 1L) ||
    any(POWER_SAMPLE_SIZES > BLOCK_SIZE)) {
  stop("Check NSET, BLOCK_SIZE, and POWER_SAMPLE_SIZES.")
}
reported_c_index <- match(REPORTED_C, C_VALUES)
if (is.na(reported_c_index)) {
  stop("REPORTED_C must occur in C_VALUES.")
}

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
    stop("PI file is missing fields: ", paste(missing_fields, collapse = ", "))
  }

  if (!isTRUE(all.equal(as.numeric(record$radius_shift), radius_shift))) {
    stop("Radius-shift metadata do not match the requested file.")
  }
  if (as.integer(record$nset) != NSET ||
      as.integer(record$sample_size) != BLOCK_SIZE ||
      !identical(as.integer(record$diagram_indices), seq_len(BLOCK_SIZE)) ||
      as.integer(record$pi_resolution) != PI_RESOLUTION) {
    stop("PI cache metadata do not match the expected 500 x 50 layout.")
  }

  images <- record$persistence_images
  if (length(images) != NSET * BLOCK_SIZE) {
    stop("Expected ", NSET * BLOCK_SIZE, " persistence images.")
  }
  dimensions_ok <- vapply(
    images,
    function(image) {
      is.matrix(image) &&
        identical(dim(image), c(PI_RESOLUTION, PI_RESOLUTION))
    },
    logical(1)
  )
  if (!all(dimensions_ok)) {
    stop("At least one persistence image has invalid dimensions.")
  }
  images
}

# Replication numbers are 1-based. The stored images are block-major,
# with BLOCK_SIZE images in each block. Ordering remains block-major after
# taking the first n images from every requested block.
images_for_blocks <- function(images, replications, n) {
  if (length(replications) == 0L) {
    return(list())
  }
  if (anyNA(replications) || any(replications < 1L) ||
      any(replications > NSET) || n < 1L || n > BLOCK_SIZE) {
    stop("Invalid replication number or within-block sample size.")
  }
  indices <- unlist(
    lapply(replications, function(b) (b - 1L) * BLOCK_SIZE + seq_len(n)),
    use.names = FALSE
  )
  images[indices]
}

power_path <- function(n) {
  file.path(RESULT_DIR, sprintf("Torus_PI_power_n%d_checkpoint.rds", n))
}

new_power_progress <- function() {
  array(
    NA_integer_,
    dim = c(length(ALTERNATIVE_SHIFTS), NSET, length(C_VALUES)),
    dimnames = list(
      radius_shift = sprintf("%.3f", ALTERNATIVE_SHIFTS),
      replication = as.character(seq_len(NSET)),
      C = as.character(C_VALUES)
    )
  )
}

load_power_progress <- function(path) {
  expected <- new_power_progress()
  if (!file.exists(path)) {
    return(expected)
  }
  result <- readRDS(path)
  if (!is.array(result) || !identical(dim(result), dim(expected)) ||
      !identical(dimnames(result), dimnames(expected)) ||
      any(!is.na(result) & !(result %in% 0:1))) {
    stop("Invalid or incompatible power checkpoint: ", path)
  }
  result
}

save_power_progress <- function(result, path) {
  temporary <- tempfile(pattern = "power_checkpoint_", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(result, temporary)
  if (!file.rename(temporary, path)) {
    stop("Could not replace checkpoint: ", path)
  }
}

# ---------------------------------------------------------------------------
# Power: same block b on the null and shifted sides, with nested first-n
# subsets. Save after each batch of at most 25 replications.
# ---------------------------------------------------------------------------
message("Loading null persistence images.")
null_images <- load_shift_persistence_images(NULL_SHIFT)
power_tables <- vector("list", length(POWER_SAMPLE_SIZES))

for (size_index in seq_along(POWER_SAMPLE_SIZES)) {
  n <- POWER_SAMPLE_SIZES[[size_index]]
  path <- power_path(n)
  power_rejections <- load_power_progress(path)
  started <- Sys.time()

  for (shift_index in seq_along(ALTERNATIVE_SHIFTS)) {
    radius_shift <- ALTERNATIVE_SHIFTS[[shift_index]]
    shifted_images <- load_shift_persistence_images(radius_shift)
    message(sprintf("Power: n=m=%d, radius_shift=%.3f", n, radius_shift))

    for (batch_start in seq.int(1L, NSET, by = CHECKPOINT_BATCH)) {
      batch_end <- min(batch_start + CHECKPOINT_BATCH - 1L, NSET)
      replications <- seq.int(batch_start, batch_end)
      missing <- replications[vapply(
        replications,
        function(b) anyNA(power_rejections[shift_index, b, ]),
        logical(1)
      )]

      if (length(missing) > 0L) {
        null_subset <- images_for_blocks(null_images, missing, n)
        shifted_subset <- images_for_blocks(shifted_images, missing, n)
        batch_result <- ts_main_power(
          onepi = list(null_subset),
          twopi = list(shifted_subset),
          sig = radius_shift,
          npc = n,
          nset = length(missing),
          range = PI_RANGE,
          res = PI_RESOLUTION,
          alpha = ALPHA
        )

        expected_shape <- c(1L, length(missing), length(C_VALUES))
        if (!identical(dim(batch_result), expected_shape)) {
          stop("Unexpected ts_main_power() output shape for n=", n)
        }
        values <- matrix(
          as.integer(batch_result),
          nrow = length(missing),
          ncol = length(C_VALUES)
        )
        if (anyNA(values) || any(!(values %in% 0:1))) {
          stop("The PI test returned a missing or invalid rejection value.")
        }

        power_rejections[shift_index, missing, ] <- values
        save_power_progress(power_rejections, path)
      }

      if (batch_end %% 50L == 0L || batch_end == NSET) {
        message(sprintf(
          "n=m=%d, shift=%.3f: %d/%d, elapsed=%.1f min",
          n, radius_shift, batch_end, NSET,
          as.numeric(difftime(Sys.time(), started, units = "mins"))
        ))
      }
    }

    rm(shifted_images)
    invisible(gc())
  }

  if (anyNA(power_rejections)) {
    stop("Power experiment did not finish for n=m=", n)
  }
  save_power_progress(power_rejections, path)

  power_by_c <- apply(power_rejections, c(1L, 3L), mean)
  power_at_reported_c <- power_by_c[, reported_c_index]
  power_tables[[size_index]] <- data.frame(
    radius_shift = ALTERNATIVE_SHIFTS,
    empirical_power = as.numeric(power_at_reported_c),
    C = REPORTED_C,
    n = n,
    m = n,
    replications = NSET
  )

  saveRDS(
    setNames(
      as.numeric(power_at_reported_c),
      sprintf("%.3f", ALTERNATIVE_SHIFTS)
    ),
    file.path(RESULT_DIR, sprintf("Torus_PI_power_n%d_result.rds", n))
  )
}

power_table <- do.call(rbind, power_tables)
row.names(power_table) <- NULL
print(power_table)
saveRDS(
  power_table,
  file.path(RESULT_DIR, "Torus_PI_power_by_sample_size.rds")
)

