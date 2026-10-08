# Build persistence images for the torus radius-shift experiment.
#
# The input RDS has the following structure:
#   torus_pd[[radius shift]][[replication]][[point cloud]]
#
# To match the Python experiment exactly, only point clouds 1:50 from
# each of the 500 replication blocks are used.  Persistence images are
# saved separately for each radius shift to keep peak memory manageable.

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
NPC_STORED <- 100L
SAMPLE_SIZE <- 50L
NSET <- 500L

PI_RESOLUTION <- 40L
PI_RANGE_X <- c(0, 2)
PI_RANGE_Y <- c(0, 2)
PI_BANDWIDTH <- 0.2
WEIGHT_POWER <- 2 

PD_RDS_PATH <- file.path(
  "..", "data", "Rdata",
  "npc=100_nset=500_torus_radius_combined.rds"
)
PI_OUTPUT_DIR <- file.path("PI_Torus_data", "torus_radius_shift")
PROGRESS_EVERY <- 25L

dir.create(PI_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

pi_file_for_shift <- function(radius_shift) {
  file.path(
    PI_OUTPUT_DIR,
    sprintf("torus_pi_linear_shift_%.3f.rds", radius_shift)
  )
}

# The R code that generated the data used
# unname(as.numeric(torus_diag$diagram)).  Since R flattened the diagram
# column by column, the inverse operation is t(matrix(x, nrow = 3,
# byrow = TRUE)).
decode_h1_birth_persistence <- function(flat_diagram) {
  flat_diagram <- as.numeric(flat_diagram)

  if (length(flat_diagram) == 0L) {
    return(data.frame(birth = numeric(0), death = numeric(0)))
  }

  if (length(flat_diagram) %% 3L != 0L) {
    stop("A serialized persistence diagram has length not divisible by 3.")
  }

  diagram <- t(matrix(flat_diagram, nrow = 3L, byrow = TRUE))
  colnames(diagram) <- c("dimension", "birth", "death")

  h1 <- diagram[diagram[, "dimension"] == 1, , drop = FALSE]
  if (nrow(h1) == 0L) {
    return(data.frame(birth = numeric(0), death = numeric(0)))
  }

  birth_persistence <- data.frame(
    birth = h1[, "birth"],
    death = h1[, "death"] - h1[, "birth"]
  )

  if (any(!is.finite(as.matrix(birth_persistence)))) {
    stop("An H1 persistence diagram contains a non-finite value.")
  }
  if (any(birth_persistence$death < 0)) {
    stop("An H1 persistence point has negative persistence.")
  }

  birth_persistence
}

make_linear_persistence_image <- function(flat_diagram) {
  pd <- decode_h1_birth_persistence(flat_diagram)

  if (nrow(pd) == 0L) {
    return(matrix(0, nrow = PI_RESOLUTION, ncol = PI_RESOLUTION))
  }

  weights <- polyweight(pd, n = WEIGHT_POWER)
  persistence_image <- pers.image(
    pd = pd,
    rangex = PI_RANGE_X,
    rangey = PI_RANGE_Y,
    wgt = weights,
    nbins = PI_RESOLUTION,
    h = PI_BANDWIDTH
  )

  if (
    !identical(dim(persistence_image), c(PI_RESOLUTION, PI_RESOLUTION)) ||
      any(!is.finite(persistence_image))
  ) {
    stop("Persistence-image construction returned an invalid matrix.")
  }

  persistence_image
}

if (!file.exists(PD_RDS_PATH)) {
  stop("Persistence-diagram RDS not found: ", normalizePath(
    PD_RDS_PATH,
    mustWork = FALSE
  ))
}

message("Loading persistence diagrams from: ", normalizePath(PD_RDS_PATH))
torus_pd <- readRDS(PD_RDS_PATH)

if (length(torus_pd) != length(RADIUS_SHIFTS)) {
  stop(
    "Expected ", length(RADIUS_SHIFTS),
    " radius-shift groups, but found ", length(torus_pd), "."
  )
}

for (shift_index in seq_along(RADIUS_SHIFTS)) {
  radius_shift <- RADIUS_SHIFTS[[shift_index]]
  simulation_blocks <- torus_pd[[shift_index]]

  if (length(simulation_blocks) != NSET) {
    stop(
      sprintf(
        "radius_shift=%.3f: expected %d blocks, but found %d.",
        radius_shift, NSET, length(simulation_blocks)
      )
    )
  }

  # Flat, replication-major order:
  # block 1 diagrams 1:50, block 2 diagrams 1:50, ..., block 500.
  # This is the ordering expected by ts_main_power().
  persistence_images <- vector("list", NSET * SAMPLE_SIZE)

  for (replication in seq_len(NSET)) {
    raw_diagrams <- simulation_blocks[[replication]]

    if (length(raw_diagrams) != NPC_STORED) {
      stop(
        sprintf(
          paste0(
            "radius_shift=%.3f, replication=%d: expected %d stored ",
            "diagrams, but found %d."
          ),
          radius_shift, replication, NPC_STORED, length(raw_diagrams)
        )
      )
    }

    for (diagram_index in seq_len(SAMPLE_SIZE)) {
      flat_index <- (replication - 1L) * SAMPLE_SIZE + diagram_index
      persistence_images[[flat_index]] <- make_linear_persistence_image(
        raw_diagrams[[diagram_index]]
      )
    }

    if (replication %% PROGRESS_EVERY == 0L || replication == NSET) {
      message(sprintf(
        "radius_shift=%.3f: %d/%d blocks completed",
        radius_shift, replication, NSET
      ))
    }
  }

  output_record <- list(
    radius_shift = radius_shift,
    persistence_images = persistence_images,
    nset = NSET,
    npc_stored = NPC_STORED,
    sample_size = SAMPLE_SIZE,
    diagram_indices = seq_len(SAMPLE_SIZE),
    pi_resolution = PI_RESOLUTION,
    pi_range_x = PI_RANGE_X,
    pi_range_y = PI_RANGE_Y,
    pi_bandwidth = PI_BANDWIDTH,
    weight = "persistence^(1/2)"
  )

  output_path <- pi_file_for_shift(radius_shift)
  saveRDS(output_record, output_path)
  message("Saved: ", normalizePath(output_path))

  rm(persistence_images, output_record, simulation_blocks)
  invisible(gc())
}

message("All radius-shift persistence-image files are complete.")
