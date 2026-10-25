# Simulation study of the paper; usage and options in README.md.

#------ Packages and paths ----------------------------------------------------
repos <- getOption("repos")
if (is.null(repos) || isTRUE(repos["CRAN"] == "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(Matrix, ggplot2, spdep, here)
source(here("script", "functions.R"))
install_inla()

#------ Command-line arguments ------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1) args[1] else "pilot"

#------ Run -------------------------------------------------------------------
if (mode == "check") {
  check_engine()
} else if (mode == "report") {
  preset <- if (length(args) >= 2) args[2] else "full"
  make_report(readRDS(here("output", preset, "config.rds")))
} else {
  only <- if (length(args) >= 2 && nzchar(args[2])) {
    strsplit(args[2], ",")[[1]]
  } else NULL
  cfg <- sim_config(mode, out_dir = here("output", mode))
  if (length(args) >= 3 && args[3] == "brazil") {
    cfg$graphs <- c(cfg$graphs, "brazil")
  }
  ensure_dir(cfg$out_dir)
  writeLines(utils::capture.output(utils::sessionInfo()),
             file.path(cfg$out_dir, "session_info.txt"))
  saveRDS(cfg, file.path(cfg$out_dir, "config.rds"))
  run_experiments(cfg, only)
  make_report(cfg)
  msg("finished")
}
