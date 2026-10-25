# Functions for the simulation study and the application of the paper.

# ---- Simulation settings: pilot (minutes) or full (paper design) -------------
sim_config <- function(preset = c("pilot", "full"), out_dir = "output") {
  preset <- match.arg(preset)
  cfg <- list(
    preset = preset, out_dir = out_dir, seed = 20261012L,
    cores = max(1L, parallel::detectCores() - 1L),
    beta = 0.2, alpha = 0.05,
    budgets = c(0, 0.05, 0.1, 0.2, 0.5, 1, 2),
    kappas = c(-1, -0.5, 0, 0.5, 1),
    c0_grid = c(0, 0.05, 0.1, 0.2, 0.5, 1, 2),
    c0_fine = c(0, 0.025, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1, 1.5, 2),
    levels = c(0.5, 2, 10, 50), shares = c(0.2, 0.5, 0.8),
    decays = c(1, 2), level_main = 10, share_main = 0.5, decay_main = 2,
    trunc_fixed = c(5L, 20L, 50L), trunc_grid = c(0L, 5L, 10L, 20L, 40L),
    phi_bym2 = 0.5, phi_grid = c(0.1, 0.5, 0.9),
    graphs = c("lattice", "scotland", "delaunay"),
    lattice_side = 20L, delaunay_n = 1000L,
    scotland_graph = NULL, use_inla = FALSE,
    brazil_rdata = "Municipios_Brasil.RData",
    brazil_exposure = "tx_analfa", brazil_impute = "complete",
    brazil_join = FALSE, sigma_grid = exp(seq(log(0.02), log(20),
                                                      length.out = 60)),
    benchmark_mult = c(0.25, 0.5, 1, 2, 4)
  )
  if (preset == "pilot") {
    cfg <- utils::modifyList(cfg, list(
      designs = 3L, reps = 40L,
      designs_full = 3L,
      reps_full = c(lattice = 3L, scotland = 3L, delaunay = 0L, brazil = 0L),
      reps_alt = 20L, graphs = c("lattice", "scotland"),
      levels = c(2, 10), shares = c(0.2, 0.8), decays = 2,
      kappas = c(-1, 0, 1), c0_grid = c(0, 0.1, 0.5, 2),
      designs_brazil = 2L, reps_brazil = 20L))
  } else {
    cfg <- utils::modifyList(cfg, list(
      designs = 50L, reps = 100L,
      designs_full = 30L,
      reps_full = c(lattice = 50L, scotland = 50L, delaunay = 10L,
                    brazil = 0L),
      reps_alt = 10L, designs_brazil = 50L, reps_brazil = 100L))
  }
  cfg
}

# ---- Brazilian municipalities: data, covariates and graph --------------------
numeric_cols <- c(
  "numero_de_homicidios", "pop_total", "pop_homens", "prop_agua_canal",
  "prop_dest_lixo", "prop_esg_rede", "iie", "tx_analfa", "prop_preto_pardo",
  "razao_rend_branca_pp", "prop_forca_ocup", "tx_estab_mil_hab",
  "tx_creche_10mil_hab", "repasse_bf_pcapita", "valor_med_repas_por_fam")

raw_covariates <- c(
  "tx_analfa", "prop_preto_pardo", "razao_rend_branca_pp", "prop_agua_canal",
  "prop_dest_lixo", "prop_esg_rede", "iie", "prop_forca_ocup",
  "tx_estab_mil_hab", "repasse_bf_pcapita", "tx_creche_10mil_hab",
  "valor_med_repas_por_fam", "prop_homens")

covariates <- c(
  "log_population_density", "tx_analfa", "prop_preto_pardo",
  "razao_rend_branca_pp", "infrastructure_index", "prop_forca_ocup",
  "log_tx_estab", "log_repasse_bf_pcapita", "log_tx_creche")

# Outcome of the application: homicides, unless the option sp.cause holds
# another cause (see script/Application_mortality.R).
homicide_cause <- list(noun = "homicide", Noun = "Homicides", dir = "brazil",
                       file = NULL)
cause <- function() getOption("sp.cause", homicide_cause)

covariate_labels <- c(
  log_population_density = "Population density (log)",
  tx_analfa = "Illiteracy rate (%)",
  prop_preto_pardo = "Black/Brown population (%)",
  razao_rend_branca_pp = "White to Black/Brown income ratio",
  infrastructure_index = "Infrastructure index",
  prop_forca_ocup = "Occupied labour force (%)",
  log_tx_estab = "Health establishments (log)",
  log_repasse_bf_pcapita = "Bolsa Família transfer (log)",
  log_tx_creche = "Childcare facilities (log)")

covariates_raw <- c("tx_analfa", "prop_preto_pardo", "razao_rend_branca_pp",
                 "prop_agua_canal", "prop_dest_lixo", "prop_esg_rede",
                 "prop_forca_ocup", "tx_estab_mil_hab", "repasse_bf_pcapita",
                 "tx_creche_10mil_hab")

coerce_numeric <- function(x) {
  out <- suppressWarnings(as.numeric(as.character(x)))
  out[trimws(as.character(x)) %in% c("", "NA", "NaN", "NULL", ".")] <- NA
  out
}

# Data file: path as given, or under the project root or its data/ folder.
locate_file <- function(path) {
  tried <- c(path, here::here(path), here::here("data", path))
  found <- tried[file.exists(tried)]
  if (!length(found)) {
    stop("Data file not found: ", path, "\n  searched: ",
         paste(tried, collapse = ", "))
  }
  normalizePath(found[1])
}

read_municipalities <- function(path, map_crs = 4674L) {
  path <- locate_file(path)
  env <- new.env(parent = baseenv())
  load(path, envir = env)
  muni <- env$muni_sf
  for (v in numeric_cols) muni[[v]] <- coerce_numeric(muni[[v]])
  muni <- sf::st_transform(sf::st_make_valid(muni), map_crs)
  muni[!sf::st_is_empty(muni), ]
}

region_labels <- c("Centro Oeste" = "Centre-West", Nordeste = "Northeast",
                   Norte = "North", Sudeste = "Southeast", Sul = "South")

region_factor <- function(x) {
  x <- as.character(x)
  hit <- x %in% names(region_labels)
  x[hit] <- region_labels[x[hit]]
  droplevels(factor(x, levels = unique(c(unname(region_labels), x))))
}

select_columns <- function(muni) {
  out <- dplyr::transmute(
    muni,
    municipality_id = as.character(code_muni),
    municipality_name = as.character(name_muni),
    state = as.factor(abbrev_state), region = region_factor(name_region),
    deaths = as.integer(round(count)), expected_ext, pop_out, q_ill,
    q_nonspec, population = pop_total, pop_homens = pop_homens,
    prop_agua_canal, prop_dest_lixo, prop_esg_rede, iie, tx_analfa,
    prop_preto_pardo, razao_rend_branca_pp, prop_forca_ocup,
    tx_estab_mil_hab, tx_creche_10mil_hab, repasse_bf_pcapita,
    valor_med_repas_por_fam)
  dplyr::filter(out, is.finite(deaths), deaths >= 0,
                is.finite(population), population > 0)
}

add_derived <- function(dat, area_crs = 5880L) {
  area <- as.numeric(sf::st_area(sf::st_transform(dat, area_crs))) / 1e6
  dat$area_km2 <- area
  dat$population_density <- dat$population / area
  dat$log_population_density <- log1p(dat$population_density)
  dat$prop_homens <- 100 * dat$pop_homens / dat$population
  den <- ifelse(is.finite(dat$pop_out), dat$pop_out, dat$population)
  dat$rate_100k <- 1e5 * dat$deaths / den
  dat$size_class <- cut(dat$population, c(0, 29999, 300000, Inf),
                  labels = c("< 30,000", "30,000 to 300,000", "> 300,000"))
  crude <- dat$population * sum(dat$deaths) / sum(dat$population)
  dat$expected <- if (all(is.finite(dat$expected_ext))) dat$expected_ext else
    crude
  dat$smr <- dat$deaths / dat$expected
  dat
}

# Counts of the active cause; those of another cause (and its expected
# counts under national age and sex rates) come from a file.
attach_cause <- function(muni) {
  muni$count <- muni$numero_de_homicidios
  muni$expected_ext <- NA_real_
  muni$pop_out <- muni$pop_total
  muni$q_ill <- muni$q_nonspec <- NA_real_
  file <- cause()$file
  if (is.null(file)) return(muni)
  ext <- readRDS(locate_file(file))
  hit <- match(as.character(muni$code_muni), ext$code_muni)
  col <- function(name) ext[[name]][hit]
  muni$count <- col(paste0("deaths_", cause()$key))
  muni$expected_ext <- col(paste0("expected_", cause()$key))
  muni$pop_out <- col(paste0("pop_", cause()$key))
  muni$q_ill <- col("q_ill")
  muni$q_nonspec <- col("q_nonspec")
  muni
}

neighbours_from_adj <- function(adj) {
  a <- as.matrix(adj)
  lapply(seq_len(nrow(a)), function(i) which(a[i, ] > 0))
}

spatial_median_impute <- function(x, nb) {
  x <- as.numeric(x)
  fallback <- stats::median(x, na.rm = TRUE)
  for (i in which(!is.finite(x))) {
    v <- x[nb[[i]]]
    v <- v[is.finite(v)]
    x[i] <- if (length(v)) stats::median(v) else fallback
  }
  x
}

infrastructure_index <- function(dat) {
  cols <- c("prop_agua_canal", "prop_dest_lixo", "prop_esg_rede")
  comps <- sf::st_drop_geometry(dat)[, cols]
  pc <- stats::prcomp(scale(comps), center = FALSE, scale. = FALSE)$x[, 1]
  if (stats::cor(pc, rowMeans(comps)) < 0) -pc else pc
}

missing_table <- function(dat) {
  d <- sf::st_drop_geometry(dat)[covariates_raw]
  n <- vapply(d, function(v) sum(!is.finite(v)), numeric(1))
  data.frame(Variable = covariates_raw, Missing = n,
             Percent = sprintf("%.2f", 100 * n / nrow(d)))
}

complete_rows <- function(dat) {
  d <- sf::st_drop_geometry(dat)[covariates_raw]
  Reduce(`&`, lapply(d, is.finite))
}

fill_missing <- function(dat, nb, impute) {
  for (v in raw_covariates) {
    if (impute == "spatial") {
      dat[[v]] <- spatial_median_impute(dat[[v]], nb)
    } else if (impute == "global") {
      dat[[v]][!is.finite(dat[[v]])] <- stats::median(dat[[v]], na.rm = TRUE)
    }
  }
  dat
}

add_covariates <- function(dat, nb, impute) {
  dat <- fill_missing(dat, nb, impute)
  dat$infrastructure_index <- infrastructure_index(dat)
  dat$log_tx_estab <- log1p(dat$tx_estab_mil_hab)
  dat$log_repasse_bf_pcapita <- log1p(dat$repasse_bf_pcapita)
  dat$log_tx_creche <- log1p(dat$tx_creche_10mil_hab)
  for (v in covariates) dat[[paste0(v, "_z")]] <- standardize(dat[[v]])
  dat$id_spatial <- seq_len(nrow(dat))
  dat
}

# Municipalities kept at each step of the selection.
sample_flow <- function(n_file, n_valid, n_kept, n_quality = NULL,
                        max_q = Inf) {
  steps <- c("Brazilian municipalities",
             paste("With valid", cause()$noun, "count and population"))
  n <- c(n_file, n_valid)
  if (!is.null(n_quality)) {
    steps <- c(steps, sprintf("Ill-defined causes at most %g%% of deaths",
                              100 * max_q))
    n <- c(n, n_quality)
  }
  data.frame(Step = c(steps, "With every covariate (analysis sample)"),
             Municipalities = c(n, n_kept))
}

prepare_brazil <- function(path, queen = TRUE, impute = "complete",
                           join = FALSE, max_q = Inf) {
  impute <- match.arg(impute, c("complete", "global", "spatial"))
  raw <- attach_cause(read_municipalities(path))
  dat <- add_derived(select_columns(raw))
  n_valid <- nrow(dat)
  if (is.finite(max_q)) dat <- dat[which(dat$q_ill <= max_q), ]
  n_quality <- if (is.finite(max_q)) nrow(dat)
  adj <- polygon_adjacency(sf::st_transform(dat, 5880L), queen, join)
  miss <- missing_table(dat)
  keep <- if (impute == "complete") complete_rows(dat) else
    rep(TRUE, nrow(dat))
  flow <- sample_flow(nrow(raw), n_valid, sum(keep), n_quality, max_q)
  dat <- dat[keep, ]
  adj <- adj[keep, keep]
  nb <- if (impute == "spatial") neighbours_from_adj(adj) else NULL
  list(data = add_covariates(dat, nb, impute), adj = adj, missing = miss,
       flow = flow, impute = impute, join = join)
}

# Intercept, states, other covariates and component indicators (full rank).
brazil_controls <- function(dat, exposure = NULL, components = NULL) {
  keep <- setdiff(covariates, exposure)
  form <- stats::reformulate(c("state", paste0(keep, "_z")))
  xmat <- stats::model.matrix(form, data = sf::st_drop_geometry(dat))
  if (is.null(components)) return(xmat)
  colnames(components) <- sprintf("comp_%d", seq_len(ncol(components)))
  both <- cbind(xmat, components)
  q <- qr(both, tol = 1e-8)
  both[, sort(q$pivot[seq_len(q$rank)]), drop = FALSE]
}

# Cached, because the eigendecomposition is slow.
brazil_spectrum <- function(adj, cache = NULL) {
  if (!is.null(cache) && file.exists(cache)) return(readRDS(cache))
  msg("eigendecomposition of the Brazilian Laplacian")
  spec <- graph_spectrum(adj)
  if (!is.null(cache)) {
    ensure_dir(dirname(cache))
    saveRDS(spec, cache)
  }
  spec
}

spectrum_cache <- function(dir, prep) {
  file.path(dir, sprintf("brazil_spec_v3_%s_%s.rds", prep$impute,
                         if (prep$join) "joined" else "components"))
}

brazil_spec <- function(cfg) {
  prep <- prepare_brazil(cfg$brazil_rdata, impute = cfg$brazil_impute,
                         join = cfg$brazil_join)
  dat <- prep$data
  spec <- brazil_spectrum(prep$adj, spectrum_cache(cfg$out_dir, prep))
  cent <- sf::st_coordinates(suppressWarnings(sf::st_point_on_surface(dat)))
  expo <- cfg$brazil_exposure
  c(list(adj = prep$adj, coords = cent, name = "brazil", label = "brazil"),
    spec,
    list(control = brazil_controls(dat, expo, spec$null_basis),
         cache = new.env(),
         base = list(o = log(dat$expected),
                     x = dat[[paste0(expo, "_z")]])))
}

# ---- Figures of the application ----------------------------------------------
budget_label <- expression("Radius " * c / c[ref])

covariate_factor <- function(v) {
  factor(v, levels = covariates, labels = covariate_labels)
}

pseudo_log_x <- function(breaks = c(0, 0.25, 1, 4, 8)) {
  ggplot2::scale_x_continuous(
    trans = scales::pseudo_log_trans(sigma = 0.1, base = 10),
    breaks = breaks)
}

facet_covariates <- function() {
  ggplot2::facet_wrap(~ covariate, scales = "free_y",
                      labeller = ggplot2::label_wrap_gen(24))
}

# Dashed and dotted lines: field energy at the BYM2 and REML scales.
fig_app_curves <- function(res, curves, reference, out_dir) {
  curves$covariate <- covariate_factor(curves$covariate)
  curves[c("est", "lo", "hi")] <- exp(curves[c("est", "lo", "hi")])
  marks <- do.call(rbind, lapply(res, function(r) {
    data.frame(covariate = covariate_factor(r$covariate),
               x = c(r$energy_fit, r$energy_reml) / r$ref,
               scale = c("fit", "reml"))
  }))
  p <- ggplot2::ggplot(curves, ggplot2::aes(mult, est)) +
    ggplot2::geom_hline(yintercept = 1, colour = "grey40") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.25,
                         fill = paper_colours[1]) +
    ggplot2::geom_line(colour = paper_colours[1]) +
    ggplot2::geom_vline(data = marks, colour = paper_colours[2],
                        ggplot2::aes(xintercept = x, linetype = scale)) +
    ggplot2::scale_linetype_manual(
      values = c("dashed", "dotted"), name = NULL,
      labels = expression(hat(sigma)[BYM2], hat(sigma)[REML])) +
    pseudo_log_x() + facet_covariates() +
    ggplot2::labs(x = budget_label,
                  y = expression(exp(beta) ~ "per standard deviation")) +
    theme_paper(9)
  if (!is.null(reference)) p <- add_reference(p, reference, max(curves$mult))
  save_figure(p, "app_curves", list(out_dir = out_dir), 10, 7)
}

# Rate ratios per standard deviation, by covariate and model.
fig_app_effects <- function(irr, out_dir) {
  irr$covariate <- factor(covariate_factor(irr$covariate),
                          levels = rev(covariate_labels))
  irr$method <- factor(irr$method, levels = unique(irr$method))
  p <- ggplot2::ggplot(irr, ggplot2::aes(irr, covariate, colour = method,
                                         shape = method)) +
    ggplot2::geom_vline(xintercept = 1, colour = "grey40") +
    ggplot2::geom_pointrange(ggplot2::aes(xmin = lo, xmax = hi), size = 0.3,
                             position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(values = paper_colours[1:3], name = NULL) +
    ggplot2::scale_shape_discrete(name = NULL) +
    ggplot2::labs(x = expression(exp(beta) ~ "per standard deviation"),
                  y = NULL) +
    theme_paper(9) + ggplot2::theme(legend.position = "bottom")
  save_figure(p, "app_effects", list(out_dir = out_dir), 7, 4.5)
}

add_reference <- function(p, ref, xmax) {
  ref$covariate <- covariate_factor(ref$covariate)
  methods <- unique(ref$method)
  ref$x <- xmax * (1.15 + 0.2 * (match(ref$method, methods) - 1))
  p + ggplot2::geom_pointrange(
    data = ref, inherit.aes = FALSE, size = 0.3,
    ggplot2::aes(x = x, y = irr, ymin = lo, ymax = hi, shape = method)) +
    ggplot2::scale_shape_discrete(name = NULL)
}

# Information fraction and worst-case risk gain of the fitted spatial weights
# over the non-spatial ones. Ticks: benchmark radii; dashed line: break-even
# radius.
fig_app_risk <- function(res, out_dir) {
  mult <- exp(seq(log(0.05), log(16), length.out = 80))
  curve <- do.call(rbind, lapply(res, function(r) {
    ct <- mult * r$ref / sqrt(r$phi)
    gain <- vapply(ct, function(cc) {
      risk_gain(r$geo, r$coords$ns, r$coords$sp, cc)
    }, numeric(1))
    data.frame(covariate = r$covariate, mult = rep(mult, 2),
               value = c(nsif(r$geo, ct), gain),
               measure = rep(c("NSIF", "G"), each = length(mult)))
  }))
  star <- do.call(rbind, lapply(res, function(r) {
    data.frame(covariate = r$covariate, x = r$cstar / r$ref)
  }))
  ticks <- do.call(rbind, lapply(res, function(r) {
    data.frame(covariate = r$covariate, x = r$bench$budget / r$ref)
  }))
  curve$covariate <- covariate_factor(curve$covariate)
  star$covariate <- covariate_factor(star$covariate)
  ticks$covariate <- covariate_factor(ticks$covariate)
  curve$measure <- factor(curve$measure, levels = c("NSIF", "G"))
  p <- ggplot2::ggplot(curve, ggplot2::aes(mult, value, colour = measure)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
    ggplot2::geom_line() +
    ggplot2::geom_vline(data = star, ggplot2::aes(xintercept = x),
                        linetype = "dashed", colour = "grey30",
                        na.rm = TRUE) +
    ggplot2::geom_rug(data = ticks, ggplot2::aes(x = x), inherit.aes = FALSE,
                      sides = "b") +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(values = paper_colours[1:2], name = NULL,
                                 labels = expression(NSIF(c), G(c))) +
    ggplot2::facet_wrap(~ covariate, labeller = ggplot2::label_wrap_gen(24)) +
    ggplot2::labs(x = budget_label, y = NULL) +
    theme_paper(9) + ggplot2::theme(legend.position = "bottom")
  save_figure(p, "app_risk", list(out_dir = out_dir), 10, 7)
}

# Decision frontier for one exposure: worst-case risk of the non-spatial,
# BYM2 and minimax weights relative to V0, and the variance and squared
# bias bound of the first two. Dashed line: break-even radius.
fig_app_frontier <- function(res, out_dir, covariate = names(res)[1]) {
  r <- res[[covariate]]
  mult <- exp(seq(log(0.01), log(16), length.out = 120))
  ct <- mult * r$ref / sqrt(r$phi)
  v0 <- sum(r$coords$ns^2)
  wc <- function(coord) vapply(ct, function(cc) {
    worst_case_risk(r$geo, coord, cc)
  }, numeric(1)) / v0
  risk <- data.frame(
    mult = rep(mult, 3), value = c(wc(r$coords$ns), wc(r$coords$sp),
                                   minimax_risk_curve(r$geo, ct) / v0),
    weights = rep(c("Non-spatial", "BYM2 weights", "Minimax"),
                  each = length(mult)))
  comp <- function(coord, label) {
    data.frame(mult = mult, variance = sum(coord^2) / v0,
               bias = ct^2 * sum(r$geo$s * coord^2) / v0, weights = label)
  }
  parts <- rbind(comp(r$coords$ns, "Non-spatial"),
                 comp(r$coords$sp, "BYM2 weights"))
  parts <- stats::reshape(parts, direction = "long",
                          varying = c("variance", "bias"), v.names = "value",
                          timevar = "part", times = c("Variance", "Bias"),
                          idvar = c("mult", "weights"))
  risk$panel <- "Worst-case risk"
  parts$panel <- "Variance and squared bias bound"
  d <- rbind(risk[c("mult", "value", "weights", "panel")],
             transform(parts, weights = paste(weights, part))[
               c("mult", "value", "weights", "panel")])
  p <- ggplot2::ggplot(d, ggplot2::aes(mult, value, colour = weights)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = r$cstar / r$ref, linetype = "dashed",
                        colour = "grey30") +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    scale_colour_paper(name = NULL) +
    ggplot2::facet_wrap(~ panel, scales = "free_y") +
    ggplot2::labs(x = budget_label, y = expression("Relative to " * V[0])) +
    theme_paper(9) + ggplot2::theme(legend.position = "bottom")
  save_figure(p, "app_frontier", list(out_dir = out_dir), 10, 4)
}

fig_app_bench <- function(bench, out_dir) {
  bench$covariate <- covariate_factor(bench$covariate)
  bench$benchmark <- covariate_factor(bench$benchmark)
  bench$ratio <- bench$ccrit / bench$budget
  bench$label <- ifelse(is.finite(bench$ratio), sprintf("%.1f", bench$ratio),
                        "")
  p <- ggplot2::ggplot(bench, ggplot2::aes(benchmark, covariate,
                       fill = log10(pmax(ratio, 1e-3)))) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 2.6) +
    ggplot2::scale_fill_gradient2(
      low = paper_colours[1], mid = "white", high = paper_colours[2],
      midpoint = 0, name = expression(log[10](c[crit] / c^(k)))) +
    ggplot2::labs(x = "Benchmark covariate", y = "Exposure") +
    theme_paper(9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))
  save_figure(p, "app_benchmark", list(out_dir = out_dir), 8, 5.5)
}

# Maps: one blue-orange palette, legend and annotations inside the frame.
rr_labels <- c("<0.25", "0.25-0.5", "0.5-0.8", "0.8-1.25", "1.25-2", "2-4",
               ">=4")

rr_class <- function(x) {
  cut(x, c(0, 0.25, 0.5, 0.8, 1.25, 2, 4, Inf), labels = rr_labels,
      right = FALSE, include.lowest = TRUE)
}

diverging_colours <- function(n) {
  grDevices::colorRampPalette(c(paper_colours[1], "white",
                                paper_colours[2]))(n)
}

class_scale <- function(name) {
  ggplot2::scale_fill_manual(
    values = diverging_colours(length(rr_labels)), drop = FALSE, name = name,
    aesthetics = c("colour", "fill"))
}

continuous_scale <- function(name, mid, limits = NULL, oob = scales::censor) {
  ggplot2::scale_fill_gradient2(
    low = paper_colours[1], mid = "white", high = paper_colours[2],
    midpoint = mid, limits = limits, name = name, oob = oob,
    aesthetics = c("colour", "fill"),
    guide = ggplot2::guide_colourbar(barheight = grid::unit(2.2, "cm"),
                                     barwidth = grid::unit(0.3, "cm")))
}

map_annotations <- function() {
  if (!requireNamespace("ggspatial", quietly = TRUE)) return(NULL)
  list(ggspatial::annotation_scale(location = "br", width_hint = 0.25,
                                   text_cex = 0.6, height = grid::unit(0.12,
                                                                       "cm")),
       ggspatial::annotation_north_arrow(
         location = "tr", which_north = "true", height = grid::unit(0.8, "cm"),
         width = grid::unit(0.8, "cm"),
         style = ggspatial::north_arrow_minimal(text_size = 6)))
}

map_layer <- function(dat, borders, var, scale) {
  ggplot2::ggplot(dat) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[[var]], colour = .data[[var]]),
                     linewidth = 0.05) +
    ggplot2::geom_sf(data = borders, fill = NA, colour = "grey25",
                     linewidth = 0.15) +
    scale + map_annotations() + ggplot2::coord_sf(expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = c(0.02, 0.02),
                   legend.justification = c(0, 0),
                   legend.title = ggplot2::element_text(size = 7),
                   legend.text = ggplot2::element_text(size = 6),
                   legend.key.size = grid::unit(0.3, "cm"),
                   plot.margin = ggplot2::margin(0, 0, 0, 0))
}

# Height that makes the panel fit the map, so that no margin is left.
map_height <- function(dat, width) {
  bb <- sf::st_bbox(dat)
  width * unname((bb["ymax"] - bb["ymin"]) / (bb["xmax"] - bb["xmin"]))
}

# P(RR > 1) from the posterior marginals of the linear predictor, which
# excludes the offset; the normal approximation is the fallback.
exceed_prob <- function(fit) {
  mg <- fit$marginals.linear.predictor
  if (is.null(mg)) {
    lp <- fit$summary.linear.predictor
    return(stats::pnorm(lp$mean / lp$sd))
  }
  vapply(mg, function(m) 1 - INLA::inla.pmarginal(0, m), numeric(1))
}

# Posterior mean rate ratio (the INLA fitted values exclude the offset E),
# probability that it exceeds one and spatial effect on the log scale.
map_data <- function(dat, fit) {
  dat$rr <- fit$summary.fitted.values$mean
  dat$exceed <- exceed_prob(fit)
  dat$spatial <- fit$summary.random$id_spatial$mean[seq_len(nrow(dat))]
  dat$smr_class <- rr_class(dat$smr)
  dat$rr_class <- rr_class(dat$rr)
  dat
}

# Equal-area projection for Brazil (SIRGAS 2000 / Brazil Polyconic).
map_crs <- 5880

# State outlines. The union is taken on a 1 m grid, so that slivers
# between municipal polygons do not leave interior lines, and holes are
# dropped.
state_borders <- function(dat) {
  geom <- sf::st_make_valid(sf::st_geometry(dat))
  sf::st_precision(geom) <- 1
  states <- sort(unique(dat$state))
  merged <- lapply(states, function(s) sf::st_union(geom[dat$state == s])[[1]])
  merged <- sf::st_cast(sf::st_sfc(merged, crs = sf::st_crs(dat)), "POLYGON")
  sf::st_sf(geometry = sf::st_sfc(
    lapply(merged, function(p) sf::st_polygon(p[1])), crs = sf::st_crs(dat)))
}

smr_label <- function() {
  if (is.null(cause()$file)) "Crude SMR" else "Standardised SMR"
}

fig_app_maps <- function(dat, fit, out_dir, width = 4) {
  dat <- sf::st_transform(map_data(dat, fit), map_crs)
  borders <- state_borders(dat)
  panels <- list(
    app_map_smr = list("smr_class", class_scale(smr_label())),
    app_map_rr = list("rr_class", class_scale(expression(hat(RR)))),
    app_map_exceed = list("exceed", continuous_scale(expression(P(RR > 1)), 0.5,
                                                     c(0, 1))),
    app_map_spatial = list("spatial", continuous_scale(expression(hat(u)), 0)))
  for (name in names(panels)) {
    p <- map_layer(dat, borders, panels[[name]][[1]], panels[[name]][[2]])
    save_figure(p, name, list(out_dir = out_dir), width,
                map_height(dat, width))
  }
}

# Context maps, same palette as above, centred at the median. The share of
# ill-defined causes is capped at 20%, the threshold of the sensitivity run.
fig_app_context <- function(dat, exposure, out_dir, width = 4) {
  dat <- sf::st_transform(dat, map_crs)
  borders <- state_borders(dat)
  dat$exposure_map <- dat[[exposure]]
  panels <- list(app_map_exposure = list(
    "exposure_map", continuous_scale(covariate_labels[[exposure]],
                                     stats::median(dat$exposure_map))))
  if (any(is.finite(dat$q_ill))) {
    dat$ill_map <- 100 * dat$q_ill
    panels$app_map_quality <- list(
      "ill_map", continuous_scale("Ill-defined causes (%)",
                                  stats::median(dat$ill_map, na.rm = TRUE),
                                  c(0, 20), scales::squish))
  }
  for (name in names(panels)) {
    p <- map_layer(dat, borders, panels[[name]][[1]], panels[[name]][[2]])
    save_figure(p, name, list(out_dir = out_dir), width,
                map_height(dat, width))
  }
}

# Standardised mortality ratio against the share of ill-defined causes.
fig_app_quality <- function(dat, out_dir) {
  d <- sf::st_drop_geometry(dat)
  d <- d[is.finite(d$q_ill) & d$deaths > 0, ]
  if (!nrow(d)) return(invisible(NULL))
  p <- ggplot2::ggplot(d, ggplot2::aes(100 * q_ill, smr)) +
    ggplot2::geom_point(ggplot2::aes(size = expected), alpha = 0.25,
                        colour = paper_colours[1]) +
    ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                         colour = paper_colours[2], linewidth = 0.7) +
    ggplot2::geom_vline(xintercept = 20, linetype = "dashed",
                        colour = "grey40") +
    ggplot2::scale_y_log10() +
    ggplot2::scale_size_area(max_size = 3, name = "Expected deaths") +
    ggplot2::labs(x = "Ill-defined causes among all deaths (%)",
                  y = smr_label()) +
    theme_paper(9) + ggplot2::theme(legend.position = "bottom")
  save_figure(p, "app_quality", list(out_dir = out_dir), 6, 4)
}

# ---- Bayesian fits of the application (INLA) ---------------------------------
sp_hyper <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                 phi = list(prior = "pc", param = c(0.5, 0.5)))

besag_hyper <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))

write_graph <- function(adj, path) {
  if (file.exists(path)) file.remove(path)
  INLA::inla.write.graph(adj, filename = path, mode = "ascii")
  INLA::inla.read.graph(path)
}

# Spatial+: each covariate minus its fitted Besag smooth. The outcome model
# keeps the full-rank spatial term, which can absorb part of the residual
# field, so this coefficient is a sensitivity check only.
spatial_plus_columns <- function(dat, graph, threads) {
  form <- x ~ f(id_spatial, model = "besag", graph = graph,
                scale.model = TRUE, constr = TRUE, hyper = besag_hyper)
  for (v in paste0(covariates, "_z")) {
    fit <- INLA::inla(form, family = "gaussian", num.threads = threads,
                      data = data.frame(x = dat[[v]],
                                        id_spatial = dat$id_spatial))
    smooth <- fit$summary.random$id_spatial$mean
    dat[[paste0(v, "_spplus")]] <- standardize(dat[[v]] - smooth)
  }
  dat
}

# Component indicators not already spanned by the intercept and the state
# effects, so that the fixed part matches the controls of the sensitivity
# analysis.
add_component_terms <- function(dat, adj) {
  comp <- graph_components(as.matrix(adj))
  if (max(comp) == 1L) return(list(data = dat, terms = character()))
  base <- stats::model.matrix(~ state, data = dat)
  ind <- outer(comp, seq_len(max(comp))[-1L], "==") * 1
  colnames(ind) <- sprintf("comp_%d", seq_len(max(comp))[-1L])
  rank_of <- function(m) qr(m, tol = 1e-8)$rank
  keep <- vapply(seq_len(ncol(ind)), function(j) {
    rank_of(cbind(base, ind[, seq_len(j), drop = FALSE])) >
      rank_of(cbind(base, ind[, seq_len(j - 1L), drop = FALSE]))
  }, logical(1))
  ind <- ind[, keep, drop = FALSE]
  list(data = cbind(dat, as.data.frame(ind)), terms = colnames(ind))
}

# Linear constraints that make the random effect orthogonal to the fixed part.
rsr_constraint <- function(dat, comps = character()) {
  design <- stats::model.matrix(
    stats::reformulate(c("state", paste0(covariates, "_z"), comps)),
    data = dat)[, -1, drop = FALSE]
  qrx <- qr(design, tol = 1e-10)
  q <- qr.Q(qrx)[, seq_len(qrx$rank), drop = FALSE]
  list(A = cbind(matrix(0, qrx$rank, nrow(dat)), t(q)),
       e = rep(0, qrx$rank))
}

# The formula is evaluated by INLA in this frame, where `graph` and `extra`
# are defined; `sp_hyper` is found in the enclosing environment.
model_formula <- function(method, extra, graph, comps = character()) {
  suffix <- if (method == "SpatialPlus") "_z_spplus" else "_z"
  term <- paste0("f(id_spatial, model = 'bym2', graph = graph, ",
                 "scale.model = TRUE, constr = TRUE, hyper = sp_hyper",
                 if (method == "RSR") ", extraconstr = extra", ")")
  rhs <- paste(c(paste0(covariates, suffix), comps, term), collapse = " + ")
  stats::as.formula(paste("deaths ~ 1 + state +", rhs),
                    env = environment())
}

fit_bym2 <- function(dat, graph, family, method, extra, threads,
                     comps = character()) {
  msg("INLA: ", family, ", ", method)
  INLA::inla(model_formula(method, extra, graph, comps), family = family,
             data = dat, E = dat$expected, num.threads = threads,
             control.inla = list(control.vb = list(enable = FALSE)),
             control.predictor = list(compute = TRUE),
             control.compute = list(waic = TRUE, dic = TRUE,
                                  return.marginals.predictor = TRUE))
}

try_bym2 <- function(...) {
  tryCatch(fit_bym2(...), error = function(e) {
    msg("fit failed: ", conditionMessage(e))
    NULL
  })
}

# Unrestricted fit and WAIC per likelihood; failed fits get an infinite WAIC.
screen_families <- function(dat, graph, candidates, threads,
                            comps = character()) {
  fits <- lapply(stats::setNames(candidates, candidates), function(f) {
    try_bym2(dat, graph, f, "Unrestricted", NULL, threads, comps)
  })
  waic <- vapply(fits, function(f) if (is.null(f)) Inf else f$waic$waic, 1)
  list(fits = fits, waic = waic)
}

# Posterior means of the BYM2 standard deviation and mixing parameter and,
# for the negative binomial, of its size.
bym2_hyper <- function(fit) {
  hp <- fit$marginals.hyperpar
  mean_of <- function(name, fun) INLA::inla.emarginal(fun, hp[[name]])
  out <- c(sigma = mean_of("Precision for id_spatial",
                           function(p) 1 / sqrt(p)),
           phi = mean_of("Phi for id_spatial", function(p) p))
  size <- grep("^size", names(hp), value = TRUE)
  if (length(size)) out["theta"] <- mean_of(size[1], function(p) p)
  out
}

# Rate ratios of the covariates, one row per covariate and method.
irr_rows <- function(fit, method) {
  fx <- fit$summary.fixed
  keep <- grepl("_z", rownames(fx))
  data.frame(
    covariate = sub("_z.*$", "", rownames(fx)[keep]), method = method,
    irr = exp(fx$mean[keep]), lo = exp(fx[["0.025quant"]][keep]),
    hi = exp(fx[["0.975quant"]][keep]), row.names = NULL)
}

# The sensitivity analysis needs a Poisson or negative binomial family
# (`families`); the zero-inflated ones (`extra`) enter the comparison and, if
# preferred by WAIC, the tables and maps.
run_inla_models <- function(dat, adj, out_dir, threads,
                            families = c("poisson", "nbinomial"),
                            extra = c("zeroinflatedpoisson1",
                                      "zeroinflatednbinomial1")) {
  graph <- write_graph(adj, file.path(out_dir, "brazil.graph"))
  dat <- sf::st_drop_geometry(spatial_plus_columns(dat, graph, threads))
  ct <- add_component_terms(dat, adj)
  scr <- screen_families(ct$data, graph, c(families, extra), threads,
                         ct$terms)
  main <- names(which.min(scr$waic))
  base <- families[which.min(scr$waic[families])]
  fit <- function(method, constraint = NULL) {
    try_bym2(ct$data, graph, main, method, constraint, threads, ct$terms)
  }
  fits <- list(Unrestricted = scr$fits[[main]],
               RSR = fit("RSR", rsr_constraint(ct$data, ct$terms)),
               `Spatial+-ICAR` = fit("SpatialPlus"))
  fits <- Filter(Negate(is.null), fits)
  list(family = base, main_family = main, waic = scr$waic, fits = fits,
       ref_fit = scr$fits[[base]],
       irr = do.call(rbind, Map(irr_rows, fits, names(fits))),
       irr_ref = irr_rows(scr$fits[[base]], "Unrestricted"),
       hyper = bym2_hyper(scr$fits[[base]]),
       checks = check_all(scr$fits, ct$data, families))
}

# Predictive checks of every candidate fit that has marginals.
check_all <- function(fits, dat, families) {
  rows <- lapply(families, function(f) {
    if (is.null(fits[[f]])) return(NULL)
    predictive_checks(fits[[f]], dat$deaths, dat$expected, f)
  })
  do.call(rbind, rows)
}

likelihood_labels <- c(poisson = "Poisson", nbinomial = "Negative binomial",
                       zeroinflatedpoisson1 = "Zero-inflated Poisson",
                       zeroinflatednbinomial1 =
                         "Zero-inflated negative binomial")

# Predictive checks of a Poisson or negative binomial fit: observed and
# expected number of zeros and of counts at or above two high thresholds,
# averaging the probabilities over the posterior of the linear predictor,
# and the Pearson dispersion at the fitted means.
count_cdf <- function(family, size) {
  if (family == "poisson") {
    function(k, mu) stats::ppois(k, mu)
  } else function(k, mu) stats::pnbinom(k, size = size, mu = mu)
}

predictive_checks <- function(fit, y, e, family, probs = c(0.95, 0.99)) {
  hp <- fit$summary.hyperpar
  size <- if (family == "nbinomial") hp[grep("^size", rownames(hp))[1], "mean"]
  cdf <- count_cdf(family, size)
  expect <- function(fun) {
    sum(vapply(seq_along(y), function(i) {
      INLA::inla.emarginal(function(x) fun(e[i] * exp(x)),
                           fit$marginals.linear.predictor[[i]])
    }, numeric(1)))
  }
  ks <- unique(ceiling(stats::quantile(y[y > 0], probs, names = FALSE)))
  mu <- e * exp(fit$summary.linear.predictor$mean)
  out <- data.frame(
    statistic = c("Zeros", sprintf("Counts >= %d", ks), "Pearson dispersion"),
    observed = c(sum(y == 0), vapply(ks, function(k) sum(y >= k), 1),
                 sum((y - mu)^2 / (mu + mu^2 / ifelse(is.null(size), Inf,
                                                      size))) / length(y)),
    expected = c(expect(function(m) cdf(0, m)),
                 vapply(ks, function(k) expect(function(m) {
                   1 - cdf(k - 1, m)
                 }), 1), 1))
  out$family <- likelihood_labels[family]
  out
}

waic_table <- function(waic) {
  data.frame(Likelihood = unname(likelihood_labels[names(waic)]),
             WAIC = ifelse(is.finite(waic), sprintf("%.1f", waic), "--"))
}

irr_table <- function(irr) {
  irr$cell <- sprintf("%.3f (%.3f, %.3f)", irr$irr, irr$lo, irr$hi)
  wide <- stats::reshape(irr[c("covariate", "method", "cell")],
                         idvar = "covariate", timevar = "method",
                         direction = "wide")
  names(wide) <- sub("cell.", "", names(wide), fixed = TRUE)
  wide$covariate <- covariate_labels[wide$covariate]
  names(wide)[1] <- "Covariate"
  wide[match(covariate_labels, wide$Covariate), ]
}

check_table_ppc <- function(checks) {
  data.frame(Likelihood = checks$family, Statistic = checks$statistic,
             Observed = sprintf("%.1f", checks$observed),
             Expected = sprintf("%.1f", checks$expected))
}

# ---- Tables, figures and CSV files of the sensitivity analysis ---------------
write_sensitivity <- function(res, inla_out, cfg, alpha = 0.05) {
  out_dir <- cfg$out_dir
  curves <- curve_table(res)
  bench <- bench_table(res)
  utils::write.csv(curves, file.path(out_dir, "sensitivity_curves.csv"),
                   row.names = FALSE)
  utils::write.csv(bench, file.path(out_dir, "sensitivity_benchmarks.csv"),
                   row.names = FALSE)
  utils::write.csv(cost_table(res), file.path(out_dir,
                   "sensitivity_cost.csv"), row.names = FALSE)
  irr_ref <- if (is.null(inla_out)) NULL else inla_out$irr_ref
  write_table(primary_table(res, irr_ref, alpha), "table_sensitivity", cfg,
              paste("Rate ratios per standard deviation from the non-spatial",
                    "fit and the BYM2 fit by INLA, information fraction at",
                    "the median",
                    "benchmark radius $c_{\\mathrm{ref}}$, break-even radius",
                    "$c^{\\ast}$ of BYM2 against the non-spatial fit, number",
                    "of benchmark radii above $c^{\\ast}$, and critical",
                    "radius of the non-spatial fit."),
              "tab:sensitivity", size = wide_size)
  write_table(grid_table(res, alpha = alpha), "table_grid", cfg,
              paste("Rate ratios (95% intervals allowing for bias) at",
                    "multiples of $c_{\\mathrm{ref}}$ and risk gain $G$ of",
                    "the BYM2 weights over the non-spatial ones at",
                    "$c_{\\mathrm{ref}}$."), "tab:grid", size = wide_size)
  if (!is.null(irr_ref)) {
    write_table(check_table(res, irr_ref), "table_check", cfg,
                paste("Rate ratios from BYM2 (INLA) and from the one-step",
                      "minimax estimate at the fitted scale, and the",
                      "difference of their logarithms."), "tab:check")
  }
  lg <- length_grid(res, alpha = alpha)
  utils::write.csv(lg, file.path(out_dir, "interval_length.csv"),
                   row.names = FALSE)
  write_table(length_table(lg), "table_length", cfg,
              paste("Ratio of the radius of the shortest interval to the",
                    "minimax radius, efficiency of the shortest and of the",
                    "minimax interval against the lower bound on the",
                    "expected length at zero field, and extra length of the",
                    "minimax interval (%), all at $c_{\\mathrm{ref}}$."),
              "tab:length")
  fig_app_length(lg, out_dir)
  fig_app_risk(res, out_dir)
  main <- if (is.null(cfg$exposure)) names(res)[1] else cfg$exposure
  fig_app_frontier(res, out_dir, main)
  fig_app_spectrum(res, out_dir, main)
  utils::write.csv(spectrum_table(res), file.path(out_dir,
                   "exposure_spectrum.csv"), row.names = FALSE)
  fig_app_curves(res, curves, if (is.null(inla_out)) NULL else inla_out$irr,
                 out_dir)
  fig_app_bench(bench, out_dir)
}

# ---- Spectrum of the reduced exposure ----------------------------------------
# Share of the squared norm of the reduced exposure in the components of the
# class with the largest prior variance (largest eigenvalues of Sigma).
spectrum_share <- function(r) {
  e <- r$geo$xi[order(r$geo$s, decreasing = TRUE)]^2
  cumsum(e) / sum(e)
}

spectrum_table <- function(res, tops = c(0.01, 0.05, 0.1, 0.25)) {
  do.call(rbind, lapply(res, function(r) {
    cum <- spectrum_share(r)
    share <- cum[ceiling(tops * length(cum))]
    names(share) <- sprintf("top_%g_pct", 100 * tops)
    data.frame(covariate = r$covariate, as.list(share),
               nsif_ref = nsif(r$geo, r$ref / sqrt(r$phi)))
  }))
}

fig_app_spectrum <- function(res, out_dir, main) {
  d <- do.call(rbind, lapply(res, function(r) {
    cum <- spectrum_share(r)
    k <- unique(round(exp(seq(0, log(length(cum)), length.out = 200))))
    data.frame(covariate = r$covariate, frac = k / length(cum),
               share = cum[k], main = r$covariate == main)
  }))
  p <- ggplot2::ggplot(d, ggplot2::aes(frac, share, group = covariate)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         colour = "grey50") +
    ggplot2::geom_line(data = d[!d$main, ], colour = "grey70") +
    ggplot2::geom_line(data = d[d$main, ], colour = paper_colours[1],
                       linewidth = 1) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "Fraction of components, largest prior variance first",
                  y = "Cumulative share of the squared norm") +
    theme_paper(9)
  save_figure(p, "app_spectrum", list(out_dir = out_dir), 6, 4)
}

# ---- Length of the bias-aware interval ---------------------------------------
# Norm and bias per unit radius of the weights a_{c'}, in noise-normalised
# units; the row c' = 0 holds the non-spatial weights and the last row
# numerically approximates the limit c' -> Inf.
weight_frontier <- function(geo, c, span = 30, m = 601) {
  grid <- c(0, c * exp(seq(-span, span, length.out = m)))
  fr <- t(vapply(grid, function(cp) {
    co <- weights_coord(geo, cp)
    c(weights_norm(co), weights_bias(geo, co))
  }, numeric(2)))
  cbind(radius = grid, norm = fr[, 1], bias = fr[, 2])
}

# Minimum over c' of f(norm, bias): grid search, then a local refinement
# around the best grid point. Returns the value and the radius c'.
frontier_min <- function(geo, fr, f) {
  val <- f(fr[, "norm"], fr[, "bias"])
  k <- which.min(val)
  out <- list(value = val[k], radius = fr[k, "radius"])
  if (k > 1 && k < nrow(fr)) {
    step <- diff(log(fr[2:3, "radius"]))
    at <- function(l) {
      co <- weights_coord(geo, exp(l))
      f(weights_norm(co), weights_bias(geo, co))
    }
    o <- stats::optimize(at, log(out$radius) + c(-step, step))
    if (o$objective < out$value) out <- list(value = o$objective,
                                             radius = exp(o$minimum))
  }
  out
}

# Lower bound on the expected length of a confidence set that covers the
# class of radius c, at zero field (Donoho 1994; Armstrong and Kolesar 2018).
length_bound <- function(geo, fr, c, alpha = 0.05) {
  z <- stats::qnorm(1 - alpha)
  modulus <- function(d) {
    frontier_min(geo, fr, function(n, b) n * d + 2 * c * b)$value
  }
  integrand <- function(v) {
    vapply(v, function(x) modulus(2 * (z - x)), numeric(1)) * stats::dnorm(v)
  }
  (1 - alpha) * stats::integrate(integrand, -Inf, z)$value / stats::pnorm(z)
}

# Shortest interval of the form CI_a(c), the minimax one (c' = c) and the
# efficiency of both against the bound; c is in noise-normalised units.
length_study <- function(geo, c, alpha = 0.05) {
  fr <- weight_frontier(geo, c)
  half <- function(n, b) n * cva(c * b / n, alpha)
  co <- weights_coord(geo, c)
  hm <- half(weights_norm(co), weights_bias(geo, co))
  best <- frontier_min(geo, fr, half)
  hb <- min(best$value, hm)
  bound <- length_bound(geo, fr, c, alpha) / 2
  data.frame(ratio = best$radius / c, eff_best = bound / hb,
             eff_minimax = bound / hm, loss = hm / hb - 1)
}

length_grid <- function(res, multiples = c(0.25, 0.5, 1, 2, 4), alpha = 0.05) {
  rows <- lapply(res, function(r) {
    do.call(rbind, lapply(multiples, function(m) {
      cbind(covariate = r$covariate, mult = m,
            length_study(r$geo, m * r$ref / sqrt(r$phi), alpha))
    }))
  })
  do.call(rbind, rows)
}

length_table <- function(lg) {
  d <- lg[lg$mult == 1, ]
  out <- data.frame(unname(covariate_labels[d$covariate]),
                    sprintf("%.2f", d$ratio), sprintf("%.3f", d$eff_best),
                    sprintf("%.3f", d$eff_minimax),
                    sprintf("%.1f", 100 * d$loss))
  names(out) <- c("Covariate", "$c'/c$", "Shortest", "Minimax",
                  "Extra length (%)")
  out
}

fig_app_length <- function(lg, out_dir) {
  long <- rbind(
    data.frame(lg[c("covariate", "mult")], value = lg$eff_best,
               series = "Shortest interval", panel = "Efficiency"),
    data.frame(lg[c("covariate", "mult")], value = lg$eff_minimax,
               series = "Minimax weights", panel = "Efficiency"),
    data.frame(lg[c("covariate", "mult")], value = 100 * lg$loss,
               series = "Minimax weights", panel = "Extra length (%)"))
  p <- ggplot2::ggplot(long, ggplot2::aes(mult, value, colour = series,
                                          group = interaction(covariate,
                                                              series))) +
    ggplot2::geom_line(alpha = 0.6) + ggplot2::geom_point(size = 0.8) +
    ggplot2::scale_x_log10() + ggplot2::facet_wrap(~ panel, scales = "free_y") +
    scale_colour_paper(name = NULL) +
    ggplot2::labs(x = "Radius c as a multiple of the reference radius",
                  y = NULL) + theme_paper()
  save_figure(p, "interval_length", list(out_dir = out_dir), 7, 3.6)
}

# ---- Sensitivity analysis of the application ---------------------------------
app_multiples <- c(0, 0.125, 0.25, 0.5, 1, 2, 4, 8)

risk_gain <- function(geo, coord0, coord_s, ct) {
  1 - worst_case_risk(geo, coord_s, ct) / worst_case_risk(geo, coord0, ct)
}

# Estimate, interval, NSIF, critical radius and gain of the spatial weights
# `coord_s` over the non-spatial ones, at raw radii c; phi is the dispersion.
curve_details <- function(geo, beta_ns, ehat, phi, budgets, alpha = 0.05,
                          coord_s = NULL) {
  coord0 <- weights_coord(geo, 0)
  rows <- lapply(budgets, function(cc) {
    ct <- cc / sqrt(phi)
    coord <- weights_coord(geo, ct)
    est <- beta_ns + sum(coord * ehat)
    sd <- sqrt(phi) * weights_norm(coord)
    bias <- weights_bias(geo, coord)
    half <- bias_aware_halfwidth(sd, bias, cc, alpha)
    gain <- if (is.null(coord_s)) NA_real_ else
      risk_gain(geo, coord0, coord_s, ct)
    data.frame(budget = cc, est = est, sd = sd, lo = est - half,
               hi = est + half, nsif = nsif(geo, ct), gain = gain,
               ccrit = critical_budget(est, sd, bias, alpha))
  })
  do.call(rbind, rows)
}

# Field that covariate k would add if left out: its effect times the part of
# it that the other controls cannot absorb.
benchmark_field <- function(ctrl, k, gamma, w, v) {
  others <- ctrl[, colnames(ctrl) != k, drop = FALSE]
  gamma * weighted_residual(others, w, v)
}

# Radius of each other covariate, c^(k) = (delta' Q delta)^(1/2) with
# delta = U'z^(k), with the interval and critical radius there.
benchmarks <- function(dat, spec, covariate, ctrl, fit, geo, ehat, phi,
                       alpha = 0.05) {
  rows <- lapply(setdiff(covariates, covariate), function(k) {
    name <- paste0(k, "_z")
    gamma <- fit$coef[match(name, colnames(ctrl))]
    zk <- benchmark_field(ctrl, name, gamma, geo$w, dat[[name]])
    delta <- drop(crossprod(spec$vecs, zk))
    ck <- sqrt(sum(delta^2 / geo$qinv))
    coord <- weights_coord(geo, ck / sqrt(phi))
    bias <- weights_bias(geo, coord)
    est <- fit$coef[fit$p] + sum(coord * ehat)
    sd <- sqrt(phi) * weights_norm(coord)
    half <- bias_aware_halfwidth(sd, bias, ck, alpha)
    data.frame(benchmark = k, budget = ck,
               est = est, lo = est - half, hi = est + half,
               ccrit = critical_budget(est, sd, bias, alpha))
  })
  do.call(rbind, rows)
}

# Everything that depends on the dispersion: benchmarks, break-even radius
# against the non-spatial fit, and the curve over the grid.
phi_summary <- function(core, geo, dat, spec, ctrl, multiples, phi, alpha) {
  bench <- benchmarks(dat, spec, core$covariate, ctrl, core$fit, geo,
                      core$ehat, phi, alpha)
  ref <- stats::median(bench$budget)
  cstar <- break_even(geo, core$coords$ns, core$coords$sp) * sqrt(phi)
  curve <- curve_details(geo, core$beta, core$ehat, phi, multiples * ref,
                         alpha, core$coords$sp)
  curve$mult <- multiples
  list(phi = phi, ref = ref, cstar = cstar, bench = bench, curve = curve)
}

# `sigma` is the scale of the spatial fit in noise-normalised units (the
# REML scale when NULL).
sensitivity_one <- function(dat, spec, covariate, qinv, multiples,
                            alpha = 0.05, theta = Inf, sigma = NULL) {
  x <- dat[[paste0(covariate, "_z")]]
  ctrl <- brazil_controls(dat, covariate, spec$null_basis)
  fit <- fit_ns(dat$deaths, log(dat$expected), cbind(ctrl, x), theta)
  geo <- reduced_geometry(fit$mu, ctrl, x, spec$vecs, qinv, theta)
  ehat <- residual_coords(geo, fit$mu, dat$deaths)
  reml <- sigma_hat(geo, ehat, fit$phi)
  sigma <- if (is.null(sigma)) reml else sigma
  core <- list(covariate = covariate, fit = fit, beta = fit$coef[fit$p],
               ehat = ehat, sigma = sigma,
               coords = list(ns = weights_coord(geo, 0),
                             sp = weights_coord(geo, sigma)))
  out <- c(core, phi_summary(core, geo, dat, spec, ctrl, multiples, fit$phi,
                             alpha))
  root_m <- sqrt(ncol(spec$vecs))
  out$energy_fit <- sigma * sqrt(fit$phi) * root_m
  out$energy_reml <- reml * sqrt(fit$phi) * root_m
  out$geo <- geo[c("s", "xi", "xt")]
  invisible(gc())
  out
}

run_sensitivity <- function(dat, spec, qinv, multiples = app_multiples,
                            vars = covariates, alpha = 0.05, theta = Inf,
                            sigma = NULL) {
  lapply(stats::setNames(vars, vars), function(v) {
    msg("sensitivity: ", v)
    invisible(gc(reset = TRUE))
    secs <- system.time(r <- sensitivity_one(dat, spec, v, qinv, multiples,
                                             alpha, theta,
                                             sigma))[["elapsed"]]
    r$cost <- c(seconds = secs, peak_mb = sum(gc()[, 6]))
    r
  })
}

# Same analysis under other shapes of the class (REML scale for c*).
run_shapes <- function(dat, spec, phi_grid = c(0.1, 0.5, 0.9), alpha = 0.05,
                       theta = Inf) {
  shapes <- c(list(ICAR = shape_qinv(spec, "icar")),
              lapply(stats::setNames(phi_grid, sprintf("BYM2 (%g)", phi_grid)),
                     function(p) shape_qinv(spec, "bym2", p)))
  rows <- lapply(names(shapes), function(nm) {
    res <- run_sensitivity(dat, spec, shapes[[nm]], alpha = alpha,
                           theta = theta)
    do.call(rbind, lapply(res, function(r) data.frame(
      shape = nm, covariate = r$covariate,
      nsif = nsif(r$geo, r$ref / sqrt(r$phi)),
      cstar = r$cstar / r$ref, ccrit = r$curve$ccrit[1] / r$ref)))
  })
  do.call(rbind, rows)
}

cost_table <- function(res) {
  do.call(rbind, lapply(res, function(r) {
    data.frame(covariate = r$covariate, seconds = r$cost[["seconds"]],
               peak_mb = r$cost[["peak_mb"]])
  }))
}

curve_table <- function(res) {
  do.call(rbind, lapply(res, function(r) cbind(covariate = r$covariate,
                                               r$curve)))
}

bench_table <- function(res) {
  do.call(rbind, lapply(res, function(r) cbind(covariate = r$covariate,
                                               r$bench)))
}

# ---- Tables of the application -----------------------------------------------
descriptive_labels <- function() c(
  deaths = cause()$Noun,
  rate_100k = paste(sub("s$", "", cause()$Noun), "rate (per 100,000)"),
  smr = "Standardised mortality ratio",
  q_ill = "Ill-defined causes among all deaths (proportion)",
  q_nonspec = "Nonspecific codes among circulatory deaths (proportion)",
  population = "Population",
  population_density = "Population density (per km$^2$)",
  tx_analfa = "Illiteracy rate (%)",
  prop_preto_pardo = "Black/Brown population (%)",
  razao_rend_branca_pp = "White to Black/Brown income ratio",
  prop_agua_canal = "Piped water coverage (proportion)",
  prop_dest_lixo = "Waste collection coverage (proportion)",
  prop_esg_rede = "Sewer network coverage (proportion)",
  prop_forca_ocup = "Occupied labour force (%)",
  tx_estab_mil_hab = "Health establishments (per 1,000)",
  repasse_bf_pcapita = "Bolsa Família transfer (per capita)",
  tx_creche_10mil_hab = "Childcare facilities (per 10,000)")

wide_size <- "\\scriptsize\\setlength{\\tabcolsep}{2pt}"

# Median (first, third quartile), with fewer decimals for large values.
quartile_text <- function(v) {
  q <- stats::quantile(v, c(0.5, 0.25, 0.75), na.rm = TRUE)
  digits <- if (all(v == round(v)) || all(abs(q) >= 100)) 0L else
    if (all(abs(q) >= 10)) 1L else 2L
  num <- formatC(q, format = "f", digits = digits, big.mark = ",")
  sprintf("%s (%s, %s)", num[1], num[2], num[3])
}

descriptive_table <- function(dat) {
  d <- sf::st_drop_geometry(dat)
  groups <- c(list(Overall = seq_len(nrow(d))),
              split(seq_len(nrow(d)), d$region))
  labels <- descriptive_labels()
  labels <- labels[vapply(names(labels), function(v) any(is.finite(d[[v]])),
                          TRUE)]
  cols <- lapply(groups, function(i) {
    vapply(names(labels), function(v) quartile_text(d[[v]][i]), "")
  })
  out <- data.frame(Variable = unname(labels), cols, check.names = FALSE)
  out[nrow(out) + 1L, ] <- c("Municipalities", vapply(groups, length, 1L))
  out
}

# Variance inflation factor of each covariate given the others and the states.
vif_table <- function(dat) {
  d <- sf::st_drop_geometry(dat)
  d <- cbind(d[paste0(covariates, "_z")], state = d$state)
  vif <- vapply(paste0(covariates, "_z"), function(v) {
    r2 <- summary(stats::lm(stats::reformulate(".", v), data = d))$r.squared
    1 / (1 - r2)
  }, 1)
  data.frame(Covariate = unname(covariate_labels[covariates]),
             VIF = sprintf("%.2f", vif))
}

# Counts by population size class, and agreement with the national rates.
size_table <- function(dat) {
  d <- sf::st_drop_geometry(dat)
  rows <- lapply(split(d, d$size_class), function(g) {
    data.frame(Municipalities = nrow(g), Deaths = sum(g$deaths),
               `No deaths` = sum(g$deaths == 0),
               `Expected below 5` = sum(g$expected < 5),
               Dispersion = sprintf("%.2f", sum((g$deaths - g$expected)^2 /
                                                  g$expected) / nrow(g)),
               check.names = FALSE)
  })
  out <- do.call(rbind, rows)
  cbind(`Population` = rownames(out), out, row.names = NULL)
}

irr_text <- function(est, lo, hi) {
  sprintf("%.3f (%.3f, %.3f)", exp(est), exp(lo), exp(hi))
}

# Non-spatial rate ratio with its model-based interval.
ns_text <- function(r, alpha) {
  z <- stats::qnorm(1 - alpha / 2)
  se <- sqrt(r$fit$phi * r$fit$bread[r$fit$p, r$fit$p])
  irr_text(r$beta, r$beta - z * se, r$beta + z * se)
}

bym2_text <- function(irr, covariate) {
  if (is.null(irr)) return("--")
  i <- irr$covariate == covariate & irr$method == "Unrestricted"
  sprintf("%.3f (%.3f, %.3f)", irr$irr[i], irr$lo[i], irr$hi[i])
}

# Number of benchmarks above the break-even radius, as "k/K".
above_text <- function(r) {
  if (is.na(r$cstar)) return("--")
  sprintf("%d/%d", sum(r$bench$budget > r$cstar), nrow(r$bench))
}

ratio_text <- function(x) ifelse(is.na(x), "--", sprintf("%.2f", x))

cref_tex <- "c_{\\mathrm{ref}}"

primary_table <- function(res, irr = NULL, alpha = 0.05) {
  out <- do.call(rbind, lapply(res, function(r) {
    data.frame(
      Covariate = covariate_labels[[r$covariate]],
      `Non-spatial` = ns_text(r, alpha),
      BYM2 = bym2_text(irr, r$covariate),
      nsif = sprintf("%.3f", nsif(r$geo, r$ref / sqrt(r$phi))),
      cstar = ratio_text(r$cstar / r$ref), above = above_text(r),
      ccrit = ratio_text(r$curve$ccrit[r$curve$mult == 0] / r$ref),
      check.names = FALSE, row.names = NULL)
  }))
  names(out)[4:7] <- c(sprintf("$\\mathrm{NSIF}(%s)$", cref_tex),
                       sprintf("$c^{\\ast}/%s$", cref_tex),
                       "Benchmarks above $c^{\\ast}$",
                       sprintf("$c_{\\mathrm{crit}}^{\\mathrm{NS}}/%s$",
                               cref_tex))
  out
}

# Main exposure in its original scale, per `step` points.
percent_covariates <- c("tx_analfa", "prop_preto_pardo", "prop_forca_ocup")

unit_table <- function(res, irr, exposure, sd_x, step = 5, alpha = 0.05) {
  r <- res[[exposure]]
  k <- step / sd_x
  z <- stats::qnorm(1 - alpha / 2)
  se <- sqrt(r$fit$phi * r$fit$bread[r$fit$p, r$fit$p])
  cell <- function(est, lo, hi) {
    sprintf("%.3f (%.3f, %.3f)", exp(k * est), exp(k * lo), exp(k * hi))
  }
  rows <- list(c("Non-spatial", cell(r$beta, r$beta - z * se,
                                      r$beta + z * se)))
  if (!is.null(irr)) {
    i <- irr$covariate == exposure & irr$method == "Unrestricted"
    rows <- c(rows, list(c("BYM2 (INLA)", cell(log(irr$irr[i]),
                                                log(irr$lo[i]),
                                                log(irr$hi[i])))))
  }
  mult <- c(0.5, 1, 2)
  cd <- r$curve[match(mult, r$curve$mult), ]
  for (i in seq_along(mult)) {
    label <- sprintf("Bias-aware, radius %g times the reference", mult[i])
    rows <- c(rows, list(c(label, cell(cd$est[i], cd$lo[i], cd$hi[i]))))
  }
  out <- as.data.frame(do.call(rbind, rows))
  names(out) <- c("Estimate", sprintf("Rate ratio per %g points", step))
  out
}

# Estimates at multiples of the reference radius.
grid_table <- function(res, multiples = c(0.5, 1, 2), alpha = 0.05) {
  out <- do.call(rbind, lapply(res, function(r) {
    cd <- r$curve[match(multiples, r$curve$mult), ]
    cells <- lapply(seq_along(multiples), function(i) {
      irr_text(cd$est[i], cd$lo[i], cd$hi[i])
    })
    data.frame(Covariate = covariate_labels[[r$covariate]],
               ns_text(r, alpha), do.call(cbind, cells),
               gain = sprintf("%.2f", cd$gain[match(1, multiples)]),
               check.names = FALSE, row.names = NULL)
  }))
  names(out) <- c("Covariate", "Non-spatial",
                  sprintf("$c=%g\\,%s$", multiples, cref_tex),
                  sprintf("$G(%s)$", cref_tex))
  out
}

# BYM2 fitted by INLA against the one-step minimax estimate at its scale.
check_table <- function(res, irr) {
  rows <- lapply(res, function(r) {
    i <- irr$covariate == r$covariate & irr$method == "Unrestricted"
    est <- r$beta + sum(r$coords$sp * r$ehat)
    data.frame(Covariate = covariate_labels[[r$covariate]],
               `BYM2 (INLA)` = sprintf("%.3f", irr$irr[i]),
               `One-step` = sprintf("%.3f", exp(est)),
               Difference = sprintf("%.3f", log(irr$irr[i]) - est),
               check.names = FALSE)
  })
  do.call(rbind, rows)
}

shape_table <- function(shapes) {
  out <- data.frame(unname(covariate_labels[shapes$covariate]), shapes$shape,
                    sprintf("%.3f", shapes$nsif), ratio_text(shapes$cstar),
                    ratio_text(shapes$ccrit))
  names(out) <- c("Covariate", "Shape",
                  sprintf("$\\mathrm{NSIF}(%s)$", cref_tex),
                  sprintf("$c^{\\ast}/%s$", cref_tex),
                  sprintf("$c_{\\mathrm{crit}}^{\\mathrm{NS}}/%s$", cref_tex))
  out
}

# ---- Data generation: exposures, unmeasured fields and counts ----------------
draw_offsets <- function(n, level, sd = 0.5) {
  e <- stats::rnorm(n, sd = sd)
  log(level) + e - stats::median(e)
}

draw_exposure_parts <- function(spec, decay) {
  sp <- spec$vecs %*% (stats::rnorm(ncol(spec$vecs)) * spec$lambda^(-decay / 2))
  list(smooth = standardize(drop(sp)), noise = stats::rnorm(spec$n))
}

# Remove from v its projection on w and scale to unit root mean square.
orthogonal_part <- function(v, w) {
  w0 <- w - mean(w)
  v0 <- v - mean(v)
  standardize(v0 - sum(v0 * w0) / sum(w0^2) * w0)
}

# Exposure x = sqrt(s) x_sp + sqrt(1 - s) noise. With a field z and rho != 0
# the spatial part has sample correlation exactly rho with z.
combine_exposure <- function(parts, share, z = NULL, rho = 0) {
  smooth <- parts$smooth
  has_field <- !is.null(z) && stats::sd(z) > 1e-10
  if (has_field && rho != 0) {
    smooth <- rho * standardize(z) +
      sqrt(1 - rho^2) * orthogonal_part(smooth, z)
  }
  noise <- orthogonal_part(parts$noise, smooth)
  standardize(sqrt(share) * smooth + sqrt(1 - share) * noise)
}

design_geometry <- function(spec, qinv, o, x, beta, ctrl, theta = Inf) {
  reduced_geometry(exp(o + beta * x), ctrl, x, spec$vecs, qinv, theta)
}

# Weights of the estimator that a confounder is aligned against.
target_weights <- function(geo, target, c0) {
  kind <- sub(":.*$", "", target)
  arg <- suppressWarnings(as.numeric(sub("^[^:]*:?", "", target)))
  switch(kind,
    ns = geo$xt / sum(geo$xt^2),
    mm = weights_vector(geo, arg),
    `c0` = weights_vector(geo, c0),
    trunc = truncation_weights(geo, as.integer(arg)),
    stop("unknown target: ", target))
}

roughness_energy <- function(delta, qinv) sum(delta^2 / qinv)

# Aligned: field in the class of radius c0 at alignment kappa with weights a.
z_aligned <- function(spec, qinv, geo, a, c0, kappa) {
  va <- sqrt(qinv) * drop(crossprod(geo$tm, a))
  va <- va / sqrt(sum(va^2))
  vp <- stats::rnorm(length(va))
  vp <- vp - sum(vp * va) * va
  vp <- vp / sqrt(sum(vp^2))
  v <- kappa * va + sqrt(1 - kappa^2) * vp
  drop(spec$vecs %*% (c0 * sqrt(qinv) * v))
}

# Smooth: smooth random field of roughness energy c0^2.
z_smooth <- function(spec, qinv, c0) {
  delta <- sqrt(qinv) * stats::rnorm(length(qinv))
  delta <- delta * c0 / sqrt(roughness_energy(delta, qinv))
  drop(spec$vecs %*% delta)
}

# Radius of a field in the class of shape qinv: the smallest radius whose
# class contains it.
field_radius <- function(spec, qinv, z) {
  sqrt(sum(drop(crossprod(spec$vecs, z))^2 / qinv))
}

# True radius of a realized field: its radius in the fitted shape, which
# equals c0 only for the families generated with that energy.
true_budget <- function(spec, st, qinv, z) field_radius(spec, qinv, z)

# Orthogonal: smooth field with the aligned component removed, so that the
# reduced confounder is exactly orthogonal to x (spatial structure, no bias).
z_orthogonal <- function(spec, qinv, geo, c0) {
  delta <- sqrt(qinv) * stats::rnorm(length(qinv))
  dir <- drop(crossprod(geo$tm, geo$xt))
  delta <- delta - sum(dir * delta) / sum(dir^2) * dir
  delta <- delta * c0 / sqrt(roughness_energy(delta, qinv))
  drop(spec$vecs %*% delta)
}

# Band-limited: field in the span of the K smoothest eigenvectors,
# with ||u|| = c0.
z_band <- function(spec, geo, c0, n_vec) {
  delta <- numeric(ncol(spec$vecs))
  delta[seq_len(n_vec)] <- stats::rnorm(n_vec)
  delta <- delta * c0 / sqrt(sum((geo$tm %*% delta)^2))
  drop(spec$vecs %*% delta)
}

matern_cov <- function(d, nu, range) {
  r <- d / range
  if (nu == 0.5) return(exp(-r))
  if (nu == 1.5) return((1 + sqrt(3) * r) * exp(-sqrt(3) * r))
  stop("nu must be 0.5 or 1.5")
}

# Matern: Matern field on the centroids, rescaled to roughness energy c0^2.
z_matern <- function(spec, qinv, coords, c0, nu, range_frac) {
  if (is.null(coords) || nrow(coords) > 2000L) {
    stop("Matern fields need coordinates and n <= 2000")
  }
  d <- as.matrix(stats::dist(coords))
  k <- matern_cov(d, nu, range_frac * max(d))
  f <- drop(crossprod(chol(k + diag(1e-8, nrow(k))), stats::rnorm(nrow(k))))
  delta <- drop(crossprod(spec$vecs, f))
  delta <- delta * c0 / sqrt(roughness_energy(delta, qinv))
  drop(spec$vecs %*% delta)
}

# Seed of a localized field: a node of a component with at least `size`
# nodes (of the largest component if none is large enough).
local_seed <- function(spec, size) {
  cs <- table(spec$comp)
  ok <- as.integer(names(cs)[cs >= size])
  if (!length(ok)) ok <- as.integer(names(cs)[which.max(cs)])
  pool <- which(spec$comp %in% ok)
  pool[sample.int(length(pool), 1L)]
}

# Local: localized field of log-rate amplitude `amp`. "hotspot" raises the
# `size` nodes closest to a seed; "region" raises the nodes within the
# median graph distance of the seed (a graph-distance neighbourhood).
# Only nodes reachable from the seed are used.
z_local <- function(spec, adj, type, size, amp) {
  seed <- local_seed(spec, if (type == "hotspot") size else 1L)
  dist <- graph_distance(adj, seed)
  reach <- which(is.finite(dist))
  z <- numeric(spec$n)
  if (type == "hotspot") {
    z[reach[order(dist[reach])][seq_len(min(size, length(reach)))]] <- amp
  } else {
    z[reach[dist[reach] <= stats::median(dist[reach])]] <- amp
  }
  z <- z - mean(z)
  if (ncol(spec$null_basis)) {
    z <- z - spec$null_basis %*% crossprod(spec$null_basis, z)
  }
  drop(z)
}

z_mixture <- function(spec, qgen, st) {
  smooth <- z_smooth(spec, qgen, st$c0)
  local <- z_local(spec, spec$adj, st$kind, st$size, 1)
  e_smooth <- roughness_energy(crossprod(spec$vecs, smooth), qgen)
  e_local <- roughness_energy(crossprod(spec$vecs, local), qgen)
  smooth + local * sqrt(e_smooth / e_local)
}

make_z <- function(spec, geo, qfit, qgen, st) {
  switch(st$fam,
    aligned = z_aligned(spec, qfit, geo,
                   target_weights(geo, st$target, st$c0), st$c0, st$kappa),
    smooth = , independent = z_smooth(spec, qgen, st$c0),
    band = z_band(spec, geo, st$c0, st$n_vec),
    orthogonal = z_orthogonal(spec, qfit, geo, st$c0),
    matern = z_matern(spec, qgen, spec$coords, st$c0, st$nu, st$range),
    local = z_local(spec, spec$adj, st$kind, st$size, st$amp),
    mixture = z_mixture(spec, qgen, st),
    stop("unknown family: ", st$fam))
}

# Reduced confounder of a log-rate field z: vector in R^n and its
# coordinates in the eigenbasis of Sigma.
reduced_vector <- function(geo, z) {
  v <- geo$w * z
  drop(v - geo$q1 %*% crossprod(geo$q1, v))
}

reduced_confounder <- function(geo, z) {
  drop(crossprod(geo$E, reduced_vector(geo, z)))
}

draw_counts <- function(mu, family = "poisson", theta = 5) {
  if (family == "poisson") return(stats::rpois(length(mu), mu))
  stats::rnbinom(length(mu), size = theta, mu = mu)
}

# ---- One-step minimax adjustment and its intervals ---------------------------
residual_coords <- function(geo, mu_hat, y) {
  r <- geo$w * (y - mu_hat) / mu_hat
  r <- r - geo$q1 %*% crossprod(geo$q1, r)
  drop(crossprod(geo$E, r))
}

curve_row <- function(geo, ct, beta_ns, ehat, phi, alpha, ucoord, leak) {
  d <- 1 / (1 + ct^2 * geo$s)
  f <- sum(geo$xi^2 * d)
  coord <- geo$xi * d / f
  nrm <- sqrt(sum(coord^2))
  bias <- sqrt(sum(geo$s * coord^2))
  sd <- sqrt(phi) * nrm
  half <- sd * cva(ct * bias / nrm, alpha)
  est <- beta_ns + sum(coord * ehat)
  row <- stat_row(est, sd, est - half, est + half, nsif = f / sum(geo$xi^2),
                  pred_bias = if (is.null(ucoord)) NA else sum(coord * ucoord),
                  pred_sd = sd)
  if (!is.na(leak)) {
    # leak is measured on the field itself, noise-normalised here
    wide <- sd * cva(ct * bias / nrm + leak / sqrt(phi), alpha)
    row[c("lo_leak", "hi_leak")] <- c(est - wide, est + wide)
  }
  row
}

# Rows for a vector of normalised radii (radius / sqrt(phi)).
adjust_curve <- function(geo, beta_ns, ehat, phi, budgets, labels,
                         alpha = 0.05, ucoord = NULL, leak = NULL) {
  if (is.null(leak)) leak <- rep(NA_real_, length(budgets))
  rows <- lapply(seq_along(budgets), function(i) {
    curve_row(geo, budgets[i], beta_ns, ehat, phi, alpha, ucoord, leak[i])
  })
  out <- do.call(rbind, rows)
  rownames(out) <- labels
  out
}

# Reduced REML estimate of sigma in noise-normalised coordinates.
sigma_hat <- function(geo, ehat, phi) {
  reml_sigma(ehat / sqrt(phi), geo$xi / sqrt(phi), geo$s)
}

# Same interval with the bias term ignored (radius set to zero).
naive_row <- function(row, alpha = 0.05) {
  half <- stats::qnorm(1 - alpha / 2) * row["sd"]
  row[c("lo", "hi")] <- c(row["est"] - half, row["est"] + half)
  row[c("lo_leak", "hi_leak")] <- NA
  row
}

# ---- Penalized spatial fits: ICAR, BYM2, RSR and Spatial+-ICAR ---------------
backsolve_chol <- function(info_mat, g) {
  chol_fac <- chol(info_mat)
  backsolve(chol_fac, backsolve(chol_fac, g, transpose = TRUE))
}

line_search <- function(obj, th, step, cur) {
  for (k in 0:25) {
    cand <- th + step * 0.5^k
    val <- obj(cand)
    if (is.finite(val) && val >= cur - 1e-10) return(list(th = cand, val = val))
  }
  list(th = th, val = cur)
}

# Penalized iteratively reweighted least squares; the first np columns of
# Z are unpenalized and the rest carry the precision matrix P.
pirls <- function(y, o, zmat, np, prec_mat, start, maxit = 40, tol = 1e-9) {
  idx <- np + seq_len(ncol(zmat) - np)
  obj <- function(th) {
    eta <- o + drop(zmat %*% th)
    sum(y * eta - exp(eta)) - 0.5 * sum(th[idx] * drop(prec_mat %*% th[idx]))
  }
  info <- function(th) {
    mu <- exp(pmin(o + drop(zmat %*% th), 30))
    info_mat <- crossprod(zmat, zmat * mu)
    info_mat[idx, idx] <- info_mat[idx, idx] + prec_mat
    list(mu = mu, info_mat = info_mat)
  }
  th <- start
  cur <- obj(th)
  for (it in seq_len(maxit)) {
    s <- info(th)
    g <- crossprod(zmat, y - s$mu)
    g[idx] <- g[idx] - prec_mat %*% th[idx]
    new <- line_search(obj, th, backsolve_chol(s$info_mat, drop(g)), cur)
    done <- abs(new$val - cur) < tol * (abs(cur) + 1)
    th <- new$th
    cur <- new$val
    if (done) break
  }
  c(list(theta = th, obj = cur, idx = idx), info(th))
}

# Laplace approximation of the log marginal likelihood. The term
# 0.5 * log|prec| does not depend on sigma but does change with the shape
# (e.g. the BYM2 mixing parameter), so it is kept for comparisons.
laplace_logml <- function(fit, m, sigma, logdet_prec) {
  chol_fac <- chol(fit$info_mat)
  fit$obj + 0.5 * logdet_prec - m * log(sigma) - sum(log(diag(chol_fac)))
}

# Posterior summary for the coefficient of column `xcol` of Z.
penalized_row <- function(fit, xcol, sigma, alpha = 0.05) {
  chol_fac <- chol(fit$info_mat)
  cov_mat <- chol2inv(chol_fac)
  wz <- crossprod(fit$zmat, fit$zmat * fit$mu)
  edf <- sum(cov_mat * wz)
  stat_row(fit$theta[xcol], sqrt(cov_mat[xcol, xcol]), extra = sigma,
           edf = edf, alpha = alpha)
}

# Fit with sigma estimated by the Laplace marginal likelihood. `prec` is
# the (scaled) precision of the penalized coefficients at sigma = 1.
fit_penalized <- function(y, o, xfix, basis, prec, alpha = 0.05) {
  np <- ncol(xfix); m <- ncol(basis)
  zmat <- cbind(xfix, basis)
  logdet <- as.numeric(determinant(prec, logarithm = TRUE)$modulus)
  start <- c(fit_ns(y, o, xfix)$coef, rep(0, m))
  state <- new.env()
  state$th <- start
  run <- function(ls) {
    f <- pirls(y, o, zmat, np, prec / exp(2 * ls), state$th)
    state$th <- f$theta
    f$zmat <- zmat
    f
  }
  nll <- function(ls) -laplace_logml(run(ls), m, exp(ls), logdet)
  opt <- stats::optimize(nll, log(c(0.02, 20)), tol = 1e-3)
  fit <- run(opt$minimum)
  list(row = penalized_row(fit, np, exp(opt$minimum), alpha),
       logml = -opt$objective)
}

# `k` restricts the random effect to the k smoothest eigenvectors.
icar_fit <- function(y, o, xfix, spec, qinv, alpha = 0.05,
                     k = length(qinv)) {
  use <- seq_len(k)
  fit_penalized(y, o, xfix, spec$vecs[, use, drop = FALSE],
                diag(1 / qinv[use], k), alpha)
}

# BYM2: the mixing parameter is chosen on a grid by marginal likelihood.
bym2_fit <- function(y, o, xfix, spec, phis, alpha = 0.05,
                     k = length(spec$g)) {
  fits <- lapply(phis, function(ph) {
    icar_fit(y, o, xfix, spec, shape_qinv(spec, "bym2", ph), alpha, k)
  })
  best <- which.max(vapply(fits, `[[`, 0, "logml"))
  fits[[best]]$row
}

# Dimensions n/10 and n/4 of the Moran basis, capped at 200.
rsr_dims <- function(n, p) {
  sort(unique(as.integer(pmin(200, round(c(n / 10, n / 4)), n - p - 1))))
}

# Restricted spatial regression: Moran basis orthogonal to the fixed part.
rsr_basis <- function(spec, xfix, q) {
  qx <- qr.Q(qr(xfix))
  proj <- function(basis) basis - qx %*% crossprod(qx, basis)
  moran <- proj(t(proj(as.matrix(spec$adj))))
  ee <- eigen(moran, symmetric = TRUE)
  npos <- sum(ee$values > 1e-8 * max(abs(ee$values)))
  ee$vectors[, seq_len(min(q, npos)), drop = FALSE]
}

rsr_fit <- function(y, o, xfix, spec, q, alpha = 0.05) {
  basis <- rsr_basis(spec, xfix, q)
  lap <- Matrix::Diagonal(x = Matrix::rowSums(spec$adj)) - spec$adj
  prec <- as.matrix(Matrix::crossprod(basis, spec$g_node * lap %*% basis))
  fit_penalized(y, o, xfix, basis, prec, alpha)$row
}

# Spatial+ (Dupont et al., 2022) keeps the spatial term of the outcome model
# in a low-rank basis, as the thin plate splines of the original do. With
# the full basis the term can also absorb the residual of x, and the
# penalty then shifts the coefficient away from beta.
splus_dim <- function(n, m) {
  as.integer(min(50L, max(10L, round(n / 10)), m - 1L))
}

# Residual of x after a Gaussian ICAR smoother on the k smoothest
# eigenvectors, with variances fitted by maximum likelihood.
smooth_residual <- function(x, spec, qinv, k) {
  use <- seq_len(k)
  xu <- drop(crossprod(spec$vecs[, use, drop = FALSE], x - mean(x)))
  nll <- function(par) {
    v <- exp(par[1]) * qinv[use] + exp(par[2])
    sum(log(v) + xu^2 / v)
  }
  par <- stats::optim(c(0, 0), nll)$par
  v1 <- exp(par[1]) * qinv[use]
  x - drop(spec$vecs[, use, drop = FALSE] %*% (v1 / (v1 + exp(par[2])) * xu))
}

spatial_plus_fit <- function(y, o, ctrl, x, spec, qinv, phis, alpha = 0.05) {
  k <- splus_dim(spec$n, length(qinv))
  xfix <- cbind(ctrl, smooth_residual(x, spec, qinv, k))
  bym2_fit(y, o, xfix, spec, phis, alpha, k)
}

# All penalized competitors, in a fixed row order.
spatial_rows <- function(y, o, ctrl, x, spec, cfg, alpha = 0.05) {
  xfix <- cbind(ctrl, x)
  qs <- rsr_dims(spec$n, ncol(xfix))
  rows <- list(
    ICAR = icar_fit(y, o, xfix, spec, shape_qinv(spec, "icar"), alpha)$row,
    BYM2 = bym2_fit(y, o, xfix, spec, cfg$phi_grid, alpha),
    `Spatial+-ICAR` = spatial_plus_fit(y, o, ctrl, x, spec,
                                  shape_qinv(spec, "icar"),
                                  cfg$phi_grid, alpha))
  for (q in qs) rows[[sprintf("RSR(q=%d)", q)]] <- rsr_fit(y, o, xfix, spec, q)
  do.call(rbind, rows)
}

inla_rows <- function(y, o, x, graph_file) {
  dat <- data.frame(y = y, x = x, id = seq_along(y))
  one <- function(model) {
    hyper <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
    if (model == "bym2") {
      hyper$phi <- list(prior = "pc", param = c(0.5, 2 / 3))
    }
    form <- y ~ 1 + x + f(id, model = model, graph = graph_file,
                          scale.model = TRUE, constr = TRUE, hyper = hyper)
    r <- INLA::inla(form, family = "poisson", data = dat, offset = o,
                    num.threads = 1)
    s <- r$summary.fixed["x", ]
    stat_row(s[["mean"]], s[["sd"]], s[["0.025quant"]], s[["0.975quant"]])
  }
  out <- rbind(`ICAR-INLA` = one("besag"), `BYM2-INLA` = one("bym2"))
  out
}

# ---- Non-spatial and truncated-spectrum Poisson fits -------------------------
stat_names <- c("est", "sd", "lo", "hi", "lo_leak", "hi_leak", "extra",
           "pred_bias", "pred_sd", "nsif", "edf")

stat_row <- function(est, sd, lo = NA, hi = NA, ..., alpha = 0.05) {
  row <- stats::setNames(rep(NA_real_, length(stat_names)), stat_names)
  z <- stats::qnorm(1 - alpha / 2)
  row[c("est", "sd")] <- c(est, sd)
  row["lo"] <- if (is.na(lo)) est - z * sd else lo
  row["hi"] <- if (is.na(hi)) est + z * sd else hi
  extra <- list(...)
  row[names(extra)] <- unlist(extra)
  row
}

# Fit with the Poisson family, or the negative binomial with known size
# `theta` (the weights W = mu / (1 + mu / theta) enter the covariances).
# With dispersion = "fixed" the scale is phi = 1 (strict Poisson, or the
# negative binomial at the known size); with "estimated" it is the Pearson
# estimate (quasi-Poisson).
fit_ns <- function(y, o, xmat, theta = Inf, dispersion = "estimated") {
  fam <- if (is.finite(theta)) MASS::negative.binomial(theta) else
    stats::poisson()
  fit <- stats::glm.fit(xmat, y, offset = o, family = fam,
                        control = stats::glm.control(maxit = 50))
  mu <- fit$fitted.values
  n <- length(y); p <- ncol(xmat)
  wt <- reduced_weights(mu, theta)
  bread <- solve(crossprod(xmat, xmat * wt))
  meat <- crossprod(xmat * (wt * (y - mu) / mu))
  pearson <- sum((y - mu)^2 / (mu + mu^2 / theta))
  list(coef = unname(fit$coefficients), mu = mu, p = p, dev = fit$deviance,
       bread = bread,
       phi = if (dispersion == "fixed") 1 else pearson / (n - p),
       robust = bread %*% meat %*% bread * n / (n - p))
}

ns_rows <- function(fit, alpha = 0.05) {
  j <- fit$p
  est <- fit$coef[j]
  se_model <- sqrt(fit$phi * fit$bread[j, j])
  se_robust <- sqrt(fit$robust[j, j])
  rbind(`NS` = stat_row(est, se_model, alpha = alpha),
        `NS-robust` = stat_row(est, se_robust, alpha = alpha))
}

fit_truncated <- function(y, o, ctrl, x, vecs, n_vec, theta = Inf,
                          dispersion = "estimated") {
  xmat <- cbind(ctrl, x, vecs[, seq_len(n_vec), drop = FALSE])
  fit <- fit_ns(y, o, xmat, theta, dispersion)
  j <- ncol(ctrl) + 1L
  list(est = fit$coef[j], sd = sqrt(fit$phi * fit$bread[j, j]),
       dev = fit$dev, p = fit$p, n_vec = n_vec)
}

gcv_score <- function(fit, n) n * fit$dev / (n - fit$p)^2

# Rows for fixed K and for K chosen by generalized cross-validation.
truncation_rows <- function(y, o, ctrl, x, vecs, fixed, grid, alpha = 0.05,
                            theta = Inf, dispersion = "estimated") {
  ks <- sort(unique(c(fixed, grid)))
  fits <- lapply(ks, function(k) {
    fit_truncated(y, o, ctrl, x, vecs, k, theta, dispersion)
  })
  pick <- function(f) stat_row(f$est, f$sd, extra = f$n_vec, edf = f$p,
                               alpha = alpha)
  rows <- lapply(fits[match(fixed, ks)], pick)
  best <- fits[[which.min(vapply(fits, gcv_score, 0, n = length(y)))]]
  out <- do.call(rbind, c(rows, list(pick(best))))
  rownames(out) <- c(sprintf("Trunc(K=%d)", fixed), "Trunc(GCV)")
  out
}

# ---- Reduced geometry, minimax weights and risks -----------------------------
# Reduced weights W = mu / (1 + mu / theta); theta = Inf is the Poisson case.
reduced_weights <- function(mu, theta = Inf) mu / (1 + mu / theta)

reduced_geometry <- function(mu, ctrl, x, vecs, qinv, theta = Inf) {
  w <- sqrt(reduced_weights(mu, theta))
  q1 <- qr.Q(qr(w * ctrl))
  proj <- function(v) v - q1 %*% crossprod(q1, v)
  xt <- drop(proj(w * x))
  tm <- proj(w * vecs)
  half <- sweep(tm, 2, sqrt(qinv), "*")
  ee <- eigen(tcrossprod(half), symmetric = TRUE)
  list(w = w, q1 = q1, xt = xt, tm = tm, s = pmax(ee$values, 0),
       E = ee$vectors, xi = drop(crossprod(ee$vectors, xt)), qinv = qinv)
}

# Weighted residual of v on the columns of Z (weights w^2 = reduced
# weights); aliased columns are dropped.
weighted_residual <- function(zmat, w, v) {
  coef <- qr.coef(qr(w * zmat), w * v)
  coef[is.na(coef)] <- 0
  drop(v - zmat %*% coef)
}

info_f <- function(geo, c) sum(geo$xi^2 / (1 + c^2 * geo$s))

nsif <- function(geo, c) {
  vapply(c, function(cc) info_f(geo, cc), numeric(1)) / sum(geo$xt^2)
}

# Coordinates of the minimax weights in the eigenbasis of Sigma.
weights_coord <- function(geo, c) {
  geo$xi / (1 + c^2 * geo$s) / info_f(geo, c)
}

weights_vector <- function(geo, c) drop(geo$E %*% weights_coord(geo, c))

# Squared norm and worst-case bias per unit radius, from coordinates.
weights_norm <- function(coord) sqrt(sum(coord^2))

weights_bias <- function(geo, coord) sqrt(sum(geo$s * coord^2))

# Coordinates of an arbitrary weight vector given in R^n.
to_coord <- function(geo, a) drop(crossprod(geo$E, a))

# Alignment of a confounder (coordinates u, radius cu) with weights coord,
# in [0, 1]; NA when the bound c * b(a) vanishes.
alignment <- function(geo, coord, u, cu) {
  den <- cu * weights_bias(geo, coord)
  if (!is.finite(den) || den <= 0) return(NA_real_)
  abs(sum(coord * u)) / den
}

minimax_risk <- function(geo, c) 1 / info_f(geo, c)

minimax_risk_curve <- function(geo, c) {
  vapply(c, function(cc) minimax_risk(geo, cc), numeric(1))
}

worst_case_risk <- function(geo, coord, c) {
  sum(coord^2) + c^2 * sum(geo$s * coord^2)
}

# Break-even radius of a2 against a1; NA when a2 never wins.
break_even <- function(geo, coord1, coord2) {
  gain <- sum(coord2^2) - sum(coord1^2)
  b1 <- sum(geo$s * coord1^2)
  b2 <- sum(geo$s * coord2^2)
  if (gain <= 0 || b1 <= b2) return(NA_real_)
  sqrt(gain / (b1 - b2))
}

# Interval half-width and critical radius for given weights.
bias_aware_halfwidth <- function(sd, bias, c, alpha = 0.05) {
  sd * cva(c * bias / sd, alpha)
}

critical_budget <- function(est, sd, bias, alpha = 0.05) {
  if (bias <= 0) return(NA_real_)
  sd / bias * cva_inverse(abs(est) / sd, alpha)
}

# Distance from u (eigen coordinates) to the ellipsoid of radius c.
leak_radius <- function(u, s, budget) {
  if (budget <= 0) return(sqrt(sum(u^2)))
  pos <- s > 1e-12
  outside <- sqrt(sum(u[!pos]^2))
  if (sum(u[pos]^2 / s[pos]) <= budget^2) return(outside)
  g <- function(l) sum(u[pos]^2 * s[pos] / (s[pos] + l)^2) - budget^2
  l <- stats::uniroot(g, c(0, 1e8), tol = 1e-10)$root
  inside <- u[pos] * l / (s[pos] + l)
  sqrt(sum(inside^2) + outside^2)
}

# Restricted likelihood of the reduced model in the eigenbasis (minus
# twice its logarithm, up to a constant); the coordinates y and xi are
# already divided by the noise scale.
reml_nll <- function(ls, y, xi, s) {
  d <- 1 + exp(2 * ls) * s
  f <- sum(xi^2 / d)
  sum(log(d)) + log(f) + sum(y^2 / d) - sum(xi * y / d)^2 / f
}

reml_sigma <- function(y, xi, s, lower = 1e-3, upper = 50) {
  opt <- stats::optimize(reml_nll, log(c(lower, upper)), y = y, xi = xi,
                         s = s)
  exp(opt$minimum)
}

# Truncated-spectrum weights: x with the first K columns of tm projected
# out.
truncation_weights <- function(geo, n_vec) {
  if (n_vec < 1L) return(geo$xt / sum(geo$xt^2))
  qk <- qr.Q(qr(geo$tm[, seq_len(n_vec), drop = FALSE]))
  r <- geo$xt - qk %*% crossprod(qk, geo$xt)
  drop(r) / sum(r^2)
}

# ---- Graphs, their spectra and the prior shapes ------------------------------
lattice_graph <- function(k) {
  id <- matrix(seq_len(k * k), k, k)
  right <- cbind(as.vector(id[, -k]), as.vector(id[, -1]))
  down <- cbind(as.vector(id[-k, ]), as.vector(id[-1, ]))
  edges <- rbind(right, down)
  adj <- Matrix::sparseMatrix(edges[, 1], edges[, 2], x = 1,
                              dims = c(k * k, k * k), symmetric = TRUE)
  coords <- cbind(rep(seq_len(k), each = k), rep(seq_len(k), times = k))
  list(adj = adj, coords = coords, name = sprintf("lattice%dx%d", k, k))
}

delaunay_graph <- function(n, seed = 1L) {
  set.seed(seed)
  coords <- cbind(stats::runif(n), stats::runif(n))
  nb <- spdep::tri2nb(coords)
  adj <- as_binary(spdep::nb2mat(nb, style = "B", zero.policy = TRUE))
  list(adj = adj, coords = coords, name = sprintf("delaunay%d", n))
}

# Adjacency from an INLA graph file, e.g. the Scottish districts.
file_graph <- function(path, name) {
  g <- INLA::inla.read.graph(path)
  adj <- as_binary(INLA::inla.graph2matrix(g))
  list(adj = adj, coords = graph_coordinates(adj), name = name)
}

# Planar coordinates for a graph without any: classical multidimensional
# scaling of the shortest-path distances (unreachable pairs are set to
# twice the largest finite distance).
graph_coordinates <- function(adj) {
  d <- vapply(seq_len(nrow(adj)), function(i) graph_distance(adj, i),
              numeric(nrow(adj)))
  d[is.infinite(d)] <- 2 * max(d[is.finite(d)])
  stats::cmdscale(d, k = 2)
}

# Neighbour list from polygons: queen contiguity, no self links. Islands
# stay isolated components unless `join` links them to the nearest unit.
polygon_adjacency <- function(geom, queen = TRUE, join = FALSE) {
  nb <- drop_self_links(spdep::poly2nb(geom, queen = queen))
  if (join) nb <- join_islands(nb, geom)
  nb <- drop_self_links(spdep::make.sym.nb(nb))
  as_binary(spdep::nb2mat(nb, style = "B", zero.policy = TRUE))
}

as_binary <- function(m) {
  m <- as.matrix(m)
  m <- (m + t(m)) > 0
  diag(m) <- FALSE
  Matrix::Matrix(m * 1, sparse = TRUE)
}

drop_self_links <- function(nb) {
  for (i in seq_along(nb)) {
    keep <- nb[[i]][nb[[i]] != i & nb[[i]] > 0L]
    nb[[i]] <- if (length(keep)) as.integer(sort(unique(keep))) else 0L
  }
  nb
}

join_islands <- function(nb, geom) {
  islands <- which(spdep::card(nb) == 0L)
  if (!length(islands)) return(nb)
  pts <- sf::st_coordinates(suppressWarnings(sf::st_point_on_surface(geom)))
  for (i in islands) {
    d <- sqrt(rowSums((pts - matrix(pts[i, ], nrow(pts), 2, byrow = TRUE))^2))
    d[i] <- Inf
    j <- which.min(d)
    nb[[i]] <- sort(unique(c(nb[[i]][nb[[i]] > 0L], j)))
    nb[[j]] <- sort(unique(c(nb[[j]][nb[[j]] > 0L], i)))
  }
  nb
}

# Connected components of an adjacency matrix, as integer labels.
graph_components <- function(adj) {
  n <- nrow(adj)
  label <- integer(n)
  k <- 0L
  for (i in seq_len(n)) {
    if (label[i]) next
    k <- k + 1L
    label[is.finite(graph_distance(adj, i))] <- k
  }
  label
}

# Eigenpairs of the Laplacian of one component, with the generalized
# variance g of its intrinsic prior: the geometric mean of the diagonal of
# the pseudo-inverse. A single node has no non-null direction (g = 1).
component_spectrum <- function(lap, tol) {
  if (nrow(lap) == 1L) return(NULL)
  ee <- eigen(lap, symmetric = TRUE)
  keep <- ee$values > tol * max(ee$values)
  lambda <- ee$values[keep]
  basis <- ee$vectors[, keep, drop = FALSE]
  g <- exp(mean(log(drop((basis^2) %*% (1 / lambda)))))
  list(lambda = lambda, basis = basis, g = g)
}

# Eigenpairs of the Laplacian on the complement of its null space, sorted
# from the smoothest direction up. Each component is scaled on its own, as
# the intrinsic prior of a disconnected graph is the product of the
# priors of its components. Returned: `g` per column of `U` and per node
# (`g_node`), and the component indicators in `null_basis`, which belong
# in the control matrix.
graph_spectrum <- function(adj, tol = 1e-8) {
  adj <- as.matrix(adj)
  n <- nrow(adj)
  lap <- diag(rowSums(adj)) - adj
  comp <- graph_components(adj)
  g_comp <- rep(1, max(comp))
  pieces <- lapply(seq_len(max(comp)), function(k) {
    idx <- which(comp == k)
    sp <- component_spectrum(lap[idx, idx, drop = FALSE], tol)
    if (is.null(sp)) return(NULL)
    g_comp[k] <<- sp$g
    vecs <- matrix(0, n, length(sp$lambda))
    vecs[idx, ] <- sp$basis
    list(vecs = vecs, lambda = sp$lambda, g = rep(sp$g, length(sp$lambda)))
  })
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  stack <- function(field) lapply(pieces, `[[`, field)
  lambda <- unlist(stack("lambda"))
  ord <- order(lambda)
  list(vecs = do.call(cbind, stack("vecs"))[, ord, drop = FALSE],
       lambda = lambda[ord], g = unlist(stack("g"))[ord],
       g_node = g_comp[comp], comp = comp, n = n,
       null_basis = component_indicators(comp))
}

# One orthonormal indicator column per component.
component_indicators <- function(comp) {
  ind <- outer(comp, seq_len(max(comp)), "==") * 1
  sweep(ind, 2, sqrt(colSums(ind)), "/")
}

# Size, components and isolated nodes of a graph spectrum; m = n - G is the
# dimension of the class and sqrt(m) converts a radius c to c / sqrt(m).
graph_summary <- function(spec) {
  size <- tabulate(spec$comp)
  data.frame(graph = spec$label, n = spec$n, components = length(size),
             isolated = sum(size == 1L), m = ncol(spec$vecs),
             sqrt_m = sqrt(ncol(spec$vecs)))
}

# Inverse prior weights q_k^{-1} of the shape Q in the paper.
shape_qinv <- function(spec, shape = "icar", phi = 0.5) {
  base <- 1 / (spec$g * spec$lambda)  # g is one value per column
  switch(shape,
    icar = base,
    bym2 = phi * base + (1 - phi),
    stop("unknown shape: ", shape))
}

# Graph distances from a seed node by breadth-first search.
graph_distance <- function(adj, seed) {
  n <- nrow(adj)
  dist <- rep(Inf, n)
  dist[seed] <- 0
  frontier <- seed
  step <- 0
  while (length(frontier)) {
    step <- step + 1
    reach <- which(Matrix::colSums(adj[frontier, , drop = FALSE]) > 0)
    frontier <- reach[is.infinite(dist[reach])]
    dist[frontier] <- step
  }
  dist
}

scotland_graph <- function(cfg) {
  path <- cfg$scotland_graph
  if (is.null(path) && requireNamespace("INLA", quietly = TRUE)) {
    path <- system.file("demodata", "scotland.graph", package = "INLA")
  }
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  file_graph(path, "scotland")
}

build_spec <- function(name, cfg) {
  g <- switch(name,
    lattice = lattice_graph(cfg$lattice_side),
    delaunay = delaunay_graph(cfg$delaunay_n, seed = 1L),
    scotland = scotland_graph(cfg),
    stop("unknown graph: ", name))
  if (is.null(g)) {
    msg("graph '", name, "' unavailable, skipped")
    return(NULL)
  }
  c(g, graph_spectrum(g$adj), list(label = name))
}

# ---- Shape, leakage, estimated radius, benchmark and timing ------------------
fig_shape <- function(res, cfg) {
  res$pair <- paste(res$model, "fit /", res$truth, "truth")
  p <- ggplot2::ggplot(res, ggplot2::aes(bound, regret, colour = pair)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey30") +
    ggplot2::geom_point(alpha = 0.5, size = 0.9) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::facet_wrap(~ graph) + scale_colour_paper(name = NULL) +
    ggplot2::labs(x = "Regret bound from the two shape constants",
                  y = "Regret of the fitted shape") + theme_paper()
  save_figure(p, "shape_misspecification", cfg, 7, 3.8)
}

leakage_join <- function(res) {
  key <- c("graph", "mode", "level", "share", "decay", "design", "setting")
  leak <- res$info[grepl("^leak:[0-9]", res$info$name), c(key, "name", "value")]
  leak$method <- paste0("MM(c=", sub("leak:", "", leak$name), ")")
  merge(res$metrics, leak[c(key, "method", "value")], by = c(key, "method"))
}

fig_leakage <- function(res, cfg) {
  d <- leakage_join(res)
  d$value <- pmax(d$value, 1e-12)
  d <- d[d$method %in% c("MM(c=0)", "MM(c=0.5)", "MM(c=2)"), ]
  long <- rbind(
    data.frame(d[c("graph", "mode", "method", "value")], cover = d$cover,
               interval = "Nominal"),
    data.frame(d[c("graph", "mode", "method", "value")], cover = d$cover_leak,
               interval = "Widened by r"))
  p <- ggplot2::ggplot(long, ggplot2::aes(value, cover, colour = interval)) +
    ggplot2::geom_hline(yintercept = 0.95, colour = "grey30") +
    ggplot2::geom_point(alpha = 0.35, size = 0.7) +
    ggplot2::scale_x_log10() + ggplot2::facet_grid(mode ~ method) +
    scale_colour_paper(name = NULL) +
    ggplot2::labs(x = "Distance r from the field to the class (floor 1e-12)",
                  y = "Coverage") + theme_paper()
  save_figure(p, "leakage", cfg, 8, 4.8)
}

fig_estimated_budget <- function(res, cfg) {
  m <- only_positive(res$metrics)
  m <- m[m$method %in% c("MM(c0)", "MM(sigma)", "MM(sigma)-no-bias"), ]
  cell <- pool(m, c("graph", "fam", "kappa", "rho", "method", "c0"))
  cell$scenario <- ifelse(cell$fam == "aligned",
                          paste("Aligned, kappa =", cell$kappa), "Smooth")
  p <- ggplot2::ggplot(cell, ggplot2::aes(c0, cover,
                       colour = plain_method(method))) +
    band_layer(stats::median(cell$cover_se)) +
    ggplot2::geom_hline(yintercept = 0.95, colour = "grey30") +
    ggplot2::geom_line() + ggplot2::geom_point(size = 1) +
    ggplot2::scale_x_log10() + ggplot2::facet_grid(scenario ~ graph) +
    scale_colour_paper(name = NULL) +
    ggplot2::labs(x = "True radius c0", y = "Coverage") + theme_paper()
  save_figure(p, "estimated_budget_coverage", cfg, 7.5, 5.5)
  fig_budget_sigma(cell, cfg)
}

fig_budget_sigma <- function(cell, cfg) {
  s <- cell[cell$method == "MM(sigma)", ]
  p <- ggplot2::ggplot(s, ggplot2::aes(c0, sigma, colour = scenario)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey30",
                         linetype = "dashed") +
    ggplot2::geom_line() + ggplot2::geom_point(size = 1) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::facet_wrap(~ graph) + scale_colour_paper(name = NULL) +
    ggplot2::labs(x = "True radius c0",
                  y = "Median restricted-likelihood estimate of sigma") +
    theme_paper()
  save_figure(p, "estimated_budget_sigma", cfg, 7, 3.8)
}

fig_budget_profile <- function(res, cfg) {
  res$panel <- sprintf("%s, c0 = %g", res$fam, res$c0)
  p <- ggplot2::ggplot(res, ggplot2::aes(sigma, nll, group = rep)) +
    ggplot2::geom_line(alpha = 0.6, colour = paper_colours[1]) +
    ggplot2::scale_x_log10() + ggplot2::coord_cartesian(ylim = c(0, 15)) +
    ggplot2::facet_grid(graph ~ panel) +
    ggplot2::labs(x = "sigma", y = "Restricted deviance above its minimum") +
    theme_paper()
  save_figure(p, "estimated_budget_profile", cfg, 8, 4.5)
}

fig_benchmark <- function(res, cfg) {
  cell <- stats::aggregate(p_robust ~ graph + beta + c0 + mult, res, mean)
  p <- ggplot2::ggplot(cell, ggplot2::aes(mult, p_robust,
                       colour = factor(beta))) +
    ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed",
                        colour = "grey40") +
    ggplot2::geom_line() + ggplot2::geom_point(size = 1) +
    ggplot2::scale_x_log10(breaks = unique(cell$mult)) +
    ggplot2::facet_grid(graph ~ c0, labeller = ggplot2::label_both) +
    scale_colour_paper(name = "True effect") +
    ggplot2::labs(x = "Assumed confounding / benchmark",
                  y = "Share declared robust") + theme_paper()
  save_figure(p, "benchmark", cfg, 8, 4.5)
}

table_timing <- function(res, cfg) {
  out <- res
  out[-1] <- lapply(out[-1], function(v) fmt(v, 2))
  out$n <- as.character(res$n)
  names(out) <- c("Graph", "$n$", "Geometry (s)", "NS fit (s)", "Curve (s)",
                  "REML (s)", "Penalized fit (s)", "Geometry (MB)")
  write_table(out, "table_timing", cfg, "Computing time of the steps of the
    procedure and memory of the geometry.", "tab:sim_timing")
}

# ---- Plot theme, pooling over designs and LaTeX tables -----------------------
# Okabe-Ito palette, safe for colour-blind readers.
paper_colours <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
             "#56B4E9", "#000000", "#F0E442")

theme_paper <- function(base = 10) {
  ggplot2::theme_bw(base_size = base) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = "grey94"),
                   legend.position = "bottom",
                   legend.key.size = ggplot2::unit(0.9, "lines"))
}

scale_colour_paper <- function(...) {
  ggplot2::scale_colour_manual(values = paper_colours, ...)
}

save_figure <- function(plot, name, cfg, width = 7, height = 4.5) {
  dir <- ensure_dir(file.path(cfg$out_dir, "figures"))
  dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(file.path(dir, paste0(name, ".pdf")), plot,
                  width = width, height = height, device = dev)
  ggplot2::ggsave(file.path(dir, paste0(name, ".png")), plot,
                  width = width, height = height, dpi = 300, bg = "white")
  invisible(plot)
}

read_result <- function(cfg, name) {
  path <- file.path(cfg$out_dir, sprintf("results_%s.rds", name))
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

# Concatenate result lists (metrics and info); NULL entries are ignored.
merge_results <- function(...) {
  parts <- Filter(Negate(is.null), list(...))
  if (!length(parts)) return(NULL)
  list(metrics = rbind_list(lapply(parts, `[[`, "metrics")),
       info = rbind_list(lapply(parts, `[[`, "info")))
}

method_family <- function(m) {
  ifelse(grepl("^MM", m), "Proposed",
  ifelse(grepl("^NS", m), "Non-spatial",
  ifelse(grepl("^Trunc", m), "Truncated spectrum", "Spatial model")))
}

plain_method <- function(m) gsub("sigma", "sigma-hat", m, fixed = TRUE)

# Monte Carlo error of the coverage averaged over designs: the spread of
# the design-level coverages (which carries both the between-design and the
# within-design variation); binomial when there is a single design.
design_se <- function(cover, reps, pooled) {
  if (length(cover) > 1L) return(stats::sd(cover) / sqrt(length(cover)))
  sqrt(pooled * (1 - pooled) / reps)
}

pool_block <- function(d) {
  w <- d$n_rep
  wm <- function(v) stats::weighted.mean(v, w, na.rm = TRUE)
  cover <- wm(d$cover)
  reps <- sum(w)
  se <- design_se(d$cover, reps, cover)
  data.frame(
    designs = nrow(d), reps = reps, bias = wm(d$bias),
    rmse = sqrt(wm(d$rmse^2)), emp_sd = sqrt(wm(d$emp_sd^2)),
    mean_sd = wm(d$mean_sd), cover = cover,
    cover_se = se,
    cover_min = min(d$cover, na.rm = TRUE),
    cover_med = stats::median(d$cover, na.rm = TRUE),
    cover_q10 = unname(stats::quantile(d$cover, 0.1, na.rm = TRUE)),
    cover_q05 = unname(stats::quantile(d$cover, 0.05, na.rm = TRUE)),
    len = wm(d$len), cover_leak = wm(d$cover_leak),
    len_leak = wm(d$len_leak), reject = wm(d$reject),
    pred_mse = wm(d$pred_bias^2 + d$pred_sd^2), nsif = wm(d$nsif),
    pred_bias = wm(d$pred_bias), pred_sd = wm(d$pred_sd),
    sigma = wm(d$extra), edf = wm(d$edf))
}

# Pool the per-design rows within the groups defined by `by`.
pool <- function(df, by) {
  key <- do.call(paste, c(df[by], sep = "\r"))
  rows <- lapply(split(df, key), function(d) {
    cbind(d[1, by, drop = FALSE], pool_block(d), row.names = NULL)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Band of +-1.96 Monte Carlo standard errors around the nominal level,
# drawn on the y or x axis.
band_layer <- function(se, axis = "y", target = 0.95) {
  half <- 1.96 * se
  lim <- c(target - half, target + half)
  if (axis == "x") {
    return(ggplot2::annotate("rect", xmin = lim[1], xmax = lim[2],
                             ymin = -Inf, ymax = Inf, alpha = 0.18,
                             fill = "grey40"))
  }
  ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = lim[1],
                    ymax = lim[2], alpha = 0.18, fill = "grey40")
}

only_positive <- function(df, col = "c0") df[df[[col]] > 0, ]

# Escape LaTeX specials outside math, and typeset sigma-hat.
tex_text <- function(x) {
  one <- function(s) {
    pat <- "\\$[^$]*\\$"
    math <- regmatches(s, gregexpr(pat, s))[[1]]
    plain <- regmatches(s, gregexpr(pat, s), invert = TRUE)[[1]]
    plain <- gsub("([%&#_])", "\\\\\\1", plain)
    plain <- gsub("sigma", "$\\\\hat\\\\sigma$", plain)
    out <- plain[1]
    for (i in seq_along(math)) out <- paste0(out, math[i], plain[i + 1L])
    out
  }
  vapply(as.character(x), one, character(1), USE.NAMES = FALSE)
}

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "--", formatC(x, digits = digits, format = "f"))
}

bold <- function(x, cond) {
  ifelse(cond & !is.na(cond), paste0("\\textbf{", x, "}"), x)
}

latex_table <- function(df, caption, label, note = NULL, size = "\\small") {
  cols <- paste0("l", strrep("r", ncol(df) - 1))
  head <- paste(tex_text(names(df)), collapse = " & ")
  body <- apply(df, 1, function(r) paste(r, collapse = " & "))
  c("\\begin{table}[ht]", paste0("\\centering", size),
    sprintf("\\caption{%s}", caption), sprintf("\\label{%s}", label),
    sprintf("\\begin{tabular}{@{}%s@{}}", cols), "\\toprule",
    paste0(head, " \\\\"), "\\midrule", paste0(body, " \\\\"),
    "\\bottomrule", "\\end{tabular}",
    if (!is.null(note)) sprintf("\\par\\smallskip\\footnotesize %s", note),
    "\\end{table}")
}

# Single-spaced caption with a literal percent sign escaped.
caption_text <- function(x) {
  gsub("(?<!\\\\)%", "\\\\%", gsub("\\s+", " ", x), perl = TRUE)
}

write_table <- function(df, stem, cfg, caption, label, note = NULL,
                        size = "\\footnotesize\\setlength{\\tabcolsep}{3pt}") {
  dir <- ensure_dir(file.path(cfg$out_dir, "tables"))
  utils::write.csv(df, file.path(dir, paste0(stem, ".csv")),
                   row.names = FALSE)
  df[] <- lapply(df, tex_text)
  writeLines(latex_table(df, caption_text(caption), label, note, size),
             file.path(dir, paste0(stem, ".tex")))
  invisible(df)
}

# Row of `df` with the smallest `var` within each group of `by`.
min_by <- function(df, by, var) {
  key <- do.call(paste, c(df[by], sep = "\r"))
  rows <- lapply(split(df, key), function(d) d[which.min(d[[var]]), ])
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# ---- Coverage across confounder families -------------------------------------
family_labels <- c(aligned = "Aligned", smooth = "Smooth",
                   band = "Band-limited", matern = "Matern",
                   local = "Localized", independent = "Independent",
                   orthogonal = "Orthogonal", mixture = "Mixture")

base_shown <- c("NS", "NS-robust", "Trunc(K=20)", "Trunc(GCV)", "ICAR",
                "BYM2", "Spatial+-ICAR")

proposed_shown <- c("MM(c0)", "MM(sigma)", "MM(sigma)-no-bias")

# Methods to display, in a fixed order, among those present.
shown_methods <- function(present) {
  c(intersect(base_shown, present), sort(present[grepl("^RSR", present)]),
    intersect(proposed_shown, present))
}

# Parameters that define a setting, with a short text for the csv files.
setting_cols <- c("c0", "kappa", "target", "rho", "n_vec", "nu", "range",
                  "kind", "size", "amp")

setting_text <- function(st) {
  one <- function(i) {
    r <- st[i, ]
    switch(r$fam,
      aligned = sprintf("target %s, kappa %g, c0 %g", r$target, r$kappa,
                        r$c0),
      smooth = sprintf("rho %g, c0 %g", r$rho, r$c0),
      band = sprintf("K %d, c0 %g", r$n_vec, r$c0),
      matern = sprintf("nu %g, range %g, c0 %g", r$nu, r$range, r$c0),
      local = if (r$kind == "hotspot") {
        sprintf("hotspot of %d units, amplitude %g", r$size, r$amp)
      } else sprintf("region, amplitude %g", r$amp),
      mixture = sprintf("c0 %g", r$c0),
      sprintf("c0 %g", r$c0))
  }
  vapply(seq_len(nrow(st)), one, "")
}

family_cells <- function(m) {
  m <- m[m$method %in% shown_methods(unique(m$method)), ]
  cell <- pool(m, c("graph", "fam", "setting", "method"))
  cell$family_label <- family_labels[cell$fam]
  keys <- c("graph", "fam", "setting")
  params <- unique(m[c(keys, setting_cols)])
  params$description <- setting_text(params)
  merge(cell, params, by = keys, sort = FALSE)
}

fig_coverage_dots <- function(cell, name, cfg, height = 7) {
  summ <- do.call(rbind, lapply(split(cell, list(cell$graph, cell$fam,
                                    cell$method), drop = TRUE), function(d) {
    data.frame(graph = d$graph[1], fam = d$family_label[1],
               method = d$method[1], med = median(d$cover),
               lo = min(d$cover), hi = max(d$cover))
  }))
  summ$method <- factor(summ$method,
                        levels = rev(shown_methods(unique(cell$method))))
  summ$class <- method_family(as.character(summ$method))
  p <- ggplot2::ggplot(summ, ggplot2::aes(med, method, colour = class)) +
    band_layer(stats::median(cell$cover_se), "x") +
    ggplot2::geom_vline(xintercept = 0.95, colour = "grey30") +
    ggplot2::geom_pointrange(ggplot2::aes(xmin = lo, xmax = hi), size = 0.25) +
    ggplot2::facet_grid(graph ~ fam) + scale_colour_paper(name = NULL) +
    ggplot2::labs(x = "Coverage: median and range over settings",
                  y = NULL) +
    ggplot2::coord_cartesian(xlim = c(0, 1)) + theme_paper(9)
  save_figure(p, name, cfg, 11, height)
}

fig_length <- function(m, name, cfg) {
  keep <- shown_methods(unique(m$method))
  cell <- pool(m[m$method %in% keep, ], c("graph", "fam", "setting", "method"))
  ref <- cell[cell$method == "NS", c("graph", "fam", "setting", "len")]
  d <- merge(cell, ref, by = c("graph", "fam", "setting"),
             suffixes = c("", "_ns"))
  d$rel <- d$len / d$len_ns
  d$method <- factor(d$method, levels = rev(keep))
  p <- ggplot2::ggplot(d, ggplot2::aes(rel, method, colour = graph)) +
    ggplot2::geom_vline(xintercept = 1, colour = "grey30") +
    ggplot2::geom_jitter(height = 0.15, alpha = 0.5, size = 0.7) +
    ggplot2::scale_x_log10() + ggplot2::facet_wrap(~ fam, nrow = 2) +
    scale_colour_paper(name = "Graph") +
    ggplot2::labs(x = "Interval length relative to NS", y = NULL) +
    theme_paper(9)
  save_figure(p, name, cfg, 10, 6)
}

# 5% quantile of a vector, the summary of the settings of a family.
q05 <- function(v) unname(stats::quantile(v, 0.05, na.rm = TRUE))

# Heat map of the coverage of every method in every family, as the
# difference from 0.95; `stat` summarises the settings of a family.
fig_coverage_heatmap <- function(cell, name, cfg, stat = q05) {
  cell$method <- factor(cell$method,
                        levels = rev(shown_methods(unique(cell$method))))
  d <- stats::aggregate(cover ~ graph + family_label + method, cell, stat)
  d$gap <- d$cover - 0.95
  p <- ggplot2::ggplot(d, ggplot2::aes(family_label, method, fill = gap)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", cover)),
                       size = 2.4) +
    ggplot2::scale_fill_gradient2(low = "#D55E00", mid = "white",
                                  high = "#0072B2", midpoint = 0,
                                  limits = c(-0.95, 0.05), oob = scales::squish,
                                  name = "Coverage - 0.95") +
    ggplot2::facet_wrap(~ graph, nrow = 1) +
    ggplot2::labs(x = NULL, y = NULL) + theme_paper(9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  save_figure(p, name, cfg, 11, 4.2)
}

# Coverage against interval length relative to the non-spatial interval:
# what the protection costs. Points are medians over the settings.
fig_coverage_width <- function(cell, name, cfg) {
  ref <- cell[cell$method == "NS", c("graph", "fam", "setting", "len")]
  d <- merge(cell, ref, by = c("graph", "fam", "setting"),
             suffixes = c("", "_ns"))
  d$rel <- d$len / d$len_ns
  d <- stats::aggregate(cbind(cover, rel) ~ graph + family_label + method,
                        d, stats::median)
  d$class <- method_family(as.character(d$method))
  p <- ggplot2::ggplot(d, ggplot2::aes(rel, cover, colour = class,
                                       shape = class)) +
    ggplot2::geom_hline(yintercept = 0.95, linetype = "dashed",
                        colour = "grey30") +
    ggplot2::geom_point(alpha = 0.8, size = 1.4) +
    ggplot2::scale_x_log10() + scale_colour_paper(name = NULL) +
    ggplot2::scale_shape_discrete(name = NULL) +
    ggplot2::facet_grid(graph ~ family_label) +
    ggplot2::labs(x = "Interval length relative to NS", y = "Coverage") +
    theme_paper(9)
  save_figure(p, name, cfg, 11, 6)
}

# Table of one summary of the pooled coverage of the settings of a family:
# the 5% quantile in the text tables, the minimum in the appendix.
table_family <- function(cell, graph, stem, cfg, stat = "q05") {
  d <- cell[cell$graph == graph, ]
  fun <- if (stat == "min") min else q05
  worst <- stats::aggregate(cover ~ method + fam, d, fun)
  wide <- stats::reshape(worst, idvar = "method", timevar = "fam",
                         direction = "wide")
  wide <- wide[match(shown_methods(wide$method), wide$method), ]
  vals <- as.matrix(wide[-1])
  low <- apply(vals, 1, min, na.rm = TRUE)
  out <- data.frame(method = wide$method, lapply(seq_len(ncol(vals)),
    function(j) bold(fmt(vals[, j]), vals[, j] == low)))
  names(out) <- c("Method", sub("cover.", "", colnames(vals), fixed = TRUE))
  what <- if (stat == "min") "Minimum" else "5% quantile"
  write_table(out, stem, cfg, sprintf("%s of the coverage over the settings
    of each family, %s graph; the lowest entry of each row is in bold.",
    what, graph),
    paste0("tab:", stem), paste0("Nominal level 0.95. Replicates per cell: ",
           paste(unique(range(d$reps)), collapse = " to "), "."))
}

# Coverage summary of every cell: pooled coverage with its Monte Carlo
# standard error and Wald interval, and the median, 10% and 5% quantiles
# and minimum of the design-level coverages. No cell is classified as
# nominal or not; the full table goes to a csv file.
coverage_summary <- function(cell) {
  keep <- c("graph", "fam", "setting", "description", setting_cols,
            "method", "designs", "reps", "cover",
            "cover_se", "cover_med", "cover_q10", "cover_q05", "cover_min")
  out <- cell[keep]
  out$lo <- out$cover - 1.96 * out$cover_se
  out$hi <- out$cover + 1.96 * out$cover_se
  out
}

report_coverage <- function(res, tag, cfg) {
  m <- res$metrics
  cell <- family_cells(m)
  utils::write.csv(coverage_summary(cell),
                   file.path(ensure_dir(cfg$out_dir),
                             sprintf("%s.csv", tag)),
                   row.names = FALSE)
  fig_coverage_dots(cell, paste0(tag, "_dots"), cfg)
  fig_length(m, paste0(tag, "_length"), cfg)
  fig_coverage_heatmap(cell, paste0(tag, "_heatmap"), cfg)
  fig_coverage_width(cell, paste0(tag, "_width"), cfg)
  for (g in unique(cell$graph)) {
    table_family(cell, g, sprintf("table_%s_%s", tag, g), cfg)
    table_family(cell, g, sprintf("table_%s_%s_min", tag, g), cfg, "min")
  }
}

# ---- Report driver -----------------------------------------------------------
# The quasi-Poisson run shares the baseline geometry: give it its own label.
relabel_quasi <- function(res) {
  if (is.null(res)) return(NULL)
  for (nm in c("metrics", "info")) {
    res[[nm]]$geometry <- paste(res[[nm]]$geometry, "quasi")
  }
  res
}

# Every report is attempted; the run fails at the end if any of them did.
make_report <- function(cfg) {
  failed <- character()
  safe <- function(label, expr) {
    tryCatch(expr, error = function(e) {
      msg("report ", label, " failed: ", conditionMessage(e))
      failed <<- c(failed, label)
    })
  }
  r <- function(n) read_result(cfg, n)
  lin <- merge_results(r("count_linearization"),
                       r("count_linearization_true"),
                       r("count_linearization_fitted"),
                       relabel_quasi(r("count_linearization_quasi")))
  gauss <- r("gaussian_exact"); cost <- r("adjustment_cost")
  even <- r("break_even")
  if (!is.null(gauss)) safe("gaussian_exact", {
    d <- gaussian_data(gauss)
    fig_gaussian_mse(d, cfg); fig_gaussian_coverage(d, cfg)
    fig_gaussian_alignment(d, cfg); table_gaussian(d, cfg)
  })
  if (!is.null(lin)) safe("count_linearization", fig_linearization(lin, cfg))
  if (!is.null(cost)) safe("adjustment_cost", fig_cost(cost, cfg))
  if (!is.null(even)) {
    safe("break_even", fig_break_even(even, cfg))
    safe("break_even random", fig_break_even_random(even, cfg))
  }
  for (tag in c("coverage", "coverage_brazil", "coverage_nb",
                "coverage_spatial_exposure")) {
    if (!is.null(r(tag))) safe(tag, report_coverage(r(tag), tag, cfg))
  }
  optional <- list(
    shape_misspecification = fig_shape, leakage = fig_leakage,
    estimated_budget = fig_estimated_budget,
    estimated_budget_profile = fig_budget_profile,
    benchmark = fig_benchmark, timing = table_timing)
  for (nm in names(optional)) {
    if (!is.null(r(nm))) safe(nm, optional[[nm]](r(nm), cfg))
  }
  if (length(failed)) stop("reports failed: ", paste(failed, collapse = ", "))
  invisible(NULL)
}

# ---- Gaussian check, count linearization, adjustment cost, break-even --------
gaussian_data <- function(res) {
  m <- res$metrics
  keep <- (m$method == "MM(c0)" & m$target == "c0") |
    (m$method == "NS" & m$target == "ns")
  m[keep, ]
}

fig_gaussian_mse <- function(m, cfg) {
  m$emp <- m$rmse^2
  m$pred <- m$pred_bias^2 + m$pred_sd^2
  p <- ggplot2::ggplot(m, ggplot2::aes(pred, emp, colour = graph)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey30") +
    ggplot2::geom_point(alpha = 0.35, size = 0.8) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::facet_wrap(~ method) + scale_colour_paper(name = "Graph") +
    ggplot2::labs(x = "Predicted mean squared error",
                  y = "Monte Carlo mean squared error") + theme_paper()
  save_figure(p, "gaussian_mse", cfg, 7, 3.8)
}

fig_gaussian_coverage <- function(m, cfg) {
  cell <- pool(only_positive(m), c("graph", "method", "c0", "kappa", "share"))
  worst <- min_by(cell, c("graph", "method", "c0", "share"), "cover")
  p <- ggplot2::ggplot(worst, ggplot2::aes(c0, cover, colour = method,
                       linetype = factor(share))) +
    band_layer(stats::median(cell$cover_se)) +
    ggplot2::geom_hline(yintercept = 0.95, colour = "grey30") +
    ggplot2::geom_line() + ggplot2::geom_point(size = 1) +
    ggplot2::scale_x_log10() + ggplot2::facet_wrap(~ graph) +
    scale_colour_paper(name = "Estimator") +
    ggplot2::scale_linetype_discrete(name = "Spatial share of x") +
    ggplot2::labs(x = "Radius c0", y = "Worst coverage over alignment") +
    theme_paper()
  save_figure(p, "gaussian_coverage", cfg, 7, 3.8)
}

fig_gaussian_alignment <- function(m, cfg) {
  cell <- pool(only_positive(m[m$method == "MM(c0)", ]),
               c("graph", "c0", "kappa"))
  p <- ggplot2::ggplot(cell, ggplot2::aes(kappa, cover,
                       colour = factor(c0), group = c0)) +
    band_layer(stats::median(cell$cover_se)) +
    ggplot2::geom_hline(yintercept = 0.95, colour = "grey30") +
    ggplot2::geom_line() + ggplot2::geom_point(size = 1) +
    ggplot2::facet_wrap(~ graph) + scale_colour_paper(name = "c0") +
    ggplot2::labs(x = "Alignment kappa", y = "Coverage") + theme_paper()
  save_figure(p, "gaussian_alignment", cfg, 7, 3.8)
}

table_gaussian <- function(m, cfg) {
  m <- only_positive(m[m$method == "MM(c0)", ])
  cell <- pool(m, c("graph", "c0", "kappa"))
  cell$ratio <- cell$rmse^2 / cell$pred_mse
  rows <- lapply(split(cell, list(cell$graph, cell$c0), drop = TRUE),
                 function(d) {
    data.frame(graph = d$graph[1], c0 = d$c0[1],
               max_dev = max(abs(d$ratio - 1)), min_cover = min(d$cover),
               max_cover = max(d$cover), nsif = mean(d$nsif))
  })
  out <- do.call(rbind, rows)
  out[] <- lapply(out, function(v) if (is.numeric(v)) fmt(v) else v)
  names(out) <- c("Graph", "$c_0$", "Max |MSE ratio - 1|", "Min coverage",
                  "Max coverage", "NSIF")
  write_table(out, "table_gaussian", cfg, "Gaussian reduced model:
    agreement with the exact risk and coverage of the bias-aware interval
    at the worst alignment (minimum over kappa).", "tab:sim_gaussian")
}

fig_linearization <- function(res, cfg) {
  m <- only_positive(res$metrics[res$metrics$method == "MM(c0)", ])
  cell <- pool(m, c("graph", "geometry", "level", "c0", "kappa"))
  cell$gap <- abs(cell$bias - cell$pred_bias) / cell$pred_sd
  cell$sdr <- cell$emp_sd / cell$mean_sd
  keys <- c("graph", "geometry", "level", "c0")
  collapse <- function(var, fun, label) {
    a <- stats::aggregate(cell[[var]], cell[keys], fun)
    names(a)[length(keys) + 1L] <- "value"
    cbind(a, measure = label)
  }
  all <- rbind(collapse("gap", max, "Bias gap (SD units)"),
               collapse("sdr", mean, "SD: empirical / estimated"),
               collapse("cover", min, "Worst coverage"))
  p <- ggplot2::ggplot(all, ggplot2::aes(c0, value, colour = factor(level),
                       linetype = geometry)) +
    ggplot2::geom_line() + ggplot2::geom_point(size = 1) +
    ggplot2::scale_x_log10() +
    ggplot2::facet_grid(measure ~ graph, scales = "free_y") +
    scale_colour_paper(name = "Median mean count") +
    ggplot2::scale_linetype_discrete(name = "Geometry") +
    ggplot2::labs(x = "Radius c0", y = NULL) + theme_paper()
  save_figure(p, "count_linearization", cfg, 7, 6)
}

cost_pairs <- function(m) {
  key <- c("graph", "mode", "share", "decay", "design", "setting")
  ns <- m[m$method == "NS", c(key, "emp_sd", "pred_sd")]
  mm <- m[grepl("^MM\\(c=", m$method), ]
  mm$budget <- as.numeric(sub("MM\\(c=(.*)\\)", "\\1", mm$method))
  keep <- c(key, "budget", "emp_sd", "pred_sd", "nsif", "fam", "c0")
  out <- merge(mm[mm$budget > 0, keep], ns, by = key,
               suffixes = c("_mm", "_ns"))
  out$mc <- (out$emp_sd_mm / out$emp_sd_ns)^2
  out$exact <- (out$pred_sd_mm / out$pred_sd_ns)^2
  out$bound <- 1 / out$nsif
  out
}

# Gaussian runs use the exact ratio of model-based variances, count runs
# the Monte Carlo ratio. The shaded region is the one the bound allows.
fig_cost <- function(res, cfg) {
  d <- cost_pairs(res$metrics)
  d$ratio <- ifelse(d$mode == "gaussian", d$exact, d$mc)
  grid <- seq(min(d$bound), max(d$bound), length.out = 100)
  region <- data.frame(bound = grid, low = 1, high = grid)
  p <- ggplot2::ggplot(d, ggplot2::aes(bound, ratio,
                       colour = factor(decay), shape = factor(budget))) +
    ggplot2::geom_ribbon(data = region, inherit.aes = FALSE,
                         ggplot2::aes(bound, ymin = low, ymax = high),
                         fill = "grey90") +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey30") +
    ggplot2::geom_hline(yintercept = 1, colour = "grey30") +
    ggplot2::geom_point(alpha = 0.3, size = 1) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::facet_grid(mode ~ graph) + scale_colour_paper(name = "Decay") +
    ggplot2::scale_shape_discrete(name = "Radius c") +
    ggplot2::labs(x = "Upper bound 1 / NSIF(c)",
                  y = "Var(adjusted) / Var(NS)") +
    theme_paper()
  save_figure(p, "cost_vs_nsif", cfg, 7, 4.5)
  table_cost(d, cfg)
}

# The exact ratio uses the model-based variances of the Gaussian run; the
# Monte Carlo ratio uses the empirical ones and is noisy at finite R.
table_cost <- function(d, cfg) {
  d <- d[d$budget %in% c(0.1, 0.5, 2) & d$mode == "gaussian", ]
  rows <- lapply(split(d, list(d$graph, d$budget), drop = TRUE),
                 function(g) {
    data.frame(graph = g$graph[1], c = g$budget[1],
               bound = stats::median(g$bound),
               exact = stats::median(g$exact),
               exact_max = max(g$exact / g$bound),
               mc = stats::median(g$mc),
               mc_med = stats::median(g$mc / g$bound),
               mc_q95 = stats::quantile(g$mc / g$bound, 0.95,
                                        names = FALSE))
  })
  out <- do.call(rbind, rows)
  out[] <- lapply(out, function(v) if (is.numeric(v)) fmt(v) else v)
  names(out) <- c("Graph", "$c$", "$1/\\mathrm{NSIF}$", "Exact",
                  "Exact / bound, max", "Monte Carlo",
                  "MC / bound, median", "MC / bound, 95th pct.")
  write_table(out, "table_cost", cfg, "Variance of the adjusted estimator
    relative to the non-spatial one when the field carries no bias
    (independent and orthogonal fields pooled), Gaussian reduced model.
    The exact ratio never exceeds the bound 1/NSIF; the Monte Carlo ratio
    scatters around it with the sampling error of two variances.",
    "tab:sim_cost")
}

own_target <- c(NS = "ns", `MM(c=0.1)` = "mm:0.1", `MM(c=0.5)` = "mm:0.5",
                `MM(c=2)` = "mm:2", `Trunc(K=5)` = "trunc:5",
                `Trunc(K=20)` = "trunc:20")

# Worst-case MSE of each estimator, on fields aligned against itself.
break_even_worst <- function(m) {
  m <- m[m$fam == "aligned" & m$method %in% names(own_target), ]
  m <- m[own_target[m$method] == m$target, ]
  m$mse <- m$rmse^2
  m
}

fig_break_even <- function(res, cfg) {
  w <- break_even_worst(res$metrics)
  cell <- pool(w, c("graph", "mode", "method", "c0"))
  cell$mse <- cell$rmse^2
  star <- break_even_cstar(res$info)
  p <- ggplot2::ggplot(only_positive(cell), ggplot2::aes(c0, mse,
                       colour = method)) +
    ggplot2::geom_line() + ggplot2::geom_point(size = 0.9) +
    ggplot2::geom_vline(data = star, ggplot2::aes(xintercept = cstar,
                        colour = method), linetype = "dotted") +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::facet_grid(mode ~ graph) + scale_colour_paper(name = "Estimator") +
    ggplot2::labs(x = "Radius c0", y = "Worst-case mean squared error") +
    theme_paper()
  save_figure(p, "break_even_worst_case", cfg, 7.5, 5)
  fig_break_even_gap(w, star, cfg)
  table_break_even(w, res$info, cfg)
}

# Paired difference of the worst-case mean squared errors, relative to the
# non-spatial one: below zero the adjustment pays off. The band is the
# Monte Carlo interval of the mean over designs.
fig_break_even_gap <- function(w, star, cfg) {
  key <- c("graph", "mode", "design", "c0")
  ns <- w[w$method == "NS", c(key, "mse")]
  d <- merge(w[grepl("^MM", w$method), ], ns, by = key,
             suffixes = c("", "_ns"))
  d$gap <- (d$mse - d$mse_ns) / d$mse_ns
  by <- c("graph", "mode", "method", "c0")
  cell <- stats::aggregate(gap ~ graph + mode + method + c0, d,
                           function(v) c(mean = mean(v),
                                         se = stats::sd(v) / sqrt(length(v))))
  cell <- do.call(data.frame, cell)
  names(cell) <- c(by, "gap", "se")
  p <- ggplot2::ggplot(only_positive(cell), ggplot2::aes(c0, gap,
                       colour = method, fill = method)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey30") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = gap - 1.96 * se,
                                      ymax = gap + 1.96 * se),
                         alpha = 0.2, colour = NA) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(data = star, ggplot2::aes(xintercept = cstar,
                        colour = method), linetype = "dotted") +
    ggplot2::scale_x_log10() +
    ggplot2::facet_grid(mode ~ graph) + scale_colour_paper(name = "Estimator") +
    ggplot2::scale_fill_manual(values = paper_colours, name = "Estimator") +
    ggplot2::labs(x = "Radius c0",
                  y = "Relative difference in worst-case MSE") +
    theme_paper()
  save_figure(p, "break_even_gap", cfg, 7.5, 5)
}

fig_break_even_random <- function(res, cfg) {
  m <- res$metrics
  m <- m[m$fam %in% c("smooth", "band") & m$method %in% names(own_target), ]
  m$scenario <- ifelse(m$fam == "smooth", "Smooth, rho = 0.6",
                       sprintf("Band-limited, K = %d", m$n_vec))
  cell <- pool(only_positive(m), c("graph", "mode", "scenario", "method",
                                   "c0"))
  cell$mse <- cell$rmse^2
  p <- ggplot2::ggplot(cell, ggplot2::aes(c0, mse, colour = method)) +
    ggplot2::geom_line() + ggplot2::geom_point(size = 0.9) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::facet_grid(scenario ~ graph + mode) +
    scale_colour_paper(name = "Estimator") +
    ggplot2::labs(x = "Radius c0", y = "Mean squared error") + theme_paper(9)
  save_figure(p, "break_even_random", cfg, 10, 5.5)
}

break_even_cstar <- function(info) {
  s <- info[grepl("^cstar:", info$name), ]
  s$method <- paste0("MM(c=", sub("cstar:", "", s$name), ")")
  s <- s[s$method %in% names(own_target) & is.finite(s$value), ]
  out <- stats::aggregate(value ~ graph + mode + method, s, median)
  names(out)[4] <- "cstar"
  out
}

# Crossing of two worst-case risks, linear in c0^2.
crossing <- function(c0, diff) {
  o <- order(c0)
  c0 <- c0[o]; diff <- diff[o]
  i <- which(diff[-length(diff)] < 0 & diff[-1] >= 0)[1]
  if (is.na(i)) return(NA_real_)
  t <- -diff[i] / (diff[i + 1] - diff[i])
  sqrt(c0[i]^2 + t * (c0[i + 1]^2 - c0[i]^2))
}

table_break_even <- function(w, info, cfg) {
  ns <- w[w$method == "NS", c("graph", "mode", "design", "c0", "mse")]
  mm <- w[grepl("^MM", w$method), ]
  d <- merge(mm, ns, by = c("graph", "mode", "design", "c0"),
             suffixes = c("", "_ns"))
  d$diff <- d$mse_ns - d$mse
  rows <- lapply(split(d, list(d$graph, d$mode, d$method, d$design),
                       drop = TRUE), function(g) {
    cbind(g[1, c("graph", "mode", "method", "design")],
          empirical = crossing(g$c0, g$diff))
  })
  emp <- stats::aggregate(empirical ~ graph + mode + method,
                          do.call(rbind, rows), median)
  out <- merge(emp, break_even_cstar(info), by = c("graph", "mode", "method"))
  out$ratio <- out$empirical / out$cstar
  out[] <- lapply(out, function(v) if (is.numeric(v)) fmt(v) else v)
  names(out) <- c("Graph", "Model", "Estimator", "Empirical crossing",
                  "$c^*$", "Ratio")
  write_table(out, "table_break_even", cfg, "Radius at which the worst-case mean
    squared error of the adjusted estimator falls below that of the
    non-spatial fit: Monte Carlo crossing against the break-even radius
    (medians over designs).", "tab:sim_break_even")
}

# ---- Shape, radius-profile, benchmark and timing experiments ---------------
shape_names <- c("icar", "bym2")

shape_regret_design <- function(job, cfg, specs) {
  spec <- specs[[job$graph]]
  dsg <- make_design(spec, cfg, job$design, job$level, job$share, job$decay)
  x <- combine_exposure(dsg$parts, dsg$share)
  qs <- lapply(stats::setNames(shape_names, shape_names), shape_qinv,
               spec = spec, phi = cfg$phi_bym2)
  geos <- lapply(qs, function(q) {
    design_geometry(spec, q, dsg$o, x, cfg$beta, control_matrix(spec))
  })
  rows <- list()
  for (mod in shape_names) for (tru in shape_names) {
    rows[[length(rows) + 1L]] <- regret_rows(geos[[mod]], geos[[tru]],
                                             qs[[mod]], qs[[tru]], cfg,
                                             mod, tru)
  }
  cbind(graph = job$graph, design = job$design, do.call(rbind, rows))
}

regret_rows <- function(geo_m, geo_u, q_m, q_u, cfg, mod, tru) {
  ratio <- q_u / q_m
  k1 <- min(ratio); k2 <- max(ratio)
  grid <- expand.grid(c_true = positive(cfg$c0_grid), mult = c(0.5, 1, 2))
  out <- lapply(seq_len(nrow(grid)), function(i) {
    c_true <- grid$c_true[i]
    c_m <- c_true * grid$mult[i]
    a <- weights_vector(geo_m, c_m)
    wc <- worst_case_risk(geo_u, to_coord(geo_u, a), c_true)
    nu <- max(c_m / (sqrt(k2) * c_true), sqrt(k2) * c_true / c_m)
    data.frame(model = mod, truth = tru, c_true = c_true, c_model = c_m,
               regret = wc / minimax_risk(geo_u, c_true),
               bound = nu^2 * k2 / k1, kappa1 = k1, kappa2 = k2)
  })
  do.call(rbind, out)
}

run_shape <- function(cfg, specs) {
  jobs <- make_jobs(cfg, specs, "shape_misspecification", NULL,
                    setdiff(names(specs), "brazil"))
  rbind_list(par_lapply(jobs, shape_regret_design, cores = cfg$cores,
                        cfg = cfg, specs = specs))
}

profile_examples <- function(job, cfg, specs, reps = 5L) {
  spec <- specs[[job$graph]]
  dsg <- make_design(spec, cfg, job$design, job$level, job$share, job$decay)
  cache <- new.env()
  out <- lapply(seq_len(nrow(job$settings)), function(k) {
    st <- job$settings[k, ]
    set.seed(task_seed(cfg$seed, dsg$id, st$id, 3))
    prep <- prepare_setting(spec, dsg, st, cfg, cache, weight_theta(job))
    ctx <- make_context(spec, dsg, st, prep, cfg, job)
    cbind(graph = job$graph, setting = st$id,
          profile_reps(ctx, reps, cfg$sigma_grid))
  })
  rbind_list(out)
}

profile_reps <- function(ctx, reps, grid) {
  geo <- ctx$geo
  rows <- lapply(seq_len(reps), function(r) {
    eps <- stats::rnorm(length(geo$w))
    eps <- eps - geo$q1 %*% crossprod(geo$q1, eps)
    yt <- geo$xt * ctx$beta + ctx$uvec + drop(eps)
    ehat <- drop(crossprod(geo$E, yt - geo$xt * sum(geo$xt * yt) /
                             sum(geo$xt^2)))
    nll <- vapply(log(grid), reml_nll, 0, y = ehat, xi = geo$xi, s = geo$s)
    data.frame(rep = r, sigma = grid, nll = nll - min(nll))
  })
  do.call(rbind, rows)
}

run_profiles <- function(cfg, specs) {
  st <- stack_settings(make_settings("aligned", c0 = c(0.5, 2), kappa = 1,
                                     target = "c0"),
                       make_settings("smooth", c0 = c(0.5, 2), rho = 0.6))
  jobs <- make_jobs(cfg, specs, "estimated_budget_profile", st,
                    setdiff(names(specs), "brazil"), mode = "gaussian",
                    designs = 1L)
  res <- rbind_list(lapply(jobs, profile_examples, cfg = cfg, specs = specs))
  res <- cbind(res, st[match(res$setting, st$id), c("fam", "c0")])
  save_result(res, "estimated_budget_profile", cfg)
}

decision_rep <- function(ctx, cfg) {
  y <- draw_counts(ctx$mu0)
  ctrl <- ctx$ctrl
  xv <- standardize(ctx$z)
  withv <- fit_ns(y, ctx$o, cbind(ctrl, ctx$x, xv))
  gamma <- withv$coef[ncol(ctrl) + 2L]
  w <- sqrt(reduced_weights(withv$mu))
  zv <- gamma * weighted_residual(ctrl, w, xv)
  delta <- drop(crossprod(ctx$spec$vecs, zv))
  cj <- sqrt(sum(delta^2 / ctx$qfit))
  fit <- fit_ns(y, ctx$o, cbind(ctrl, ctx$x))
  geo <- ctx$geo
  coord <- weights_coord(geo, cj / sqrt(fit$phi))
  ehat <- residual_coords(geo, fit$mu, y)
  est <- fit$coef[fit$p] + sum(coord * ehat)
  sd <- sqrt(fit$phi) * weights_norm(coord)
  bias <- weights_bias(geo, coord)
  u <- reduced_confounder(geo, drop(ctx$spec$vecs %*% delta))
  c(est = est, sd = sd, cj = cj, ccrit = critical_budget(est, sd, bias),
    align = alignment(geo, coord, u, cj))
}

decision_setting <- function(spec, dsg, st, cfg, job, cache) {
  set.seed(task_seed(cfg$seed, dsg$id, st$id, 1))
  prep <- prepare_setting(spec, dsg, st, cfg, cache, weight_theta(job))
  ctx <- make_context(spec, dsg, st, prep, cfg, job)
  set.seed(task_seed(cfg$seed, dsg$id, st$id, 2))
  res <- do.call(rbind, lapply(seq_len(job$reps), function(r) {
    decision_rep(ctx, cfg)
  }))
  robust <- vapply(cfg$benchmark_mult, function(m) {
    mean(res[, "ccrit"] >= m * res[, "cj"], na.rm = TRUE)
  }, numeric(1))
  data.frame(mult = cfg$benchmark_mult, p_robust = robust,
             cj = mean(res[, "cj"]), ccrit = stats::median(res[, "ccrit"]),
             align = mean(res[, "align"], na.rm = TRUE), n_rep = nrow(res))
}

run_decision_job <- function(job, cfg, specs) {
  cfg$beta <- job$beta
  spec <- specs[[job$graph]]
  dsg <- make_design(spec, cfg, job$design, job$level, job$share, job$decay)
  cache <- new.env()
  rows <- lapply(seq_len(nrow(job$settings)), function(k) {
    st <- job$settings[k, ]
    cbind(graph = job$graph, design = dsg$id, beta = job$beta,
          st[setdiff(names(st), "id")],
          decision_setting(spec, dsg, st, cfg, job, cache))
  })
  rbind_list(rows)
}

run_benchmark <- function(cfg, specs) {
  st <- stack_settings(make_settings("smooth", c0 = positive(cfg$c0_grid),
                                     rho = 0.6))
  jobs <- list()
  for (b in c(0, cfg$beta)) {
    jobs <- c(jobs, lapply(make_jobs(cfg, specs, "benchmark", st,
                                     setdiff(names(specs), "brazil")),
                           function(j) c(j, list(beta = b))))
  }
  rbind_list(par_lapply(jobs, run_decision_job, cores = cfg$cores,
                        cfg = cfg, specs = specs))
}

time_steps <- function(spec, cfg) {
  n <- spec$n
  set.seed(cfg$seed)
  o <- draw_offsets(n, cfg$level_main)
  x <- standardize(stats::rnorm(n))
  ctrl <- control_matrix(spec)
  qinv <- shape_qinv(spec, "icar")
  tm <- function(expr) unname(system.time(expr)[["elapsed"]])
  geo <- NULL
  t_geo <- tm(geo <- reduced_geometry(exp(o + cfg$beta * x), ctrl, x, spec$vecs,
                                      qinv))
  y <- draw_counts(exp(o + cfg$beta * x))
  fit <- NULL
  t_ns <- tm(fit <- fit_ns(y, o, cbind(ctrl, x)))
  ehat <- residual_coords(geo, fit$mu, y)
  t_curve <- tm(adjust_curve(geo, fit$coef[fit$p], ehat, fit$phi,
                             cfg$budgets, cfg$budgets))
  t_sigma <- tm(sigma_hat(geo, ehat, fit$phi))
  t_pen <- if (n <= 1000) tm(icar_fit(y, o, cbind(ctrl, x), spec, qinv)) else NA
  data.frame(graph = spec$label, n = n, geometry = t_geo, ns_fit = t_ns,
             curve = t_curve, sigma = t_sigma, penalized = t_pen,
             geometry_mb = as.numeric(utils::object.size(geo)) / 2^20)
}

# ---- Designs, contexts and replicates of one setting -------------------------
task_seed <- function(...) {
  v <- as.numeric(c(...))
  w <- c(1, 7919, 104729, 1299709, 15485863)[seq_along(v)]
  as.integer(sum(v * w) %% 2147483000)
}

make_design <- function(spec, cfg, id, level, share, decay, base = NULL) {
  set.seed(task_seed(cfg$seed, id, round(100 * decay), round(100 * share)))
  fixed <- !is.null(base$x)
  list(id = id, level = level, share = share, decay = decay,
       o = if (is.null(base$o)) draw_offsets(spec$n, level) else base$o,
       parts = if (fixed) NULL else draw_exposure_parts(spec, decay),
       x_fixed = base$x)
}

# Geometry depends on x only; it is shared by settings with the same x.
cached_geometry <- function(cache, key, build) {
  if (is.null(cache[[key]])) cache[[key]] <- build()
  cache[[key]]
}

prepare_setting <- function(spec, dsg, st, cfg, cache, theta = Inf) {
  qfit <- shape_qinv(spec, st$fit_shape, cfg$phi_bym2)
  qgen <- shape_qinv(spec, st$gen_shape, cfg$phi_bym2)
  needs_geo <- st$fam %in% c("aligned", "band", "orthogonal")
  z <- if (needs_geo) NULL else make_z(spec, NULL, qfit, qgen, st)
  x <- if (is.null(dsg$x_fixed)) {
    combine_exposure(dsg$parts, dsg$share, z, st$rho)
  } else dsg$x_fixed
  shared <- st$rho == 0 || !is.null(dsg$x_fixed)
  key <- if (shared) st$fit_shape else paste("setting", st$id)
  geo <- cached_geometry(cache, key, function() {
    design_geometry(spec, qfit, dsg$o, x, cfg$beta, control_matrix(spec),
                    theta)
  })
  if (is.null(z)) z <- make_z(spec, geo, qfit, qgen, st)
  list(x = x, z = z, geo = geo, qfit = qfit)
}

# Control columns: graph components, or the user's own (Brazil).
control_matrix <- function(spec) {
  if (is.null(spec$control)) spec$null_basis else spec$control
}

leak_radii <- function(ucoord, s, budgets) {
  vapply(budgets, function(b) leak_radius(ucoord, s, b), numeric(1))
}

# Size used in the reduced weights: the true size when the job asks for the
# negative binomial procedure, otherwise Inf (Poisson / quasi-Poisson).
weight_theta <- function(job) {
  if (isTRUE(job$nb_weights)) job$theta else Inf
}

make_context <- function(spec, dsg, st, prep, cfg, job) {
  geo <- prep$geo
  uvec <- reduced_vector(geo, prep$z)
  ucoord <- drop(crossprod(geo$E, uvec))
  c0 <- true_budget(spec, st, prep$qfit, prep$z)
  ks <- cfg$trunc_fixed[cfg$trunc_fixed <= spec$n / 4]
  grid <- cfg$trunc_grid[cfg$trunc_grid <= spec$n / 4]
  list(spec = spec, geo = geo, x = prep$x, z = prep$z, qfit = prep$qfit,
       o = dsg$o, ctrl = control_matrix(spec), st = st, cfg = cfg, job = job,
       beta = cfg$beta, alpha = cfg$alpha, uvec = uvec, ucoord = ucoord,
       theta_w = weight_theta(job),
       dispersion = if (is.null(job$dispersion)) "estimated" else
         job$dispersion,
       mu0 = exp(dsg$o + cfg$beta * prep$x + prep$z),
       budgets = cfg$budgets, c0 = c0, ks = ks, kgrid = grid,
       n_full = if (st$id %in% job$full_ids && dsg$id <= cfg$designs_full)
         job$n_full else 0L,
       leak = leak_radii(ucoord, geo$s, c(cfg$budgets, c0)),
       trunc = trunc_setup(geo, sort(unique(c(ks, grid)))))
}

trunc_setup <- function(geo, ks) {
  lapply(stats::setNames(ks, ks), function(k) {
    basis <- cbind(geo$xt, geo$tm[, seq_len(k), drop = FALSE])
    list(k = k, a = truncation_weights(geo, k), q = qr.Q(qr(basis)))
  })
}

trunc_gaussian_rows <- function(ctx, yt) {
  n <- length(yt) - ncol(ctx$geo$q1)
  one <- function(t) {
    rss <- sum(yt^2) - sum(crossprod(t$q, yt)^2)
    list(est = sum(t$a * yt), sd = sqrt(sum(t$a^2)), k = t$k,
         gcv = n * rss / (n - t$k - 1)^2)
  }
  fits <- lapply(ctx$trunc, one)
  pick <- function(f) stat_row(f$est, f$sd, extra = f$k, alpha = ctx$alpha)
  fixed <- lapply(fits[as.character(ctx$ks)], pick)
  cand <- fits[as.character(ctx$kgrid)]
  best <- cand[[which.min(vapply(cand, `[[`, 0, "gcv"))]]
  out <- do.call(rbind, c(fixed, list(pick(best))))
  rownames(out) <- c(sprintf("Trunc(K=%d)", ctx$ks), "Trunc(GCV)")
  out
}

curve_rows <- function(ctx, geo, beta_ns, ehat, phi, ucoord) {
  ct <- c(ctx$budgets, ctx$c0) / sqrt(phi)
  labels <- c(sprintf("MM(c=%g)", ctx$budgets), "MM(c0)")
  rows <- adjust_curve(geo, beta_ns, ehat, phi, ct, labels, ctx$alpha,
                       ucoord, ctx$leak)
  sig <- sigma_hat(geo, ehat, phi)
  srow <- adjust_curve(geo, beta_ns, ehat, phi, sig, "MM(sigma)",
                       ctx$alpha, ucoord)
  srow[, "extra"] <- sig
  rbind(rows, srow, `MM(sigma)-no-bias` = naive_row(srow[1, ], ctx$alpha))
}

rep_gaussian <- function(ctx) {
  geo <- ctx$geo
  eps <- stats::rnorm(length(geo$w))
  eps <- eps - geo$q1 %*% crossprod(geo$q1, eps)
  yt <- geo$xt * ctx$beta + ctx$uvec + drop(eps)
  nx <- sum(geo$xt^2)
  b_ns <- sum(geo$xt * yt) / nx
  ehat <- drop(crossprod(geo$E, yt - geo$xt * b_ns))
  ns <- stat_row(b_ns, 1 / sqrt(nx), alpha = ctx$alpha,
                 pred_bias = sum(geo$xt * ctx$uvec) / nx,
                 pred_sd = 1 / sqrt(nx))
  rbind(NS = ns, curve_rows(ctx, geo, b_ns, ehat, 1, ctx$ucoord),
        trunc_gaussian_rows(ctx, yt))
}

mean_geometry <- function(ctx, mu_hat) {
  reduced_geometry(mu_hat, ctx$ctrl, ctx$x, ctx$spec$vecs, ctx$qfit,
                   ctx$theta_w)
}

# Geometry of one count replicate: at the mean without the confounder
# ("baseline"), at the true mean that contains it ("true"), or at the
# fitted mean ("fitted"). The confounder's reduced field follows.
count_geometry <- function(ctx, fit) {
  geo <- switch(ctx$job$geometry,
    baseline = ctx$geo,
    true = mean_geometry(ctx, ctx$mu0),
    fitted = mean_geometry(ctx, fit$mu),
    stop("unknown geometry: ", ctx$job$geometry))
  ucoord <- if (ctx$job$geometry == "baseline") ctx$ucoord else {
    reduced_confounder(geo, ctx$z)
  }
  list(geo = geo, ucoord = ucoord)
}

rep_counts <- function(ctx, with_full) {
  y <- draw_counts(ctx$mu0, ctx$job$family, ctx$job$theta)
  fit <- fit_ns(y, ctx$o, cbind(ctx$ctrl, ctx$x), ctx$theta_w,
                ctx$dispersion)
  gu <- count_geometry(ctx, fit)
  geo <- gu$geo
  ucoord <- gu$ucoord
  ehat <- residual_coords(geo, fit$mu, y)
  rows <- rbind(ns_rows(fit, ctx$alpha),
                curve_rows(ctx, geo, fit$coef[fit$p], ehat, fit$phi, ucoord),
                truncation_rows(y, ctx$o, ctx$ctrl, ctx$x, ctx$spec$vecs,
                                ctx$ks, ctx$kgrid, ctx$alpha, ctx$theta_w,
                                ctx$dispersion))
  if (with_full) rows <- rbind(rows, try_rows(full_rows(ctx, y)))
  rows
}

full_rows <- function(ctx, y) {
  rows <- spatial_rows(y, ctx$o, ctx$ctrl, ctx$x, ctx$spec, ctx$cfg,
                       ctx$alpha)
  if (!is.null(ctx$spec$graph_file) && isTRUE(ctx$cfg$use_inla)) {
    rows <- rbind(rows, inla_rows(y, ctx$o, ctx$x, ctx$spec$graph_file))
  }
  rows
}

full_names <- function(ctx) {
  nm <- c("ICAR", "BYM2", "Spatial+-ICAR",
          sprintf("RSR(q=%d)", rsr_dims(ctx$spec$n, ncol(ctx$ctrl) + 1L)))
  if (!is.null(ctx$spec$graph_file) && isTRUE(ctx$cfg$use_inla)) {
    nm <- c(nm, "ICAR-INLA", "BYM2-INLA")
  }
  nm
}

fast_names <- function(ctx, mode) {
  lin <- c(sprintf("MM(c=%g)", ctx$budgets), "MM(c0)", "MM(sigma)",
           "MM(sigma)-no-bias")
  tr <- c(sprintf("Trunc(K=%d)", ctx$ks), "Trunc(GCV)")
  c(if (mode == "counts") c("NS", "NS-robust") else "NS", lin, tr)
}

# Arrange the rows of one replicate in a fixed order, NA when missing.
align_rows <- function(rows, names) {
  out <- matrix(NA_real_, length(names), length(stat_names),
                dimnames = list(names, stat_names))
  if (!is.null(rows)) {
    keep <- intersect(names, rownames(rows))
    out[keep, ] <- rows[keep, stat_names, drop = FALSE]
  }
  out
}

one_replicate <- function(r, ctx, names) {
  rows <- if (ctx$job$mode == "gaussian") {
    try_rows(rep_gaussian(ctx))
  } else {
    try_rows(rep_counts(ctx, r <= ctx$n_full))
  }
  align_rows(rows, names)
}

try_rows <- function(expr) tryCatch(expr, error = function(e) NULL)

mn <- function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)

summarise_method <- function(a, beta) {
  ok <- !is.na(a["est", ])
  g <- function(k) a[k, ok]
  hit <- function(lo, hi) mn(g(lo) <= beta & beta <= g(hi))
  est <- g("est")
  data.frame(
    n_rep = sum(ok), bias = mean(est) - beta, emp_sd = stats::sd(est),
    rmse = sqrt(mean((est - beta)^2)), mean_sd = mn(g("sd")),
    cover = hit("lo", "hi"), len = mn(g("hi") - g("lo")),
    cover_leak = hit("lo_leak", "hi_leak"),
    len_leak = mn(g("hi_leak") - g("lo_leak")),
    reject = mn(g("lo") > 0 | g("hi") < 0), extra = stats::median(g("extra")),
    pred_bias = mn(g("pred_bias")), pred_sd = mn(g("pred_sd")),
    nsif = mn(g("nsif")), edf = mn(g("edf")))
}

summarise_reps <- function(arr, beta) {
  out <- lapply(seq_len(dim(arr)[1]), function(m) {
    a <- matrix(arr[m, , ], nrow = length(stat_names),
                dimnames = list(stat_names, NULL))
    cbind(method = dimnames(arr)[[1]][m], summarise_method(a, beta))
  })
  do.call(rbind, out)
}

# Quantities that depend on the geometry and the confounder, not on y.
setting_info <- function(ctx) {
  geo <- ctx$geo
  b <- ctx$budgets
  ns <- geo$xi / sum(geo$xi^2)
  star <- vapply(b[b > 0], function(cc) {
    break_even(geo, ns, weights_coord(geo, cc))
  }, numeric(1))
  vals <- c(
    stats::setNames(nsif(geo, b), sprintf("nsif:%g", b)),
    stats::setNames(star, sprintf("cstar:%g", b[b > 0])),
    stats::setNames(vapply(b, minimax_risk, 0, geo = geo),
                    sprintf("minimax:%g", b)),
    stats::setNames(ctx$leak, c(sprintf("leak:%g", b), "leak:c0")),
    xnorm = sqrt(sum(geo$xi^2)), u_norm = sqrt(sum(ctx$ucoord^2)),
    bias_ns = sum(ns * ctx$ucoord))
  data.frame(name = names(vals), value = unname(vals))
}

# ---- Running tasks and experiments -------------------------------------------
run_setting <- function(spec, dsg, st, cfg, job, cache) {
  set.seed(task_seed(cfg$seed, dsg$id, st$id, 1))
  prep <- prepare_setting(spec, dsg, st, cfg, cache, weight_theta(job))
  ctx <- make_context(spec, dsg, st, prep, cfg, job)
  names <- c(fast_names(ctx, job$mode),
             if (job$mode == "counts" && ctx$n_full > 0) full_names(ctx))
  set.seed(task_seed(cfg$seed, dsg$id, st$id, 2))
  reps <- lapply(seq_len(job$reps), one_replicate, ctx = ctx, names = names)
  arr <- array(unlist(reps), dimnames = list(names, stat_names, NULL),
               dim = c(length(names), length(stat_names), job$reps))
  list(metrics = summarise_reps(arr, cfg$beta), info = setting_info(ctx))
}

label_columns <- function(df, job, dsg, st) {
  cbind(exp = job$tag, graph = job$graph, mode = job$mode,
        geometry = job$geometry, family = job$family,
        dispersion = if (is.null(job$dispersion)) "estimated" else
          job$dispersion,
        nb_weights = isTRUE(job$nb_weights), level = job$level,
        share = job$share, decay = job$decay, design = dsg$id,
        st[setdiff(names(st), "id")], setting = st$id, df,
        row.names = NULL, stringsAsFactors = FALSE)
}

# One task is one setting of one job. Seeds depend only on the design and
# setting ids, so the split into tasks does not change any result.
job_tasks <- function(job) {
  lapply(seq_len(nrow(job$settings)), function(k) list(job = job, row = k))
}

run_task <- function(task, cfg, specs) {
  job <- task$job
  st <- job$settings[task$row, ]
  spec <- specs[[job$graph]]
  dsg <- make_design(spec, cfg, job$design, job$level, job$share,
                     job$decay, spec$base)
  cache <- if (is.null(spec$cache)) new.env() else spec$cache
  res <- tryCatch(run_setting(spec, dsg, st, cfg, job, cache),
    error = function(e) {
      stop(sprintf("setting %d (%s) of %s on %s failed: %s", st$id, st$fam,
                   job$tag, job$graph, conditionMessage(e)), call. = FALSE)
    })
  list(metrics = label_columns(res$metrics, job, dsg, st),
       info = label_columns(res$info, job, dsg, st))
}

run_jobs <- function(jobs, cfg, specs) {
  tasks <- unlist(lapply(jobs, job_tasks), recursive = FALSE)
  msg(length(tasks), " tasks on ", cfg$cores, " cores")
  res <- par_lapply(tasks, run_task, cfg = cfg, specs = specs,
                    cores = cfg$cores)
  list(metrics = rbind_list(lapply(res, `[[`, "metrics")),
       info = rbind_list(lapply(res, `[[`, "info")))
}

build_specs <- function(cfg) {
  syn <- setdiff(cfg$graphs, "brazil")
  specs <- lapply(stats::setNames(syn, syn), build_spec, cfg = cfg)
  specs <- specs[!vapply(specs, is.null, logical(1))]
  if ("brazil" %in% cfg$graphs) specs$brazil <- brazil_spec(cfg)
  specs
}

# The Brazilian geometry is the same for every design: compute it once,
# before the jobs are forked, and share it through the spec's cache.
warm_brazil_cache <- function(spec, cfg) {
  dsg <- make_design(spec, cfg, 0L, 1, 1, 1, spec$base)
  spec$cache[["icar"]] <- design_geometry(
    spec, shape_qinv(spec, "icar"), dsg$o, dsg$x_fixed, cfg$beta,
    control_matrix(spec))
  invisible(NULL)
}

save_result <- function(res, name, cfg) {
  ensure_dir(cfg$out_dir)
  saveRDS(res, file.path(cfg$out_dir, sprintf("results_%s.rds", name)))
  invisible(res)
}

save_graph_summary <- function(specs, cfg) {
  out <- rbind_list(lapply(specs, graph_summary))
  utils::write.csv(out, file.path(ensure_dir(cfg$out_dir), "graphs.csv"),
                   row.names = FALSE)
}

run_experiments <- function(cfg, which = NULL) {
  specs <- build_specs(cfg)
  save_graph_summary(specs, cfg)
  jobs <- experiment_jobs(cfg, specs)
  if (!is.null(which)) {
    base <- sub("_(fitted|true|quasi|nb|spatial_exposure)$", "", names(jobs))
    jobs <- jobs[names(jobs) %in% which | base %in% which]
  }
  if ("coverage_brazil" %in% names(jobs)) warm_brazil_cache(specs$brazil, cfg)
  for (nm in names(jobs)) {
    msg("running ", nm)
    save_result(run_jobs(jobs[[nm]], cfg, specs), nm, cfg)
  }
  extras <- list(shape_misspecification = run_shape, benchmark = run_benchmark)
  for (nm in names(extras)) {
    if (is.null(which) || nm %in% which) {
      msg("running ", nm)
      save_result(extras[[nm]](cfg, specs), nm, cfg)
    }
  }
  if (is.null(which) || "timing" %in% which) {
    save_result(rbind_list(lapply(specs, time_steps, cfg = cfg)), "timing", cfg)
  }
  if (is.null(which) || "estimated_budget" %in% which) run_profiles(cfg, specs)
  invisible(NULL)
}

# ---- Settings, jobs and the experiment lists ---------------------------------
setting_defaults <- list(
  fam = "aligned", c0 = 0, kappa = 0, target = "c0", rho = 0, n_vec = 20L,
  nu = 0.5, range = 0.1, kind = "hotspot", size = 5L, amp = 0.25,
  fit_shape = "icar", gen_shape = "icar")

# Full factorial of the given values, other columns at their defaults.
make_settings <- function(fam, ...) {
  grid <- expand.grid(c(list(fam = fam), list(...)),
                      stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  for (nm in setdiff(names(setting_defaults), names(grid))) {
    grid[[nm]] <- setting_defaults[[nm]]
  }
  grid[names(setting_defaults)]
}

stack_settings <- function(...) {
  out <- rbind_list(list(...))
  out$id <- seq_len(nrow(out))
  out
}

make_jobs <- function(cfg, specs, tag, settings, graphs = names(specs),
                      mode = "counts", geometry = "baseline",
                      family = "poisson", theta = NA_real_,
                      dispersion = NULL,
                      nb_weights = FALSE, levels = cfg$level_main,
                      shares = cfg$share_main, decays = cfg$decay_main,
                      reps = cfg$reps, designs = cfg$designs,
                      full_ids = integer(0)) {
  if (is.null(dispersion)) {
    dispersion <- if (family == "poisson" || nb_weights) "fixed" else
      "estimated"
  }
  grid <- expand.grid(graph = graphs, level = levels, share = shares,
                      decay = decays, design = seq_len(designs),
                      stringsAsFactors = FALSE)
  lapply(seq_len(nrow(grid)), function(i) {
    g <- as.list(grid[i, ])
    c(g, list(tag = tag, settings = settings, mode = mode,
              geometry = geometry, family = family, theta = theta,
              dispersion = dispersion,
              nb_weights = nb_weights, reps = reps,
              n_full = if (length(full_ids)) cfg$reps_full[[g$graph]] else 0L,
              full_ids = full_ids))
  })
}

positive <- function(v) v[v > 0]

settings_gaussian <- function(cfg) {
  stack_settings(make_settings("aligned", c0 = positive(cfg$c0_grid),
                               kappa = cfg$kappas,
                               target = c("c0", "ns")))
}

settings_linearization <- function(cfg) {
  stack_settings(make_settings("aligned", c0 = positive(cfg$c0_grid),
                               kappa = c(-1, 0, 1), target = "c0"))
}

settings_cost <- function(cfg) {
  stack_settings(make_settings("independent", c0 = c(0, 0.5)),
                 make_settings("orthogonal", c0 = 0.5))
}

settings_break_even <- function(cfg) {
  targets <- c("ns", "mm:0.1", "mm:0.5", "mm:2", "trunc:5", "trunc:20")
  stack_settings(
    make_settings("aligned", c0 = cfg$c0_fine, kappa = 1, target = targets),
    make_settings("smooth", c0 = cfg$c0_fine, rho = 0.6),
    make_settings("band", c0 = cfg$c0_fine, n_vec = c(5L, 20L)))
}

settings_coverage <- function(cfg) {
  c3 <- c(0.5, 1, 2)
  stack_settings(
    make_settings("aligned", c0 = c3, kappa = c(-1, 0, 1),
                  target = c("ns", "c0")),
    make_settings("smooth", c0 = c3, rho = c(0, 0.3, 0.6, 0.9)),
    make_settings("band", c0 = c3, n_vec = c(5L, 20L, 50L)),
    make_settings("matern", c0 = c3, nu = c(0.5, 1.5), range = c(0.1, 0.3)),
    make_settings("local", kind = "hotspot", size = c(1L, 5L, 25L),
                  amp = c(0.25, 0.5)),
    make_settings("local", kind = "region", amp = c(0.25, 0.5)),
    make_settings("independent", c0 = c(0, 1)),
    make_settings("orthogonal", c0 = c(0.5, 1)),
    make_settings("mixture", c0 = c3))
}

# Settings that also receive the (slow) penalized competitors.
full_subset <- function(st) {
  keep <- (st$fam == "aligned" & st$c0 == 1 & st$target == "ns") |
    (st$fam == "smooth" & st$c0 == 1 & st$rho %in% c(0, 0.6)) |
    (st$fam == "band" & st$c0 == 1 & st$n_vec == 20) |
    (st$fam == "matern" & st$c0 == 1) |
    (st$fam == "local" & st$amp == 0.5) |
    (st$fam %in% c("independent", "orthogonal", "mixture") & st$c0 == 1) |
    (st$fam == "independent" & st$c0 == 0)
  st$id[keep]
}

settings_nb <- function(cfg) {
  stack_settings(
    make_settings("aligned", c0 = 1, kappa = c(-1, 0, 1), target = "ns"),
    make_settings("smooth", c0 = 1, rho = 0.6),
    make_settings("orthogonal", c0 = 1), make_settings("independent", c0 = 0))
}

settings_spatial_exposure <- function(cfg) {
  stack_settings(make_settings("aligned", c0 = c(0.5, 1), kappa = c(0, 1),
                               target = "ns"),
                 make_settings("smooth", c0 = c(0.5, 1), rho = 0.9))
}

settings_leakage <- function(cfg) {
  stack_settings(
    make_settings("local", kind = "hotspot", size = c(1L, 5L, 25L),
                  amp = c(0.1, 0.25, 0.5)),
    make_settings("local", kind = "region", amp = c(0.1, 0.25, 0.5)))
}

settings_budget <- function(cfg) {
  pos <- positive(cfg$c0_grid)
  stack_settings(make_settings("aligned", c0 = pos, kappa = c(0, 1),
                               target = "c0"),
                 make_settings("smooth", c0 = pos, rho = 0.6))
}

settings_brazil <- function(cfg) {
  c3 <- c(0.5, 1, 2)
  stack_settings(
    make_settings("aligned", c0 = c3, kappa = c(-1, 0, 1),
                  target = c("ns", "c0")),
    make_settings("smooth", c0 = c3),
    make_settings("band", c0 = c3, n_vec = 20L),
    make_settings("local", kind = "hotspot", size = c(1L, 25L), amp = 0.25),
    make_settings("local", kind = "region", amp = 0.25))
}

# Negative binomial data, fitted with Poisson (quasi) weights and with the
# negative binomial weights at the true size.
nb_jobs <- function(cfg, specs, graphs) {
  grid <- expand.grid(theta = c(2, 10), nb = c(FALSE, TRUE))
  do.call(c, lapply(seq_len(nrow(grid)), function(i) {
    make_jobs(cfg, specs, "coverage_nb", settings_nb(cfg), graphs,
              family = "nb", theta = grid$theta[i], nb_weights = grid$nb[i])
  }))
}

# The count linearization is run at the mean without the confounder
# (baseline), at the true mean, with the Pearson dispersion, and at the
# fitted mean.
linearization_jobs <- function(cfg, specs, syn) {
  one <- function(...) {
    make_jobs(cfg, specs, "count_linearization",
              settings_linearization(cfg), syn, levels = cfg$levels, ...)
  }
  list(count_linearization = one(),
       count_linearization_true = one(geometry = "true",
                                      reps = cfg$reps_alt),
       count_linearization_quasi = one(dispersion = "estimated",
                                       reps = cfg$reps_alt),
       count_linearization_fitted = one(geometry = "fitted",
                                        reps = cfg$reps_alt))
}

# Gaussian and count runs of the same settings.
both_modes <- function(cfg, specs, tag, settings, syn, ...) {
  c(make_jobs(cfg, specs, tag, settings, syn, mode = "gaussian", ...),
    make_jobs(cfg, specs, tag, settings, syn, ...))
}

experiment_jobs <- function(cfg, specs) {
  syn <- setdiff(names(specs), "brazil")
  st_cov <- settings_coverage(cfg)
  jobs <- c(
    list(gaussian_exact = make_jobs(
           cfg, specs, "gaussian_exact", settings_gaussian(cfg), syn,
           mode = "gaussian", shares = c(cfg$share_main, 0.8))),
    linearization_jobs(cfg, specs, syn),
    list(
      adjustment_cost = c(
        make_jobs(cfg, specs, "adjustment_cost", settings_cost(cfg), syn,
                  mode = "gaussian", shares = cfg$shares, decays = c(1, 2)),
        make_jobs(cfg, specs, "adjustment_cost", settings_cost(cfg), syn,
                  shares = cfg$shares, decays = c(1, 2))),
      break_even = both_modes(cfg, specs, "break_even",
                              settings_break_even(cfg), syn),
      coverage = make_jobs(cfg, specs, "coverage", st_cov, syn,
                           full_ids = full_subset(st_cov)),
      coverage_nb = nb_jobs(cfg, specs, syn),
      coverage_spatial_exposure = make_jobs(
        cfg, specs, "coverage_spatial_exposure",
        settings_spatial_exposure(cfg), syn, shares = c(0.95, 1)),
      leakage = both_modes(cfg, specs, "leakage", settings_leakage(cfg),
                           syn),
      estimated_budget = make_jobs(cfg, specs, "estimated_budget",
                                   settings_budget(cfg), syn)))
  if ("brazil" %in% names(specs)) {
    jobs$coverage_brazil <- make_jobs(
      cfg, specs, "coverage_brazil", settings_brazil(cfg), "brazil",
      designs = cfg$designs_brazil, reps = cfg$reps_brazil)
  }
  jobs
}

# ---- General helpers ---------------------------------------------------------
msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

# INLA is not on CRAN: install it from its own repository when missing.
install_inla <- function() {
  if (requireNamespace("INLA", quietly = TRUE)) return(invisible(TRUE))
  msg("INLA not found, installing from inla.r-inla-download.org")
  repos <- getOption("repos")
  if (is.null(repos) || identical(unname(repos["CRAN"]), "@CRAN@")) {
    repos <- c(CRAN = "https://cloud.r-project.org")
  }
  repos <- c(repos, INLA = "https://inla.r-inla-download.org/R/stable")
  try(utils::install.packages("INLA", repos = repos, dependencies = TRUE),
      silent = TRUE)
  ok <- requireNamespace("INLA", quietly = TRUE)
  if (!ok) msg("INLA could not be installed; Bayesian fits will be skipped")
  invisible(ok)
}

par_lapply <- function(x, fun, cores = 1L, ...) {
  if (cores > 1L && .Platform$OS.type == "unix") {
    out <- parallel::mclapply(x, fun, ..., mc.cores = cores,
                              mc.preschedule = FALSE)
    failed <- vapply(out, inherits, logical(1), what = "try-error")
    if (any(failed)) stop(as.character(out[[which(failed)[1]]]))
    return(out)
  }
  lapply(x, fun, ...)
}

rbind_list <- function(x) do.call(rbind, x[!vapply(x, is.null, logical(1))])

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

standardize <- function(x) {
  x <- x - mean(x)
  x / sqrt(mean(x^2))
}

# Noncentral chi-square form of the (1 - alpha) quantile of |N(t, 1)|.
cva <- function(t, alpha = 0.05) {
  sqrt(stats::qchisq(1 - alpha, df = 1, ncp = t^2))
}

# Inverse of cva on [z_alpha, Inf): the t with cva(t) = h.
cva_inverse <- function(h, alpha = 0.05) {
  z <- stats::qnorm(1 - alpha / 2)
  one <- function(hh) {
    if (hh <= z) return(0)
    root <- function(t) stats::pnorm(hh - t) - stats::pnorm(-hh - t) -
      (1 - alpha)
    stats::uniroot(root, c(0, hh + 10), tol = 1e-10)$root
  }
  vapply(h, one, numeric(1))
}

# ---- Identity and sanity checks of the engine (simulations.R check) ----------
check_engine <- function() {
  failures <- 0L
  check <- function(name, ok) {
    cat(sprintf("%-62s %s\n", name, if (isTRUE(ok)) "ok" else "FAILED"))
    if (!isTRUE(ok)) failures <<- failures + 1L
  }
  close_to <- function(a, b, tol = 1e-8) max(abs(a - b)) < tol

  cfg <- sim_config("pilot")
  g <- lattice_graph(10)
  spec <- c(g, graph_spectrum(g$adj), list(label = "lattice"))
  n <- spec$n
  set.seed(11)
  o <- draw_offsets(n, 10)
  x <- standardize(stats::rnorm(n))
  qinv <- shape_qinv(spec, "icar")
  geo <- design_geometry(spec, qinv, o, x, cfg$beta, spec$null_basis)

  # ---- Graph and prior ----------------------------------------------------
  lap <- diag(rowSums(as.matrix(g$adj))) - as.matrix(g$adj)
  check("spectrum reconstructs the Laplacian",
        close_to(spec$vecs %*% diag(spec$lambda) %*% t(spec$vecs), lap))
  marg <- rowSums(sweep(spec$vecs^2, 2, 1 / (spec$g * spec$lambda), "*"))
  check("scaled ICAR has unit geometric mean marginal variance",
        close_to(exp(mean(log(marg))), 1))

  # Disconnected graph: a path, a triangle and an isolated node.
  adj_d <- Matrix::Matrix(0, 8, 8, sparse = TRUE)
  for (e in list(c(1, 2), c(2, 3), c(3, 4), c(5, 6), c(6, 7), c(5, 7)))
    adj_d[e[1], e[2]] <- adj_d[e[2], e[1]] <- 1
  sd <- graph_spectrum(adj_d)
  check("disconnected graph: three components, U orthonormal",
        max(sd$comp) == 3 && ncol(sd$vecs) == 5 &&
          close_to(crossprod(sd$vecs), diag(5)))
  marg_d <- rowSums(sweep(sd$vecs^2, 2, 1 / (sd$g * sd$lambda), "*"))
  check("each component has unit geometric mean marginal variance",
        close_to(exp(mean(log(marg_d[1:4]))), 1) &&
          close_to(exp(mean(log(marg_d[5:7]))), 1))

  # ---- Minimax weights ----------------------------------------------------
  c_grid <- c(0, 0.1, 0.5, 2)
  coords <- lapply(c_grid, function(cc) weights_coord(geo, cc))
  check("weights are unbiased for beta: a'x = 1",
        all(vapply(coords, function(a) abs(sum(a * geo$xi) - 1), 0) < 1e-10))
  check("worst-case risk of a_c equals 1 / f(c)",
        close_to(vapply(seq_along(c_grid), function(i) {
          worst_case_risk(geo, coords[[i]], c_grid[i])
        }, 0), vapply(c_grid, minimax_risk, 0, geo = geo)))
  check("NSIF(0) = 1 and NSIF decreases in c",
        abs(nsif(geo, 0) - 1) < 1e-12 && all(diff(nsif(geo, c_grid)) <= 0))
  check("minimax risk is nondecreasing in c",
        all(diff(vapply(c_grid, minimax_risk, 0, geo = geo)) >= -1e-12))

  perturb_ok <- all(vapply(1:20, function(i) {
    d <- stats::rnorm(n); d <- d - sum(d * geo$xi) * geo$xi / sum(geo$xi^2)
    a <- coords[[3]] + 0.1 * d
    worst_case_risk(geo, a, 0.5) >= minimax_risk(geo, 0.5) - 1e-10
  }, logical(1)))
  check("no feasible perturbation beats the minimax weights", perturb_ok)

  ns <- geo$xi / sum(geo$xi^2)
  cstar <- break_even(geo, ns, coords[[3]])
  check("worst-case risks cross at the break-even radius",
        !is.na(cstar) && abs(worst_case_risk(geo, ns, cstar) -
          worst_case_risk(geo, coords[[3]], cstar)) < 1e-8)

  # ---- Intervals ----------------------------------------------------------
  check("cva(0) is the normal quantile",
        abs(cva(0) - stats::qnorm(0.975)) < 1e-10)
  t <- c(0.3, 1, 2.5)
  check("cva_inverse inverts cva", close_to(cva_inverse(cva(t)), t, 1e-6))
  check("cva grows with the bias, slower than one-for-one",
        all(diff(cva(c(0, 1, 2, 3))) > 0 & diff(cva(c(0, 1, 2, 3))) < 1))

  # ---- Truncation and leakage ---------------------------------------------
  a_k <- truncation_weights(geo, 5L)
  check("truncated weights annihilate the first K directions",
        max(abs(crossprod(geo$tm[, 1:5], a_k))) < 1e-8 &&
          abs(sum(a_k * geo$xt) - 1) < 1e-10)
  u <- stats::rnorm(n) * sqrt(geo$s)
  check("leak radius is zero inside the class",
        leak_radius(u, geo$s, 1e3) < 1e-10)
  r <- leak_radius(u, geo$s, 0.5)
  pos <- geo$s > 1e-12
  feasible <- replicate(200, {
    v <- stats::rnorm(n) * sqrt(geo$s); v[!pos] <- 0
    v * 0.5 / sqrt(sum(v[pos]^2 / geo$s[pos]))
  })
  dist <- apply(feasible, 2, function(v) sqrt(sum((u - v)^2)))
  check("leak radius is a lower bound of distances to feasible points",
        r <= min(dist) + 1e-8)

  # ---- Exposure, orthogonal field, alignment, weights ---------------------
  parts <- draw_exposure_parts(spec, 2)
  zz <- z_smooth(spec, qinv, 1)
  check("realized correlation of x and z equals rho exactly",
        abs(stats::cor(combine_exposure(parts, 1, zz, 0.6), zz) - 0.6) < 1e-10)
  z_o <- z_orthogonal(spec, qinv, geo, 1)
  check("orthogonal field leaves x'u at zero",
        abs(sum(geo$xt * reduced_vector(geo, z_o))) < 1e-8)
  z_a <- z_aligned(spec, qinv, geo, weights_vector(geo, 0.5), 1, 0.5)
  check("alignment of an aligned field equals kappa",
        abs(alignment(geo, coords[[3]], reduced_confounder(geo, z_a), 1) -
              0.5) < 1e-8)
  check("reduced weights: Poisson and negative binomial",
        close_to(reduced_weights(c(1, 4), Inf), c(1, 4)) &&
          close_to(reduced_weights(c(1, 4), 2), c(1 / 1.5, 4 / 3)))
  y_nb <- draw_counts(exp(o + cfg$beta * x), "nb", 5)
  nb_fit <- fit_ns(y_nb, o, cbind(spec$null_basis, x), theta = 5)
  check("negative binomial fit returns finite estimates",
        is.finite(nb_fit$coef[2]) && nb_fit$phi > 0)
  cells <- c(0.93, 0.95, 0.97, 0.94)
  check("coverage error uses the spread between designs",
        abs(design_se(cells, 400, mean(cells)) - stats::sd(cells) / 2) < 1e-12)

  # ---- Gaussian replicate: one-step form equals the weighted sum ----------
  st <- cbind(make_settings("aligned", c0 = 1, kappa = 1, target = "c0"),
              id = 1L)
  dsg <- make_design(spec, cfg, 1L, 10, 0.5, 2)
  job <- list(mode = "gaussian", geometry = "baseline", family = "poisson",
              theta = NA, reps = 10L, n_full = 0L, full_ids = integer(0))
  cache <- new.env()
  prep <- prepare_setting(spec, dsg, st, cfg, cache)
  ctx <- make_context(spec, dsg, st, prep, cfg, job)
  yt <- geo$xt * 0.2 + ctx$uvec
  b_ns <- sum(ctx$geo$xt * yt) / sum(ctx$geo$xt^2)
  ehat <- drop(crossprod(ctx$geo$E, yt - ctx$geo$xt * b_ns))
  row <- adjust_curve(ctx$geo, b_ns, ehat, 1, 0.5, "m")
  a_vec <- weights_vector(ctx$geo, 0.5)
  check("one-step estimator equals a'y in the linear model",
        abs(row[1, "est"] - sum(a_vec * yt)) < 1e-8)

  # Coverage at kappa = 1: exact in the reduced model
  reps <- 3000
  cover <- mean(replicate(reps, {
    rows <- rep_gaussian(ctx)["MM(c0)", ]
    rows["lo"] <= cfg$beta && cfg$beta <= rows["hi"]
  }))
  se <- sqrt(0.95 * 0.05 / reps)
  check("coverage at the worst alignment is at least nominal",
        cover > 0.95 - 3 * se)

  # ---- Estimation of sigma ------------------------------------------------
  sig <- 0.4
  est <- vapply(1:30, function(i) {
    u0 <- sig * sqrt(geo$s) * stats::rnorm(n)
    reml_sigma(u0 + stats::rnorm(n), geo$xi, geo$s)
  }, 0)
  check("REML estimate of sigma is centred near the truth",
        abs(median(est) - sig) < 0.25)

  # ---- Counts and penalized fit -------------------------------------------
  mu <- exp(o + cfg$beta * x)
  y <- draw_counts(mu)
  xmat <- cbind(spec$null_basis, x)
  ns_fit <- fit_ns(y, o, xmat)
  check("quasi-Poisson fit returns finite estimates",
        is.finite(ns_fit$coef[2]) && ns_fit$phi > 0)
  check("strict Poisson fit fixes the dispersion at one",
        fit_ns(y, o, xmat, dispersion = "fixed")$phi == 1)
  rb <- rsr_basis(spec, xmat, 1000L)
  check("RSR basis keeps only positive Moran eigenvalues", ncol(rb) < 1000L)
  set.seed(3)
  check("hotspot field has the requested number of raised nodes",
        sum(z_local(spec, spec$adj, "hotspot", 5L, 1) != 0) >= 5L)
  parts0 <- draw_exposure_parts(spec, 2)
  x0 <- combine_exposure(parts0, 0.5, rep(0, n), 0.6)
  check("a null field with a requested correlation leaves x finite",
        all(is.finite(x0)))
  z_loc <- z_local(spec, spec$adj, "region", 1L, 0.5)
  check("true radius of a localized field is its radius",
        isTRUE(all.equal(true_budget(spec, list(fam = "local", c0 = 0), qinv,
                                       z_loc),
                         field_radius(spec, qinv, z_loc))))
  check("planar coordinates for a graph without any",
        identical(dim(graph_coordinates(spec$adj)), c(n, 2L)))
  pen <- icar_fit(y, o, xmat, spec, qinv)$row
  check("penalized ICAR fit returns finite estimate and interval",
        all(is.finite(pen[c("est", "sd", "lo", "hi")])))
  m <- ncol(spec$vecs)
  tiny <- pirls(y, o, cbind(xmat, spec$vecs), ncol(xmat), diag(1e8, m),
                c(ns_fit$coef, rep(0, m)))
  check("penalized fit with a huge penalty reduces to the Poisson fit",
        abs(tiny$theta[2] - ns_fit$coef[2]) < 1e-4)

  # ---- Reproducibility ----------------------------------------------------
  job$mode <- "counts"
  run_once <- function() run_setting(spec, dsg, st, cfg, job, new.env())
  check("a setting is reproducible given its seed",
        isTRUE(all.equal(run_once()$metrics, run_once()$metrics)))

  if (failures > 0L) stop(failures, " check(s) failed")
  cat("all checks passed\n")
}
