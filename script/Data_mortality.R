# Municipal mortality, Brazil 2022: deaths and expected deaths by cause

#------ Packages and paths ----------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, sidrar, stringi, sf, here)
if (!requireNamespace("microdatasus", quietly = TRUE)) {
  pacman::p_load_gh("rfsaldanha/microdatasus")
}
year <- 2022
raw_dir <- here("data", "raw")
dir.create(file.path(raw_dir, "sim"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(raw_dir, "ibge"), recursive = TRUE, showWarnings = FALSE)

# Nonspecific circulatory causes among the leading garbage codes listed by
# França et al. (2020): heart failure, essential hypertension, pulmonary
# embolism, unspecified stroke and unspecified cardiomyopathy. This is not the
# complete garbage-code classification.
nonspecific_pattern <- "^(I50|I10|I26|I64|I678|I679|I694|I698|I420|I429)"

uf <- c(RO = 11, AC = 12, AM = 13, RR = 14, PA = 15, AP = 16, TO = 17,
        MA = 21, PI = 22, CE = 23, RN = 24, PB = 25, PE = 26, AL = 27,
        SE = 28, BA = 29, MG = 31, ES = 32, RJ = 33, SP = 35, PR = 41,
        SC = 42, RS = 43, MS = 50, MT = 51, GO = 52, DF = 53)

# Outcomes: ICD-10 pattern of the underlying cause and age range in years.
spec <- data.frame(
  key = c("cvd", "cvd_premature", "ihd", "stroke", "hypertensive",
          "diabetes", "neoplasm", "homicide", "ill_defined", "all"),
  pattern = c("^I[0-9]{2}", "^I[0-9]{2}", "^I2[0-5]", "^I6[0-9]", "^I1[0-5]",
              "^E1[0-4]", "^C[0-9]{2}", "^(X8[5-9]|X9[0-9]|Y0[0-9]|Y3[56])",
              "^R[0-9]{2}", "."),
  min = c(0, 30, 0, 0, 0, 0, 0, 0, 0, 0),
  max = c(Inf, 69, rep(Inf, 8)))

# Download once and keep the raw file; later runs read it back. A file that
# lacks a required column is downloaded again. A failed or empty download is
# retried and then stops the script, so no incomplete state is ever analysed.
options(timeout = 900)
cached <- function(file, fun, need = NULL, tries = 5) {
  if (file.exists(file)) {
    out <- readRDS(file)
    if (is.data.frame(out) && all(need %in% names(out))) return(out)
  }
  for (k in seq_len(tries)) {
    out <- tryCatch(fun(), error = function(e) {
      message("Attempt ", k, " failed: ", conditionMessage(e))
      NULL
    })
    if (is.data.frame(out) && nrow(out) > 0 && all(need %in% names(out))) {
      saveRDS(out, file)
      return(out)
    }
    Sys.sleep(5 * k)
  }
  stop("Download failed after ", tries, " attempts: ", basename(file))
}

age_labels <- c(paste(seq(0, 95, 5), seq(5, 100, 5) - 1, sep = "-"), "100+")
to_group <- function(age) age_labels[ifelse(age >= 100, 21, age %/% 5 + 1)]

#------ Deaths: SIM-DO, all causes, municipality of residence -------------
sim_vars <- c("CODMUNRES", "DTOBITO", "CAUSABAS", "IDADE", "SEXO", "TIPOBITO")
sim <- lapply(names(uf), function(s) {
  cached(file.path(raw_dir, "sim", paste0("sim_", s, ".rds")), function() {
    message("SIM-DO ", year, ": ", s)
    microdatasus::fetch_datasus(
      year_start = year, year_end = year, uf = s,
      information_system = "SIM-DO", vars = sim_vars)
  }, need = sim_vars)
})
sim <- bind_rows(sim)

# IDADE: 4xx is age in years, 5xx is 100 years or more, a code below 400 is a
# period shorter than a year, and 000 or 999 is an ignored age.
sim <- sim |>
  transmute(code6 = substr(as.character(CODMUNRES), 1, 6),
            year = substr(stri_pad_left(DTOBITO, 8, "0"), 5, 8),
            fetal = TIPOBITO %in% "1",
            cause = toupper(trimws(as.character(CAUSABAS))),
            sex = unname(c("1" = "M", "2" = "F")[as.character(SEXO)]),
            code_age = suppressWarnings(as.integer(IDADE))) |>
  mutate(age = case_when(code_age %in% c(0L, 999L) ~ NA_integer_,
                         code_age < 400 ~ 0L,
                         code_age >= 400 & code_age < 600 ~ code_age - 400L,
                         TRUE ~ NA_integer_)) |>
  filter(year == as.character(.env$year), !fetal)
msg <- function(...) message(paste0(...))
msg("SIM-DO records of non-fetal deaths: ", nrow(sim))

#------ Population: Census 2022, table 9514 --------------------------------
sidra_get <- function(...) {
  tibble::as_tibble(sidrar::get_sidra(..., format = 4,
                                      value_type = "character"))
}
col <- function(x, pattern) {
  nm <- names(x)[grepl(pattern, tolower(stri_trans_general(names(x),
                                                              "Latin-ASCII")))]
  if (length(nm) != 1) stop("SIDRA column not found: ", pattern)
  nm
}
clean_age <- function(x) {
  x <- tolower(stri_trans_general(trimws(x), "Latin-ASCII"))
  x <- gsub(" anos", "", gsub(" a ", "-", x, fixed = TRUE), fixed = TRUE)
  ifelse(x == "100 ou mais", "100+", x)
}

# One municipality fixes the codes of the 21 age groups.
probe <- cached(file.path(raw_dir, "ibge", "9514_probe.rds"), function() {
  sidra_get(x = 9514, variable = 93, period = as.character(year),
            geo = "City", geo.filter = list(City = 3550308),
            classific = c("c2", "c287", "c286"),
            category = list(4, "all", 113635))
})
groups <- probe |>
  transmute(label = clean_age(.data[[col(probe, "^idade$")]]),
            code = as.character(.data[[col(probe, "^idade.*codigo")]])) |>
  filter(label %in% age_labels) |>
  distinct()
if (nrow(groups) != 21) stop("Expected 21 age groups, found ", nrow(groups))

pop <- lapply(uf, function(s) {
  cached(file.path(raw_dir, "ibge", paste0("pop_", s, ".rds")), function() {
    message("Census 2022, table 9514: UF ", s)
    sidra_get(x = 9514, variable = 93, period = as.character(year),
              geo = "City", geo.filter = list(State = s),
              classific = c("c2", "c287", "c286"),
              category = list(c(4, 5), groups$code, 113635))
  })
})
pop <- bind_rows(pop)
value <- trimws(pop$Valor)
value[value == "-"] <- "0"
pop <- pop |>
  transmute(code7 = as.character(.data[[col(pop, "municipio.*codigo")]]),
            sex = ifelse(grepl("^homens", tolower(
              stri_trans_general(.data[[col(pop, "^sexo$")]],
                                 "Latin-ASCII"))), "M", "F"),
            age_group = clean_age(.data[[col(pop, "^idade$")]]),
            pop = as.numeric(value))
stopifnot(!anyNA(pop$pop), n_distinct(pop$code7) == 5570,
          nrow(pop) == 5570 * 42, sum(pop$pop) > 190e6)
pop$code6 <- substr(pop$code7, 1, 6)
pop$age_lo <- as.integer(sub("[-+].*", "", pop$age_group))

#------ Indirect standardisation by cause -----------------------------------
# National rates by sex and age group, from deaths that match a municipality
# and have a known sex and age.
sim$matched <- sim$code6 %in% pop$code6
sim$valid <- sim$matched & !is.na(sim$sex) & !is.na(sim$age)
sim$age_group <- ifelse(is.na(sim$age), NA, to_group(sim$age))

hits <- function(k) {
  ok <- grepl(spec$pattern[k], sim$cause)
  ok & (is.na(sim$age) | (sim$age >= spec$min[k] & sim$age <= spec$max[k]))
}
audit <- do.call(rbind, lapply(seq_len(nrow(spec)), function(k) {
  h <- hits(k)
  data.frame(outcome = spec$key[k], records = sum(h),
             outside_municipalities = sum(h & !sim$matched),
             missing_sex = sum(h & sim$matched & is.na(sim$sex)),
             missing_age_only = sum(h & sim$matched & !is.na(sim$sex) &
                                      is.na(sim$age)),
             kept = sum(h & sim$valid))
}))
write.csv(audit, here("data", "Mortality_audit.csv"), row.names = FALSE)
print(audit)

standardise <- function(k) {
  d <- sim[sim$valid & hits(k), ] |>
    count(code6, sex, age_group, name = "deaths")
  strata <- pop |>
    filter(age_lo >= spec$min[k], age_lo + 4 <= spec$max[k]) |>
    left_join(d, by = c("code6", "sex", "age_group")) |>
    mutate(deaths = coalesce(deaths, 0L))
  rates <- strata |>
    group_by(sex, age_group) |>
    summarise(rate = sum(deaths) / sum(pop), .groups = "drop")
  strata |>
    left_join(rates, by = c("sex", "age_group")) |>
    group_by(code_muni = code7) |>
    summarise(deaths = sum(deaths), expected = sum(pop * rate),
              pop = sum(pop), .groups = "drop") |>
    rename_with(~ paste0(.x, "_", spec$key[k]), c(deaths, expected, pop))
}
out <- Reduce(function(a, b) left_join(a, b, by = "code_muni"),
              lapply(seq_len(nrow(spec)), standardise))

#------ Quality of the registered causes ------------------------------------
by_muni <- function(keep) {
  tab <- table(factor(sim$code6[keep], levels = unique(pop$code6)))
  as.numeric(tab[match(substr(out$code_muni, 1, 6), names(tab))])
}
m <- sim$matched
deaths_all_raw <- by_muni(m)
cvd_all <- by_muni(m & grepl("^I[0-9]{2}", sim$cause))
out$deaths_all_raw <- deaths_all_raw
out$q_ill <- by_muni(m & grepl("^R[0-9]{2}", sim$cause)) / deaths_all_raw
out$q_nonspec <- by_muni(m & grepl(nonspecific_pattern, sim$cause)) / cvd_all

#------ Complete municipal table: covariates and outcomes -----------------
env <- new.env()
load(here("data", "Municipios_Brasil.RData"), envir = env)
base <- left_join(sf::st_drop_geometry(env$muni_sf), out, by = "code_muni")

# The old table counts deaths by residence from another extraction.
gap <- sum(base$deaths_all_raw != base$obitos_p_residenc, na.rm = TRUE)
msg("Municipalities where deaths differ from obitos_p_residenc: ", gap)
if (gap > 0) warning("SIM totals differ from obitos_p_residenc in ", gap,
                     " municipalities; check data/Mortality_audit.csv")

saveRDS(as.data.frame(out), here("data", "Mortality_2022.rds"))
write.csv(base, here("data", "Base_municipal_2022.csv"), row.names = FALSE,
          fileEncoding = "UTF-8")
msg("Saved data/Mortality_2022.rds and data/Base_municipal_2022.csv (",
    nrow(base), " municipalities)")
