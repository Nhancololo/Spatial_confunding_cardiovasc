# Application to municipal mortality, Brazil 2022 (run Data_mortality.R first)
# Environment variables: outcome, exposure, max_q, shapes (see README).

outcome <- Sys.getenv("outcome", "cvd")
nouns <- c(cvd = "cardiovascular death",
           cvd_premature = "premature cardiovascular death",
           ihd = "ischaemic heart disease death",
           stroke = "cerebrovascular death",
           hypertensive = "hypertensive disease death",
           diabetes = "diabetes death", neoplasm = "cancer death",
           homicide = "homicide")
noun <- nouns[[outcome]]
options(sp.cause = list(key = outcome, noun = noun,
                        Noun = paste0(toupper(substring(noun, 1, 1)),
                                      substring(noun, 2), "s"),
                        dir = paste0("brazil_", outcome),
                        file = "Mortality_2022.rds"))
Sys.setenv(shapes = Sys.getenv("shapes", "TRUE"),
           require_inla = Sys.getenv("require_inla", "TRUE"))
source(here::here("script", "Real_data_application.R"))
