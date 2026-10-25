# Simulation study and application: spatial confounding

Code for Section 8 (simulation) and the application of the paper
"When Should Spatial Confounding Be Adjusted? A Minimax Sensitivity
Analysis for Count Regression".

## Layout

    JoseBerilo.Rproj                  project root (read by `here`)
    data/Municipios_Brasil.RData      Brazilian municipalities
    script/functions.R                all functions: settings, graphs,
                                      geometry, fits, study, report,
                                      application, checks
    script/simulations.R              simulation study
    script/Real_data_application.R    Brazilian homicide data
    output/                           created by the scripts

## Use

Open `JoseBerilo.Rproj` or run the commands from any folder inside the
project. 

    Rscript script/simulations.R full                     # the design of Section 8
    Rscript script/simulations.R full coverage,leakage    # selected experiments
    Rscript script/simulations.R report full              # redraw figures, tables
    Rscript script/Real_data_application.R Municipios_Brasil.RData \
        output/brazil complete FALSE TRUE

Methods in the result files: `NS` (non-spatial), `MM(c=...)` (minimax at a
fixed radius), `MM(c0)` (minimax at the true radius), `MM(sigma)` (radius
set to the fitted prior scale) and `MM(sigma)-no-bias` (the same estimate
with an interval that ignores the bias term).

Spatial+-ICAR in the simulation uses a low-rank spatial term (the
`splus_dim` smoothest eigenvectors); with the full basis the penalized
term absorbs part of the residual of x and the coefficient is biased even
without confounding. The INLA version in the application is full rank.

`count_linearization` runs the count analysis with the geometry at the
mean without the confounder ("baseline"); the suffixes `_true` and
`_fitted` use the true and the fitted mean. Poisson jobs use the
dispersion phi = 1 (strict Poisson); `_quasi` repeats the baseline run
with the Pearson estimate of phi (quasi-Poisson). The Brazilian analysis writes `sensitivity_cost.csv`
with the time and peak memory of each covariate.


The Scottish graph is read from INLA's `demodata/scotland.graph` unless
`cfg$scotland_graph` gives another file; without it, that graph is skipped.

`redraw` is `none` (default), `maps` (maps and effects figure) or
`sensitivity` (sensitivity tables, CSV files and figures); both read
`analysis_workspace.rds` of a finished run and do not refit the models.

`redraw = sensitivity` also writes `interval_length.csv`,
`tables/table_length.tex` and `figures/interval_length.pdf`: the radius of
the shortest bias-aware interval, its efficiency and that of the minimax
interval against the lower bound on the expected length of any confidence
set, at multiples of the reference radius.
