# Project Context & AI Rules

## Project Scope

- **Objective:** Simulate animal datasets with known "true" biodiversity. Test extrapolative model performance using sampled data to predict total gamma diversity and extinction from Amazonian deforestation.
- **Data Types:** Spatial occurrence data, diversity matrices, area-loss simulations.

## Technologies & Frameworks

- **Language:** R (Base R syntax preferred; avoid heavy `dplyr`/tidy verbs).
- **Key Packages:** `data.table` (for performance), `sf` (for spatial data), `ggplot2`, `tidyr`.
- **Custom Tools:** Uses external sampling functions from `./5-SensitivityTesting/scripts/utils/samplingFunctions.R`.

## Coding & Architecture Rules

- **Style:** Follow Tidyverse layout rules (snake_case, spacing). Do not use `dplyr` chaining/pipes (`%>%` or `|>`).
- **Modularity:** Write small, functional code. For ease of line-by-line debugging, internal function variable names must mirror external global variable names.
- **Documentation:** Use `roxygen2` for functions. Include thorough inline comments explaining algorithmic steps.
- **Performance:** Use vectorized operations or `data.table`. Suggest parallelization for multi-dataset operations.
- **Defensive Design:** Include `verbose = TRUE` options and frequent status `message()` printing in functions.

## Scientific Integrity

- **Reproducibility:** Maintain consistent spatial/tabular formats. Functions must return/log their input parameters alongside results.
- **Output:** Publication-ready `ggplot2` figures with consistent, clean themes.

## Interaction Rules

- **Tone:** Concise, objective, zero conversational filler.
- **Diagnostics:** For errors, prioritize structural checks (e.g., `str()`, `dim()`, `head()`, coordinate reference system checks, variable names) before modifying logic.
- **Context:** Assume strong Base R proficiency and low-moderate domain knowledge in ecology and statistics.
