# Setup Log

## Why renv, not a global R library or conda

Previously, R packages were installed globally (shared across all
projects). This project uses `renv` instead, which creates an
isolated, project-local package library recorded in `renv.lock`.
This means:
- Packages installed for other past projects (e.g. CRISPRa scRNA-seq)
  are NOT visible here, even if already installed elsewhere on this
  Mac — this project has its own copy.
- Anyone (including future me) can run `renv::restore()` in this
  project folder and get the exact same package versions back.

Conda was deliberately not used for the R side of this project —
mixing conda-managed R packages with renv's project library is a
common source of dependency conflicts, especially for compiled
Bioconductor packages.

## System setup (Terminal, one-time, not project-specific)

1. Xcode command line tools — already installed (`xcode-select -p`
   returned a path).
2. Homebrew — already installed.
3. System libraries needed to compile some Bioconductor/CRAN packages
   on Apple Silicon:
   `brew install hdf5 gcc libxml2 pkg-config`
4. Attempted to install R via `rig` (`brew install r-lib/tap/rig`) for
   clean version management — failed due to a git/GitHub auth issue
   when Homebrew tried to clone the rig tap (root cause not fully
   resolved; a conda env may have been interfering with PATH, but
   `which git` still resolved to system git, so this wasn't fully
   confirmed).
5. Fallback: installed R directly via
   `brew install --cask r` — pulls the arm64 build straight from CRAN,
   no GitHub clone involved. Confirmed via `R --version`:
   R 4.6.1, `aarch64-apple-darwin` (native arm64, not Rosetta).

## Project environment (R Console, inside this project)

1. `install.packages("renv")` then `renv::init(bare = TRUE)` — creates
   an isolated project library.
2. CRAN packages installed: tidyverse, Seurat, SeuratObject, Matrix,
   patchwork, harmony, future, BiocManager, qs2
   - Note: originally tried `qs` (fast serialization), but it was
     removed from CRAN (2026-01-17, unresolved maintenance issues).
     Switched to its actively maintained successor, `qs2`
     (`qs_save()`/`qs_read()` instead of `qsave()`/`qread()`).
3. Bioconductor packages installed via
   `BiocManager::install(..., update = FALSE, ask = FALSE)`:
   GEOquery, SingleCellExperiment, scran, scater, glmGamPoi
4. GitHub-only package: `remotes::install_github("immunogenomics/presto")`
   — fast marker/DE tests for Seurat, not on CRAN.
5. `renv::snapshot()` — recorded exact versions to `renv.lock`.

## To reproduce this environment on another machine

```r
install.packages("renv")
renv::restore()
```