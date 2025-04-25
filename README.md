
<!-- README.md is generated from README.Rmd. Please edit that file -->

# READ ME

Source code and essential data to support Julian et al. (*Submitted*)
Planning for the future, algae bloom dynamics in water management and
ecosystem restoration efforts.

Corresponding Author: Paul Julian (<pjulian@evergladesfoundation.org>)

## File Structure

`\LOK_AlgaeEval_manu`

### `\Data\`

- `LOK.gpkg`
  - shapefiles created for project that are not readily available via
    `library(EVERSpatDat)` or other public sources
- `LOK_bath.tif`
  - bathymetry geotiff in feet, NGVD29 from USACE 50ft composite
    bathymetry data. (originally retrieve from SFWMD GIS)

Some files are too large to store on Github, however `dat_proc.R` has
the code necessary to download publicly available files or retrieve data
from publicly available sources. The following files were saved to a
local drive:

- `LOK_Q.csv`
  - Inflow and outflow daily discharge data
- `LOK_STG.csv`
  - Daily stage elevation data
- `LOK_WQ.csv`
  - Water quality data from SFWMD DBHYDRO for locations within Lake
    Okeechobee
- `LOK_WQ_str.csv`
  - Water quality data from SFWMD DBHYDRO for structure discharge
    locations specific to Lake Okeechobee
- `LOK_sonde.csv`
  - Daily average sonde data including phycocyanin sensor data within
    Lake Okeechobee
- `RSMBN_LOK.csv`
  - RSMBN output data of modeled Lake Okeechobee stage elevation across
    LOSOM and LOCAR alternatives
- `RSMBN_LOK_Q.csv`
  - RSMBN output of modeled structure discharge into and out of Lake
    Okeechobee across LOSOM and LOCAR alternatives
- `LOK_AlgaeEval_analysis.RData`
  - Complete R session results of data analysis created by running
    `\src\dat_analysis.R`

### `\src\`

- `func.R`
  - file with custom functions used in the analysis of data
- `dat_proc.R`
  - initial data retrieval processes, where possible, files were saved
    as csv or other file format
- `dat_analyze.R`
  - code used to read in data, format and perform analyses presented in
    manuscript. This does not include all data visualization or final
    visualization included in the final publications

### `\Plots\`

- Final plots used in manuscript

### `\_documents\`

- place holder *submitted*
- place holder *final accepted*
