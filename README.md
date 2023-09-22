# DrugAnalyzer

## Development Environment
- install.packages('readxl')

## Example
```
setwd("/data/projects/DrugAnalyzer/src")
source("DataHandler.R")
source("DataHandler.R")
source("DataMatcher.R")
source("DataProcessor.R")
source("DrugDataHandler.R")

drug_printer_file <- "/data/projects/DrugAnalyzer/data/autophagy.xlsx"
output_dir <- "/data/projects/DrugAnalyzer/output"

drug_screen_table <- deconvolute_multidrug(drug_printer_file, output_dir)
```