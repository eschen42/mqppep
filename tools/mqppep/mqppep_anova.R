#!/usr/bin/env Rscript
# libraries
library(optparse)
library(data.table)
library(stringr)
# bioconductor-preprocesscore
#  - libopenblas
#  - r-data.table
#  - r-rmarkdown
#  - r-ggplot2
#  - texlive-core

# ref for parameterizing Rmd document: https://stackoverflow.com/a/37940285

# parse options
option_list <- list(
  make_option(
    c("-i", "--inputFile"),
    action = "store",
    default = NA,
    type = "character",
    help = "Phosphopeptide Intensities sparse input file path"
  ),
  make_option(
    c("-a", "--alphaFile"),
    action = "store",
    default = NA,
    type = "character",
    help = paste0("List of alpha cutoff values for significance testing;",
             " path to text file having one column and no header")
  ),
  make_option(
    c("-f", "--firstDataColumn"),
    action = "store",
    default = "10",
    type = "character",
    help = "First column of intensity values"
  ),
  make_option(
    c("-m", "--imputationMethod"),
    action = "store",
    default = "random",
    type = "character",
    help = paste0("Method for missing-value imputation,",
             " one of c('group-median','median','mean','random')")
  ),
  make_option(
    c("-p", "--meanPercentile"),
    action = "store",
    default = 3,
    type = "integer",
    help = paste0("Mean percentile for randomly generated imputed values;",
              ", range [1,99]")
  ),
  make_option(
    c("-d", "--sdPercentile"),
    action = "store",
    default = 3,
    type = "double",
    help = paste0("Adjustment value for standard deviation of",
              " randomly generated imputed values; real")
  ),
  make_option(
    c("-s", "--regexSampleNames"),
    action = "store",
    default = "\\.(\\d+)[A-Z]$",
    type = "character",
    help = "Regular expression extracting sample-names"
  ),
  make_option(
    c("-g", "--regexSampleGrouping"),
    action = "store",
    default = "(\\d+)",
    type = "character",
    help = paste0("Regular expression extracting sample-group",
             " from an extracted sample-name")
  ),
  make_option(
    c("-o", "--imputedDataFile"),
    action = "store",
    default = "output_imputed.tsv",
    type = "character",
    help = "Imputed Phosphopeptide Intensities output file path"
  ),
  make_option(
    c("-n", "--imputedQNLTDataFile"),
    action = "store",
    default = "output_imp_qn_lt.tsv",
    type = "character",
    help =
      paste(
        "Imputed, Quantile-Normalized Log-Transformed Phosphopeptide",
        "Intensities output file path"
        )
  ),
  make_option(
    c("-r", "--reportFile"),
    action = "store",
    default = "QuantDataProcessingScript.html",
    type = "character",
    help = "HTML report file path"
  )
)
args <- parse_args(OptionParser(option_list = option_list))
print("args is:")
cat(str(args))

# Check parameter values

if (! file.exists(args$inputFile)) {
  stop((paste("Input file", args$inputFile, "does not exist")))
}
input_file <- args$inputFile
alpha_file <- args$alphaFile
first_data_column <- args$firstDataColumn
imputation_method <- args$imputationMethod
print(
  grepl(
    pattern = imputation_method,
    x = c("group-median", "median", "mean", "random")
    )
  )

if (
  sum(
    grepl(
      pattern = imputation_method,
      x = c("group-median", "median", "mean", "random")
      )
    ) < 1
  ) {
    print(sprintf("bad imputationMethod argument: %s", imputation_method))
    return(-1)
    }

mean_percentile <- args$meanPercentile
print("mean_percentile is:")
cat(str(mean_percentile))

sd_percentile <- args$sdPercentile
print("sd_percentile is:")
cat(str(mean_percentile))


regex_sample_names    <- gsub("^[ \t\n]*", "",
                         readChar(args$regexSampleNames,  1000)
                       )
regex_sample_names    <- gsub("[ \t\n]*$", "",
                         regex_sample_names
                       )
cat(regex_sample_names)
cat("\n")

regex_sample_grouping <- gsub("^[ \t\n]*", "",
                           readChar(args$regexSampleGrouping, 1000)
                         )
regex_sample_grouping <- gsub("[ \t\n]*$", "",
                           regex_sample_grouping
                         )
cat(regex_sample_grouping)
cat("\n")

imputed_data_file_name <- args$imputedDataFile
imp_qn_lt_data_filenm <-  args$imputedQNLTDataFile
report_file_name <- args$reportFile

print("regex_sample_names is:")
cat(str(regex_sample_names))

print("regex_sample_grouping is:")
cat(str(regex_sample_grouping))

# from: https://github.com/molgenis/molgenis-pipelines/wiki/
#   How-to-source-another_file.R-from-within-your-R-script
# Function location_of_this_script returns the location of this .R script
#   (may be needed to source other files in same dir)
location_of_this_script <- function() {
    this_file <- NULL
    # This file may be 'sourced'
    for (i in - (1:sys.nframe())) {
        if (identical(sys.function(i), base::source)) {
            this_file <- (normalizePath(sys.frame(i)$ofile))
        }
    }

    if (!is.null(this_file)) return(dirname(this_file))

    # But it may also be called from the command line
    cmd_args <- commandArgs(trailingOnly = FALSE)
    cmd_args_trailing <- commandArgs(trailingOnly = TRUE)
    cmd_args <- cmd_args[
      seq.int(
        from = 1,
        length.out = length(cmd_args) - length(cmd_args_trailing)
        )
      ]
    res <- gsub("^(?:--file=(.*)|.*)$", "\\1", cmd_args)

    # If multiple --file arguments are given, R uses the last one
    res <- tail(res[res != ""], 1)
    if (0 < length(res)) return(dirname(res))

    # Both are not the case. Maybe we are in an R GUI?
    return(NULL)
}

script_dir <-  location_of_this_script()

rmarkdown_params <- list(
    inputFile = input_file
  , alphaFile = alpha_file
  , firstDataColumn = first_data_column
  , imputationMethod = imputation_method
  , meanPercentile = mean_percentile
  , sdPercentile = sd_percentile
  , regexSampleNames = regex_sample_names
  , regexSampleGrouping = regex_sample_grouping
  , imputedDataFilename = imputed_data_file_name
  , imputedQNLTDataFile = imp_qn_lt_data_filenm
  )

print("rmarkdown_params")
str(rmarkdown_params)

# freeze the random number generator so the same results will be produced
#  from run to run
set.seed(28571)

# BUG (or "opportunity")
# To render as PDF for the time being requires installing the conda
# package `r-texlive` until this issue in `texlive-core` is resolved:
#   https://github.com/conda-forge/texlive-core-feedstock/issues/19
# This workaround is detailed in the fourth comment of:
#   https://github.com/conda-forge/texlive-core-feedstock/issues/61

library(tinytex)
tinytex::install_tinytex()
rmarkdown::render(
  input = paste(script_dir, "mqppep_anova_script.Rmd", sep = "/")
, output_format = rmarkdown::pdf_document(toc = TRUE)
, output_file = report_file_name
, params = rmarkdown_params
)
