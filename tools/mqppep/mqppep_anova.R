#!/usr/bin/env Rscript
# libraries
library(optparse)
library(data.table)
library(stringr)
#library(ggplot2)
#library(PTXQC)
#require(PTXQC)
#require(methods)
# bioconductor-preprocesscore
#  - libopenblas
#  - r-data.table
#  - r-rmarkdown
#  - r-ggplot2
#  - texlive-core

# ref for parameterizing Rmd document: https://stackoverflow.com/a/37940285

# parse options
option_list <- list(
  # <param name="inputFilename" type="data" format="tabular" label="Phosphopeptide Intensities" help="First column label 'Phosphopeptide'; sample-intensities must begin in column 10 and must have column labels to match argument regexSampleNames"/>
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
    help = "List of alpha cutoff values for significance testing; path to text file having one column and no header"
  ),
  make_option(
    c("-f", "--firstDataColumn"),
    action = "store",
    default = "10",
    type = "character",
    help = "First column of intensity values"
  ),
  make_option( # imputationMethod <- c("group-median","median","mean","random")[1]
    c("-m", "--imputationMethod"),
    action = "store",
    default = "group-median",
    type = "character",
    help = "Method for missing-value imputation, one of c('group-median','median','mean','random')"
  ),
  make_option(
    c("-p", "--meanPercentile"),
    action = "store",
    default = 3,
    type = "integer",
    help = "Mean percentile for randomly generated imputed values; range [1,99]"
  ),
  make_option(
    c("-d", "--sdPercentile"),
    action = "store",
    default = 3,
    type = "double",
    help = "Adjustment value for standard deviation of randomly generated imputed values; real"
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
    help = "Regular expression extracting sample-group from an extracted sample-name"
  ),
  # <data name="imputed_data_file" format="tabular" label="${input_file.name}.intensities_${imputation.imputation_method}-imputed_QN_LT" ></data>
  make_option(
    c("-o", "--imputedDataFile"),
    action = "store",
    default = "output_imputed.tsv",
    type = "character",
    help = "Imputed Phosphopeptide Intensities output file path"
  ),
  # <data name="report_file" format="html" label="report (download/unzip to view)" ></data>
  make_option(
    c("-r", "--reportFile"),
    action = "store",
    default = "QuantDataProcessingScript.html",
    type = "character",
    help = "HTML report file path"
  )
)
args <- parse_args(OptionParser(option_list=option_list))
# Check parameter values

if (! file.exists(args$inputFile)) {
  stop((paste("Input file", args$inputFile, "does not exist")))
}
inputFile <- args$inputFile
alphaFile <- args$alphaFile
firstDataColumn <- args$firstDataColumn
imputationMethod <- args$imputationMethod
meanPercentile <- args$meanPercentile
sdPercentile <- args$sdPercentile

regexSampleNames    <- gsub('^[ \t\n]*', ''  , readChar(args$regexSampleNames,  1000))
regexSampleNames    <- gsub('[ \t\n]*$', ''  ,               regexSampleNames        )
# regexSampleNames    <- gsub('\\\\'     , '@@',               regexSampleNames        )
# regexSampleNames    <- gsub('@@'       , '\\',               regexSampleNames        )
cat(regexSampleNames)
cat('\n')

regexSampleGrouping <- gsub('^[ \t\n]*', '', readChar(args$regexSampleGrouping, 1000))
regexSampleGrouping <- gsub('[ \t\n]*$', '',               regexSampleGrouping       )
# regexSampleGrouping <- gsub('\\\\'     , '@@',             regexSampleGrouping       )
cat(regexSampleGrouping)
cat('\n')

# regexSampleGrouping <- gsub('@@'       , '\\',             regexSampleGrouping       )
imputedDataFilename <- args$imputedDataFile
reportFileName <- args$reportFile

print("args is:")
cat(str(args))

print("regexSampleNames is:")
cat(str(regexSampleNames))

print("regexSampleGrouping is:")
cat(str(regexSampleGrouping))

# from: https://github.com/molgenis/molgenis-pipelines/wiki/How-to-source-another_file.R-from-within-your-R-script
LocationOfThisScript = function() # Function LocationOfThisScript returns the location of this .R script (may be needed to source other files in same dir)
{
    this.file = NULL
    # This file may be 'sourced'
    for (i in -(1:sys.nframe())) {
        if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
    }

    if (!is.null(this.file)) return(dirname(this.file))

    # But it may also be called from the command line
    cmd.args = commandArgs(trailingOnly = FALSE)
    cmd.args.trailing = commandArgs(trailingOnly = TRUE)
    cmd.args = cmd.args[seq.int(from=1, length.out=length(cmd.args) - length(cmd.args.trailing))]
    res = gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)

    # If multiple --file arguments are given, R uses the last one
    res = tail(res[res != ""], 1)
    if (0 < length(res)) return(dirname(res))

    # Both are not the case. Maybe we are in an R GUI?
    return(NULL)
}

script.dir <-  LocationOfThisScript()

rmarkdown_params <- list(
    inputFile = inputFile
  , alphaFile = alphaFile
  , firstDataColumn = firstDataColumn
  , imputationMethod = imputationMethod
  , meanPercentile = meanPercentile
  , sdPercentile = sdPercentile
  , regexSampleNames = regexSampleNames
  , regexSampleGrouping = regexSampleGrouping
  , imputedDataFilename = imputedDataFilename
  )

str(rmarkdown_params)

# BUG
# Must render as HTML for the time being until this issue is resolved:
#   https://github.com/conda-forge/texlive-core-feedstock/issues/19
# for reason:
#   "The following dependencies are not available in conda"
# reported here:
#   https://github.com/ami-iit/bipedal-locomotion-framework/pull/457/commits/e98ccef8c8cb63e207df36628192af6ce22feb13

# freeze the random number generator so the same results will be produced from run to run
set.seed(28571)

rmarkdown::render(
  input = paste(script.dir, "mqppep_anova_script.Rmd", sep="/")
, output_format = rmarkdown::html_document(pandoc_args = "--self-contained")
, output_file = reportFileName
, params = rmarkdown_params
)
