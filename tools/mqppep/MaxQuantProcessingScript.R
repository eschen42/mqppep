#!/usr/bin/env Rscript

# This is the implementation for the 
#   "MaxQuant Phosphopeptide Localization Probability Cutoff"
#   Galaxy tool (mqppep_lclztn_filter)
# It is adapted from the MaxQuant Processing Script written by Larry Cheng.

# libraries
library(optparse)
library(data.table)
library(stringr)
library(ggplot2)
#library(PTXQC)
#require(PTXQC)
#require(methods)

# title: "MaxQuant Processing Script"
# author: "Larry Cheng"
# date: "February 19, 2018"
#
# # MaxQuant Processing Script
# Takes MaxQuant Phospho (STY)sites.txt file as input and performs the following (in order):
# 1) Runs the Proteomics Quality Control software
# 2) Remove contaminant and reverse sequence rows
# 3) Filters rows based on localization probability
# 4) Extract the quantitative data
# 5) Sequences phosphopeptides
# 6) Merges multiply phosphorylated peptides
# 7) Filters out phosphopeptides based on enrichment
# The output file contains the phosphopeptide (first column) and the quantitative values for each sample
#
# ## Revision History
# Rev. 2022-02-10 :wrap for inclusion in Galaxy
# Rev. 2018-02-19 :break up analysis script into "MaxQuant Processing Script" and "Phosphopeptide Processing Script"
# Rev. 2017-12-12 :added PTXQC
#                  added additional plots and table outputs for quality control
#                  allowed for more than 2 samples to be grouped together (up to 26 (eg, 1A, 1B, 1C, etc))regexSampleNames <-
#                  "\\.(\\d+)[A-Z]$"
#                  converted from .r to .rmd file to knit report for quality control
# Rev. 2016-09-11 :automated the FDR cutoffs; removed the option to data impute multiple times
# Rev. 2016-09-09 :added filter to eliminate contaminant and reverse sequence rows
# Rev. 2016-09-01 :moved the collapse step from after ANOVA filter to prior to preANOVA file output
# Rev. 2016-08-22 :changed regexpression to regexSampleNames <- "\\.(\\d+)[AB]$" so that it looks at the end of string
# Rev. 2016-08-05 :Removed vestigial line (ppeptides <- ....)
# Rev. 2016-07-03 :Removed row names from the write.table() output for ANOVA and PreANOVA
# Rev. 2016-06-25 :Set default Localization Probability cutoff to 0.75
# Rev. 2016-06-23 :fixed a bug in filtering for pY enrichment by resetting the row numbers afterwards
# Rev. 2016-06-21 :test18 + standardized the regexpression in protocol


### FUNCTION DECLARATIONS begin ----------------------------------------------

# Read first line of file at filePath
# adapted from: https://stackoverflow.com/a/35761217/15509512
readFirstLine <- function(filepath) {
  con = file(filepath, "r")
  line = readLines(con, n = 1)
  close(con)
  return(line)
}

# Move columns to the end of dataframe
# - data: the dataframe
# - move: a vector of column names, each of which is an element of names(data)
movetolast <- function(data, move) {
  data[c(setdiff(names(data), move), move)]
}

# Generate phosphopeptide and build list when applied
phosphopeptide_func <- function(df) {

  #generate peptide sequence and list of phosphopositions
  phosphoprobsequence <- strsplit(as.character(df["Phospho (STY) Score diffs"]), "")[[1]]
  output <- vector()
  phosphopeptide <- ""
  counter <- 0 #keep track of position in peptide
  phosphopositions <- vector() #keep track of phosphorylation positions in peptide
  score_diff <- ""
  for (chara in phosphoprobsequence){
    #build peptide sequence
    if (!(chara == " " | chara == "(" | chara == ")" | chara =="." | chara =="-" | chara == "0" | chara == "1" | chara == "2" | chara == "3" | chara =="4" | chara == "5" | chara == "6" | chara == "7" | chara =="8" | chara =="9")) {
      phosphopeptide <- paste(phosphopeptide,chara,sep="")
      counter <- counter + 1
    }
    #generate score_diff
    if (chara == "-" | chara =="." | chara == "0" | chara == "1" | chara == "2" | chara == "3" | chara =="4" | chara == "5" | chara == "6" | chara == "7" | chara =="8" | chara =="9"){
      score_diff <- paste(score_diff,chara,sep="")
    }
    #evaluate score_diff
    if (chara == ")" ){
      score_diff <- as.numeric(score_diff)
      #only consider a phosphoresidue if score_diff > 0
      if (score_diff > 0) {
        phosphopositions <- append(phosphopositions, counter)
      }
      score_diff <- ""
    }
  }

  #generate phosphopeptide sequence (ie, peptide sequence with "p"'s)
  counter <- 1
  phosphoposition_correction1 <- -1 #used to correct phosphosposition as "p"'s are inserted into the phosphopeptide string
  phosphoposition_correction2 <- 0 #used to correct phosphosposition as "p"'s are inserted into the phosphopeptide string
  while (counter <= length(phosphopositions) ) {
    phosphopeptide <- paste(substr(phosphopeptide,0,phosphopositions[counter]+phosphoposition_correction1),"p",substr(phosphopeptide,phosphopositions[counter]+phosphoposition_correction2,nchar(phosphopeptide)),sep="")
    counter <- counter + 1
    phosphoposition_correction1 <- phosphoposition_correction1 + 1
    phosphoposition_correction2 <- phosphoposition_correction2 + 1
  }

  #building phosphopeptide list
  output <- append(output,phosphopeptide)
  return(output)
}

### FUNCTION DECLARATIONS end ------------------------------------------------


### EXTRACT ARGUMENTS begin --------------------------------------------------

# parse options
option_list <- list(
  make_option(
    c("-i", "--input"),
    action = "store",
    type = "character",
    help = "A MaxQuant Phospho (STY)Sites.txt"
  )
,  make_option(
    c("-o", "--output"),
    action = "store",
    type = "character",
    help = "path to output file"
  )
, make_option(
    c("-E", "--enrichGraph"),
    action = "store",
    type = "character",
    help = "path to enrichment graph PDF"
  )
, make_option(
    c("-F", "--enrichGraph_svg"),
    action = "store",
    type = "character",
    help = "path to enrichment graph SVG"
  )
, make_option(
    c("-L", "--locProbCutoffGraph"),
    action = "store",
    type = "character",
    help = "path to location-proability cutoff graph PDF"
  )
, make_option(
    c("-M", "--locProbCutoffGraph_svg"),
    action = "store",
    type = "character",
    help = "path to location-proability cutoff graph SVG"
  )
, make_option(
    c("-e", "--enriched"),
    action = "store",
    type = "character",
    help = "pY or pST enriched samples (ie, 'Y' or 'ST')"
  )
  # default = "^Number of Phospho [(]STY[)]$",
, make_option(
    c("-p", "--phosphoCol"),
    action = "store",
    type = "character",
    help = "PERL-compatible regular expression matching header of column having number of 'Phospho (STY)'"
  )
  # default = "^Intensity[^_]",
, make_option(
    c("-s", "--startCol"),
    action = "store",
    type = "character",
    help = "PERL-compatible regular expression matching column header having first sample intensity"
  )
  # default = 1,
, make_option(
    c("-I", "--intervalCol"),
    action = "store",
    type = "integer",
    help = "Column interval between the Intensities of samples (eg, 1 if subsequent column; 2 if every other column"
  )
  # default = 0.75,
, make_option(
    c("-l", "--localProbCutoff"),
    action = "store",
    type = "double",
    help = "Localization Probability Cutoff"
  )
  # default = "sum",
, make_option(
    c("-f", "--collapse_func"),
    action = "store",
    type = "character",
    help = "merge identical phosphopeptides by ('sum' or 'average') the intensities"
  )
  # default = "filteredData.txt",
, make_option(
    c("-r", "--filtered_data"),
    action = "store",
    type = "character",
    help = "filteredData.txt"
  )
  # default = "quantData.txt",
, make_option(
    c("-q", "--quant_data"),
    action = "store",
    type = "character",
    help = "quantData.txt"
  )
)
args <- parse_args(OptionParser(option_list=option_list))
# Check parameter values

### EXTRACT ARGUMENTS end ----------------------------------------------------


### EXTRACT PARAMETERS from arguments begin ----------------------------------

if (! file.exists(args$input)) {
  stop((paste("File", args$input, "does not exist")))
}

phosphoColPattern <- "^Number of Phospho [(][STY][STY]*[)]$"
startColPattern <- "^Intensity[^_]"
phosphoColPattern <- readFirstLine(args$phosphoCol)
startColPattern <- readFirstLine(args$startCol)

sink(getConnection(2))
#ACE print(paste("phosphoColPattern", phosphoColPattern))
#ACE print(paste("startColPattern", startColPattern))

inputFilename <- args$input
filteredFilename <- args$filtered_data
quantFilename <- args$quant_data
intervalCol <- as.integer(args$intervalCol)

firstLine <- readFirstLine(inputFilename)
columnHeaders <- unlist(strsplit(x=firstLine, split=c('\t'), fixed=TRUE))
sink(getConnection(2))
#ACE print("columnHeaders")
#ACE print(columnHeaders)
sink()


intensityHeaderCols <- grep(pattern=startColPattern, x=columnHeaders, perl=TRUE)
if ( length(intensityHeaderCols) == 0) {
    err_msg <- paste("Found no intensity columns matching pattern:", startColPattern)
    # Divert output to stderr
    sink(getConnection(2))
    print(err_msg)
    sink()
    stop(err_msg)
    }


phosphoCol <- grep(pattern=phosphoColPattern, x=columnHeaders, perl=TRUE)[1]
if (is.na(phosphoCol)) {
    err_msg <- paste("Found no 'number of phospho sites' columns matching pattern:", phosphoColPattern)
    # Divert output to stderr
    sink(getConnection(2))
    print(err_msg)
    sink()
    stop(err_msg)
    }


i_count <- 0
this_column <- 1
last_value <- intensityHeaderCols[1]
intensityCols <- c(last_value)

while ( length(intensityHeaderCols) >= intervalCol * i_count ) {
  i_count <- 1 + i_count
  this_column <- intervalCol + this_column
  if ( last_value + intervalCol != intensityHeaderCols[this_column] ) break
  last_value <- intensityHeaderCols[this_column]
  if (length(intensityHeaderCols) < intervalCol * i_count) break
  intensityCols <- c(intensityCols, intensityHeaderCols[this_column])
  }

startCol <- intensityCols[1]
numSamples <- i_count

outputfilename <- args$output
enrichGraphFilename <- args$enrichGraph
locProbCutoffGraphFilename <- args$locProbCutoffGraph
enrichGraphFilename_svg <- args$enrichGraph_svg
locProbCutoffGraphFilename_svg <- args$locProbCutoffGraph_svg

localProbCutoff <- args$localProbCutoff
enriched <- args$enriched
collapse_FUN <- args$collapse_func

### EXTRACT PARAMETERS from arguments end ------------------------------------


# Proteomics Quality Control for MaxQuant Results
#  (Bielow C et al. J Proteome Res. 2016 PMID: 26653327)
# is run by the Galaxy MaxQuant wrapper and need not be invoked here.


# Read data, filtering out contaminants, reverse sequences, and localization probability
# ---
fullData <- read.table(file = inputFilename, sep ="\t", header=T, quote="")

#Filter out contaminant rows and reverse rows
filteredData <- subset(fullData,!grepl("CON__", Proteins))
filteredData <- subset(filteredData,!grepl("_MYCOPLASMA", Proteins))
filteredData <- subset(filteredData,!grepl("CONTAMINANT_", Proteins))
filteredData <- subset(filteredData,!grepl("REV__", Protein)) #since REV__ rows are blank in the first column (Proteins)
write.table(filteredData, file = filteredFilename, sep = "\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
# ...


# Filter out data with localization probability below localProbCutoff
# ---
#Data filtered by localization probability
locProbFilteredData <- filteredData[filteredData$Localization.prob>=localProbCutoff,]
# ...


# Localization probability -- visualize locprob cutoff
# ---
locProbGraphData <- data.frame(
  group = c(paste(">",toString(localProbCutoff),sep=""), paste("<",toString(localProbCutoff),sep="")),
  value = c(nrow(locProbFilteredData)/nrow(filteredData)*100, (nrow(filteredData)-nrow(locProbFilteredData))/nrow(filteredData)*100)
)
gigi <-
  ggplot(locProbGraphData, aes(x = "", y = value, fill = group)) +
  geom_bar(width = 0.5, stat = "identity", color = "black") +
  labs(
    x = NULL
  , y = "percent"
  , title = "Phosphopeptides partitioned by localization-probability cutoff"
  ) +
  scale_fill_discrete(name = "phosphopeptide\nlocalization-\nprobability") +
  theme_minimal() +
  theme(
         legend.position = "right"
       , legend.title=element_text()
       , plot.title = element_text(hjust = 0.5)
       , plot.subtitle = element_text(hjust = 0.5)
       , plot.title.position = "plot"
       )
pdf(locProbCutoffGraphFilename)
print(gigi)
dev.off()
svg(locProbCutoffGraphFilename_svg)
print(gigi)
dev.off()
# ...


# Extract quantitative values from filtered data
# ---
quantData <- locProbFilteredData[,seq(from=startCol, by=intervalCol, length.out=numSamples)]
# ...


# Generate Phosphopeptide Sequence
#   for latest version of MaxQuant (Version 1.5.3.30)
# ---
dataTable <- data.frame(locProbFilteredData[,1:8],locProbFilteredData[,phosphoCol],locProbFilteredData[,phosphoCol+1],locProbFilteredData[,phosphoCol+2],locProbFilteredData[,phosphoCol+3],locProbFilteredData[,phosphoCol+4],locProbFilteredData[,phosphoCol+5],locProbFilteredData[,phosphoCol+6],locProbFilteredData[,phosphoCol+7],quantData)
colnames(dataTable) <- c("Proteins","Positions within proteins", "Leading proteins", "Protein", "Protein names", "Gene names", "Fasta headers", "Localization prob", "Number of Phospho (STY)", "Amino Acid", "Sequence window","Modification window", "Peptide window coverage", "Phospho (STY) Probabilities", "Phospho (STY) Score diffs", "Position in peptide", colnames(quantData))
# 'phosphopeptide_func' generates a phosphopeptide sequence for each row of data.
#   for the 'apply' function: MARGIN 1 == rows, 2 == columns, c(1,2) = both
dataTable$Phosphopeptide <- apply(X=dataTable, MARGIN=1, FUN=phosphopeptide_func)
# Move the quant data columns to the right end of the data.frame
dataTable <- movetolast(dataTable,c(colnames(quantData)))
# ...


# Write quantitative values for debugging purposes
# ---
quantWrite <- cbind( dataTable[,"Sequence window"], quantData ) 
colnames(quantWrite)[1] <- "Sequence.Window"
write.table(quantWrite, file = quantFilename, sep = "\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
# ...


# Make new data frame containing only Phosphopeptides to be mapped to quant data (merge_df)
# ---
dataTable <- setDT(dataTable, keep.rownames=TRUE) #row name will be used to map
merge_df <- data.frame(as.integer(dataTable$rn), dataTable$Phosphopeptide) #row index to merge data frames
colnames(merge_df) <- c("rn", "Phosphopeptide")
# ...


# Add Phosphopeptide column to quant columns for quality control checking
# ---
quantData_qc <- as.data.frame(quantData)
setDT(quantData_qc, keep.rownames=TRUE) #will use to match rowname to data
quantData_qc$rn <- as.integer(quantData_qc$rn)
quantData_qc <- merge(merge_df,quantData_qc, by="rn")
quantData_qc$rn <- NULL #remove rn column
# ...


# Collapse multiphosphorylated peptides
# ---
quantData_qc_collapsed <- data.table(quantData_qc, key = "Phosphopeptide")
quantData_qc_collapsed <- aggregate(. ~ Phosphopeptide,quantData_qc, FUN= collapse_FUN)
# ...


# Compute (as string) % of phosphopeptides that are multiphosphorylated (for use in next step)
# ---
pct_multiphos <- (nrow(quantData_qc) - nrow(quantData_qc_collapsed)) / (2 * nrow(quantData_qc))
pct_multiphos <- sprintf("%0.1f%s", 100 * pct_multiphos, "%")
# ...


# Compute and visualize breakdown of pY, pS, and pT before enrichment filter
# ---
pY_data <- quantData_qc_collapsed[str_detect(quantData_qc_collapsed$Phosphopeptide, "pY"),]
pS_data <- quantData_qc_collapsed[str_detect(quantData_qc_collapsed$Phosphopeptide, "pS"),]
pT_data <- quantData_qc_collapsed[str_detect(quantData_qc_collapsed$Phosphopeptide, "pT"),]

pY_num <- nrow(pY_data)
pS_num <- nrow(pS_data)
pT_num <- nrow(pT_data)

# Visualize enrichment
enrichGraphData <- data.frame(
  group = c("pY", "pS", "pT"),
  value = c(pY_num, pS_num, pT_num)
)

enrichGraphData <- enrichGraphData[enrichGraphData$value > 0,]

# Plot pie chart with legend
# start: https://stackoverflow.com/a/62522478/15509512
# refine: https://www.statology.org/ggplot-pie-chart/
# colors: https://colorbrewer2.org/#type=diverging&scheme=BrBG&n=8
slices <- enrichGraphData$value
phosphoresidue <- enrichGraphData$group
pct    <- round(100 * slices / sum(slices))
lbls   <- paste(enrichGraphData$group,"\n",pct, "%\n(", slices, ")", sep="")
slc_ctr <- c()
run_tot <- 0
for (p in pct) {
  slc_ctr <- c(slc_ctr, run_tot + p/2.0)
  run_tot <- run_tot + p
}
lbl_y  <- 100 - slc_ctr
df     <- data.frame(slices, pct, lbls, phosphoresidue = factor(phosphoresidue, levels = phosphoresidue))
gigi <- ggplot(
  df
, aes(x = 1, y = pct, fill = phosphoresidue)) +
  geom_col(position = "stack", orientation = "x") +
  geom_text(aes(x = 1, y = lbl_y, label = lbls), col = "black") +
  coord_polar(theta = "y", direction = -1) +
  labs(
    x = NULL
  , y = NULL
  , title = "Percentages (and counts) of phosphosites, by type of residue"
  , caption = sprintf("Roughly %s of peptides have multiple phosphosites.", pct_multiphos)
  ) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_classic() +
  theme( legend.position="right"
       , axis.line = element_blank()
       , axis.text = element_blank()
       , axis.ticks = element_blank()
       , plot.title = element_text(hjust = 0.5)
       , plot.subtitle = element_text(hjust = 0.5)
       , plot.caption = element_text(hjust = 0.5)
       , plot.title.position = "plot"
       ) +
  scale_fill_manual(breaks = phosphoresidue, values=c("#c7eae5", "#f6e8c3", "#dfc27d"))

pdf(enrichGraphFilename)
print(gigi)
dev.off()
svg(enrichGraphFilename_svg)
print(gigi)
dev.off()
# ...


# Filter phosphopeptides by enrichment
# --
if (enriched == "Y"){
  quantData_qc_enrichment <- quantData_qc_collapsed[str_detect(quantData_qc_collapsed$Phosphopeptide, "pY"),]
} else if ( enriched == "ST" ) {
  quantData_qc_enrichment <- quantData_qc_collapsed[str_detect(quantData_qc_collapsed$Phosphopeptide, "pS") | str_detect(quantData_qc_collapsed$Phosphopeptide, "pT"),]
} else {
  print("Error in enriched variable. Set to either 'Y' or 'ST'")
}
# ...


# Write phosphopeptides filtered by enrichment
# --
write.table(quantData_qc_enrichment, file=outputfilename, sep="\t", quote = FALSE, row.names = FALSE)
# ...
