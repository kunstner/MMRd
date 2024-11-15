# Credits -----------------------------------------------------------------

# Determine the MMR status of a patient based on their NGS data
# and a mutational signature
# @author: Vanessa Klauenberg
# @created: 2024-07-16
# @lastchange: 2024-11-14

# Libraries ---------------------------------------------------------------

library(Biostrings)
library(GenomicRanges)
library(writexl)

# Data --------------------------------------------------------------------

# load the result of 02_bedHomopolymers.r
load(file = ".../matching_homopolymers.Rdata")
# load the total length of the sequenced chromosomes
load (file = ".../total_length.Rdata")

# read the .maf file
path <- "MAF"
variants <- read.table(path, header = TRUE, sep = "\t")

# check if the chromosomes start with 'chr', and add it if not present
variants <- variants %>%
    mutate(Chromosome = ifelse(grepl("^chr", Chromosome), Chromosome, paste0("chr", Chromosome)))

# Determine MMR-Status ----------------------------------------------------

# get the unique tumor sample barcodes
tumors <- unique(variants$Tumor_Sample_Barcode)

# prepare an empty data frame for the outcomes
out <- data.frame(Tumor = character(),
                  Variants = integer(),
                  Indels = integer(),
                  SBP_Indels = integer(),
                  Overlaps = integer(),
                  Score = numeric(),
                  Result = character(),
                  stringsAsFactors = FALSE)

# loop through each tumor
for (i in tumors) {
  tumor_variants <- variants[variants$ Tumor_Sample_Barcode == i, ]
  # check if there are more than 12 variants, if not skip to the next barcode
  if (nrow(tumor_variants) < 12) {
    print(paste("Too little variants, MMR-Status of ", i,
                " cannot be assessed."))
    next
  }
  # select only the indels from the variant file, with their relevant columns
  indels <- tumor_variants[tumor_variants$Variant_Type %in% c("INS", "DEL"), ]
  columns <- c("Chromosome", "Start_Position", "End_Position",
               "Tumor_Sample_Barcode")
  indels <- indels[, columns]

  # convert the indels to GRanges object
  indel_gr <- GRanges(seqnames = indels$Chromosome,
                      ranges = IRanges(start = indels$Start_Position,
                                       end = indels$End_Position),
                      tumor_barcode = indels$Tumor_Sample_Barcode)

  # filter the single base indels
  single_base_indels <- indel_gr[width(indel_gr) == 1]

  # find the indels located in homopolymers
  overlaps <- findOverlaps(single_base_indels, matching_homopolymers)
  # count the number of indels located in homopolymers
  counter <- length(unique(queryHits(overlaps)))

  # calculate the indel density per MB
  indel_density <- counter/total_length
  indel_density <- round(indel_density,2)

  # calculate the likelihood score
  if (indel_density > 0) {
    score <- 1.5/indel_density
  } else {
    stop ("The patient with barcode ", i,
          " has 0 indels in homopolymers and is therefore most likely MMR-proficient")
  }

  # print the MMR-status
  if (indel_density >= 1.5) {
    score <- 100 - score*100
	score <- round(score,2)
    result <- "D"
    print(paste("The patient with barcode ", i,
                " is MMR-deficient, with an indel density of ",
                indel_density, " and a likelihood score of ", score))
  }else {
	score <- round(score,2)
    result <- "P"
    print(paste("The patient with barcode ", i,
                " is MMR-proficient, with an indel density of ",
                indel_density, " and a likelihood score of ", score))
  }

  # save the outcome as a data frame
  out <- rbind(out, data.frame(Tumor = i,
                               Variants = nrow(tumor_variants),
                               Indels = nrow(indels),
                               SBP_Indels = length(single_base_indels),
                               Overlaps = counter,
                               Score = score,
                               Result = result,
                               stringsAsFactors = FALSE))

}

# Save results ------------------------------------------------------------

# save all the outcomes in an excel file
write_xlsx(list("MMR-Status" = out), path = "MMRstatus.xlsx")
