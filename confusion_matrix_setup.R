library(data.table)
library(dplyr)
library(ggplot2)

# Load inversion data
inversions <- fread("/home/celia/Downloads/Mosaicatcher/variants_callset_inversions.tsv", sep = '\t') %>%
  select(-misorient, -CALLERSET_LIST, -lowConfCount, -noreadsCount, -categ) %>%
  filter(len > 400)  # Filter inversions larger than 400bp

# Load mosaicatcher data
inv_mosaicatcher <- fread("/home/celia/Downloads/Mosaicatcher/lenient_filterFALSE.tsv", sep = '\t')
inv_mosaicatcher <- inv_mosaicatcher %>%
  mutate(cell = sub("x.*", "", sub("_.*", "", sub("HG02059102.*", "HG02059x102", cell)))) %>%
  filter(!sv_call_name %in% c("dup_h1", "complex", "dup_h2", "idup_h1", "idup_h2", "del_h1", 
                              "del_h2", "dup_hom", "del_hom"))

# Define genotype categories
no_inv <- c("0|0", "0|0_lowconf")
inv_hom <- c("1|1", "1|1_lowconf", "inv_hom", "0101")
inv_het <- c("0|1", "1|0", "0110", "1001", "0|1_lowconf", "1|0_lowconf", "inv_h1", "inv_h2")

# Recode genotype variables in mosaicatcher and inversions datasets
inv_mosaicatcher <- inv_mosaicatcher %>%
  mutate(sv_call_name = case_when(
    sv_call_name %in% inv_hom ~ "inv_hom",
    sv_call_name %in% inv_het ~ "inv_het",
    TRUE ~ NA_character_
  ))

inversions <- inversions %>%
  mutate(across(contains(c("NA1", "NA20", "NA24", "HG")), ~case_when(
    . %in% no_inv ~ "no_inv",
    . %in% inv_hom ~ "inv_hom",
    . %in% inv_het ~ "inv_het",
    TRUE ~ NA_character_
  )))

# Save BED files for further analysis
fwrite(inv_mosaicatcher, "/home/celia/Downloads/Mosaicatcher/mosaic_new.bed", sep = '\t', col.names = FALSE)
fwrite(inversions, "/home/celia/Downloads/Mosaicatcher/inv_paper_400kb.bed", sep = '\t', col.names = FALSE)

# Use BEDTools for intersecting BED files
mosaic_file <- "/home/celia/Downloads/Mosaicatcher/mosaic_new.bed"
inversions_all <- "/home/celia/Downloads/Mosaicatcher/inv_paper_400kb.bed"
output <- "/home/celia/Downloads/Mosaicatcher/intersect400.bed"

system(paste0("bedtools intersect -a ", inversions_all, " -b ", mosaic_file, " -f 0.7 -r -wao > ", output))

# Load and process the intersect output
output_intersect_4 <- fread(output, sep = '\t')

output_confusion <- output_intersect_4 %>%
  select(-V4,-V5,-V7,-V8,-V9,-V10,-V11,-V12,-V13,-V14,-V62,-V64,-V65,-V66,-V68,-V69,-V70,-V71,-V72,-V73,-V74)

colnames(output_confusion) <- c('chrom', 'start', 'end', 'stack', "NA19434", "HG00096", "HG00171", "HG00268",
                                "HG00512", "HG00513", "HG00514", "HG00731", "HG00732", "HG00733", "HG00864",
                                "HG01114", "HG01352", "HG01505", "HG01573", "HG01596", "HG02011", "HG02018",
                                "HG02059", "HG02106", "HG02492", "HG02587", "HG02818", "HG03009", "HG03065",
                                "HG03125", "HG03371", "HG03486", "HG03683", "HG03732", "HG04217", "NA12329",
                                "NA12878", "NA18534", "NA18939", "NA19036", "NA19238", "NA19239", "NA19240",
                                "NA19650", "NA19983", "NA20509", "NA20847", "NA24385", "chrom_mosaic", 
                                "start_mosaic", "end_mosaic", "sample", "SV")

output_confusion_unique <- distinct(output_confusion, stack, sample, SV, .keep_all = TRUE)

#fwrite(output_confusion_unique, "/home/celia/Downloads/Mosaicatcher/celia-testdir/celia-testdir/output_confusion_400.csv", sep = '\t', quote = FALSE)

