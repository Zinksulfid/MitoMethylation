######################################
## Schober et al. The Mitomethylome ##
######################################

### The code contains the most important steps during analysis. Figure plots are ###
### derived from the data analysed herein. ###
### The code requires a folder "Data", in which all evidence files, MaxQuant PTM ###
### identification files and MethylQuant output files (named "MethQ_H_out_.*.csv") ###
### are in subfolder uniquely named. ###
### Additionally, a Drosophila mitocarta file as well as a mapping of mammalian homologs ###
### will be required. ###

library("tidyverse")
library("Biostrings")
library("Rcpi")
library("UniProt.ws")
library("msa")
library("circlize")
library("parallel")
library("doParallel")
library("ggpubr")
library("ggseqlogo")
library("org.Dm.eg.db")

#######################################################
## A: Unique sites list of all sites in MQ PTM files ##
#######################################################

PTM.site.files <- list.files(pattern = "ethyl.*\\(K.*txt", recursive = TRUE) # Identify all files in subolders that follow the MQ pattern of variable modification tables

d <- lapply(PTM.site.files, read_tsv) # Read in files
names(d) <- sub("Data/", "", PTM.site.files) # Give elements of list the names of the files and cut "Data/"

# Replace stuff in columns to facilitate combination
replace_tmp <- function(x){
  x %>% 
    select_all(funs(gsub(" ", "_", .))) %>%
    select_all(funs(gsub(".*Probabilities", "MethylProbabilities", .))) %>%
    filter(!grepl("REV_|CON_", Protein)) 
}

d <- lapply(d, replace_tmp)
#up <- UniProt.ws(taxId=7227)
load("UniProt.ws_7227")

# Reduce the number of columns
reduce_d <- function(x){
  x %>%
    dplyr::select(matches("Protein$|Localization_prob$|MethylProbabilities$|PEP$|Intensity$|Mass_error_.*$|Position_in_peptide$|Best_score_raw_file$"))# Select here which columns to include
}
d <- lapply(d, reduce_d)

# Replace old Uniprot IDs
replace_names <- function(x){
  x %>%
    mutate(Protein = replace(Protein, Protein == "E2RTZ8", "Q4V4U2")) %>%
    mutate(Protein = replace(Protein, Protein == "M9PBP0", "P84181")) %>%
    mutate(Protein = replace(Protein, Protein == "P84249", "C0HL66")) %>%
    mutate(Protein = replace(Protein, Protein == "Q9VGX3", "C0HK92")) %>%
    mutate(Protein = replace(Protein, Protein == "A0A0B4LFD2", "P29845")) %>%
    mutate(Protein = replace(Protein, Protein == "A0A0B4KHL7", "P91928")) %>%
    mutate(Protein = replace(Protein, Protein == "X2JB48", "Q26365")) %>%
    mutate(Protein = replace(Protein, Protein == "B7Z0E0", "Q7KUB0")) %>% 
    mutate(Protein = replace(Protein, Protein == "P07487", "M9PJN8")) %>% 
    mutate(Protein = replace(Protein, Protein == "P29844", "F3YDH0")) %>%
    mutate(Protein = replace(Protein, Protein == "O97125;P11146", "O97125")) %>%
    mutate(Protein = replace(Protein, Protein == "B5RIV0", "Q9V3S9")) %>%
    mutate(Protein = replace(Protein, Protein == "Q0E9E2;Q7KN97", "QQ0E9E2")) %>%
    mutate(Protein = replace(Protein, Protein == "A0A0B4K6V1;A0A0B4KHW3", "Q9VA73")) %>%
    mutate(Protein = replace(Protein, Protein == "A0A0B4K6V1", "Q9VA73")) %>%
    mutate(Protein = replace(Protein, Protein == "Q8IQW2", "Q9VWD1")) %>%
    mutate(Protein = replace(Protein, Protein == "A0A0B4LF45", "P82712")) %>%
    mutate(Protein = replace(Protein, Protein == "P82712;Q7KR10", "P82712")) %>%
    mutate(Protein = replace(Protein, Protein == "C7LA75", "P11147")) %>%
    mutate(Protein = replace(Protein, Protein == "P92192", "C0HL62")) -> x
}

d <- lapply(d, replace_names)

# Add sequence and ENSEMBL/FLYBASE ID. The uniprot.ws search cannot cope with large chunks, so the requests get limited to chunks of 200
seq_annotate <- function(x){
  
  # Function A to be able to rbind all lists together
  reduce_internal <- function(y){
    y %>%
      dplyr::select(matches("MethylProbablities$|Protein$"))# Select here which columns to include
  }
  
  x_reduced <- lapply(x, reduce_internal)
  
  # Take unique uniprot IDs
  dplyr::bind_rows(x_reduced, .id = "Origin") %>%
    distinct(Protein) -> x_bound
  
  # Split the request list to avoid crashing the pipe. Function B.
  x_fracs_list <- split(x_bound, (as.numeric(rownames(x_bound))-1) %/% 100)
  
  seq_annotate_internal <- function(w){ # Requiring some sort of nested function with reduced capabilities
    columns <- c("FLYBASE","SEQUENCE")
    kt <- "UNIPROTKB"
    mapping <- UniProt.ws::select(up, unique(sub(";.*","",w$Protein)), columns, kt)
    
    w %>%
      mutate(ENSEMBL = mapping[match(sub(";.*","",w$Protein), mapping$UNIPROTKB),]$FLYBASE) %>%
      mutate(Protein_Sequence = mapping[match(sub(";.*","",w$Protein), mapping$UNIPROTKB),]$SEQUENCE)
  }
  
  x_fracs_list <- lapply(x_fracs_list, seq_annotate_internal)
  
  # Unlist the fractions again and combine as initial x
  x_UNIPROT.library <- bind_rows(x_fracs_list)
  
  # Function C: Add annotation to list elements
  add_to_list <- function(z){
    dplyr::left_join(z, x_UNIPROT.library, by = "Protein") -> z
  }
  
  x <- lapply(x, add_to_list)
  return(x)
}

d <- seq_annotate(d)

# Annotate a unique peptide
pos_protein <- function(x){
  x$Position_in_peptide <- as.numeric(x$Position_in_peptide) # Some are encoded as character
  replace_vec <- setNames("", "\\(.{1,5}\\)")
  peptide <- str_replace_all(x$MethylProbabilities, replace_vec) # Remove probabilities and brackets from detected peptide
  names(peptide) <- x$Protein
  
  counter <- c()
  for(i in c(1:length(peptide))){
    nTermAA <- as.numeric(nchar(sub(paste(peptide[i], ".*", sep = ""), "", x$Protein_Sequence[i]))) # Loop through the FASTA sequence, cutting the peptide sequence and everything C-terminally of that of.
    methyl_position <- x$Position_in_peptide[i] + nTermAA # Calculate the length of the remaining N-terminal part and the absolute position of the PTM.
    counter <- c(counter, methyl_position)
  }
  
  x$Specific_position_in_protein <- counter
  return(x)
}

d <- lapply(d, pos_protein)

# Bind everything to one large data.frame and attach a column with the origin
d <- dplyr::bind_rows(d, .id = "Origin")
d %>%
  mutate(Origin = stringr::str_replace_all(d$Origin, "^.*_data/|Sites.txt|\\(.+\\)", "")) %>%
  mutate(ID = paste(Protein, Specific_position_in_protein, stringr::str_replace_all(Origin, "heavy |.*/|me|thyl ", ""), sep = ";")) -> d

# Reduce based on minimal PEP value
d %>% 
  group_by(ID, Origin) %>%
  top_n(n = -1, PEP) -> d

# Add MitoXplorer information
mitocarta <- read_tsv("Fly_Interactome_200515.txt")

mitocarta$SYMBOL <- mapIds(org.Dm.eg.db, 
                           keys=mitocarta$Gene.ID,
                           column="SYMBOL", 
                           keytype="ENSEMBL",
                           multiVals="first")

mitocarta$UNIPROT <- mapIds(org.Dm.eg.db, 
                            keys=mitocarta$Gene.ID,
                            column="UNIPROT", 
                            keytype="ENSEMBL",
                            multiVals="first")

read_tsv("Fly_Interactome_200515_orthologs.txt") %>%
  mutate(Gene.ID = `Fly GeneID`) %>%
  arrange(`Weighted Score`) %>% 
  distinct(FlyBaseID, .keep_all = TRUE) -> orthologs

mitocarta <- left_join(mitocarta, orthologs, by = "Gene.ID")

d %>%
  mutate(SYMBOL = mitocarta[match(ENSEMBL, mitocarta$Gene.ID),]$SYMBOL) %>%
  mutate(Mito.process = mitocarta[match(ENSEMBL, mitocarta$Gene.ID),]$Mitochondrial.Process) %>%
  mutate(Hs.homolog = mitocarta[match(ENSEMBL, mitocarta$Gene.ID),]$`Human Symbol`) -> d

# Overlap light and heavy identifications
## Subset
d %>%
  mutate(Label = str_replace(Origin, "/.*", "")) %>%
  filter(Label == "Light") -> d.Light

d %>%
  mutate(Label = str_replace(Origin, "/.*", "")) %>%
  filter(Label == "Heavy") -> d.Heavy

# Uniquely in light and with shift
d.Light %>%
  mutate(HL = ID %in% d.Heavy$ID) %>%
  mutate(HL = ifelse(HL == TRUE, "L+H", "L_unique")) -> d.Light

# Unique in heavy
d.Heavy %>%
  mutate(HL = ID %in% d.Light$ID) %>%
  mutate(HL = ifelse(HL == TRUE, "L+H", "H_unique")) %>%
  filter(HL == "H_unique") -> d.Heavy # To not get redundant hits

d <- rbind(d.Light, d.Heavy)

tbl_df(d) %>%
  dplyr::select(-c(Origin,Protein_Sequence, Label)) %>%
  filter(!is.na(Mito.process)) %>%
  arrange(Mito.process, Protein, Specific_position_in_protein) -> siteome

#######################################################
## B: Unique sites list of all sites in MQ PTM files ##
#######################################################

methQ.files <- list.files(pattern = "MethQ_._out", recursive = TRUE)

d <- lapply(methQ.files, read_csv) # Read in files

names(d) <- stringr::str_replace_all(methQ.files, "Data/|\\.csv", "") # Give elements of list the names of the files and cut "Data/"

# Replace stuff in columns to facilitate combination [[and throw MethylQuant Score 0]]
replace_tmp <- function(x){
  x %>%
    # filter(`MethylQuant Score` != 0) %>% ##### !! 
    select_all(funs(gsub(" ", "_", .))) %>%
    filter(!grepl("REV_|CON_", uniprot)) 
}

d <- lapply(d, replace_tmp)

# Reduce the number of columns
reduce_d <- function(x){
  x %>%
    dplyr::select(matches("Sequence$|Modifications$|uniprot$|mod_seq$|Int$|Mass_Difference$|Isotope_Distribution_Correlation$|\\#_Good_Elution_Profile_Correlations$|MethylQuant_Score$|MethylQuant_Confidence$|Data_File$"))# Select here which columns to include
}

d <- lapply(d, reduce_d)

# Replace old Uniprot IDs
replace_names <- function(x){
  x %>%
    mutate(uniprot = replace(uniprot, uniprot == "E2RTZ8", "Q4V4U2")) %>%
    mutate(uniprot = replace(uniprot, uniprot == "M9PBP0", "P84181")) %>%
    mutate(uniprot = replace(uniprot, uniprot == "P84249", "C0HL66")) %>%
    mutate(uniprot = replace(uniprot, uniprot == "Q9VGX3", "C0HK92")) %>%
    mutate(uniprot = replace(uniprot, uniprot == "A0A0B4LFD2", "P29845")) %>%
    mutate(uniprot = replace(uniprot, uniprot == "A0A0B4KHL7", "P91928")) %>%
    mutate(uniprot = replace(uniprot, uniprot == "X2JB48", "Q26365")) %>%
    mutate(uniprot = replace(uniprot, uniprot == "B7Z0E0", "Q7KUB0")) %>% 
    mutate(uniprot = replace(uniprot, uniprot == "P07487", "M9PJN8")) %>% 
    mutate(uniprot = replace(uniprot, uniprot == "P29844", "F3YDH0")) %>%
    mutate(uniprot = replace(uniprot, uniprot == "O97125;P11146", "O97125")) %>%
    mutate(uniprot = replace(uniprot, uniprot == "B5RIV0", "Q9V3S9")) %>%
    mutate(uniprot = replace(uniprot, uniprot == "Q0E9E2;Q7KN97", "Q0E9E2")) %>%
    mutate(uniprot = replace(uniprot, uniprot == "A0A0B4K6V1;A0A0B4KHW3", "Q9VA73")) %>%
    mutate(uniprot = replace(uniprot, uniprot == "A0A0B4K6V1", "Q9VA73")) %>%
    mutate(uniprot = replace(uniprot, uniprot == "Q8IQW2", "Q9VWD1")) %>%
    mutate(uniprot = replace(uniprot, uniprot == "A0A0B4LF45", "P82712")) %>%
    mutate(uniprot = replace(uniprot, uniprot == "P82712;Q7KR10", "P82712")) %>%
    mutate(uniprot = replace(uniprot, uniprot == "C7LA75", "P11147")) %>%
    mutate(uniprot = replace(uniprot, uniprot == "P92192", "C0HL62")) -> x
}

d <- lapply(d, replace_names)

# Add sequence and ENSEMBL/FLYBASE ID. The uniprot.ws search cannot cope with large chunks, so the requests get limited to chunks of 200
load("UniProt.ws_7227")

seq_annotate <- function(x){
  
  # Function A to be able to rbind all lists together
  reduce_internal <- function(y){
    y %>%
      dplyr::select(matches("Sequence$|uniprot$"))# Select here which columns to include
  }
  
  x_reduced <- lapply(x, reduce_internal)
  
  # Take unique uniprot IDs
  dplyr::bind_rows(x_reduced, .id = "Origin") %>%
    distinct(uniprot) -> x_bound
  
  # Split the request list to avoid crashing the pipe. Function B.
  x_fracs_list <- split(x_bound, (as.numeric(rownames(x_bound))-1) %/% 100)
  
  seq_annotate_internal <- function(w){ # Requiring some sort of nested function with reduced capabilities
    columns <- c("FLYBASE","SEQUENCE")
    kt <- "UNIPROTKB"
    mapping <- UniProt.ws::select(up, unique(sub(";.*","",w$uniprot)), columns, kt)
    
    w %>%
      mutate(ENSEMBL = mapping[match(sub(";.*","",w$uniprot), mapping$UNIPROTKB),]$FLYBASE) %>%
      mutate(Protein_Sequence = mapping[match(sub(";.*","",w$uniprot), mapping$UNIPROTKB),]$SEQUENCE)
  }
  
  x_fracs_list <- lapply(x_fracs_list, seq_annotate_internal)
  
  # Unlist the fractions again and combine as initial x
  x_UNIPROT.library <- bind_rows(x_fracs_list)
  
  # Function C: Add annotation to list elements
  add_to_list <- function(z){
    dplyr::left_join(z, x_UNIPROT.library, by = "uniprot") -> z
  }
  
  x <- lapply(x, add_to_list)
  return(x)
}

d <- seq_annotate(d)

remove_semicolon <- function(x){
  x$Modifications <- sub("^;", "", x$Modifications)
  return(x)
}

d <- lapply(d, remove_semicolon)

# Annotate a unique peptide
pos_protein <- function(x){
  
  # Multiplying peptides with more than one PTM annotated
  x %>%
    mutate(PTMs = stringr::str_count(x$Modifications, ";")) %>%
    mutate(PTMs_counter = 1) -> x
  
  rbind(x,
        x %>% 
          filter(PTMs == 1) %>% 
          mutate(PTMs_counter = 2),
        x %>% 
          filter(PTMs == 2) %>% 
          mutate(PTMs_counter = 2),
        x %>% 
          filter(PTMs == 2) %>% 
          mutate(PTMs_counter = 3)) -> x
  
  x %>%
    group_by(mod_seq) %>%
    mutate(Modifications = unlist(str_split(Modifications, ";"))[PTMs_counter]) -> x
  
  x$Position_in_peptide <- as.numeric(stringr::str_replace_all(x$Modifications, "^;|^M.+\\(Oxidation\\);|[RK]|\\(.*\\)|\\;.*", "")) # Some are encoded as character
  
  counter <- c()
  for(i in c(1:length(x$Sequence))){
    nTermAA <- as.numeric(nchar(sub(paste(x$Sequence[i], ".*", sep = ""), "", x$Protein_Sequence[i]))) # Loop through the FASTA sequence, cutting the peptide sequence and everything C-terminally of that of.
    methyl_position <- x$Position_in_peptide[i] + nTermAA # Calculate the length of the remaining N-terminal part and the absolute position of the PTM.
    counter <- c(counter, methyl_position)
  }
  
  x$Specific_position_in_protein <- counter
  return(x)
}

d <- lapply(d, pos_protein)

d <- dplyr::bind_rows(d, .id = "Origin")
d$Origin <- stringr::str_replace_all(d$Origin, "^.*_data/|_out", "") # Simplify the origin column

# Reduce based on maximal MethylQuant_Score value
d$ID <- paste(d$uniprot, d$Specific_position_in_protein, stringr::str_replace_all(d$Modifications, "^;|^M.+\\(Oxidation\\);|.*\\(|\\).*|-Heavy|me|thyl", ""), sep = ";")

tbl_df(d) %>%
  group_by(ID) %>%
  top_n(n = 1, MethylQuant_Score) %>%
  distinct(MethylQuant_Score, .keep_all = TRUE) -> d

# Add MitoXplorer information
d %>%
  mutate(SYMBOL = mitocarta[match(ENSEMBL, mitocarta$Gene.ID),]$SYMBOL) %>%
  mutate(Mito.process = mitocarta[match(ENSEMBL, mitocarta$Gene.ID),]$Mitochondrial.Process) %>%
  mutate(Hs.homolog = mitocarta[match(ENSEMBL, mitocarta$Gene.ID),]$`Human Symbol`) -> d

# Export shiftome
d %>%
  filter(!is.na(Mito.process)) %>%
  filter(!grepl("^M", Modifications)) %>%
  arrange(Mito.process, uniprot, Specific_position_in_protein) %>%
  dplyr::select(-c(Protein_Sequence, Origin)) -> shiftome

###########################################
## C: Shiftome and Mitomethylome overlap ##
###########################################

load("Shiftome")
load("Siteome")

shiftome %>%
  dplyr::select(-c(Sequence)) %>%
  dplyr::rename(Protein = uniprot, Intensity = Int, Sequence = mod_seq) %>%
  mutate(Sequence = stringr::str_replace_all(Sequence, "_", "")) %>%
  dplyr::select(-c(MethylQuant_Confidence, Modifications, Mass_Difference, Isotope_Distribution_Correlation, `#_Good_Elution_Profile_Correlations`, PTMs, PTMs_counter)) %>%
  mutate(Confidence_Measure = ifelse(MethylQuant_Score > 40, "C", ifelse(MethylQuant_Score > 20, "D", ifelse(MethylQuant_Score > 0, "H", "I")))) -> shiftome

siteome %>%
  dplyr::rename(Data_File = Best_score_raw_file) %>%
  mutate(MethylQuant_Score = shiftome[match(ID, shiftome$ID),]$MethylQuant_Score) %>%
  mutate(Confidence_Measure = ifelse(MethylQuant_Score > 40 & !is.na(MethylQuant_Score) & Localization_prob > 0.8, "A", ifelse(HL == "L+H" & Localization_prob > 0.8, "B", ifelse(Localization_prob == 1, "E", ifelse(Localization_prob > 0.8, "F", "G"))))) %>%
  dplyr::rename(Sequence = MethylProbabilities) -> siteome

mitomethylome <- full_join(shiftome, siteome) %>%
  arrange(Mito.process, Protein, Specific_position_in_protein) %>%
  dplyr::rename(siteome_score = HL, shiftome_score = MethylQuant_Score) %>%
  mutate(Readable_ID = paste(SYMBOL, "; ", sub("^.*\\;", "", ID), Specific_position_in_protein, sep = "")) %>%
  dplyr::select(Confidence_Measure, Readable_ID, Mito.process, Hs.homolog, everything())

mitomethylome %>%
  group_by(ID) %>%
  top_n(n = -1, Confidence_Measure) %>%#  Counter intuitive that A is the lower end of the scale, this is why n = -1
  distinct(ID, .keep_all = TRUE) -> mitomethylome # Counter intuitive that A is the lower end of the scale, this is why n = -1

write_tsv(mitomethylome, paste("Methylomes/Mitomethylome_", date(), ".txt", sep = ""))

mitomethylome_hc <- filter(mitomethylome, Confidence_Measure != "F" & Confidence_Measure != "G" & Confidence_Measure != "H" & Confidence_Measure != "I")

mitomethylome_hc %>%
  mutate(Amino.ID = paste(Protein, Specific_position_in_protein, sep = ";")) %>%
  distinct(Amino.ID, .keep_all = TRUE) -> cat.mitomethylome_hc

#####################
## D: Conservation ##
#####################

# Read in uniprot accession libraries
load("UniProt.ws_9606") # as up.human
load("UniProt.ws_7227") # as up.fly
up.fly <- up
rm(up)
load("UniProt.ws_10090") # as up.mouse
load("UniProt.ws_10116") # as up.rat

# Externally blast for DROME homologs in DIOPT and import the output table, settings = "Exclude low score (score > 1, unless only match score is 1)", "All" for everything else
diopt <- read_tsv("DIOPT_200522.txt")
tbl_df(diopt) %>%
  dplyr::rename(ID1 = `Species 2 Species Gene ID`, 
                ID2 = `Species 2 Gene ID`,
                Species = `Species 2`, Weight = `Weighted Score`) %>%
  filter(Species %in% c("Human", "Mouse", "Rat")) %>%
  group_by(Species, FlyBaseID) %>%
  top_n(n = 1, Weight) %>%
  distinct(Weight, .keep_all = TRUE) -> diopt_red

diopt_red %>%
  filter(Species == "Human") %>%
  mutate(ID2 = paste("HGNC:", ID2, sep = "")) -> diopt_human

diopt_red %>%
  filter(Species == "Rat") -> diopt_rat

diopt_red %>%
  filter(Species == "Mouse") %>%
  mutate(ID2 = paste("MGI:", ID2, sep = "")) -> diopt_mouse

# Find uniprot accession numbers from ID2 column and filter on the best scored 
mapping_hs <- UniProt.ws::select(up.human, unique(diopt_human$ID2), "UNIPROTKB", "HGNC")
mapping_hs_score <- UniProt.ws::select(up.human, mapping_hs$UNIPROTKB, "SCORE", "UNIPROTKB")
mapping_hs %>%
  mutate(Score = as.numeric(as.character(sub(" out of 5", "", mapping_hs_score[match(mapping_hs$UNIPROTKB, mapping_hs_score$UNIPROTKB),]$SCORE)))) %>%
  group_by(HGNC) %>%
  top_n(n = 1, Score) %>%
  distinct(Score, .keep_all = TRUE) %>%
  dplyr::rename(ID = HGNC) -> mapping_hs

mapping_rn <- UniProt.ws::select(up.rat, unique(diopt_rat$ID2), "UNIPROTKB", "RGD")
mapping_rn_score <- UniProt.ws::select(up.rat, mapping_rn$UNIPROTKB, "SCORE", "UNIPROTKB")
mapping_rn %>%
  mutate(Score = as.numeric(as.character(sub(" out of 5", "", mapping_rn_score[match(mapping_rn$UNIPROTKB, mapping_rn_score$UNIPROTKB),]$SCORE)))) %>%
  group_by(RGD) %>%
  top_n(n = 1, Score) %>%
  distinct(Score, .keep_all = TRUE) %>%
  dplyr::rename(ID = RGD) -> mapping_rn

mapping_mm <- UniProt.ws::select(up.mouse, unique(diopt_mouse$ID2), "UNIPROTKB", "MGI")
mapping_mm_score <- UniProt.ws::select(up.mouse, mapping_mm$UNIPROTKB, "SCORE", "UNIPROTKB")
mapping_mm %>%
  mutate(Score = as.numeric(as.character(sub(" out of 5", "", mapping_mm_score[match(mapping_mm$UNIPROTKB, mapping_mm_score$UNIPROTKB),]$SCORE)))) %>%
  group_by(MGI) %>%
  top_n(n = 1, Score) %>%
  distinct(Score, .keep_all = TRUE) %>%
  dplyr::rename(ID = MGI) -> mapping_mm

# Bind all together
mapping_total <- as.data.frame(rbind(mapping_hs, mapping_rn, mapping_mm))
mapping_total %>%
  dplyr::mutate(ID = str_replace(ID, ".*\\:", "")) -> mapping_total 

diopt_red %>%
  dplyr::rename(ID = ID2) %>%
  dplyr::select(-ID1) -> diopt_red

dplyr::left_join(diopt_red, mapping_total, by = "ID") %>%
  group_by(Species) %>%
  distinct(FlyBaseID, .keep_all = TRUE) %>%
  dplyr::select(FlyBaseID, Species, UNIPROTKB) %>%
  tidyr::spread(Species, UNIPROTKB) -> diopt_spread

diopt_spread$Fly <- mitomethylome_hc[match(diopt_spread$FlyBaseID, mitomethylome_hc$ENSEMBL),]$Protein
diopt_gather <- gather(diopt_spread, "Species", "Uniprot", 2:5)

diopt_list <- setNames(split(diopt_spread, seq(nrow(diopt_spread))), diopt_spread$FlyBaseID)

conservator <- function(x){
  ID <- x$FlyBaseID
  
  x %>%
    dplyr::select(-FlyBaseID) -> x
  
  x <- apply(x, 2, as.list)
  
  x <- Filter(Negate(anyNA), x) # Throws all list elements with NA
  names <- names(x) # Saving the species and their order for later
  
  make_url <- function(x){ #Enables the search with a pre-defined URL
    x <- as.character(paste("http://www.uniprot.org/uniprot/", x, ".fasta", sep = ""))
  } 
  
  x <- lapply(x, make_url) # for every element of the list
  x <- as.vector(unlist(x)) # make a vector out of the list
  
  fasta <- readAAStringSet(x) # convert to appropriate format
  alignment <- msaClustalOmega(fasta) # to make an alignment 
  
  msaConvert <- msaConvert(alignment, "seqinr::alignment") # Convert to a readible format
  msaConvert_seq <- msaConvert$seq
  msaConvert$nam %>%
    stringr::str_extract(pattern = "DROME|RAT|HUMAN|MOUSE") -> names(msaConvert_seq)
  
  # Retrieve peptide of relevance
  mitomethylome_hc[mitomethylome_hc$ENSEMBL == ID,] %>%
    dplyr::select(ID, Sequence, Position_in_peptide) -> list.input 
  mitomethylome_hc.list <- setNames(split(list.input, seq(nrow(list.input))), list.input$ID) # Store several peptides in a list
  
  absolute_position <- function(x){ # The fasta output are all equally long. Get a distance measure for fly and apply it to all others
    if(grepl(stringr::str_replace_all(x$Sequence, "\\.|\\(|\\)|[[:digit:]]|[[a-z]]", ""), msaConvert_seq["DROME"]) == TRUE){
      
      counter <- as.numeric(
        nchar(
          stringr::str_replace_all(
            msaConvert_seq["DROME"],
            paste(stringr::str_replace_all(x$Sequence, "\\.|\\(|\\)|[[:digit:]]|[[a-z]]", ""), ".*", sep = ""),
            "")))
      
    } else { # Because a the alignment might be interrupter by "-", it might occur that a string cannot be detected.
      # Strategy here is to just take the front 7AA or, if also not found, the back 7AA of a particular peptide.
      
      new_query <- paste(substring(stringr::str_replace_all(x$Sequence,"\\.|\\(|\\)|[[:digit:]]|[[a-z]]", ""), 1, 6), 
                         ".*", sep = "") # Front
      
      if(grepl(new_query, msaConvert_seq["DROME"]) == TRUE){
        
        counter <- as.numeric(
          nchar(
            stringr::str_replace_all(
              msaConvert_seq["DROME"],
              paste(stringr::str_replace_all(new_query, "\\.|\\(|\\)|[[:digit:]]|[[a-z]]", ""), ".*", sep = ""),
              "")))
        
      } else {
        new_query <- paste(stringr::str_sub(stringr::str_replace_all(x$Sequence,"\\.|\\(|\\)|[[:digit:]]|[[a-z]]", ""), -6, -1), 
                           ".*", sep = "") # Back
        
        counter <- as.numeric(
          nchar(
            stringr::str_replace_all(
              msaConvert_seq["DROME"],
              paste(stringr::str_replace_all(new_query, "\\.|\\(|\\)|[[:digit:]]|[[a-z]]", ""), ".*", sep = ""),
              "")))
      }
    }
    
    position_in_peptide <- x$Position_in_peptide
    
    msaConvert_seq_list <- as.list(msaConvert_seq, all.names = names(msaConvert_seq)) #Convert to list
    msaConvert_seq_list <- msaConvert_seq_list[names(msaConvert_seq_list) %in% "DROME" == FALSE] # Remove the fly sequence, because we know that position already
    
    residue_at_pos <- function(x){ # Get residue at the very same position as in fly
      #x <- msaConvert_seq_list$HUMAN
      substr <- substr(x, counter + position_in_peptide, counter + position_in_peptide)
      counter <- as.numeric(nchar(stringr::str_replace_all(stringr::str_sub(x, 1, counter + position_in_peptide), "-", "")))
      
      paste(substr, counter, sep = ";")
      
    }
    
    residues <- lapply(msaConvert_seq_list , residue_at_pos)
  }
  
  residues_list <- lapply(mitomethylome_hc.list, absolute_position)
  return(residues_list)
}

conserved_residues <- mclapply(diopt_list, conservator, mc.cores = getOption("mc.cores", 20L))

names <- names(unlist(conserved_residues))
conserved_residues <- data.frame(unlist(conserved_residues))

conserved_residues %>%
  mutate(Species = stringr::str_replace_all(names, ".*\\..*\\.", "")) %>%
  mutate(ID = stringr::str_replace_all(names, "FBgn.......\\.|\\.|HUMAN|RAT|MOUSE", "")) %>%
  dplyr::rename(Conserved_AA = unlist.conserved_residues.) %>%
  spread(Species, Conserved_AA) -> conserved_residues 

mitomethylome_hc <- data.frame(mitomethylome_hc)

mitomethylome_hc$conserved <- conserved_residues[match(mitomethylome_hc$ID, conserved_residues$ID),]

mitomethylome_hc$origin.FLY <- substring(stringr::str_replace_all(mitomethylome_hc$Sequence, "\\.|\\(|\\)|[[:digit:]]|[[a-z]]", ""), 
                                         mitomethylome_hc$Position_in_peptide, mitomethylome_hc$Position_in_peptide)

write.table(mitomethylome_hc, file = paste("Methylomes/Mitomethylome_conserved", date(), ".txt", sep = ""), sep = "\t", row.names = F)

#################
## E. Analysis ##
#################

evo_tracing <- data.frame(as.matrix(mitomethylome_hc[,grepl("conserved|origin", colnames(mitomethylome_hc))]))

# Rat
evo_tracing %>%
  mutate(conserved.RAT = str_replace_all(evo_tracing$conserved.RAT, ";.*", "")) %>%
  spread(origin.FLY, conserved.RAT) -> evo_tracing_rat

circo_input_rat <- data.frame(from = c(rep("Dm.K", length(table(evo_tracing_rat$K))), rep("Dm.R", length(table(evo_tracing_rat$R)))),
                              to = c(names(table(evo_tracing_rat$K)), names(table(evo_tracing_rat$R))),
                              value = c(table(evo_tracing_rat$K), table(evo_tracing_rat$R)))

circo_input_rat %>%
  filter(to != "") -> circo_input_rat

# Mouse
evo_tracing %>%
  mutate(conserved.MOUSE = str_replace_all(evo_tracing$conserved.MOUSE, ";.*", "")) %>%
  spread(origin.FLY, conserved.MOUSE) -> evo_tracing_mouse

circo_input_mouse <- data.frame(from = c(rep("Dm.K", length(table(evo_tracing_mouse$K))), rep("Dm.R", length(table(evo_tracing_mouse$R)))),
                                to = c(names(table(evo_tracing_mouse$K)), names(table(evo_tracing_mouse$R))),
                                value = c(table(evo_tracing_mouse$K), table(evo_tracing_mouse$R)))

circo_input_mouse %>%
  filter(to != "") -> circo_input_mouse

# Human
evo_tracing %>%
  mutate(conserved.HUMAN = str_replace_all(evo_tracing$conserved.HUMAN, ";.*", "")) %>%
  spread(origin.FLY, conserved.HUMAN) -> evo_tracing_human

circo_input_human <- data.frame(from = c(rep("Dm.K", length(table(evo_tracing_human$K))), rep("Dm.R", length(table(evo_tracing_human$R)))),
                                to = c(names(table(evo_tracing_human$K)), names(table(evo_tracing_human$R))),
                                value = c(table(evo_tracing_human$K), table(evo_tracing_human$R)))

circo_input_human %>%
  filter(to != "") -> circo_input_human

# Calculate mutation rate in percent
sum(circo_input_human[circo_input_human$from == "Dm.K",]$value)
sum(circo_input_human[circo_input_human$from == "Dm.R",]$value)

circo_input_human$fraction <- circo_input_human$value / sum(circo_input_human[circo_input_human$from == "Dm.K",]$value)
circo_input_human$fraction[circo_input_human$from == "Dm.R"] <- circo_input_human$value[circo_input_human$from == "Dm.R"] / sum(circo_input_human[circo_input_human$from == "Dm.R",]$value)

methylome_hs_mut.rate <- circo_input_human

#################################
## F. Background mutation rate ##
#################################

conservator <- function(x){
  ID <- x$FlyBaseID
  
  x %>%
    dplyr::select(-FlyBaseID) -> x
  
  x <- apply(x, 2, as.list)
  
  x <- Filter(Negate(anyNA), x) # Throws all list elements with NA
  names <- names(x) # Saving the species and their order for later
  
  make_url <- function(x){ #Enables the search with a pre-defined URL
    x <- as.character(paste("http://www.uniprot.org/uniprot/", x, ".fasta", sep = ""))
  } 
  
  x <- lapply(x, make_url) # for every element of the list
  x <- as.vector(unlist(x)) # make a vector out of the list
  
  fasta <- readAAStringSet(x) # convert to appropriate format
  alignment <- msaClustalOmega(fasta) # to make an alignment 
  
  msaConvert <- msaConvert(alignment, "seqinr::alignment") # Convert to a readible format
  msaConvert_seq <- msaConvert$seq
  msaConvert$nam %>%
    stringr::str_extract(pattern = "DROME|RAT|HUMAN|MOUSE") -> names(msaConvert_seq)
  
  # Retrieve positions of Rs and Ks
  location.KR <- unlist(str_locate_all(msaConvert_seq["DROME"], "K|R")) # Numerical position
  residue.KR <- unlist(sapply(location.KR, function(x) substr(msaConvert_seq["DROME"], x, x))) # Amino acid at that position
  
  location.KR <- data.frame(Residue = residue.KR, Numerical = location.KR)
  location.KR.list <- split(location.KR, seq(nrow(location.KR)))
  
  absolute_position <- function(x){ # The fasta output are all equally long. Get a distance measure for fly and apply it to all others
    position_in_protein <- x$Numerical
    
    msaConvert_seq_list <- as.list(msaConvert_seq, all.names = names(msaConvert_seq)) #Convert to list
    #msaConvert_seq_list <- msaConvert_seq_list[names(msaConvert_seq_list) %in% "DROME" == FALSE] # Remove the fly sequence, because we know that position already
    
    residue_at_pos <- function(x){ # Get residue at the very same position as in fly
      substr(x, position_in_protein, position_in_protein)
    }
    
    residues.evo <- lapply(msaConvert_seq_list , residue_at_pos)
  }
  
  residues_list <- lapply(location.KR.list, absolute_position)
  return(residues_list)
}

conserved_residues_background <- mclapply(diopt_list, conservator, mc.cores = getOption("mc.cores", 20L))

names <- names(unlist(conserved_residues_background))
conserved_residues_background <- data.frame(unlist(conserved_residues_background))

conserved_residues_background %>%
  mutate(Species = stringr::str_replace_all(names, ".*\\..*\\.", "")) %>%
  mutate(ID = stringr::str_replace_all(names, "FBgn.......\\.|\\.|HUMAN|RAT|MOUSE|DROME", "")) %>%
  dplyr::rename(Conserved_AA = unlist.conserved_residues_background.) %>%
  mutate(ENSEMBL = stringr::str_replace_all(names, "\\..*", "")) %>%
  spread(Species, Conserved_AA) -> conserved_residues_background

# Rat
evo_tracing_rat <- spread(conserved_residues_background, DROME, RAT)
circo_input_rat <- data.frame(from = c(rep("Dm.K", length(table(evo_tracing_rat$R))), rep("Dm.R", length(table(evo_tracing_rat$R)))),
                              to = c(names(table(evo_tracing_rat$K)), names(table(evo_tracing_rat$R))),
                              value = c(table(evo_tracing_rat$K), table(evo_tracing_rat$R)))
# Mouse
evo_tracing_mouse <- spread(conserved_residues_background, DROME, MOUSE)
circo_input_mouse <- data.frame(from = c(rep("Dm.K", length(table(evo_tracing_mouse$R))), rep("Dm.R", length(table(evo_tracing_mouse$R)))),
                                to = c(names(table(evo_tracing_mouse$K)), names(table(evo_tracing_mouse$R))),
                                value = c(table(evo_tracing_mouse$K), table(evo_tracing_mouse$R)))

# Human
evo_tracing_human <- spread(conserved_residues_background, DROME, HUMAN)
circo_input_human <- data.frame(from = c(rep("Dm.K", length(table(evo_tracing_human$R))), rep("Dm.R", length(table(evo_tracing_human$R)))),
                                to = c(names(table(evo_tracing_human$K)), names(table(evo_tracing_human$R))),
                                value = c(table(evo_tracing_human$K), table(evo_tracing_human$R)))

# Calculate mutation rate in percent
sum(circo_input_human[circo_input_human$from == "Dm.K",]$value)
sum(circo_input_human[circo_input_human$from == "Dm.R",]$value)

circo_input_human$fraction <- circo_input_human$value / sum(circo_input_human[circo_input_human$from == "Dm.K",]$value)
circo_input_human$fraction[circo_input_human$from == "Dm.R"] <- circo_input_human$value[circo_input_human$from == "Dm.R"] / sum(circo_input_human[circo_input_human$from == "Dm.R",]$value)

methylome.BG_hs_mut.rate <- circo_input_human

# Combine mutation rate in background and methylated residues
left_join(methylome.BG_hs_mut.rate, methylome_hs_mut.rate, by = c("from", "to")) %>%
  dplyr::rename(Background_fraction = fraction.x) %>%
  dplyr::rename(Methylome_fraction = fraction.y) %>% 
  mutate(ID = paste(from, to, sep = ";")) %>%
  dplyr::select(-c(value.x, value.y, from , to)) %>% 
  gather("type", "fraction", 1:2) -> hs_mut.rate

hs_mut.rate$fraction[is.na(hs_mut.rate$fraction)] <- 0

##########################
## G. Proteome coverage ##
##########################

# Step 1: Load all amino acid sequences of proteins in the fly mitocarta
UniProt.ws::select(up.fly, mitocarta$Gene.ID, "UNIPROTKB", "FLYBASE") -> mitocarta.uniprot

# Send and receive requests in parallel or it would take too long and crash
cl <- makeCluster(50)
registerDoParallel(cl)
fasta <- foreach(i = mitocarta.uniprot$UNIPROTKB, .packages="Rcpi") %dopar% getFASTAFromUniProt(i)

mitocarta.uniprot$fasta <- unlist(fasta)

mitocarta.uniprot %>% # Removes the headers and line breaks
  mutate(fasta = str_replace(fasta, ".*[0-9]\\\n", "")) %>%
  mutate(fasta = str_replace_all(fasta, "\\\n", "")) %>% 
  mutate(length = nchar(fasta)) %>%
  arrange(FLYBASE, length) %>%
  group_by(FLYBASE) %>%
  top_n(n = 1) %>% # Meaning, only the longest isoform per FLYBASE ID is kept to avoid bias due to several isoforms
  distinct(FLYBASE, .keep_all = TRUE) %>%
  mutate(counter_trypsin = paste(rep("0", times = nchar(fasta)), collapse = "")) %>%
  mutate(counter_chymo = paste(rep("0", times = nchar(fasta)), collapse = "")) %>%
  mutate(counter_common = paste(rep("0", times = nchar(fasta)), collapse = "")) -> mitocarta_uniprot

# Step 2: Read in all peptides in evidence files
evidence.files <- list.files(pattern = "evidence\\.txt", recursive = TRUE) # Identify all files in subolders that follow the MQ pattern of variable modification tables

## Split into tryspin and chymotrypsin
evidence.files[!grepl("Chymo.*", evidence.files)] -> evidence.files_trypsin
evidence.files[grepl("Chymo.*", evidence.files)] -> evidence.files_chymo

read_table_private <- function(x){
  read.table(x, header = T, sep = "\t")
}

d_trypsin <- lapply(evidence.files_trypsin, read_table_private) # Read in files
names(d_trypsin) <- sub("Data/", "", evidence.files_trypsin) # Give elements of list the names of the files and cut "Data/"
d_chymo <- lapply(evidence.files_chymo, read_table_private)
names(d_chymo) <- sub("Data/", "", evidence.files_chymo) # Give elements of list the names of the files and cut "Data/"

bind_sequences <- function(x){
  unique(x$Sequence)
}

all_peptides_trypsin <- unique(unlist(lapply(d_trypsin, bind_sequences)))
all_peptides_chymo <- unique(unlist(lapply(d_chymo, bind_sequences)))

#  Only use unique peptides

# Step 3: Mark amino acid positions that have been picked up

mod <- function(x,m = 10){
  t1<-floor(x/m)
  return(x-t1*m)
}

peptide_mapper <- function(file = mitocarta_uniprot, all_peptides, value, counter_column){
  
  for(h in c(1:length(all_peptides))){ # Substitute the counter of 0s with 1s once a certain peptide has been detected
    
    test_peptide <- all_peptides[h]
    
    # Counter
    if(mod(which(all_peptides == test_peptide)) == 0){
      print(round(which(all_peptides == test_peptide) / length(all_peptides), 4))
    }
    
    for(i in c(1:length(file$FLYBASE))){
      
      pos <- str_locate(pattern = as.character(test_peptide), file[i, "fasta"])
      
      if(is.na(pos[1]) == FALSE){
        x <- as.character(file[i, counter_column])
        substr(x, pos[1], pos[2]) <- paste(rep(value, times = pos[2] - pos[1] + 1), collapse = "")
        file[i, counter_column] <- x
        
        x <- as.character(file[i, "counter_common"])
        substr(x, pos[1], pos[2]) <- paste(rep(".", times = pos[2] - pos[1] + 1), collapse = "")
        file[i, "counter_common"] <- x
      }
    }
  }
  return(file)
}

mitocarta_uniprot <- peptide_mapper(file = mitocarta_uniprot, all_peptides = all_peptides_chymo, value = "2", counter_column = "counter_chymo")
mitocarta_uniprot <- peptide_mapper(file = mitocarta_uniprot, all_peptides = all_peptides_trypsin, value = "1", counter_column = "counter_trypsin")

# Step 4a: Calculate the number of .s over total: all
counter_collapsed <- paste(mitocarta_uniprot$counter_common, collapse = "")
str_count(counter_collapsed, "\\.") / nchar(counter_collapsed)

# Step 4b: Calculate the number of .s over total: Chymo
counter_collapsed <- paste(mitocarta_uniprot$counter_chymo, collapse = "")
str_count(counter_collapsed, "2") / nchar(counter_collapsed)

# Step 4c: Calculate the number of .s over total: Trypsin
counter_collapsed <- paste(mitocarta_uniprot$counter_trypsin, collapse = "")
str_count(counter_collapsed, "1") / nchar(counter_collapsed)

# Step 5: QC
## Add individual coverage values per protein
count_fraction <- function(x){
  str_count(x, "1") / nchar(x)
}

fraction_pp <- c()
for (i in c(1:length(mitocarta_uniprot$counter_common))){
  fraction_pp_tmp <- str_count(mitocarta_uniprot[i, "counter_common"], "\\.") / nchar(mitocarta_uniprot[i, "counter_common"])
  fraction_pp <- c(fraction_pp, fraction_pp_tmp)
}
mitocarta_uniprot$fraction_pp <- fraction_pp

## Split in mitochondrial categories and boxplot
mitocarta_uniprot$process <- mitocarta[match(mitocarta_uniprot$FLYBASE, mitocarta$Gene.ID),]$Mitochondrial.Process

## Determine amino acids that get preferentially lost
covered_all <- str_locate_all(pattern = "\\.", mitocarta_uniprot$counter_common)
covered_trypsin <- str_locate_all(pattern = "1", mitocarta_uniprot$counter_trypsin)
covered_chymo <- str_locate_all(pattern = "2", mitocarta_uniprot$counter_chymo)

mitocarta_uniprot$uncoveredAA_all <- mitocarta_uniprot$fasta
mitocarta_uniprot$uncoveredAA_chymo <- mitocarta_uniprot$fasta
mitocarta_uniprot$uncoveredAA_trypsin <- mitocarta_uniprot$fasta

for(h in c(1:length(covered_all))){
  for(i in c(1:(length(covered_all[[h]])/2))){
    if(length(covered_all[[h]]) != 0){ # Or a coverage of 0 will result in an error
      substr(mitocarta_uniprot$uncoveredAA_all[h], covered_all[[h]][[i]], covered_all[[h]][[i]]) <- "."
    }
  }
}

for(h in c(1:length(covered_chymo))){
  for(i in c(1:(length(covered_chymo[[h]])/2))){
    if(length(covered_chymo[[h]]) != 0){ # Or a coverage of 0 will result in an error
      substr(mitocarta_uniprot$uncoveredAA_chymo[h], covered_chymo[[h]][[i]], covered_chymo[[h]][[i]]) <- "."
    }
  }
}

for(h in c(1:length(covered_trypsin))){
  for(i in c(1:(length(covered_trypsin[[h]])/2))){
    if(length(covered_trypsin[[h]]) != 0){ # Or a coverage of 0 will result in an error
      substr(mitocarta_uniprot$uncoveredAA_trypsin[h], covered_trypsin[[h]][[i]], covered_trypsin[[h]][[i]]) <- "."
    }
  }
}

uncoveredAA_collapsed_all <- str_replace_all(paste(mitocarta_uniprot$uncoveredAA_all, collapse = ""), "\\.", "")
1-table(str_split(uncoveredAA_collapsed_all, pattern = "")) / table(str_split(paste(mitocarta_uniprot$fasta, collapse = ""), pattern = ""))

uncoveredAA_collapsed_trypsin <- str_replace_all(paste(mitocarta_uniprot$uncoveredAA_trypsin, collapse = ""), "\\.", "")
1-table(str_split(uncoveredAA_collapsed_trypsin, pattern = "")) / table(str_split(paste(mitocarta_uniprot$fasta, collapse = ""), pattern = ""))

uncoveredAA_collapsed_chymo <- str_replace_all(paste(mitocarta_uniprot$uncoveredAA_chymo, collapse = ""), "\\.", "")
1-table(str_split(uncoveredAA_collapsed_chymo, pattern = "")) / table(str_split(paste(mitocarta_uniprot$fasta, collapse = ""), pattern = ""))

########################
## H. Sequence logos ##
########################
mitomethylome_hc_logo <- mitomethylome_hc[!duplicated(sub("(Me|Di|Tri)", "", mitomethylome_hc$Readable_ID)),]
mapping <- UniProt.ws::select(up.fly, mitomethylome_hc_logo$Protein, "SEQUENCE", "UNIPROTKB")
mitomethylome_hc_logo$Sequence <- mapping[match(mitomethylome_hc_logo$Protein, mapping$UNIPROTKB),]$SEQUENCE

windows <- c()
for(i in c(1:dim(mitomethylome_hc_logo)[1])){
  #i = 56
  pos <- mitomethylome_hc_logo[i,"Specific_position_in_protein"]
  length <- nchar(mitomethylome_hc_logo[i,"Sequence"])
  
  # A series of checkpoints in case the sequence window is too small
  if(pos-10 <= 0){
    k = pos-1
    counter_Nterm = abs(pos-11)
  } else {
    k = 10
  }
  if(pos + 10 >= length){
    l = length - pos
    counter_Cterm = abs(length - pos - 10)
  } else {
    l = 10
  }
  
  window_tmp <- substr(mitomethylome_hc_logo[i,"Sequence"], pos - k, pos + l)
  
  # Add X if sequence window too small
  if(k < 10){
    add_N <- paste(rep("X", counter_Nterm), collapse = "")
    window_tmp <- paste(add_N, window_tmp, sep = "")
  }
  if(l < 10){
    add_C <- paste(rep("X", counter_Cterm), collapse = "")
    window_tmp <- paste(window_tmp, add_C, sep = "")
  }
  print(i)
  windows <- c(windows, window_tmp)
}

mitomethylome_hc_logo$window <- windows

# Plot sequence logos with ggseqlogo
seqTot <- ggseqlogo(mitomethylome_hc_logo$window)
seqK <- ggseqlogo(mitomethylome_hc_logo[mitomethylome_hc_logo$origin.FLY == "K",]$window)
seqR <- ggseqlogo(mitomethylome_hc_logo[mitomethylome_hc_logo$origin.FLY == "R",]$window)


##########################################
## I. Position along polypeptide chain ##
#########################################

mitomethylome_hc_logo$fraction_in_chain <- mitomethylome_hc_logo$Specific_position_in_protein / nchar(mitomethylome_hc_logo$Sequence)
mitomethylome_hc_logo$random <- sample(-10:10, dim(mitomethylome_hc_logo)[1], replace=T)/10