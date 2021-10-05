#!/usr/bin/env Rscript

# Check, install and load packages
install.packages("pacman")

pacman::p_load(optparse, phylotools, tidyverse, readxl, stringr)

option_list <- list(
  make_option(c("-p", "--platform"), type="character", default=NULL, 
              help="Type platform. Må være enten: Swift_FHI, Swift_MIK, Artic_Illumina, Artic_Nanopore", metavar="character"),
  make_option(c("-o", "--oppsett"), type = "character", default = NULL,
              help = "Navn på plate/oppsett. Må tilsvare mappenavn. F.eks. FHI133 for Swift_FHI, MIK172 for Swift_MIK, 646 for Artic_Illumina eller Nr134A/Nano for Artic_Nanopore", metavar = "character"),
  make_option(c("-d", "--drop"), type="character", default=NULL, 
              help="Samples to drop (e.g. due to frameshift). Must be on the format of 252116806,252116930,252116980 and match the BN Key", metavar="character"),
  make_option(c("-s", "--specimen"), type="character", default="Unknown", 
              help="Type prøvemateriale. E.g. \"Upper respiratory swab\" [default= %default]", metavar="character"),
  make_option(c("-f", "--fasta"), type="character", default="out_fasta.fasta", 
              help="Navn på fasta-fil som skal submittes. F.eks. FHI133.fasta [default= %default]", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default="FHI133.csv", 
              help="Navn på metadata-fil som skal submittes. F.eks. FHI133.csv [default= %default]", metavar="character")
); 

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

if (is.null(opt$platform)){
  print_help(opt_parser)
  stop("At least one argument must be supplied", call.=FALSE)
}

# Hvis ingen prøver skal droppes, settes til "none"
if (!is.null(opt$drop)){
  drop <- stringr::str_split(opt$drop, ",", simplify = TRUE)
} else if (is.null(opt$drop)) {
  drop <- "none"
}


setwd("/home/jonr/Desktop/")

# Trekke ut KEYs fra BN ---------------------------------------------------

# Les inn BN spørring. Husk å Refreshe og lagre den originale excel-fila først (N:/Virologi/Influensa/2021/Spørringsfiler BN/SQLSERVER_TestBN_Spørring_Entrytable.xlsx)
BN <- suppressWarnings(read_excel("/mnt/N/Virologi/Influensa/2021/Spørringsfiler BN/SQLSERVER_TestBN_Spørring_Entrytable.xlsx") %>% 
  select(KEY, REKVNR, PROVE_TATT, FYLKENAVN, MATERIALE, PROSENTDEKNING_GENOM, DEKNING_NANOPORE, SEKV_OPPSETT_NANOPORE, DEKNING_NANOPORE, SEKV_OPPSETT_SWIFT7, SEQUENCEID_NANO29, SEQUENCEID_SWIFT, COVERAGE_BREADTH_SWIFT, GISAID_PLATFORM, GISAID_EPI_ISL, GENOTYPE_SVART_I_LABWARE, COVERAGE_BREATH_EKSTERNE, SAMPLE_CATEGORY, INNSENDER) %>% 
  rename("Dekning_Artic" = PROSENTDEKNING_GENOM, 
         "Dekning_Swift" = COVERAGE_BREADTH_SWIFT,
         "Dekning_Nano" = DEKNING_NANOPORE))

# Set parameters ----------------------------------------------------------

# Set platform
platform <- opt$platform

# Plate / oppsett
oppsett <- opt$oppsett

# Add fixed columns for metadata
submitter <- "jonbra"
type <- "betacoronavirus"
passage <- "original"
host <- "Human"
gender <- "Unknown"
age <- "Unknown"
status <- "Unknown"
covv_subm_sample_id <- "Unknown"
covv_outbreak <- "Unknown"
covv_add_host_info <- "Unknown"
covv_add_location <- "Unknown"
covv_provider_sample_id <- "Unknown"
covv_last_vaccinated <- "Unknown"
covv_treatment <- "Unknown"
covv_coverage <- "Unknown"
orig_lab <- NA
orig_adr <- NA

# Add specimen
specimen <- opt$specimen


# Set originating labs and adresses ---------------------------------------
lab_lookup_table <- tribble(
  ~`Lab code`, ~`Lab`, ~`Lab address`,
  0,	"Norwegian Institute of Public Health, Department of Virology",	"P.O.Box 222 Skoyen, 0213 Oslo, Norway",
  1,	"Ostfold Hospital Trust - Kalnes, Centre for Laboratory Medicine, Section for gene technology and infection serology",	"P.O.Box 300, N-1714 Graalum, Norway",
  2,	"Akershus University Hospital, Department for Microbiology and Infectious Disease Control",	"P.O.Box 1000, N-1478 Loerenskog, Norway",
  3,	"Oslo University Hospital, Department of Microbiology",	"P.O.Box 4956 Nydalen, N-0424 Oslo, Norway",
  4,	"Furst Medical Laboratory",	"Soeren Bulls vei 25, N-1051 Oslo, Norway",
  5,	"Innlandet Hospital Trust, Division Lillehammer, Department for Medical Microbiology",	"P.O.Box 990, N-2629 Lillehammer, Norway",
  6,	"Medical Microbiology Unit, Department for Laboratory Medicine, Drammen Hospital, Vestre Viken Health Trust", "P.O.Box 800, N-3004 Drammen, Norway",
  7,	"Vestfold Hospital, Toensberg  Department of Microbiology",	"P.O.Box 2168, N-3103 Toensberg, Norway",
  8,	"Unilabs Laboratory Medicine",	"Leirvollen 19, N-3736 Skien, Norway",
  9, NA, NA,
  10,	"Hospital of Southern Norway - Kristiansand, Department of Medical Microbiology",	"P.O.Box 416 Lundsiden, N-4604 Kristiansand, Norway",
  11,	"Dept. of Medical Microbiology, Stavanger University Hospital, Helse Stavanger HF", "P.O.Box 8100, N-4068 Stavanger, Norway",
  12,	"Haukeland University Hospital, Dept. of Microbiology",	"P.O.Box 1400, N-5021 Bergen, Norway",
  13,	"Haugesund Hospital, laboratory for Medical Microbiology",	"P.O.Box 2170, N-5504 Haugesund, Norway",
  14,	"Foerde Hospital, Department of Microbiology", "P.O.Box 1000, N-6807 Foerde, Norway",
  15,	"Department of Medical Microbiology - section Molde, Molde Hospital",	"Parkveien 84, N-6407 Molde, Norway",
  16,	"Department of Medical Microbiology, St. Olavs hospital",	"P.O.box 3250 Torgarden, N-7006 Trondheim, Norway",
  17,	"Haugesund Hospital, laboratory for Medical Microbiology",	"P.O.Box 2170, N-5504 Haugesund, Norway",
  18,	"Nordland Hospital - Bodo, Laboratory Department, Molecular Biology Unit",	"P.O.Box 1480, N-8092 Bodo, Norway",
  19,	"University Hospital of Northern Norway, Department for Microbiology and Infectious Disease Control",	"P.O.Box 56, N-9038 Tromsoe, Norway",
  20, NA, NA,
  21,	NA, NA,
  22,	"Department of medical microbiology, section Aalesund, Aalesund Hospital",	"N-6026 Aalesund, Norway",
  23,	NA, NA,
  24, NA, NA,
  25,	"Telemark Hospital Trust – Skien, Dept. of Medical Microbiology",	"P.O.Box 2900 Kjørbekk, N-3710 Skien",
  26,	"Unilabs Laboratory Medicine",	"Silurveien 2 B, N-0380 Oslo, Norway",
  27,	"Oslo Helse", "Hegdehaugsveien 36, 0352 Oslo"
) %>% 
  mutate(`Lab code` = as.character(`Lab code`))


# Define lookup function to decide originating lab and andress ------------
lookup_function <- function(metadata) {
  for (row in seq_along(metadata$INNSENDER)) {
    if (metadata[row,]$INNSENDER == 0){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[1,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[1,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 1){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[2,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[2,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 2){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[3,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[3,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 3){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[4,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[4,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 4){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[5,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[5,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 5){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[6,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[6,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 6){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[7,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[7,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 7){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[8,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[8,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 8){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[9,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[9,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 9){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[10,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[10,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 10){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[11,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[11,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 11){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[12,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[12,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 12){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[13,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[13,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 13){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[14,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[14,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 14){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[15,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[15,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 15){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[16,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[16,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 16){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[17,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[17,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 17){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[18,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[18,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 18){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[19,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[19,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 19){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[20,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[20,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 20){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[21,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[21,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 21){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[22,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[22,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 22){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[23,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[23,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 23){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[24,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[24,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 24){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[25,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[25,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 25){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[26,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[26,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 26){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[27,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[27,]$`Lab address`
    } else if (metadata[row,]$INNSENDER == 27){
      metadata[row,]$covv_orig_lab <- lab_lookup_table[28,]$Lab
      metadata[row,]$covv_orig_lab_addr <- lab_lookup_table[28,]$`Lab address`
    }
  }
  return(metadata)
}

# Start script ------------------------------------------------------------

if (platform == "Swift_FHI"){
  print ("Platform is Swift FHI")
  # Swift FHI ---------------------------------------------------------------
  
  #### Trekke ut sekvenser ####

  # Trekk ut relevant oppsett og filtrer
  oppsett_details <- BN %>% 
    filter(SEKV_OPPSETT_SWIFT7 == oppsett) %>% 
    # Fjerne evt positiv controll
    filter(str_detect(KEY, "pos", negate = TRUE)) %>% 
    # Filtrer på coverage >= 97% 
    filter(Dekning_Swift >=97) %>% 
    # Fjerne de som mangler Fylkenavn
    filter(!is.na(FYLKENAVN))
  
  # Fjerne keys pga frameshift
  suppressWarnings(
  if (drop != "none"){
    oppsett_details <- oppsett_details[!(oppsett_details$KEY %in% drop),]
  }
  )
  #### Lage metadata ####
  
  # Add platform-specific columns.
  fasta_filename <- opt$fasta
  seq_tech <- "Illumina Swift Amplicon SARS-CoV-2 protocol at Norwegian Sequencing Centre"
  ass_method <- "Assembly by reference based mapping using Bowtie2 with iVar majority rules consensus"
  sub_lab <- "Norwegian Institute of Public Health, Department of Virology"
  address <- "P.O.Box 222 Skoyen, 0213 Oslo, Norway"
  authors <- "Kathrine Stene-Johansen, Kamilla Heddeland Instefjord, Hilde Elshaug, Garcia Llorente Ignacio, Jon Bråte, Engebretsen Serina Beate, Pedersen Benedikte Nevjen, Line Victoria Moen, Debech Nadia, Atiya R Ali, Marie Paulsen Madsen, Rasmus Riis Kopperud, Hilde Vollan, Karoline Bragstad, Olav Hungnes"

  metadata <- oppsett_details %>% 
    # Lage kolonne for "year"
    separate(PROVE_TATT, into = c("Year", NA, NA), sep = "-", remove = FALSE) %>% 
    # Trekke ut sifrene fra 5 og til det siste fra BN KEY
    mutate("Uniq_nr" = str_sub(KEY, start = 5, end = -1)) %>% 
    # Legge til kolonner med fast informasjon for å lage "Virus name" senere
    add_column("Separator" = "/",
               "GISAID_prefix" = "hCoV-19/",
               "Country" = "Norway/", 
               "Continent" = "Europe/") %>% 
    # Make "Virus name" column
    unite("covv_virus_name", c(GISAID_prefix, Country, Uniq_nr, Separator, Year), sep = "", remove = FALSE) %>% 
    # Lage Location-kolonne
    unite("covv_location", c(Continent, Country, FYLKENAVN), sep = "", remove = FALSE) %>% 
    # Legge til faste kolonner
    add_column("submitter" = submitter,
               "fn" = fasta_filename,
               "covv_type" = type,
               "covv_passage" = passage,
               "covv_host" = host,
               "covv_gender" = gender,
               "covv_patient_age" = age,
               "covv_patient_status" = status,
               "covv_specimen" = specimen,
               "covv_seq_technology" = seq_tech,
               "covv_assembly_method" = ass_method,
               "covv_orig_lab" = orig_lab,
               "covv_orig_lab_addr" = orig_adr,
               "covv_subm_lab" = sub_lab,
               "covv_subm_lab_addr" = address,
               "covv_authors" = authors,
               "covv_subm_sample_id" = covv_subm_sample_id,
               "covv_outbreak" = covv_outbreak,
               "covv_add_host_info" = covv_add_host_info,
               "covv_add_location" = covv_add_location,
               "covv_provider_sample_id" = covv_provider_sample_id,
               "covv_last_vaccinated" = covv_last_vaccinated,
               "covv_treatment" = covv_treatment,
               "covv_coverage" = covv_coverage) %>% 
    # Beholde endelige kolonner og rekkefølge
    select("submitter",
           "fn",
           "covv_virus_name", 
           "covv_type",
           "covv_passage",
           "covv_collection_date" = PROVE_TATT,
           "covv_location",
           "covv_host",
           "covv_gender",
           "covv_patient_age",
           "covv_patient_status",
           "covv_specimen",
           "covv_seq_technology",
           "covv_assembly_method",
           "covv_orig_lab",
           "covv_orig_lab_addr",
           "covv_subm_lab",
           "covv_subm_lab_addr",
           "covv_authors",
           "covv_subm_sample_id",
           "covv_outbreak",
           "covv_add_host_info",
           "covv_add_location",
           "covv_provider_sample_id",
           "covv_last_vaccinated",
           "covv_treatment",
           "covv_coverage",
           "INNSENDER")
  
  # Legge inn orig lab og adresse
  metadata <- lookup_function(metadata)
  
  # Remove column INNSENDER
  metadata <- metadata %>% select(-INNSENDER)
  
  # Write csv file
  write_csv(metadata, file = opt$metadata)
  
  
  #### Make fasta file ####
  
  # Search the N: disk for consensus sequences. This could take a few minutes.
  # List relevant folders to search through
  dirs_fhi <- list.dirs("/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina_NSC_FHI/2021/", 
                        recursive = FALSE)
  # Pick our the relevant oppsett
  dir <- dirs_fhi[grep(paste0(oppsett, "\\b"), dirs_fhi)]
  
  # List the files 
  filepaths <- list.files(path = dir,
                          pattern = "ivar\\.consensus\\.masked_Nremoved\\.fa$", 
                          full.names = TRUE,
                          recursive = TRUE)
  
  samples <- str_sub(gsub("SWIFT", "", gsub("_.*","", gsub(".*/","", filepaths))), start = 1, end = -1)
  
  # Find which filepaths to keep
  keep <- vector()
  for (i in seq_along(oppsett_details$KEY)){
    keep[i] <- filepaths[grep(oppsett_details$KEY[i], filepaths)]
  }
  
  # Read each fasta file and combine them to create one file
  # First create empty data frame to fill
  fastas <- data.frame(seq.name = character(),
                       seq.text = character())
  
  for (f in seq_along(keep)){
    tmp <- read.fasta(keep[f])      # read the file
    fastas <- rbind(fastas, tmp)    # append the current file
  }
  
  # Convert to tibble for easier manipulation
  fastas <- as_tibble(fastas)
  
  # Fix names to match KEY
  fastas <- fastas %>% 
    mutate(tmp = str_remove(seq.name, "_ivar_masked")) %>% 
    mutate(KEY = str_remove(tmp, "SWIFT"))
  
  # Sett Virus name som fasta header
  # Først lage en mapping mellom KEY og virus name
  KEY_virus_mapping <- oppsett_details %>% 
    # Lage kolonne for "year"
    separate(PROVE_TATT, into = c("Year", NA, NA), sep = "-", remove = FALSE) %>% 
    # Trekke ut sifrene fra 5 og til det siste fra BN KEY
    mutate("Uniq_nr" = str_sub(KEY, start = 5, end = -1)) %>% 
    # Legge til kolonner med fast informasjon for å lage "Virus name" senere
    add_column("Separator" = "/",
               "GISAID_prefix" = "hCoV-19/",
               "Country" = "Norway/", 
               "Continent" = "Europe/") %>% 
    # Make "Virus name" column
    unite("covv_virus_name", c(GISAID_prefix, Country, Uniq_nr, Separator, Year), sep = "", remove = FALSE) %>% 
    select(KEY, covv_virus_name)
  
  
  fastas <- left_join(fastas, KEY_virus_mapping, by = "KEY") %>% 
    select(covv_virus_name, seq.text) %>% 
    rename(`seq.name` = covv_virus_name)
  
  # Lagre fastafilen
  dat2fasta(fastas, outfile = fasta_filename)
} else if (platform == "Swift_MIK") {
  print ("Platform is Swift MIK")
  # Swift MIK/OUS------------------------------------------------------------
  #### Trekke ut sekvenser ####

  
  # Trekk ut relevant oppsett og filtrer
  oppsett_details <- BN %>% 
    filter(SEKV_OPPSETT_SWIFT7 == oppsett) %>% 
    # Fjerne evt positiv controll
    filter(str_detect(KEY, "pos", negate = TRUE)) %>% 
    # Filtrer på coverage >= 97% 
    filter(Dekning_Swift >=97) %>% 
    # Fjerne de som mangler Fylkenavn
    filter(!is.na(FYLKENAVN))
  
  # Fjerne keys pga frameshift
  suppressWarnings(
    if (drop != "none"){
      oppsett_details <- oppsett_details[!(oppsett_details$KEY %in% drop),]
    }
  )
  
  
  #### Lage metadata ####
  
  # Add platform-specific columns.
  fasta_filename <- opt$fasta
  seq_tech <- "Illumina Swift Amplicon SARS-CoV-2 protocol at Norwegian Sequencing Centre"
  ass_method <- "Assembly by reference based mapping using Bowtie2 with iVar majority rules consensus"
  sub_lab <- "Norwegian Institute of Public Health, Department of Virology"
  address <- "P.O.Box 222 Skoyen, 0213 Oslo, Norway"
  authors <- "Mona Holberg-Petersen, Lise Andresen, Cathrine Fladeby, Mariann Nilsen, Teodora Plamenova Ribarska, Pål Marius Bjørnstad, Gregor D. Gilfillan, Arvind Yegambaram Meenakshi Sundaram,Kathrine Stene-Johansen, Kamilla Heddeland Instefjord, Hilde Elshaug, Garcia Llorente Ignacio, Jon Bråte, Pedersen Benedikte Nevjen, Line Victoria Moen, Rasmus Riis Kopperud, Hilde Vollan, Olav Hungnes, Karoline Bragstad"
  
  
  metadata <- IDs %>% 
    # Lage kolonne for "year"
    separate(PROVE_TATT, into = c("Year", NA, NA), sep = "-", remove = FALSE) %>% 
    # Trekke ut sifrene fra 5 og til det siste fra BN KEY
    mutate("Uniq_nr" = str_sub(KEY, start = 1, end = -1)) %>%
    # Legge til kolonner med fast informasjon for å lage "Virus name" senere
    add_column("Separator" = "/",
               "GISAID_prefix" = "hCoV-19/",
               "Country" = "Norway/", 
               "Continent" = "Europe/") %>% 
    # Make "Virus name" column
    unite("covv_virus_name", c(GISAID_prefix, Country, Uniq_nr, Separator, Year), sep = "", remove = FALSE) %>% 
    # Lage Location-kolonne
    unite("covv_location", c(Continent, Country, FYLKENAVN), sep = "", remove = FALSE) %>% 
    # Legge til faste kolonner
    add_column("submitter" = submitter,
               "fn" = fasta_filename,
               "covv_type" = type,
               "covv_passage" = passage,
               "covv_host" = host,
               "covv_gender" = gender,
               "covv_patient_age" = age,
               "covv_patient_status" = status,
               "covv_specimen" = specimen,
               "covv_seq_technology" = seq_tech,
               "covv_assembly_method" = ass_method,
               "covv_orig_lab" = orig_lab,
               "covv_orig_lab_addr" = orig_adr,
               "covv_subm_lab" = sub_lab,
               "covv_subm_lab_addr" = address,
               "covv_authors" = authors,
               "covv_subm_sample_id" = covv_subm_sample_id,
               "covv_outbreak" = covv_outbreak,
               "covv_add_host_info" = covv_add_host_info,
               "covv_add_location" = covv_add_location,
               "covv_provider_sample_id" = covv_provider_sample_id,
               "covv_last_vaccinated" = covv_last_vaccinated,
               "covv_treatment" = covv_treatment,
               "covv_coverage" = covv_coverage) %>% 
    # Beholde endelige kolonner og rekkefølge
    select("submitter",
           "fn",
           "covv_virus_name", 
           "covv_type",
           "covv_passage",
           "covv_collection_date" = PROVE_TATT,
           "covv_location",
           "covv_host",
           "covv_gender",
           "covv_patient_age",
           "covv_patient_status",
           "covv_specimen",
           "covv_seq_technology",
           "covv_assembly_method",
           "covv_orig_lab",
           "covv_orig_lab_addr",
           "covv_subm_lab",
           "covv_subm_lab_addr",
           "covv_authors",
           "covv_subm_sample_id",
           "covv_outbreak",
           "covv_add_host_info",
           "covv_add_location",
           "covv_provider_sample_id",
           "covv_last_vaccinated",
           "covv_treatment",
           "covv_coverage",
           "INNSENDER")
  
  # Legge inn orig lab og adresse
  metadata <- lookup_function(metadata)
  
  # Remove column INNSENDER
  metadata <- metadata %>% select(-INNSENDER)
  
  
  # Write csv file
  write_csv(metadata, file = opt$metadata)
  
  #### Make fasta file ####
  
  # Search the N: disk for consensus sequences. This could take a few minutes.
  # List relevant folders to search through
  dirs_fhi <- list.dirs("/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina_NSC_MIK", 
                        recursive = FALSE)
  # Pick our the relevant oppsett
  dir <- dirs_fhi[grep(paste0(oppsett), dirs_fhi)]
  
  # List the files 
  filepaths <- list.files(path = dir,
                          pattern = "ivar\\.consensus\\.masked_Nremoved\\.fa$", 
                          full.names = TRUE,
                          recursive = TRUE)
  
  samples <- gsub("_.*","", gsub(".*/","", filepaths))
  
  # Find which filepaths to keep
  keep <- vector()
  for (i in seq_along(oppsett_details$SEQUENCEID_SWIFT)){
    keep[i] <- filepaths[grep(oppsett_details$SEQUENCEID_SWIFT[i], filepaths)]
  }
  
  # Read each fasta file and combine them to create one file
  # First create empty data frame to fill
  fastas <- data.frame(seq.name = character(),
                       seq.text = character())
  
  for (f in seq_along(keep)){
    tmp <- read.fasta(keep[f])      # read the file
    fastas <- rbind(fastas, tmp)    # append the current file
  }
  
  # Convert to tibble for easier manipulation
  fastas <- as_tibble(fastas)
  
  # Fix names to match SEQUENCEID_SWIFT
  fastas <- fastas %>% 
    mutate(SEQUENCEID_SWIFT = str_remove(seq.name, "_ivar_masked")) #%>% 
    #mutate(SEQUENCEID_SWIFT = str_remove(tmp, protocol))
  
  # Sett Virus name som fasta header
  # Først lage en mapping mellom SEQUENCEID_SWIFT og virus name
  KEY_virus_mapping <- oppsett_details %>% 
    # Lage kolonne for "year"
    separate(PROVE_TATT, into = c("Year", NA, NA), sep = "-", remove = FALSE) %>% 
    # Trekke ut sifrene fra 5 og til det siste fra BN KEY
    mutate("Uniq_nr" = str_sub(KEY, start = 1, end = -1)) %>% 
    # Legge til kolonner med fast informasjon for å lage "Virus name" senere
    add_column("Separator" = "/",
               "GISAID_prefix" = "hCoV-19/",
               "Country" = "Norway/", 
               "Continent" = "Europe/") %>% 
    # Make "Virus name" column
    unite("covv_virus_name", c(GISAID_prefix, Country, Uniq_nr, Separator, Year), sep = "", remove = FALSE) %>% 
    select(SEQUENCEID_SWIFT, KEY, covv_virus_name)
  
  
  fastas <- left_join(fastas, KEY_virus_mapping, by = "SEQUENCEID_SWIFT") %>% 
    select(`seq.name` = covv_virus_name, 
           seq.text)
  
  # Lagre fastafilen
  dat2fasta(fastas, outfile = fasta_filename)
} else if (platform == "Artic_Illumina") {
  print ("Platform is Artic Illumina")
  # Artic Illumina ----------------------------------------------------------
  
  ##### Trekke ut sekvenser ####
 
  # Trekk ut relevant oppsett og filtrer
  oppsett_details <- BN %>%
    filter(str_detect(SAMPLE_CATEGORY, oppsett)) %>% 
    # Fjerne evt positiv controll
    filter(str_detect(KEY, "pos", negate = TRUE)) %>% 
    # Filtrer på coverage >= 97% 
    filter(Dekning_Artic >=97) %>% 
    # Fjerne de som mangler Fylkenavn
    filter(!is.na(FYLKENAVN))
  
  # Fjerne keys pga frameshift
  suppressWarnings(
    if (drop != "none"){
      oppsett_details <- oppsett_details[!(oppsett_details$KEY %in% drop),]
    }
  )
  
  
  #### Lage metadata ####
  
  # Add platform-specific columns.
  fasta_filename <- opt$fasta
  seq_tech <- "Illumina Swift Amplicon SARS-CoV-2 protocol at Norwegian Sequencing Centre"
  ass_method <- "Assembly by reference based mapping using Bowtie2 with iVar majority rules consensus"
  sub_lab <- "Norwegian Institute of Public Health, Department of Virology"
  address <- "P.O.Box 222 Skoyen, 0213 Oslo, Norway"
  authors <- "Kathrine Stene-Johansen, Kamilla Heddeland Instefjord, Hilde Elshaug, Garcia Llorente Ignacio, Jon Bråte, Engebretsen Serina Beate, Pedersen Benedikte Nevjen, Line Victoria Moen, Debech Nadia, Atiya R Ali, Marie Paulsen Madsen, Rasmus Riis Kopperud, Hilde Vollan, Karoline Bragstad, Olav Hungnes"
  
  metadata <- oppsett_details %>% 
    # Lage kolonne for "year"
    separate(PROVE_TATT, into = c("Year", NA, NA), sep = "-", remove = FALSE) %>%
    # Trekke ut sifrene fra 5 og til det siste fra BN KEY
    mutate("Uniq_nr" = str_sub(KEY, start = 5, end = -1)) %>% 
    # Legge til kolonner med fast informasjon for å lage "Virus name" senere
    add_column("Separator" = "/",
               "GISAID_prefix" = "hCoV-19/",
               "Country" = "Norway/", 
               "Continent" = "Europe/") %>% 
    # Make "Virus name" column
    unite("covv_virus_name", c(GISAID_prefix, Country, Uniq_nr, Separator, Year), sep = "", remove = FALSE) %>% 
    # Lage Location-kolonne
    unite("covv_location", c(Continent, Country, FYLKENAVN), sep = "", remove = FALSE) %>% 
    # Legge til faste kolonner
    add_column("submitter" = submitter,
               "fn" = fasta_filename,
               "covv_type" = type,
               "covv_passage" = passage,
               "covv_host" = host,
               "covv_gender" = gender,
               "covv_patient_age" = age,
               "covv_patient_status" = status,
               "covv_specimen" = specimen,
               "covv_seq_technology" = seq_tech,
               "covv_assembly_method" = ass_method,
               "covv_orig_lab" = orig_lab,
               "covv_orig_lab_addr" = orig_adr,
               "covv_subm_lab" = sub_lab,
               "covv_subm_lab_addr" = address,
               "covv_authors" = authors,
               "covv_subm_sample_id" = covv_subm_sample_id,
               "covv_outbreak" = covv_outbreak,
               "covv_add_host_info" = covv_add_host_info,
               "covv_add_location" = covv_add_location,
               "covv_provider_sample_id" = covv_provider_sample_id,
               "covv_last_vaccinated" = covv_last_vaccinated,
               "covv_treatment" = covv_treatment,
               "covv_coverage" = covv_coverage) %>% 
    # Beholde endelige kolonner og rekkefølge
    select("submitter",
           "fn",
           "covv_virus_name", 
           "covv_type",
           "covv_passage",
           "covv_collection_date" = PROVE_TATT,
           "covv_location",
           "covv_host",
           "covv_gender",
           "covv_patient_age",
           "covv_patient_status",
           "covv_specimen",
           "covv_seq_technology",
           "covv_assembly_method",
           "covv_orig_lab",
           "covv_orig_lab_addr",
           "covv_subm_lab",
           "covv_subm_lab_addr",
           "covv_authors",
           "covv_subm_sample_id",
           "covv_outbreak",
           "covv_add_host_info",
           "covv_add_location",
           "covv_provider_sample_id",
           "covv_last_vaccinated",
           "covv_treatment",
           "covv_coverage",
           "INNSENDER")
  
  # Legge inn orig lab og adresse
  metadata <- lookup_function(metadata)
  
  # Remove column INNSENDER
  metadata <- metadata %>% select(-INNSENDER)
  
  # Write csv file
  write_csv(metadata, file = opt$metadata)
  
  #### Make fasta file ####
  
  # Search the N: disk for consensus sequences. This could take a few minutes.
  # List relevant folders to search through
  dirs_fhi <- list.dirs("/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Illumina/2021", 
                        recursive = FALSE)
  # Pick our the relevant oppsett
  dir <- dirs_fhi[grep(paste0(oppsett), dirs_fhi)]
  
  # List the files 
  filepaths <- list.files(path = dir,
                          pattern = "consensus\\.fa$", 
                          full.names = TRUE,
                          recursive = TRUE)

  # Dropper det siste tallet.
  samples <- str_sub(gsub("Artic", "", gsub("_.*","", gsub(".*/","", filepaths))), start = 1, end = -2)
  
  # Find which filepaths to keep
  keep <- vector()
  for (i in seq_along(oppsett_details$KEY)){
    keep[i] <- filepaths[grep(oppsett_details$KEY[i], filepaths)]
  }
  
  # Read each fasta file and combine them to create one file
  # First create empty data frame to fill
  fastas <- data.frame(seq.name = character(),
                       seq.text = character())
  
  for (f in seq_along(keep)){
    tmp <- read.fasta(keep[f])      # read the file
    fastas <- rbind(fastas, tmp)    # append the current file
  }
  
  # Convert to tibble for easier manipulation
  fastas <- as_tibble(fastas)
  
  # Fix names to match KEY
  fastas <- fastas %>% 
    mutate(tmp = str_remove(seq.name, "Artic")) %>%
    mutate(KEY = str_sub(tmp, start = 1, end = -2))
  
  # Sett Virus name som fasta header
  # Først lage en mapping mellom KEY og virus name
  KEY_virus_mapping <- oppsett_details %>% 
    # Lage kolonne for "year"
    separate(PROVE_TATT, into = c("Year", NA, NA), sep = "-", remove = FALSE) %>% 
    # Trekke ut sifrene fra 5 og til det siste fra BN KEY
    mutate("Uniq_nr" = str_sub(KEY, start = 5, end = -1)) %>% 
    # Legge til kolonner med fast informasjon for å lage "Virus name" senere
    add_column("Separator" = "/",
               "GISAID_prefix" = "hCoV-19/",
               "Country" = "Norway/", 
               "Continent" = "Europe/") %>% 
    # Make "Virus name" column
    unite("covv_virus_name", c(GISAID_prefix, Country, Uniq_nr, Separator, Year), sep = "", remove = FALSE) %>% 
    select(KEY, covv_virus_name)
  
  
  fastas <- left_join(fastas, KEY_virus_mapping, by = "KEY") %>% 
    select(`seq.name` = covv_virus_name, 
           seq.text)
  
  # Lagre fastafilen
  dat2fasta(fastas, outfile = fasta_filename)
} else if (platform == "Artic_Nanopore") {
  print ("Platform is Artic Nanopore")
  # Artic Nanopore ----------------------------------------------------------
  
  ##### Trekke ut sekvenser ####
  
  # Trekk ut relevant oppsett og filtrer
  oppsett_details <- BN %>%
    filter(str_detect(SEKV_OPPSETT_NANOPORE, oppsett)) %>% 
    # Fjerne evt positiv controll
    filter(str_detect(KEY, "pos", negate = TRUE)) %>% 
    # Filtrer på coverage >= 97% 
    filter(Dekning_Nano >=97) %>% 
    # Fjerne de som mangler Fylkenavn
    filter(!is.na(FYLKENAVN))
  
  # Fjerne keys pga frameshift
  suppressWarnings(
    if (drop != "none"){
      oppsett_details <- oppsett_details[!(oppsett_details$KEY %in% drop),]
    }
  )
  
  
  #### Lage metadata ####
  
  # Add platform-specific columns.
  fasta_filename <- opt$fasta
  seq_tech <- "Nanopore GridIon, Midnight protocol modified"
  ass_method <- "Assembly by reference based mapping using the Artic Nanopore protocol with medaka"
  sub_lab <- "Norwegian Institute of Public Health, Department of Virology"
  address <- "P.O.Box 222 Skoyen, 0213 Oslo, Norway"
  authors <- "Kathrine Stene-Johansen, Kamilla Heddeland Instefjord, Hilde Elshaug, Garcia Llorente Ignacio, Jon Bråte, Engebretsen Serina Beate, Pedersen Benedikte Nevjen, Line Victoria Moen, Debech Nadia, Atiya R Ali, Marie Paulsen Madsen, Rasmus Riis Kopperud, Hilde Vollan, Karoline Bragstad, Olav Hungnes"
  
  metadata <- oppsett_info %>% 
    # Lage kolonne for "year"
    separate(PROVE_TATT, into = c("Year", NA, NA), sep = "-", remove = FALSE) %>%
    # Trekke ut sifrene fra 5 og til det siste fra BN KEY
    mutate("Uniq_nr" = str_sub(KEY, start = 5, end = -1)) %>% 
    # Legge til kolonner med fast informasjon for å lage "Virus name" senere
    add_column("Separator" = "/",
               "GISAID_prefix" = "hCoV-19/",
               "Country" = "Norway/", 
               "Continent" = "Europe/") %>% 
    # Make "Virus name" column
    unite("covv_virus_name", c(GISAID_prefix, Country, Uniq_nr, Separator, Year), sep = "", remove = FALSE) %>% 
    # Lage Location-kolonne
    unite("covv_location", c(Continent, Country, FYLKENAVN), sep = "", remove = FALSE) %>% 
    # Legge til faste kolonner
    add_column("submitter" = submitter,
               "fn" = fasta_filename,
               "covv_type" = type,
               "covv_passage" = passage,
               "covv_host" = host,
               "covv_gender" = gender,
               "covv_patient_age" = age,
               "covv_patient_status" = status,
               "covv_specimen" = specimen,
               "covv_seq_technology" = seq_tech,
               "covv_assembly_method" = ass_method,
               "covv_orig_lab" = orig_lab,
               "covv_orig_lab_addr" = orig_adr,
               "covv_subm_lab" = sub_lab,
               "covv_subm_lab_addr" = address,
               "covv_authors" = authors,
               "covv_subm_sample_id" = covv_subm_sample_id,
               "covv_outbreak" = covv_outbreak,
               "covv_add_host_info" = covv_add_host_info,
               "covv_add_location" = covv_add_location,
               "covv_provider_sample_id" = covv_provider_sample_id,
               "covv_last_vaccinated" = covv_last_vaccinated,
               "covv_treatment" = covv_treatment,
               "covv_coverage" = covv_coverage) %>% 
    # Beholde endelige kolonner og rekkefølge
    select("submitter",
           "fn",
           "covv_virus_name", 
           "covv_type",
           "covv_passage",
           "covv_collection_date" = PROVE_TATT,
           "covv_location",
           "covv_host",
           "covv_gender",
           "covv_patient_age",
           "covv_patient_status",
           "covv_specimen",
           "covv_seq_technology",
           "covv_assembly_method",
           "covv_orig_lab",
           "covv_orig_lab_addr",
           "covv_subm_lab",
           "covv_subm_lab_addr",
           "covv_authors",
           "covv_subm_sample_id",
           "covv_outbreak",
           "covv_add_host_info",
           "covv_add_location",
           "covv_provider_sample_id",
           "covv_last_vaccinated",
           "covv_treatment",
           "covv_coverage",
           "INNSENDER")
  
  # Legge inn orig lab og adresse
  metadata <- lookup_function(metadata)
  
  # Remove column INNSENDER
  metadata <- metadata %>% select(-INNSENDER)
  
  # Write csv file
  write_csv(metadata, file = opt$metadata)
  
  #### Make fasta file ####
  
  # Search the N: disk for consensus sequences. This could take a few minutes.
  # List relevant folders to search through
  dirs_fhi <- list.dirs("/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/1-Nanopore/2021", 
                        recursive = FALSE)
  # Pick our the relevant oppsett
  oppsett <- gsub("Nr", "", (gsub("/Nano", "", oppsett)))
  dir <- dirs_fhi[grep(paste0(oppsett), dirs_fhi)]
  #NB! Denne må automatiseres - hvordan?
  dir <- dir[-grep("RAPID", dir)]
  
  # List the files 
  filepaths <- list.files(path = dir,
                          pattern = "consensus\\.fasta$", 
                          full.names = TRUE,
                          recursive = TRUE)
  
  # Dropper det siste tallet.
  samples <- str_sub(gsub("_.*","", gsub(".*/","", filepaths)), start = 1, end = -2)
  
  # Find which filepaths to keep
  keep <- vector()
  for (i in seq_along(oppsett_info$KEY)){
    keep[i] <- filepaths[grep(oppsett_info$KEY[i], filepaths)]
  }
  
  # Read each fasta file and combine them to create one file
  # First create empty data frame to fill
  fastas <- data.frame(seq.name = character(),
                       seq.text = character())
  
  for (f in seq_along(keep)){
    tmp <- read.fasta(keep[f])      # read the file
    fastas <- rbind(fastas, tmp)    # append the current file
  }
  
  # Convert to tibble for easier manipulation
  fastas <- as_tibble(fastas)
  
  # Fix names to match KEY
  fastas <- fastas %>% 
    # Legger inn denne først for da kan jeg senere slice stringen fra første til nest siste karakter. Mer robust
    mutate(tmp = gsub("_.*", "", seq.name)) %>% 
    mutate(KEY = str_sub(tmp, start = 1, end = -2))
  
  # Sett Virus name som fasta header
  # Først lage en mapping mellom KEY og virus name
  KEY_virus_mapping <- oppsett_info %>% 
    # Lage kolonne for "year"
    separate(PROVE_TATT, into = c("Year", NA, NA), sep = "-", remove = FALSE) %>% 
    # Trekke ut sifrene fra 5 og til det siste fra BN KEY
    mutate("Uniq_nr" = str_sub(KEY, start = 5, end = -1)) %>% 
    # Legge til kolonner med fast informasjon for å lage "Virus name" senere
    add_column("Separator" = "/",
               "GISAID_prefix" = "hCoV-19/",
               "Country" = "Norway/", 
               "Continent" = "Europe/") %>% 
    # Make "Virus name" column
    unite("covv_virus_name", c(GISAID_prefix, Country, Uniq_nr, Separator, Year), sep = "", remove = FALSE) %>% 
    select(KEY, covv_virus_name)
  
  
  fastas <- left_join(fastas, KEY_virus_mapping, by = "KEY") %>% 
    select(`seq.name` = covv_virus_name, 
           seq.text)
  
  # Lagre fastafilen
  dat2fasta(fastas, outfile = fasta_filename)
}
