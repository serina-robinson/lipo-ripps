# Install packages
pacman::p_load("tidyverse", "Biostrings", "data.table")

# RODEO v 2.2.0.
# python3 rodeo_main.py -fn 8 -j 4 -out 5000_NHLP_accessions_pull.txt
# python3 rodeo_main.py -fn 8 -j 4 -out 5849_Nif11_accessions_pull.txt

# Read in the rodeo2 output
nhlp_dat <- read_csv("data/NHLP_all_queries/5000_NHLP_nbsize_8_rodeo2/main_co_occur.csv") %>%
  janitor::clean_names()

nif11_dat <- read_csv("data/Nif11/5849_Nif11_nbsize_8_rodeo2_out/main_co_occur.csv") %>%
  janitor::clean_names()

# Read in precursors with L....[GSA][GSA] cleavage site
nhlp <- readAAStringSet("data/combined_Nif11_NHLP_cleavage_motifs/1574_NHLP_no_cat_with_cleavage_site.faa")
nif11 <- readAAStringSet("data/combined_Nif11_NHLP_cleavage_motifs/2819_Nif11_with_cleavage_site.faa")
nhlp_queries <- unique(word(names(nhlp), sep = " ", 1))
nif11_queries <- unique(word(names(nif11), sep = "\\:", 1))

# NHLP and Nif11 trimmed
nhlp_trim <- nhlp_dat %>%
  dplyr::filter(query %in% nhlp_queries)
length(unique(nhlp_trim$query)) # 1565

nif11_trim <- nif11_dat %>%
  dplyr::filter(query %in% nif11_queries)
length(unique(nif11_trim$query)) # 2163 

# Search for NHase_beta hits
nhase_beta <- nhlp_trim %>%
  dplyr::filter(grepl("NHase_beta", name1) | 
                  grepl("NHase_beta", name2) | 
                  grepl("NHase_beta", name3)) %>%
  dplyr::select(query) %>%
  pull() %>%
  unique()

# Combine everything
# Read in polytheonamide
poly <- read_csv("data/characterized_nps/Entotheonella_polytheonamide_nb.csv")
aero <- read_csv("data/characterized_nps/Microvirgula_aeronamide_nb.csv") %>%
  janitor::clean_names()
 
dat_comb <- poly %>% # add polytheonamide cluster
  bind_rows(aero) %>%  # add aeronamide cluster
  bind_rows(nhlp_trim) %>%
  dplyr::filter(!query %in% nhase_beta) %>% # remove NHase beta
  bind_rows(nif11_trim) 

# Pull genus list to merge with Genome Taxonomy Database
dat_tax <- dat_comb %>%
  mutate(genus_1 = case_when(grepl("Candidatus|^cf\\.|^aff\\.|uncultured|^bacterium", genus_species) ~ word(genus_species, sep = " ", 2),
                             grepl("\\]", genus_species) ~ gsub("\\[", "", word(genus_species, sep = "\\]", 1)),
                             TRUE ~ word(genus_species, sep = " ", 1))) %>%
  select(genus_1) %>%
  pull() %>%
  unique()

# Remove taxonomy exceptions
exceptions <- dat_tax[!grepl("[A-Z]", dat_tax)] # remove if it isn't a Latin genus name
dat_tax_fin <- dat_tax[!grepl(paste0(paste0("^",exceptions), collapse = "|"), dat_tax)]

# GTDB taxonomy 
unparsed_tax <- fread("data/taxonomy/gtdb_bac120_tax.tsv", data.table = F,
                      fill = TRUE, header = F, sep = "") %>%
  dplyr::mutate(unparse_with_GCF = word(V1, sep = ";s__", 1)) %>%
  dplyr::mutate(unparse = word(unparse_with_GCF, sep = "d__Bacteria;", 2)) %>%
  dplyr::select(unparse) %>%
  distinct() %>%
  pull() 

gtdb_tax <- lapply(1:length(dat_tax_fin), function(x) {
  unparsed_tax[grep(dat_tax_fin[x], unparsed_tax)[1]]})

newdf <- data.frame(dat_tax_fin, unlist(gtdb_tax))
colnames(newdf) <- c("genus_key", "gtdb_tax")

# Merge taxonomy with neighborhoods
dat_match <- dat_comb %>%
  mutate(genus_1 = case_when(grepl("Candidatus|^cf\\.|^aff\\.|uncultured|^bacterium", genus_species) ~ word(genus_species, sep = " ", 2),
                             grepl("\\]", genus_species) ~ gsub("\\[", "", word(genus_species, sep = "\\]", 1)),
                             TRUE ~ word(genus_species, sep = " ", 1))) %>%
  dplyr::left_join(., newdf, by = c("genus_1" = "genus_key")) %>%
  dplyr::mutate(phylum = word(word(word(gtdb_tax, sep = "p__", 2), sep = ";", 1), sep = "_", 1)) %>%
  dplyr::mutate(class = word(word(word(gtdb_tax, sep = "c__", 2), sep = ";", 1), sep = "_", 1)) %>%
  dplyr::mutate(order = word(word(word(gtdb_tax, sep = "o__", 2), sep = ";", 1), sep = "_", 1)) %>%
  dplyr::mutate(family = word(word(word(gtdb_tax, sep = "f__", 2), sep = ";", 1), sep = "_", 1))

# Check how many NAs
dat_na <- dat_match %>%
  filter(is.na(gtdb_tax))

# Now get counts for PFAMs 
dat_trim <- dat_match %>%
  dplyr::filter(!is.na(pfam_id1)) %>%
  arrange(query) %>%
  dplyr::mutate(idvar = paste0(nucleotide_acc, "_", protein_acc)) %>%
  group_by(idvar) %>%
  dplyr::slice(1)

dat_pfam_tally <- dat_trim %>%
  select(-query) %>%
  distinct(.keep_all = T) %>%
  mutate(pfam_long = paste0(name1, "_", pfam_id1)) %>% 
  group_by(nucleotide_acc) %>%
  dplyr::add_count(pfam_long) %>%
  mutate(genvar = paste0(genus_species, '_', nucleotide_acc))

top_pfams <- sort(table(dat_pfam_tally$pfam_long), decreasing = T)
pfams_to_keep <- names(top_pfams)[top_pfams > 9]
length(pfams_to_keep) 

# WHICH KEEP - DO NOT INCLUDE IF QUERY APPEARS MORE THAN ONCE IN PROTEIN ACC
quers <- unique(dat_match$query) #3724 queries 

# Count how many BGCs
dat_quers <- dat_match %>%
  dplyr::filter(protein_acc %in% quers) %>%
  mutate(unique_idvar = paste0(nucleotide_acc, "_", protein_acc)) %>%
  group_by(unique_idvar) %>%
  dplyr::slice(1) # only one combination of nucleotide acc and query 
# in the same neighborhood, to prevent double-counting

dat_pfam_tally_trim <- dat_match %>%
  dplyr::filter(query %in% dat_quers$query) %>%
  dplyr::filter(!pfam_id1 %in% c('PF02979', 'PF14407', 'PF07862')) %>%
  dplyr::filter(!pfam_id2 %in% c('PF02979', 'PF14407', 'PF07862')) %>%
  #  "PF02211", # NHase_beta should already be removed 
  #  "PF02979", "PF14407", # NHLP precursors
  #  "PF07862", # Nif11 precursor
  dplyr::filter(!is.na(pfam_id1)) %>%
  dplyr::mutate(unid = paste0(nucleotide_acc, "_", protein_acc, "_", start)) %>%
  distinct(unid, .keep_all = T) %>%
  mutate(pfam_long = paste0(pfam_id1, "_", name1))
length(unique(dat_pfam_tally_trim$query)) # 2345, excluding neighborhoods of 1

precs <- dat_match %>%
  dplyr::filter(query %in% dat_pfam_tally_trim$query) %>%
  dplyr::filter(protein_acc %in% dat_pfam_tally_trim$query) %>%
  distinct(, .keep_all = T) %>%
  mutate(pfam_long = paste0(pfam_id1, "_", name1))

length(unique(precs$query)) # 2347 unique queries

# Condense into clans
clans <- read_tsv("data/pfam34_clans.tsv")

dat_pfam_clan <- dat_pfam_tally_trim %>%
  dplyr::bind_rows(precs) %>%
  distinct(.keep_all = T) %>%
  left_join(., clans, by = c("pfam_id1" = "pfamA_acc")) %>%
  dplyr::mutate(clan_keep = case_when(clan_acc == "NULL" ~ pfam_long,
                                      TRUE ~ paste0(clan_acc, "_", clan_id))) %>%
  dplyr::mutate(genvar = paste0(query, "_", nucleotide_acc, "_", genus_species))# %>%

length(unique(dat_pfam_clan$nucleotide_acc)) # 1922 genomes
length(unique(dat_pfam_clan$genvar)) # 2347 BGCs

write_csv(dat_pfam_clan, "data/combined_Nif11_NHLP_cleavage_motifs/1922_full_Nif11_NHLP_genome_nbs_with_clans_precursors.csv")
