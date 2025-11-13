#!/usr/bin/env Rscript
# Extended reference workflow for adding tobacco to trnL database
# Based on code/Extend reference.Rmd with taxonomy building from Taxa names.Rmd
# Builds both trnLGH.fasta and trnLGH_taxonomy.fasta

library(Biostrings)
library(here)
library(rentrez)
library(tidyverse)
library(metacoder)

# Set working directory
setwd('/home/articulatus/git_repos/food-dbs')

# Source helper functions
source(here('code', 'functions', 'find_primer_pair.R'))
source(here('code', 'functions', 'query_ncbi.R'))

cat('Starting extended reference workflow for tobacco...\n')

# ===== Read in current reference =====
cat('\nReading current trnL reference...\n')
current <-
    here('data',
         'outputs',
         'dada2-compatible',
         'trnL',
         'trnLGH.fasta') %>%
    readDNAStringSet()

cat('Current sequences:', length(current), '\n')

# Organize as dataframe
current.df <-
    data.frame(
        name = names(current),
        seq = as.character(current)
    )

row.names(current.df) <- NULL

# Separate accession from species name
current.df <-
    separate(current.df,
             name,
             sep = ' ',
             into = c('accession', 'species'),
             extra = 'merge',
             remove = FALSE) %>%
    mutate(accession = gsub('^>', '', accession))  # Strip leading >

# Combine species and sequence to make unique identifier
current.df <- mutate(current.df,
                    ID = paste(species, seq))

# ===== QC on current reference =====
cat('\nQC on current reference:\n')
cat('Duplicates:', any(duplicated(current.df$ID)), '\n')
cat('Degenerate nucleotides:', length(grep('[AGCT]*[^AGCT]+', current.df$seq)), '\n')

# ===== Pull sequence data =====
cat('\nQuerying NCBI for additions...\n')

fs <-
    here('data',
         'processed',
         'sql-compatible') %>%
    list.files(pattern = 'additions.csv',
               full.names = TRUE)

cat('Found addition files:', basename(fs), '\n')

add <-
    sapply(fs, read_csv, simplify = FALSE) %>%
    bind_rows()

cat('Addition entries to process:', nrow(add), '\n')
cat('Columns:', colnames(add), '\n')

# ===== Search for sequences =====
cat('\nSearching NCBI for sequences...\n')

# Extract organisms from accession_num column (format: "trnL Nicotiana tabacum")
organisms <- add[[1]] %>%
    str_replace('trnL ', '') %>%
    unique()

cat('Organisms:', paste(organisms, collapse=', '), '\n')

# Query using query_ncbi function
records <- query_ncbi(marker = 'trnL',
                     organisms = organisms)

cat('Found:', length(records), 'sequences\n')

if (length(records) == 0) {
    cat('No sequences found!\n')
    quit(status = 1)
}

# ===== Trim primers =====
cat('\nTrimming to trnL primers...\n')

trnLG <- DNAString('GGGCAATCCTGAGCCAA')
trnLH <- DNAString('CCATTGAGTCTCTGCACCTATC')

records.trim <- find_primer_pair(records,
                                trnLG,
                                trnLH)

cat('After primer trimming:', length(records.trim), 'sequences\n')

if (length(records.trim) == 0) {
    cat('No sequences retained after trimming!\n')
    quit(status = 1)
}

# ===== Prepare for merging =====
cat('\nPreparing sequences for merging...\n')

# Format as dataframe
add.df <-
    data.frame(name = names(records.trim),
               seq = as.character(records.trim)) %>%
    separate(name,
             into = c('accession', 'full_name'),
             extra = 'merge',
             sep = '\\s',
             remove = FALSE) %>%
    mutate(accession = gsub('^>', '', accession))  # Strip leading >

# Extract proper species name (Genus species [subsp./var./f. infraspecific])
add.df <- add.df %>%
    mutate(species = sapply(full_name, function(x) {
        parts <- str_split(x, '\\s+')[[1]]

        # Always take first two words (Genus species)
        if (length(parts) >= 2) {
            result <- paste(parts[1], parts[2])

            # Check if third word is a rank indicator
            if (length(parts) >= 4 && parts[3] %in% c('subsp.', 'var.', 'f.')) {
                result <- paste(result, parts[3], parts[4])
            }

            return(result)
        } else {
            return(x)  # Fallback if parsing fails
        }
    })) %>%
    mutate(ID = paste(species, seq))

# ===== QC on new sequences =====
cat('\nQC on new sequences:\n')

# Check duplicates within
dups_within <- any(duplicated(add.df$ID))
cat('Duplicates within new sequences:', dups_within, '\n')

if (dups_within) {
    before <- nrow(add.df)
    add.df <- add.df[!duplicated(add.df$ID), ]
    cat('  Removed:', before - nrow(add.df), 'duplicates\n')
}

# Check duplicates with existing reference
dups_existing <- any(add.df$ID %in% current.df$ID)
cat('Already in reference:', dups_existing, '\n')

if (dups_existing) {
    before <- nrow(add.df)
    add.df <- add.df[!(add.df$ID %in% current.df$ID), ]
    cat('  Removed:', before - nrow(add.df), 'duplicates\n')
}

# Check for degenerate nucleotides
degen <- any(grepl(pattern = '[AGCT]*[^AGCT]+', add.df$seq))
cat('Degenerate nucleotides:', degen, '\n')

if (degen) {
    bad_seqs <- add.df %>% filter(grepl(pattern = '[AGCT]*[^AGCT]+', seq))
    cat('  Removing:\n')
    for (i in 1:nrow(bad_seqs)) {
        cat('   ', bad_seqs$accession[i], bad_seqs$species[i], '\n')
    }
    add.df <- filter(add.df, !grepl(pattern = '[AGCT]*[^AGCT]+', seq))
}

cat('New sequences after QC:', nrow(add.df), '\n')

if (nrow(add.df) == 0) {
    cat('No sequences passed QC!\n')
    quit(status = 1)
}

# ===== Merge references =====
cat('\nMerging with existing reference...\n')

# Select needed columns
current.simple <- current.df %>% select(accession, species, seq)
add.simple <- add.df %>% select(accession, species, seq)

# Collapse duplicates: for same seq, keep first accession (first alphabetically by species)
add.simple <- add.simple %>%
    arrange(species, accession) %>%
    distinct(seq, .keep_all = TRUE)

cat('New sequences after deduplication:', nrow(add.simple), '\n')

update.df <- bind_rows(current.simple, add.simple)

cat('Total sequences after merge:', nrow(update.df), '\n')

# Sort alphabetically
update.df <- arrange(update.df,
                    species,
                    accession)

# ===== Save full updated reference =====
cat('\nSaving updated full reference...\n')

# Convert to DNAStringSet
trnL <- update.df$seq
names(trnL) <- paste(update.df$accession, update.df$species)
trnL <- DNAStringSet(trnL)

output_path <- here('data',
                   'outputs',
                   'dada2-compatible',
                   'trnL',
                   'trnLGH.fasta')

writeXStringSet(trnL, output_path)

cat('Saved full reference to:', output_path, '\n')

# ===== Build taxonomy version using metacoder =====
cat('\nBuilding taxonomy version...\n')

# Lookup taxonomy for deduplicated new sequences
accs <- gsub('^>', '', add.simple$accession)
cat('Looking up taxonomy for', length(accs), 'unique accessions...\n')

taxmap <- metacoder::lookup_tax_data(accs, type = 'seq_id')

# Get taxonomy table using metacoder
cat('Extracting taxonomy table...\n')

taxonomy <- metacoder::taxonomy_table(taxmap,
                                     use_ranks = c('domain',
                                                  'phylum',
                                                  'class',
                                                  'order',
                                                  'family',
                                                  'genus',
                                                  'species',
                                                  'subspecies',
                                                  'varietas',
                                                  'forma'),
                                     add_id_col = TRUE)

# Join with query data
taxtab <- data.frame(acc = taxmap$data$query_data,
                     taxon_id = names(taxmap$data$query_data))

taxtab <- left_join(taxtab, taxonomy, by = 'taxon_id')

# Build taxonomy strings with unite()
taxtab <- taxtab %>%
    unite(col = 'taxonomy_string',
          domain:forma,
          sep = ';',
          remove = FALSE)

cat('Taxonomy strings built\n')

# ===== Append to existing full taxonomy file =====
cat('\nAppending to full taxonomy file...\n')

tax_file <- here('data',
                'outputs',
                'dada2-compatible',
                'trnL',
                'trnLGH_taxonomy.fasta')

# Build new taxonomy entries and append
for (i in 1:nrow(add.simple)) {
    acc <- gsub('^>', '', add.simple$accession[i])
    seq <- add.simple$seq[i]

    tax_row <- taxtab %>% filter(acc == !!acc)
    if (nrow(tax_row) > 0) {
        tax_string <- tax_row$taxonomy_string[1]

        # Append directly to file
        cat('>', tax_string, '\n', sep = '', file = tax_file, append = TRUE)
        cat(seq, '\n', sep = '', file = tax_file, append = TRUE)

        cat('  Added:', tax_string, '\n')
    }
}

cat('Saved full taxonomy version to:', tax_file, '\n')

# ===== Save tobacco-only references (optional) =====
# Uncomment to generate tobacco-only reference files for review/distribution
# cat('\nSaving tobacco-only references...\n')
#
# # Tobacco-only main reference
# trnL_tobacco <- add.simple$seq
# names(trnL_tobacco) <- paste(add.simple$accession, add.simple$species)
# trnL_tobacco <- DNAStringSet(trnL_tobacco)
#
# tobacco_output_path <- here('data',
#                             'outputs',
#                             'dada2-compatible',
#                             'trnL',
#                             'trnLGH_tobacco.fasta')
#
# writeXStringSet(trnL_tobacco, tobacco_output_path)
# cat('Saved tobacco-only to:', tobacco_output_path, '\n')
#
# # Tobacco-only taxonomy reference
# tobacco_tax_file <- here('data',
#                         'outputs',
#                         'dada2-compatible',
#                         'trnL',
#                         'trnLGH_tobacco_taxonomy.fasta')
#
# for (i in 1:nrow(add.simple)) {
#     acc <- gsub('^>', '', add.simple$accession[i])
#     seq <- add.simple$seq[i]
#
#     tax_row <- taxtab %>% filter(acc == !!acc)
#     if (nrow(tax_row) > 0) {
#         tax_string <- tax_row$taxonomy_string[1]
#
#         cat('>', tax_string, '\n', sep = '', file = tobacco_tax_file, append = TRUE)
#         cat(seq, '\n', sep = '', file = tobacco_tax_file, append = TRUE)
#     }
# }
#
# cat('Saved tobacco-only taxonomy to:', tobacco_tax_file, '\n')

# ===== Summary =====
cat('\n================================\n')
cat('Workflow complete!\n')
cat('Summary:\n')
cat('  Previous sequences:', nrow(current.simple), '\n')
cat('  New sequences added:', nrow(add.simple), '\n')
cat('  Total sequences (full reference):', length(trnL), '\n')
cat('\nOutput files:\n')
cat('  Full reference:', output_path, '\n')
cat('  Full taxonomy:', tax_file, '\n')
cat('================================\n')
