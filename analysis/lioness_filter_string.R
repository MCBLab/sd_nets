library(dplyr)
library(STRINGdb)
library(httr)
library(jsonlite)
library(readr)
library(vroom)

# 1. Load Significant Edges
input_file <- "results/lioness_significant_edges.csv"
output_file <- "results/lioness_significant_edges_string.csv"

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

cat("Loading significant edges from", input_file, "...\n")
sig_edges <- vroom(input_file, show_col_types = FALSE)

# Ensure we have gene columns
if (!all(c("gene1", "gene2") %in% colnames(sig_edges))) {
  if ("edge" %in% colnames(sig_edges)) {
    cat("Reconstructing gene columns from edge names...\n")
    parts <- strsplit(sig_edges$edge, "_")
    sig_edges$gene1 <- sapply(parts, `[`, 1)
    sig_edges$gene2 <- sapply(parts, `[`, 2)
  } else {
    stop("Input file must contain 'gene1' and 'gene2' columns or an 'edge' column (format: GENE1_GENE2).")
  }
}

# 2. Map Genes to STRING IDs
unique_genes <- unique(c(sig_edges$gene1, sig_edges$gene2))
cat("Found", length(unique_genes), "unique genes involved in significant edges.\n")

# Initialize STRING connection
# Using settings inspired by lioness_clustering.qmd
cat("Initializing STRING db connection...\n")
string_db <- STRINGdb$new(
  version = "12.0", 
  species = 9606, 
  score_threshold = 400, 
  network_type = "full"
)

cat("Mapping genes to STRING IDs...\n")

# Using the API method from lioness_clustering.qmd for robustness with large lists
payload <- list(
  identifiers = paste(unique_genes, collapse = "\r"), # one per line
  species = 9606,
  limit = 1
)

# Call STRING API to map identifiers
res <- POST("https://string-db.org/api/json/get_string_ids",
            body = payload,
            encode = "form")

if (status_code(res) != 200) {
  stop("Failed to map genes to STRING IDs. HTTP Status: ", status_code(res))
}

ids <- fromJSON(content(res, as = "text", encoding = "UTF-8"))

if (nrow(ids) == 0) {
  stop("No genes could be mapped to STRING IDs.")
}

cat("Successfully mapped", nrow(ids), "genes to STRING identifiers.\n")

# 3. Retrieve Interactions
cat("Retrieving interactions from STRING for mapped genes...\n")
interactions <- string_db$get_interactions(ids$stringId)

cat("Found", nrow(interactions), "known interactions in STRING among these genes.\n")

# 4. Filter LIONESS Edges
# We need to match our edges (Gene A <-> Gene B) with STRING edges (ID A <-> ID B)
# First, map STRING interactions back to our original gene IDs

# Create a mapping dictionary: stringId -> queryItem (original ID)
id_map <- ids %>%
  select(stringId, queryItem) %>%
  distinct()

# Join interactions with the map to get original IDs
interactions_mapped <- interactions %>%
  inner_join(id_map, by = c("from" = "stringId")) %>%
  rename(gene1_orig = queryItem) %>%
  inner_join(id_map, by = c("to" = "stringId")) %>%
  rename(gene2_orig = queryItem)

# Helper function to create direction-independent keys for edges
# (e.g., "A_B" is same as "B_A")
make_key <- function(a, b) {
  ifelse(a < b, paste(a, b, sep = "_"), paste(b, a, sep = "_"))
}

cat("Filtering edges based on STRING validation...\n")

# Add keys to both dataframes
sig_edges <- sig_edges %>%
  mutate(key = make_key(gene1, gene2))

interactions_mapped <- interactions_mapped %>%
  mutate(key = make_key(gene1_orig, gene2_orig))

# Filter: Keep LIONESS edges that have a matching key in STRING interactions
filtered_edges <- sig_edges %>%
  filter(key %in% interactions_mapped$key) %>%
  select(-key) # Remove the helper key column

# 5. Save Results
cat("Filtering complete.\n")
cat("Original significant edges:", nrow(sig_edges), "\n")
cat("Edges validated by STRING:", nrow(filtered_edges), "\n")

write_csv(filtered_edges, output_file)
cat("Results saved to", output_file, "\n")
