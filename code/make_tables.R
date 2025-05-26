
# Load required packages
library(tidyverse)
library(writexl)
library(here)

set.seed(123)

# Load R objects
f <- list.files(
    here("products", "result_files"), full.names = TRUE, pattern = "\\.rds"
)
obj <- lapply(f, readRDS)
names(obj) <- gsub("\\.rds", "", basename(f))

# Sup. Table S1: Overrepresented GO terms in fams with HSPs and HDPs ----
st1 <- obj$SEA_highly_similar_dissimilar_fams

write_tsv(
    st1, file = here("products", "tables", "Sup_Table_S1.tsv")
)


# Combine tables in a single .xlsx file ----
sheets <- list(
    S1 = st1
)
write_xlsx(
    sheets, path = here("products", "tables", "Sup_Tables.xlsx")
)
