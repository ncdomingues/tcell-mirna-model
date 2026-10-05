library(readr)

load_rules_data <- function(data_dir = file.path("data")) {
  rules_df <- read_tsv(file.path(data_dir, "rules.tsv"), col_types = cols(
    node = col_character(), value = col_integer(),
    formula = col_character(), description = col_character()
  ))
  aliases_df <- read_tsv(file.path(data_dir, "aliases.tsv"), col_types = cols(
    alias = col_character(), canonical = col_character()
  ))
  aliases <- as.list(setNames(aliases_df$canonical, aliases_df$alias))
  list(rules_df = rules_df, aliases = aliases)
}
