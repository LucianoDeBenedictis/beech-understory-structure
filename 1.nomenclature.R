# setup -------------------------------------------------------------------

library(tidyverse)
library(readxl)


# taxonomic harmonization -------------------------------------------------

#with TNRS
speciesharm <- read_csv("data/tnrs_result.csv") |> 
  select(c("Name_submitted", "Accepted_species", "Accepted_family"))
speciesharm[146,2] <- "Cyclamen repandum" #not accepted yet
speciesharm[146,3] <- "Primulaceae"
speciesharm[147,2] <- "Stellaria nemorum" #error because of subsp.

specieslist <- read_excel("data/Species_list.xlsx")
specieslist <- specieslist |> 
  left_join(speciesharm, by = join_by("Species name"=="Name_submitted")) |> 
  rename(species = "Accepted_species") |> 
  rename(code = Code) |>
  rename(family = "Accepted_family") |> 
  filter(!is.na(species)) |> 
  select(!"Species name")

write_csv(specieslist, "data/specieslookup_harm.csv")
