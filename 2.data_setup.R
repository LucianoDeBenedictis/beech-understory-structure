# setup -------------------------------------------------------------------

library(tidyverse)
library(readxl)


# data import -------------------------------------------------------------

read_me <- list.files(path="data/transects", full.names = TRUE)
read_me <- read_me[order(as.numeric(gsub("[^0-9]+", "", read_me)))]

#lists of T
read_me2 <- read_me |> 
  basename() |> 
  tools::file_path_sans_ext()

#species names
specieslookup <- read_csv("data/specieslookup_harm.csv", show_col_types = F)


## transect data -----------------------------------------------------------

dat <- vector(mode = "list", length = length(read_me2))
for(i in seq_along(read_me2)){
  df <- read_csv(read_me[i], trim_ws = T, show_col_types = F) |> #header cut off already
    select(!2) |>
    rename(X=1) |> 
    mutate(across(where(is.character), toupper)) |>
    mutate(across(where(is.character), function(x) gsub("\\s+", "", x))) |>
    mutate(across(where(is.character), function(x) gsub("(PL)","", x, fixed = TRUE))) |> 
    pivot_longer(cols = !1, names_to = NULL, values_to = "code", values_drop_na = T) |> 
    filter(!code %in% c("MOSS", "TRUNK", "FERN", "D", "R","B", "ROOT", "")) |> 
    left_join(specieslookup, by = "code")
  dat[[i]] <- df
  names(dat)[i] <- read_me2[i]
}
rm(df, i, specieslookup)

#remove sites without structural and canopy data
dat[c("T15", "T16", "T17", "T32")] <- NULL

#check records without species names
NAs <- dat |> 
  lapply(function (x) filter(x, is.na(species))) |> 
  bind_rows() 
NAs |> 
  distinct(code) #all "sp."

#remove few genus-level records
dat <- dat |>  
  lapply(function (x) filter(x, !is.na(species)))

rm(NAs)

## species frequencies -----------------------------------------------------

#by plot
speciesfreqs_plot <- lapply(dat, function(x) count(x, code, species, family, sort = T)) |> 
  bind_rows(.id = "Plot")

#overall
speciesfreqs_all <- speciesfreqs_plot |>
  count(code, species, family, wt = n, sort = T) 

#clean dat
dat <- map(dat, ~select(.x, 1, 3))

# trait DB ----------------------------------------------------------------

##aboveground traits Diaz et al. 2022-----
gspff <- read_excel("data/GSPFF_enhanced/Species_mean_traits.xlsx", sheet = 1,
                    col_types = "text") |> 
  select(2, 16, 18, 20, 22, 24, 31) |> 
  rename(species=`Species name standardized against TPL`) |> 
  mutate(across(!species, as.numeric)) |> 
  mutate(presentgspff = T) |> 
  janitor::clean_names()

#synonyms
GSPFF_syno <- read_csv("data/synonyms_forGSPFF.txt", show_col_types = F) |> 
  filter(!is.na(`DB name`)) |> 
  rename(species = `DB name`)
GSPFF_syno <- setNames(GSPFF_syno$`Our name`, GSPFF_syno$species)

gspff <- gspff |>
  mutate(species = recode(species, !!!GSPFF_syno))
rm(GSPFF_syno)

##clonal from CLO-PLA---------
clopla <- read.csv("data/CLO-PLA-traits.txt", sep = "\t") |>
  mutate(Species_name = gsub("\"","", .data$Species_name)) |> 
  slice(-113) |> #allium urs. subsp., multiple match
  mutate(Species_name = word(.data$Species_name, 1, 2)) |>
  rename(species=Species_name) |> 
  mutate(presentclo = T) |> 
  select(species, presentclo, clonal, BBRsize, persistence, offspring, spread) |> 
  mutate(clonal = as.logical(clonal))

#synonyms
clopla_syno <- read_csv("data/synonyms_forCLOPLA.txt", show_col_types = F) |> 
  filter(!is.na(`DB name`)) |> 
  rename(species = `DB name`)
clopla_syno <- setNames(clopla_syno$`Our name`, clopla_syno$species)

clopla <- clopla |>
  mutate(species = recode(species, !!!clopla_syno))
rm(clopla_syno)

##join aboveground----
speciestraits <- left_join(speciesfreqs_all, gspff) |> 
  mutate(presentgspff = if_else(is.na(.data$presentgspff), F, T)) |> 
  relocate(presentgspff, .after=4)

#which species have height > 1.3 m?
tall <- speciestraits |> 
  filter(!plant_height_m < 1.3) |> 
  pull(species)

#how many records?
speciesfreqs_all |> 
  filter(species %in% tall)
rm(tall)

#remove seedlings/juveniles (above 1.3 m)
speciestraits <- speciestraits |> 
  filter(plant_height_m < 1.3 | is.na(plant_height_m))

##join clonal----
speciestraits <- left_join(speciestraits, clopla) |> 
  mutate(presentclo = if_else(is.na(.data$presentclo), F, T)) |> 
  relocate(presentclo, clonal, .after=5)

##join measured LS----
#Viola reichembachiana and Euphorbia amygdaloides have measured LS

measured <- read_csv('data/measured_LS.csv') |> 
  select(2,3) |> 
  rename(LS = 2) |> 
  group_by(Species) |> 
  summarise(LS = mean(LS))

speciestraits <- speciestraits |> 
  left_join(measured, by = join_by(species == Species)) |> 
  mutate(spread = if_else(!is.na(spread), spread, LS)) |> 
  select(!LS)

#two measured species that are actually clonal

speciestraits[speciestraits$species %in% measured$Species, "clonal"] <- T

rm(gspff,clopla, speciesfreqs_plot, speciesfreqs_all, measured)

##set species known as clonal if missing----

clonals <- c("Geranium nodosum", "Arisarum proboscideum", "Clinopodium nepeta")
speciestraits[speciestraits$species %in% clonals, "clonal"] <- T
rm(clonals)

##and species known as not clonal----

speciestraits <- speciestraits |> 
  mutate(clonal = if_else(is.na(clonal), F, clonal))

##set clonal traits at 0 if not clonal----

speciestraits <- speciestraits |> 
  mutate(across(persistence:spread, ~if_else(condition = !clonal & is.na(.x), true = 0, false = .x)))

# missingness -------------------------------------------------------------

# incomplete cases

sum(!complete.cases(speciestraits)) #26/83 overall
sum(!complete.cases(speciestraits[8:13])) #22 for gspff
sum(!complete.cases(speciestraits[14:17])) #9 for clo-pla

# missing species gspff
speciestraits |> 
  filter(presentgspff==F) |> 
  pull(species)

# missing species clo-pla
speciestraits|> 
  filter(presentclo==F) |> 
  pull(species)

## empty cells for present species----

#gspff
speciestraits |> 
  filter(if_any(8:13, ~is.na(.) & presentgspff == T)) |> 
  select(leaf_area_mm2:ssd_combined_mg_mm3) |> 
  print(n =50)

#clonals
#these are the two species with measured LS and the three species added as clonal
speciestraits |> 
  filter(if_any(14:17, ~is.na(.)) & clonal == T) |> 
  select(species, presentclo, clonal, BBRsize:spread)

#overall
speciestraits |> 
  group_by(presentgspff, presentclo) |> 
  summarise(across(where(is.numeric) & !n, ~ sum(is.na(.)))) |> 
  glimpse()

## percent completeness----

onlytraits <- speciestraits |> 
  select(where(is.numeric)&!n)

#overall
100-sum(is.na(onlytraits))/prod(dim(onlytraits))*100

#gspff
100-sum(is.na(onlytraits[1:6]))/prod(dim(onlytraits[1:6]))*100

#clo-pla
100-sum(is.na(onlytraits[7:10]))/prod(dim(onlytraits[7:10]))*100

#by trait
onlytraits |> 
  summarise(across(everything(), ~ 100-sum(is.na(.x))/n()*100)) |> 
  glimpse()

# imputation --------------------------------------------------------------
library(V.PhyloMaker2)
library(funspace)

## make phylogenetic tree----

for_phylo <- speciestraits |> 
  select(species, family) |> 
  mutate(genus = word(species), .after = 1)
for_phylo[for_phylo$genus == "Adoxa", "family"] <- "Adoxaceae"  #old family used by phylo

tree <- phylo.maker(for_phylo)
table(tree$species.list$status) #13 bind, 70 prune

#check the phylogenetic tree
plot(tree$scenario.3, type = "phylogram", cex = 0.7)

rm(for_phylo)

#add status to traits
speciestraits <- tree$species.list |> 
  select(species, status) |> 
  right_join(speciestraits) |> 
  relocate(status, .after = presentclo)

## imputation --------------------------------------------------------------

#dataframe for imputation
imputation <- speciestraits |> 
  mutate(species = str_replace_all(species, " ", "_")) |> 
  select(!2:7) |> 
  mutate(clonal = as.factor(clonal)) |> 
  column_to_rownames(var = "species")

#impute missing traits
imputed <- impute(as.data.frame(imputation), phylo = tree$scenario.3)

imputed <- imputed$imputed |> 
  rownames_to_column("species") |> 
  mutate(species = str_replace_all(species, "_", " "))

#check imputed values

if (sum(!complete.cases(imputed)) != 0){
warning("There are still missing values after imputation!")
}


#add imputed values

speciestraits <- speciestraits |> 
  select(species:status) |> 
  left_join(imputed)

rm(imputation, tree, imputed)


# save data ---------------------------------------------------------------

saveRDS(speciestraits, file = "traits_imputation.rds")
saveRDS(dat, file = "transect_data.rds")
