##############################################################################################################
#'  Data Preparations of ALKIS, Farmsubsidy & Futtermittel data for String Matching
#'  Oktober 2021
#  Author:
##############################################################################################################

##############################################################################################################
# 0 - load libraries 
##############################################################################################################

library(dplyr)
library(stringr)
library(readxl)
library(xlsx)
library(tidyverse)

##############################################################################################################
# I - set working directory and source files 
##############################################################################################################
setwd("C:/Users/IAMO/Documents/work_data/ownership_paper")

##############################################################################################################
# II - define functions & lists
##############################################################################################################

# function to categorize the farmdata into either Private (=1) or Not-Private (=2)
categorize_farmsubsidy_data <- function(input_df, df_column, name_new_owner_cat) {
  
  # 1. split category based on comma
  possible_private_persons <- as.data.frame(input_df[grepl(",", input_df[[df_column]]),])
  possible_companies <- as.data.frame(input_df[!grepl(",", input_df[[df_column]]),])
  
  # 2. split based on names
  additional_companies <- possible_private_persons %>% 
    filter(str_detect(new_names, str_c(not_private_list, collapse = "|")))
  
  # combine categorization 
  companies <- rbind(possible_companies, additional_companies )
  
  # create new col and assign a value (1 = private; 2 = not private)
  input_df[[name_new_owner_cat]][input_df[[df_column]] %in% companies[[df_column]]] <- 2
  input_df[[name_new_owner_cat]][is.na(input_df[[name_new_owner_cat]])] <- 1
  
  return(input_df)
}

# function to bring strings in the same order as ALKIS string names 
clean_up_strings_function <- function(input_df, df_column) {
  
  input_df[[df_column]] <- gsub('[.]', '', input_df[[df_column]])
  input_df[[df_column]] <- gsub('landwirtschaftl ', 'landwirtschaftliche', input_df[[df_column]])
  input_df[[df_column]] <- gsub('landw[.]', 'landwirtschaftliche', input_df[[df_column]])
  input_df[[df_column]] <- gsub('landw[.]', 'landwirtschaftliche', input_df[[df_column]])
  input_df[[df_column]] <- gsub('landwirtsch[.]', 'landwirtschaftliche', input_df[[df_column]])
  input_df[[df_column]] <- gsub('lwb', 'landwirtschaftsbetrieb', input_df[[df_column]])
  input_df[[df_column]] <- gsub('apg', 'agrarproduktivgenossenschaft', input_df[[df_column]])
  input_df[[df_column]] <- gsub('leg ', 'landwirtschaftsgesellschaft', input_df[[df_column]])
  input_df[[df_column]] <- gsub('ag ', 'agrargenossenschaft', input_df[[df_column]])
  input_df[[df_column]] <- gsub('landwirtschaftsb ', 'landwirtschaftsbetrieb', input_df[[df_column]])
  input_df[[df_column]] <- gsub('agrarges ', 'agrargesellschaft', input_df[[df_column]])
  input_df[[df_column]] <- gsub('agrargen[.]', 'agrargenossenschaft', input_df[[df_column]])
  input_df[[df_column]] <- gsub('gesellschaftmbh', 'gmbh', input_df[[df_column]])
  input_df[[df_column]] <- gsub('gemeinnuetzigegmbh', 'ggmbh', input_df[[df_column]])
  input_df[[df_column]] <- gsub('dienstl[.]', 'dienstleistungs', input_df[[df_column]])
  input_df[[df_column]] <- gsub('[.]', '', input_df[[df_column]])
  input_df[[df_column]] <- gsub('ä', 'ae', input_df[[df_column]])
  input_df[[df_column]] <- gsub('ä', 'ae', input_df[[df_column]])
  input_df[[df_column]] <- gsub('ö', 'oe', input_df[[df_column]])
  input_df[[df_column]] <- gsub('ü', 'ue', input_df[[df_column]])
  input_df[[df_column]] <- gsub('-', '', input_df[[df_column]])
  input_df[[df_column]] <- gsub('-', '', input_df[[df_column]])
  input_df[[df_column]] <- gsub('´', '', input_df[[df_column]])
  input_df[[df_column]] <- gsub('`', '', input_df[[df_column]])
  input_df[[df_column]] <- gsub("'", '', input_df[[df_column]])
  input_df[[df_column]] <- gsub("'", '', input_df[[df_column]])
  input_df[[df_column]] <- gsub("&", '', input_df[[df_column]])
  input_df[[df_column]] <- gsub(' ', '', input_df[[df_column]])
  input_df[[df_column]] <- gsub('ß', 'ss', input_df[[df_column]])
  input_df[[df_column]] <- gsub('/', '', input_df[[df_column]])
  input_df[[df_column]] <- gsub(',', '', input_df[[df_column]])
  input_df[[df_column]] <- gsub('["]', '', input_df[[df_column]])
  
  return(input_df)
}

clean_up_adresses_function <- function(input_df, df_column) {
  input_df[[df_column]] <- gsub('[.]', '', input_df[[df_column]])
  input_df[[df_column]] <- gsub('straße', 'str', input_df[[df_column]])
  input_df[[df_column]] <- gsub('ä', 'ae', input_df[[df_column]])
  input_df[[df_column]] <- gsub('ä', 'ae', input_df[[df_column]])
  input_df[[df_column]] <- gsub('ö', 'oe', input_df[[df_column]])
  input_df[[df_column]] <- gsub('ü', 'ue', input_df[[df_column]])
  input_df[[df_column]] <- gsub('-', '', input_df[[df_column]])
  input_df[[df_column]] <- gsub('´', '', input_df[[df_column]])
  input_df[[df_column]] <- gsub('`', '', input_df[[df_column]])
  input_df[[df_column]] <- gsub("'", '', input_df[[df_column]])
  input_df[[df_column]] <- gsub("'", '', input_df[[df_column]])
  input_df[[df_column]] <- gsub("&", '', input_df[[df_column]])
  input_df[[df_column]] <- gsub('ß', 'ss', input_df[[df_column]])
  input_df[[df_column]] <- gsub('/', '', input_df[[df_column]])
  input_df[[df_column]] <- gsub(' ', '', input_df[[df_column]])
  return(input_df)
}


# list with company hints 
not_private_list <- c('gmbh', 'mbh', ' ag ', ' ag$', ' gbr ', ' gbr','buergerlichen rechts', ' gbr$', 'genossenschaft', 'aktiengesellschaft', 'agrargenossenschaft', 
                      'agrargesellschaft', ' fonds ', 'stiftung*', ' eg ', ' eg$', ' e g' , ' kg ', ' kg$', ' ug ', ' ug$', 'vkg', ' lpg', ' ohg', 'beschraenkt',
                      'aktien', 'betrieb', 'baumschule', 'gbr', 'gutsverwaltung', "lwb ")


##############################################################################################################
# III - set input file paths 
##############################################################################################################
farmsubsidy_2019_bb_path <- "00_data/tables/subsidies/subs_BB_2019.csv"
farmsubsidy_2020_ge_path <- "00_data/tables/subsidies/de_2020.csv"
futtermittel_2020_bb_path <- "00_data/tables/subsidies/fodder_companies_edit_utf_8.csv"
postcode_list_ge_path <- "00_data/tables/subsidies/zuordnung_plz_ort.csv"
alkis_2020_path <- "00_data/tables/subsidies/unique_owners_and_addresses_stretched_level3.csv"

##############################################################################################################
#' 1. process farmsubsidy data 2019
##############################################################################################################

# categorize unique entries -----------------------------------------------

farmsubsidy_2019_bb <- read.csv(farmsubsidy_2019_bb_path)

# get unique recipient names and postcodes
farmsubs_2019_unique_owner <- farmsubsidy_2019_bb %>% 
  distinct(recipient_name, recipient_postcode) %>% 
  filter(! recipient_name == "")

# create new col and bring all characters to lower case
farmsubs_2019_unique_owner$new_names <- farmsubs_2019_unique_owner$recipient_name
farmsubs_2019_unique_owner$new_names <- tolower(farmsubs_2019_unique_owner$new_names)

# add 0 to plz
farmsubs_2019_unique_owner$recipient_postcode <- ifelse((nchar(farmsubs_2019_unique_owner$recipient_postcode) == 4), paste0("0", farmsubs_2019_unique_owner$recipient_postcode), as.integer(farmsubs_2019_unique_owner$recipient_postcode))

# apply categorize function 
farmsubs_2019_unique_owner <- categorize_farmsubsidy_data(farmsubs_2019_unique_owner, "recipient_name", "owner_cat")

# bring strings in shape for alkis matching -------------------------------

farmsubs_2019_unique_owner <- clean_up_strings_function(farmsubs_2019_unique_owner, "new_names")

# filter only others
farmsubs_2019_unique_owner_cat_others <- farmsubs_2019_unique_owner %>% 
  filter(owner_cat==2)

# write out dataset 
write.table(farmsubs_2019_unique_owner, "07_owner_name_cleaning/string_match_alkis_farmer_information/farmsubsidy_2019_all_categories.csv", fileEncoding = "UTF-8",row.names=FALSE , sep = ";")
write.table(farmsubs_2019_unique_owner_cat_others, "07_owner_name_cleaning/string_match_alkis_farmer_information/farmsubsidy_2019_only_cat_others.csv", fileEncoding = "UTF-8",row.names=FALSE , sep = ";")

# remove variables 
rm(farmsubs_2019_unique_owner, farmsubs_2019_unique_owner_cat_others,farmsubsidy_2019_bb)

##############################################################################################################
#' 2. process farmsubsidy data 2020
##############################################################################################################

#  subset farmsubsidy to Brandenburg via PLZ ---------------------------------

# load postcodes 
postcodes_bb <- read.csv(postcode_list_ge_path) %>% 
  filter(bundesland == "Brandenburg")

# add leading zero to postcodes with only 4 digits
postcodes_bb$plz <- ifelse((nchar(postcodes_bb$plz) == 4), paste0("0", postcodes_bb$plz), as.integer(postcodes_bb$plz))

# load data
farmsubsidy_2020_ge <- read.csv(farmsubsidy_2020_ge_path)

# delete DE- from postcodes from subsidies
farmsubsidy_2020_ge$recipient_postcode <- str_remove_all(farmsubsidy_2020_ge$recipient_postcode, "DE-")

# get only subsidies in Brandenburg based on plz list brandenburg 
farmsubsidy_2020_bb <- farmsubsidy_2020_ge %>% 
  filter(recipient_postcode %in% postcodes_bb$plz)

# delete variables
rm(farmsubsidy_2020_ge, postcodes_bb)

# categorize unique entries -----------------------------------------------

# get unique recipient names and postcodes
farmsubs_2020_unique_owner <- farmsubsidy_2020_bb %>% 
  distinct(recipient_name, recipient_postcode) %>% 
  filter(! recipient_name == "")

# create new col and bring all characters to lower case
farmsubs_2020_unique_owner$new_names <- farmsubs_2020_unique_owner$recipient_name
farmsubs_2020_unique_owner$new_names <- tolower(farmsubs_2020_unique_owner$new_names)

farmsubs_2020_unique_owner <- categorize_farmsubsidy_data(farmsubs_2020_unique_owner, "recipient_name", "owner_cat")


# bring strings in shape for alkis matching -------------------------------

farmsubs_2020_unique_owner <- clean_up_strings_function(farmsubs_2020_unique_owner, "new_names")

# filter only others
farmsubs_2020_unique_owner_cat_others <- farmsubs_2020_unique_owner %>% 
  filter(owner_cat==2)

# write out dataset 
write.table(farmsubs_2020_unique_owner, "07_owner_name_cleaning/string_match_alkis_farmer_information/farmsubsidy_2020_all_categories.csv", fileEncoding = "UTF-8",row.names=FALSE , sep = ";")
write.table(farmsubs_2020_unique_owner_cat_others, "07_owner_name_cleaning/string_match_alkis_farmer_information/farmsubsidy_2020_only_cat_others.csv", fileEncoding = "UTF-8",row.names=FALSE , sep = ";")

rm(farmsubs_2020_unique_owner, farmsubs_2020_unique_owner_cat_others,farmsubsidy_2020_bb)

##############################################################################################################
#' 3. process Futtermittel data 2020
##############################################################################################################

# categorize unique entries -----------------------------------------------

# load data 
futtermittel_2020_bb <- read.csv(futtermittel_2020_bb_path)

# separate string based on comma
sep_df <- futtermittel_2020_bb %>% 
  separate(Name_oder_Firmenbezeichnung, c("A", "B", "C", "D"), sep = ",", extra = "merge", remove = F )

# strings without comma = companies
fst_save_comp <- sep_df %>% 
  filter(is.na(B))

scd_save_comp <-  sep_df %>% filter(grepl(' ', A)  & grepl(' ', B) )

save_companies <- rbind(fst_save_comp, scd_save_comp) %>% 
  select(! c(A, B, C, D))

futtermittel_2020_bb$owner_cat[futtermittel_2020_bb$Name_oder_Firmenbezeichnung %in% save_companies$Name_oder_Firmenbezeichnung] <- 2
single_persons <- futtermittel_2020_bb %>% 
  filter(is.na(owner_cat))

single_persons$new_names <- tolower(single_persons$Name_oder_Firmenbezeichnung)

additional_companies_futtermittel <- single_persons %>% 
  filter(str_detect(new_names, str_c(not_private_list, collapse = "|"))) %>% 
  select(! c(owner_cat, new_names))

companies_futtermittel <- rbind(save_companies, additional_companies_futtermittel)
futtermittel_2020_bb$owner_cat[futtermittel_2020_bb$Name_oder_Firmenbezeichnung %in% companies_futtermittel$Name_oder_Firmenbezeichnung] <- 2
futtermittel_2020_bb$owner_cat[is.na(futtermittel_2020_bb$owner_cat)] <- 1

# delete variables
rm(single_persons, fst_save_comp, scd_save_comp, save_companies, additional_companies_futtermittel, companies_futtermittel, sep_df)


# bring adresses into the same form as alkis from clemens -----------------

# add leading zero to postcodes with only 4 digits
futtermittel_2020_bb$PLZ <- ifelse((nchar(futtermittel_2020_bb$PLZ) == 4), paste0("0", futtermittel_2020_bb$PLZ), as.integer(futtermittel_2020_bb$PLZ))
futtermittel_2020_bb <- futtermittel_2020_bb %>% 
  filter(!is.na(PLZ))

# remove characters before first comma
futtermittel_2020_bb$Ortsteil_Straße_Nr_edit <- lapply(futtermittel_2020_bb$Ortsteil_Straße_Nr, function(y) gsub(".*,","", y))
futtermittel_2020_bb$clean_address <- paste0(futtermittel_2020_bb$Ortsteil_Straße_Nr_edit, futtermittel_2020_bb$PLZ)
futtermittel_2020_bb$clean_address <- paste(futtermittel_2020_bb$clean_address, futtermittel_2020_bb$Stadt_Großgemeinde, sep = " ")

# bring adress into shape
futtermittel_2020_bb$clean_address <- tolower(futtermittel_2020_bb$clean_address)
futtermittel_2020_bb <- clean_up_adresses_function(futtermittel_2020_bb, "clean_address")

# create new col and bring all characters to lower case
futtermittel_2020_bb$new_names <- futtermittel_2020_bb$Name_oder_Firmenbezeichnung
futtermittel_2020_bb$new_names <- tolower(futtermittel_2020_bb$new_names)

# remove firstnames and add adressse for merge with alkis
futtermittel_2020_bb$new_names <- ifelse(futtermittel_2020_bb$owner_cat == 1, gsub(",.*$", "", futtermittel_2020_bb$new_names), futtermittel_2020_bb$new_names)
futtermittel_2020_bb$new_names <- ifelse(futtermittel_2020_bb$owner_cat == 1, paste0(futtermittel_2020_bb$new_names, futtermittel_2020_bb$clean_address), futtermittel_2020_bb$new_names)


# clean names and bring it into alkis shape -------------------------------

# clean names and bring it into alkis shape
futtermittel_2020_bb <- clean_up_strings_function(futtermittel_2020_bb, "new_names")

# bring into correct shape 
futtermittel_2020_bb <- apply(futtermittel_2020_bb,2,as.character)

# write out dataset 
write.table(futtermittel_2020_bb, "07_owner_name_cleaning/string_match_alkis_farmer_information/futtermittel_2020_all_categories.csv", fileEncoding = "UTF-8",row.names=FALSE, sep = ";" )

# delete variable 
rm(futtermittel_2020_bb)

##############################################################################################################
#' 4. prepare ALKIS 
##############################################################################################################

# open alkis 
alkis_2020 <- read.csv(alkis_2020_path, sep = ";")

# create new col with postcode 
alkis_2020$postcode <- sub(".*\\b(\\d{5})\\b.*", "\\1", alkis_2020$clean_address)

alkis_2020 <- alkis_2020 %>% 
  filter(! clean_address == "")

alkis_2020_unique_owner <- alkis_2020 %>% 
  distinct(owner_merge, postcode, level1 , .keep_all=T) %>% 
  select(owner_merge, postcode, owner_clean, level1)

# create new col that can be adjusted
alkis_2020_unique_owner$new_names <- alkis_2020_unique_owner$owner_merge

# bring adresses in shape 
alkis_2020_unique_owner$new_names <- tolower(alkis_2020_unique_owner$new_names)
alkis_2020_unique_owner <- clean_up_adresses_function(alkis_2020_unique_owner, "new_names")

# adjust alkis unique names 
alkis_2020_unique_owner <- clean_up_strings_function(alkis_2020_unique_owner, "new_names")

write.table(alkis_2020_unique_owner, "07_owner_name_cleaning/string_match_alkis_farmer_information/alkis_2020_unique_owner_clean.csv", fileEncoding = "UTF-8", row.names = F, sep = ";")

rm(alkis_2020, alkis_2020_unique_owner)

##############################################################################################################
# bring futtermittel and subsidy in shape and bind together
##############################################################################################################
# load data
alkis <- read.csv("07_owner_name_cleaning/string_match_alkis_farmer_information/alkis_2020_unique_owner_clean.csv",sep = ";")
futtermittel <- read.csv("07_owner_name_cleaning/string_match_alkis_farmer_information/futtermittel_2020_all_categories.csv",sep = ";")
farmsubsidy_2019 <- read.csv("07_owner_name_cleaning/string_match_alkis_farmer_information/farmsubsidy_2019_all_categories.csv",sep = ";")
farmsubsidy_2020 <- read.csv("07_owner_name_cleaning/string_match_alkis_farmer_information/farmsubsidy_2020_all_categories.csv",sep = ";")


# bring in the same shape 
futtermittel_edit <- futtermittel %>% 
  select(Name_oder_Firmenbezeichnung, PLZ, new_names, owner_cat) %>% 
  rename(original_names = Name_oder_Firmenbezeichnung,
         postcode = PLZ) %>% 
  mutate(origin = "futtermittel")

farmsubsidy_2019_edit <- farmsubsidy_2019 %>% 
  select(recipient_name,recipient_postcode, new_names, owner_cat) %>% 
  rename(original_names = recipient_name,
         postcode = recipient_postcode) %>% 
  mutate(origin = "farmsubsidy_2019")

farmsubsidy_2020_edit <- farmsubsidy_2020 %>% 
  select(recipient_name,recipient_postcode, new_names, owner_cat) %>% 
  rename(original_names = recipient_name,
         postcode = recipient_postcode) %>% 
  mutate(origin = "farmsubsidy_2020")

# bind dfs together
binded_df <- rbind(futtermittel_edit, farmsubsidy_2019_edit, farmsubsidy_2020_edit )
binded_df$postcode <- ifelse((nchar(binded_df$postcode) == 4), paste0("0", binded_df$postcode), as.integer(binded_df$postcode))

res <- binded_df[!duplicated(binded_df[c("postcode","new_names")]),] 
oth <- binded_df[duplicated(binded_df[c("postcode","new_names")]),] 

write.table(res, "07_owner_name_cleaning/string_match_alkis_farmer_information/binded_output_all_dfs.csv", fileEncoding = "UTF-8", row.names = F, sep=";")