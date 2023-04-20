##### Summarize the outbreak data.frame to generate tables for the supplementary materials 

library(dplyr)

outbreak_df <- readr::read_csv("/Users/kaiyuezou/Desktop/vimc_write_up/outbreak_df_full.csv")
outbreak_df <- outbreak_df[!is.na(outbreak_df$admin0), ]
outbreak_df$admin0 <- countrycode::countrycode(outbreak_df$admin0, "iso3c", "country.name")
all_country <- unique(outbreak_df$admin0)
noquote(all_country)
table(outbreak_df$admin0)

sum(!is.na(outbreak_df[outbreak_df$admin0 == "ago", ]$attack_rate))
sum(!is.na(outbreak_df[outbreak_df$admin0 == "moz", ]$attack_rate))
sum(!is.na(outbreak_df[outbreak_df$admin0 == "nam", ]$attack_rate))
#The countries that only have one line of data are: ago, moz, and nam. However, they do have valid data

## The countries that have valid admin 2 level data 
admin2_country <- unique(outbreak_df[ !is.na(outbreak_df$admin2) & !is.na(outbreak_df$attack_rate), ]$admin0)
noquote(admin2_country)

## Generate a summarized table
new_admin1_only <- outbreak_df %>% 
  filter(!is.na(attack_rate) & spatial_scale == "admin1") %>% 
  group_by(admin0, admin1) %>% 
  summarise(num_outbreak = length(unique(outbreak_number)))
names(new_admin1_only) <- c("Country Name", "Admin 1 District", "Number of Outbreaks")
  
new_admin2_only <- outbreak_df %>% 
  filter(!is.na(attack_rate) & spatial_scale == "admin2") %>% 
  group_by(admin0, admin2) %>% 
  summarise(num_outbreak = length(unique(outbreak_number)))
names(new_admin2_only) <- c("Country Name", "Admin 2 District", "Number of Outbreaks")

new_admin3_only <- outbreak_df %>% 
  filter(!is.na(attack_rate) & spatial_scale == "admin3") %>% 
  group_by(admin0) %>% 
  summarise(num_outbreak = length(unique(outbreak_number)))
names(new_admin3_only) <- c("Country Name", "Number of Outbreaks")

readr::write_csv(new_admin1_only, "/Users/kaiyuezou/Desktop/vimc_write_up/outbreak_admin1_table.csv")
readr::write_csv(new_admin2_only, "/Users/kaiyuezou/Desktop/vimc_write_up/outbreak_admin2_table.csv")
readr::write_csv(new_admin3_only, "/Users/kaiyuezou/Desktop/vimc_write_up/outbreak_admin3_table.csv")
