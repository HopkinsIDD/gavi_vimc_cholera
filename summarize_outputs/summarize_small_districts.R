## 11/10
## to look into the small areas issue raised by fasterize() function and get a summary table of these districts
library(dplyr)
library(GADMTools)
library(sf)
library(mapview)
library(raster)
library(exactextractr)
df_prob <- read.csv("data/df_small_targets.csv")
pop <- raster("data/ppp_2020_5km_Aggregated.tif")
inc <- raster("data/afro_2010-2016_lambda_5k_mean.tif")
case <- pop * inc
ISOs <- unique(df_prob$ISO)
districts <- unique(df_prob$NAME_2)

# get area of the problematic district
shp2 <- gadm_sf_loadCountries(ISOs, level = 2)$sf
shp2 <- shp2 %>% filter(NAME_2 %in% districts)
shp2$area_km2 <- round(as.numeric(st_area(shp2) / 1000000),1)

# get pop(2020), incidence rate (from baseline raster of mean suspected incidence rate) of these districts
pop_admin2 <- exact_extract(pop, shp2, "sum")
shp2$pop <- pop_admin2
case_admin2 <- exact_extract(case, shp2, "sum")
shp2$case <- case_admin2
shp2$incid <- shp2$case / shp2$pop

df_demo <- st_drop_geometry(shp2) %>%
  dplyr::select(ISO, NAME_1, NAME_2, area_km2, pop, case, incid)

# then look at the scenarios these districts were targeted
# (first, note that all these were admin2 targeting)
df_confirmation <- df_prob %>%
  group_by(ISO, NAME_1, NAME_2, confirmation_lens) %>%
  count() %>%
  pivot_wider(names_from = "confirmation_lens", values_from = "n") %>%
  mutate(`global-estimate` = case_when(!is.na(`global-estimate`) ~ "X"),
         `no-estimate` = case_when(!is.na(`no-estimate`) ~ "X"),
         `district-estimate` = case_when(!is.na(`district-estimate`) ~ "X"))
  
df_threshold <- df_prob %>%
  group_by(ISO, NAME_1, NAME_2, threshold) %>%
  count() %>%
  pivot_wider(names_from = "threshold", values_from = "n") %>%
  mutate(`1e-04` = case_when(!is.na(`1e-04`) ~ "X"),
         `2e-04` = case_when(!is.na(`2e-04`) ~ "X"),
         `0.001` = case_when(!is.na(`0.001`) ~ "X")) 
  # rename(`10/10,000 per year` = `0.001`,
  #        `2/10,000 per year` = `2e-04`,
  #        `1/10,000 per year` = `1e-04`)

df_small_targets_summary <- df_confirmation %>%
  left_join(df_threshold, by = c("ISO", "NAME_1", "NAME_2")) %>%
  left_join(df_demo, by = c("ISO", "NAME_1", "NAME_2")) %>%
  dplyr::select(1:3, 11, 13, 12, 10, 4:9)

write.csv(df_small_targets_summary, "data/df_small_targets_summary.csv", row.names = FALSE)









  


