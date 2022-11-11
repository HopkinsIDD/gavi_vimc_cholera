# this is the sketch board to make plots for model illustration in the paper
library(raster)
library(sf)
library(GADMTools)
library(exactextractr)
library(broom)
library(ggplot2)
library(dplyr)
setwd("/home/hanmeng/VIMC/gavi_vimc_cholera")

# read in population data and incidence data
inc_world <- raster("input_data/incidence/afro_2010-2016_lambda_5k_mean.tif")
pop_world <- raster("input_data/worldpop/ppp_2020_5km_Aggregated.tif")


# using COD as a example country
country <- "CMR"

# read in COD shapefiles 
shp0 <- gadm_sf_loadCountries(country, level = 0)$sf
shp1 <- gadm_sf_loadCountries(country, level = 1)$sf
shp2 <- gadm_sf_loadCountries(country, level = 2)$sf

# crop the pop and inc raster
cropped <- raster::crop(pop_world, shp0, snap = "out")
pop <- raster::mask(cropped, shp0, updatevalue = NA)
cropped <- raster::crop(inc_world, shp0, snap = "out")
inc <- raster::mask(cropped, shp0, updatevalue = NA)

# calculate incidence at admin1 and admin2
case <- inc * pop
case_admin1 <- exact_extract(case, shp1, 'sum')
pop_admin1 <- exact_extract(pop, shp1, "sum")
incid1 <- case_admin1 / pop_admin1

case_admin2 <- exact_extract(case, shp2, 'sum')
pop_admin2 <- exact_extract(pop, shp2, "sum")
incid2 <- case_admin2 / pop_admin2

# join incidence with shapefiles
shp1 <- cbind(shp1, incid1) 
shp2 <- cbind(shp2, incid2) 


# plot incidence (display threshold difference)
map <- 
  ggplot() + 
  geom_sf(data = shp1, aes(fill = log(incid1)), lwd = 0.1, color = "grey") +
  geom_sf(data = shp1 %>% filter(incid1 > quantile(incid1, 0.40)), aes(fill = log(incid1)), lwd = 1, color = "black") +
  scale_fill_gradient(
    low = "mistyrose", high = "red3") +
  theme_void()

ggsave(plot = map, filename = "model_illustration/plots/maps.pdf")


# plot diagnostic surveillance strategy (no, global, district)
omicron_dataset <- read.csv("input_data/confirmation_rate/parameters.csv")
cr <- rbeta(nrow(shp2), shape1 = omicron_dataset$shape1, shape2 = omicron_dataset$shape2)
shp2 <- cbind(shp2, cr)

plot_cr <- 
  ggplot() + 
  geom_sf(data = shp2, aes(fill = 0.49), color = NA) +
  scale_fill_gradient(
    low = "white", high = "steelblue") +
  theme_void()

ggsave(plot = plot_cr, filename = "model_illustration/plots/cr.pdf")

# no-estimate
plot_cr <- 
  ggplot() + 
  geom_sf(data = shp2, aes(fill = 0.49), color = NA) +
  scale_fill_gradient(
    low = "grey", high = "grey") +
  theme_void()

ggsave(plot = plot_cr, filename = "model_illustration/plots/cr.pdf")


# plot targeting scale
plot_admin <- 
  ggplot() + 
  geom_sf(data = shp1, aes(fill = 0.1)) +
  scale_fill_gradient(
    low = "white", high = "white") +
  theme_void()

ggsave(plot = plot_admin, filename = "model_illustration/plots/plot_admin.pdf")


#### added on 9/13: extracting baseline incidence rate (clinical) at all admin2 units
countries <- c("AGO", "BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "COG", "ETH", "GHA", "GIN", "GNB", "KEN", "LBR", "MDG", "MLI", "MOZ", "MRT", "MWI", "NAM", "NER", "NGA", "RWA", "SEN", "SLE", "SOM", "SSD", "TCD", "TGO", "TZA", "UGA", "ZAF", "ZMB", "ZWE")

# load all the admin2 units
shp0 <- gadm_sf_loadCountries(countries, level = 0)$sf
shp2 <- gadm_sf_loadCountries(countries, level = 2)$sf

# calculate case numbers at country and admin2 levels 
case_world <- pop_world * inc_world
case_admin0 <- exact_extract(case_world, shp0, 'sum')
pop_admin0 <- exact_extract(pop_world, shp0, "sum")
incid0 <- case_admin0 / pop_admin0

case_admin2 <- exact_extract(case_world, shp2, 'sum')
pop_admin2 <- exact_extract(pop_world, shp2, "sum")
incid2 <- case_admin2 / pop_admin2

df_incid0 <- cbind(shp0, baseline_incidence = incid0, baseline_pop = pop_admin0) %>% st_drop_geometry()
df_incid2 <- cbind(shp2, baseline_incidence = incid2, baseline_pop = pop_admin2) %>% st_drop_geometry()

write.csv(df_incid0, "df_incid0.csv")
write.csv(df_incid2, "df_incid2.csv")


#### (added on 9/19) update model illustration
# plot incidence 
map <- 
  ggplot() + 
  geom_sf(data = shp2, aes(fill = log(incid2) * 0.5), lwd = 0.1, color = "grey") +
  scale_fill_gradient(
    low = "mistyrose", high = "lightsalmon1") +
  theme_void()

ggsave(plot = map, filename = "model_illustration/plots/incid_half.pdf")


## (10/18 added summarze output)
## make plot of doses administered in SSA 
# set up
library(dplyr)
library(stringr)
library(GADMTools)
library(ggplot2)
source("packages/ocvImpact/R/utils_sum_output_median.R")
output_final_path <- "output_final/202110gavi-3"
setwd(output_final_path)

# read in doses administered by country
df_dose <- read.csv("eff_table/tp_eff_ac_byISO_median.csv")
df_dose <- read.csv("eff_table/tp_eff_ac_byISO_mean.csv")

# format a bit 
df_dose <- df_dose %>%
  mutate(doses = tp_cumu_median / 1000000 * 2) %>%
  dplyr::select(ISO, threshold, confirmation_lens, admin_level, efficiency_median, doses) %>%
  rename(efficiency = efficiency_median)



# focus only on admin2 level targeting and district-estimate confirmation
df_map <- df_dose %>% filter(admin_level == "admin2" & confirmation_lens == "district-estimate")

# load shapefiles of all African countries
all_cc <- c("DZA", "AGO", "BEN", "BWA", "BFA", "BDI", "CMR", "CPV", "CAF", "TCD", "COM", "COG", "COD",
            "CIV", "DJI", "EGY", "GNQ", "ERI", "ETH", "GAB", "GMB", "GHA", "GIN", "GNB", "KEN", "LSO", 
            "LBR", "LBY", "MDG", "MLI", "MWI", "MRT", "MAR", "MOZ", "NAM", "NER", "NGA", "REU", "RWA", 
            "STP", "SEN", "SYC", "SLE", "SOM", "ZAF", "SSD", "SDN", "SWZ", "TZA", "TGO", "TUN", "UGA",
            "ESH", "ZMB", "ZWE")
# modeled countries
cc <- c(unique(df_map$ISO), "SEN", "ZAF")

# add back those without vaccine targeting going on
df_temp <- as.data.frame(matrix(NA, ncol = 2, nrow = 3 * length(all_cc)))
colnames(df_temp) <- c("ISO", "threshold")
df_temp$ISO <- rep(all_cc, 3)
df_temp$threshold <- c(rep(0.001, length(all_cc)), rep(2e-04, length(all_cc)), rep(1e-04, length(all_cc)))
df_map <- df_temp %>% left_join(df_map, by = c("ISO", "threshold"))
# replace all the NA in the dose column (of modeled countries) with 0 (no OCV targeting going on)
df_map <- df_map %>%
  mutate(doses = ifelse(is.na(doses) & ISO %in% cc, 0, doses))

# merge with shapefile of all African countries
shp <- gadm_sf_loadCountries(all_cc, level = 0)$sf
shp_dose <- shp %>% left_join(df_map, by = "ISO")

# make the map
map <- 
   shp_dose %>% 
   #filter(threshold == 2e-04) %>%
   filter(threshold == 1e-04 | threshold == 0.001) %>%
   mutate(threshold = case_when(threshold == 0.001 ~ " \nThreshold: 10/10,000 per year \n",
                                threshold == 2e-04 ~ " \nThreshold: 2/10,000 per year \n",
                                threshold == 1e-04 ~ " \nThreshold: 1/10,000 per year \n")) %>%
   ggplot() + 
   geom_sf(aes(fill = doses), lwd = 0.1) +
   scale_fill_gradient(low = "white", high = "steelblue", na.value = "lightgrey") +
   theme_void() +
   facet_wrap(~ threshold, strip.position = "top") +
   theme(legend.position = "bottom", 
         legend.text=element_text(size = 10),
         legend.title=element_text(size=14)) +
   #theme(strip.background =element_rect(fill="lightgrey"))+
   theme(strip.text = element_text(size=15, lineheight=2)) + 
   theme(plot.margin = unit(c(-20,0,-20,0), "cm")) +
   labs(fill = "Doses administered 2022-2030 (million)", x=NULL, y =NULL, title =NULL)

saveRDS(map,"plots/dose_maps_0.001_2e-04.rds")

# supplement map
map <- 
   shp_dose %>% 
   filter(threshold == 2e-04) %>%
   #filter(threshold == 1e-04 | threshold == 0.001) %>%
   mutate(threshold = case_when(threshold == 0.001 ~ " \nThreshold: 10/10,000 per year \n",
                                threshold == 2e-04 ~ " \nThreshold: 2/10,000 per year \n",
                                threshold == 1e-04 ~ " \nThreshold: 1/10,000 per year \n")) %>%
   ggplot() + 
   geom_sf(aes(fill = doses), lwd = 0.1) +
   scale_fill_gradient(low = "white", high = "steelblue", na.value = "lightgrey") +
   theme_void() +
   facet_wrap(~ threshold, strip.position = "top") +
   theme(legend.position = "bottom", 
         legend.text=element_text(size = 10),
         legend.title=element_text(size=14)) +
   #theme(strip.background =element_rect(fill="lightgrey"))+
   theme(strip.text = element_text(size=15, lineheight=2)) + 
   theme(plot.margin = unit(c(-20,0,-20,0), "cm")) +
   labs(fill = "Doses administered 2022-2030 (million)", x=NULL, y =NULL, title =NULL)

saveRDS(map,"plots/dose_maps_no-district_2e-04.rds")



