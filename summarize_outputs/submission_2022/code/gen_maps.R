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
library(dplyr, lib = r_lib)
library(stringr)
library(tidyr)
library(tidyverse)
library(sf, lib=r_lib)
library(sp, lib=r_lib)
library(classInt, lib=r_lib)
library(rgdal, lib=r_lib)
library(GADMTools, lib = r_lib)
library(raster, lib=r_lib)
library(exactextractr, lib=r_lib)
library(s2, lib=r_lib)
library(ggplot2)
output_final_path <- "output_final/202110gavi-3"
setwd(output_final_path)

# read in doses administered by country
df_dose <- read.csv("intermediate_table/tp_eff_ac_byISO_medianCI.csv") %>%
  dplyr::select(ISO, threshold, confirmation_lens, admin_level, tp_cumu_median)


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
#df_map <- df_map %>%
#  mutate(doses = ifelse(is.na(tp_cumu_median) & ISO %in% cc, 0, tp_cumu_median))

# merge with shapefile of all African countries
shp <- gadm_sf_loadCountries(all_cc, level = 0)$sf
shp_dose <- shp %>% left_join(df_map, by = "ISO")

# manually add untargeted indicator
untargeted_ISO_1e_04 <- (df_map %>% filter(threshold == 1e-04 & tp_cumu_median == 0))$ISO
untargeted_ISO_0.001 <- (df_map %>% filter(threshold == 0.001 & tp_cumu_median == 0))$ISO
untargeted_ISO_2e_04 <- (df_map %>% filter(threshold == 2e-04 & tp_cumu_median == 0))$ISO


shp_dose <- shp_dose %>%
    mutate(ind_untargeted = case_when(threshold == 1e-04 & ISO %in% untargeted_ISO_1e_04 ~ 1,
                                      threshold == 0.001 & ISO %in% untargeted_ISO_0.001 ~ 1,
                                      threshold == 2e-04 & ISO %in% untargeted_ISO_2e_04 ~ 1,
                                      TRUE ~ 0))

shp_dose <- 
   shp_dose %>% 
   mutate(threshold = case_when(threshold == 0.001 ~ " \nThreshold: 10/10,000 per year \n",
                                threshold == 2e-04 ~ " \nThreshold: 2/10,000 per year \n",
                                threshold == 1e-04 ~ " \nThreshold: 1/10,000 per year \n"))

# make two figures for the manuscript
df_dose_1 <- shp_dose %>% filter(threshold == " \nThreshold: 1/10,000 per year \n" | threshold == " \nThreshold: 10/10,000 per year \n")
df_dose_2 <- shp_dose %>% filter(threshold == " \nThreshold: 2/10,000 per year \n")



map1 <-                                
   ggplot() + 
   geom_sf(data = df_dose_1, aes(fill = 2* tp_cumu_median / 1000000), lwd = 0.1) +
   scale_fill_gradient(low = "white", high = "steelblue", na.value = "grey31") +
   geom_sf(data = df_dose_1 %>% filter(ind_untargeted == 1), fill = "lightgrey", lwd = 0.1) +
   theme_void() +
   facet_wrap(~ threshold, strip.position = "top") +
   theme(legend.position = "bottom", 
         legend.text=element_text(size = 10),
         legend.title=element_text(size=14)) +
   #theme(strip.background =element_rect(fill="lightgrey"))+
   theme(strip.text = element_text(size=15, lineheight=2)) + 
   theme(plot.margin = unit(c(-20,0,-20,0), "cm")) +
   labs(fill = "Total doses administered (million)", x=NULL, y =NULL, title =NULL)

ggsave(plot = map1, filename = "/data/aazman1/hxu70/gavi-modeling/gavi_vimc_cholera/output_final/202110gavi-3/plots/map1.pdf")
saveRDS(map1,"plots/doses_maps_0.001_1e-04.rds")

# supplement map
map2 <-                                
   ggplot() + 
   geom_sf(data = df_dose_2, aes(fill = 2* tp_cumu_median / 1000000), lwd = 0.1) +
   scale_fill_gradient(low = "white", high = "steelblue", na.value = "grey31") +
   geom_sf(data = df_dose_2 %>% filter(ind_untargeted == 1), fill = "lightgrey", lwd = 0.1) +
   theme_void() +
   facet_wrap(~ threshold, strip.position = "top") +
   theme(legend.position = "bottom", 
         legend.text=element_text(size = 10),
         legend.title=element_text(size=14)) +
   #theme(strip.background =element_rect(fill="lightgrey"))+
   theme(strip.text = element_text(size=15, lineheight=2)) + 
   theme(plot.margin = unit(c(-20,0,-20,0), "cm")) +
   labs(fill = "Total doses administered (million)", x=NULL, y =NULL, title =NULL)

saveRDS(map2,"plots/doses_maps_2e-04.rds")



