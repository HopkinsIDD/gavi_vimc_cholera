### This is the script that will generate the diagnostic report for the stochastic output

### Setup
library(dplyr)
library(readr)
library(ggplot2)
setting_num <- "all"

### Read in the data
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
camp <- read_csv(paste0(current_working_dir, "/", list.files(current_working_dir, pattern = "^stochastic-burden.*campaign-default.csv$")))
nova <- read_csv(paste0(current_working_dir, "/", list.files(current_working_dir, pattern = "^stochastic-burden.*no-vaccination.csv$")))

### Create year only data (aggregated across ages)
camp_clean <- camp %>% 
  group_by(country, run_id, year) %>% 
  summarise(cases_total = sum(cases, na.rm = TRUE), 
            deaths_total = sum(deaths, na.rm = TRUE), 
            dalys_total = sum(dalys, na.rm = TRUE))
camp_clean$run_id <- as.character(camp_clean$run_id)
  
camp_clean_mean <- camp_clean %>% 
  group_by(country, year) %>% 
  summarise(run_id = "mean",
            cases_total_tmp = mean(cases_total, na.rm = TRUE), 
            deaths_total_tmp = mean(deaths_total, na.rm = TRUE), 
            dalys_total_tmp = mean(dalys_total, na.rm = TRUE)) %>% 
  rename(cases_total = cases_total_tmp,
         deaths_total = deaths_total_tmp,
         dalys_total = dalys_total_tmp) %>% 
  select(names(camp_clean))

camp_clean_min <- camp_clean %>% 
  group_by(country, year) %>% 
  summarise(run_id = "q025",
            cases_total_tmp = quantile(cases_total, 0.025, na.rm = TRUE), 
            deaths_total_tmp = quantile(deaths_total, 0.025, na.rm = TRUE), 
            dalys_total_tmp = quantile(dalys_total, 0.025, na.rm = TRUE)) %>% 
  rename(cases_total = cases_total_tmp,
         deaths_total = deaths_total_tmp,
         dalys_total = dalys_total_tmp) %>% 
  select(names(camp_clean))

camp_clean_max <- camp_clean %>% 
  group_by(country, year) %>% 
  summarise(run_id = "q975",
            cases_total_tmp = quantile(cases_total, 0.975, na.rm = TRUE), 
            deaths_total_tmp = quantile(deaths_total, 0.975, na.rm = TRUE), 
            dalys_total_tmp = quantile(dalys_total, 0.975, na.rm = TRUE)) %>% 
  rename(cases_total = cases_total_tmp,
         deaths_total = deaths_total_tmp,
         dalys_total = dalys_total_tmp) %>% 
  select(names(camp_clean))

camp_clean <- rbind(camp_clean, camp_clean_mean, camp_clean_min, camp_clean_max)

### Plot and save
#camp_clean <- camp_clean[camp_clean$country == "BEN", ] #tmp
#pdf(paste0(current_working_dir, "/camp_stochastic_setting", setting_num, ".pdf"), width = 9*3, height = 9*3)
pdf(paste0(current_working_dir, "/camp_stochastic_setting_", setting_num, ".pdf"), width = 6*3, height = 30*3)
camp_clean %>%
  ggplot( aes(x=year, y=(cases_total), group=run_id)) +
  geom_line(color = "lightblue") + 
  geom_line(data = camp_clean[camp_clean$run_id == "mean", ], 
            aes(x=year, y=(cases_total), group=run_id), 
            color = "red", size = 2) +
  #geom_line(data = camp_clean[camp_clean$run_id == "1", ], 
  #          aes(x=year, y=(cases_total), group=run_id), 
  #          color = "pink", size = 2) +
  geom_line(data = camp_clean[camp_clean$run_id == "q025", ], 
            aes(x=year, y=(cases_total), group=run_id), 
            color = "black", size = 1) +
  geom_line(data = camp_clean[camp_clean$run_id == "q975", ], 
            aes(x=year, y=(cases_total), group=run_id), 
            color = "black", size = 1) +
  facet_wrap( ~ country, ncol = 3, scales = "free") +
  theme_minimal()
camp_clean %>%
  ggplot( aes(x=year, y=deaths_total, group=run_id)) +
  geom_line(color = "lightblue") + 
  geom_line(data = camp_clean[camp_clean$run_id == "mean", ], 
            aes(x=year, y=deaths_total, group=run_id), 
            color = "red", size = 2) +
  geom_line(data = camp_clean[camp_clean$run_id == "q025", ], 
            aes(x=year, y=deaths_total, group=run_id), 
            color = "black", size = 1) +
  geom_line(data = camp_clean[camp_clean$run_id == "q975", ], 
            aes(x=year, y=deaths_total, group=run_id), 
            color = "black", size = 1) +
  facet_wrap( ~ country, ncol = 3, scales = "free") +
  theme_minimal()
camp_clean %>%
  ggplot( aes(x=year, y=dalys_total, group=run_id)) +
  geom_line(color = "lightblue") + 
  geom_line(data = camp_clean[camp_clean$run_id == "mean", ], 
            aes(x=year, y=dalys_total, group=run_id), 
            color = "red", size = 2) +
  geom_line(data = camp_clean[camp_clean$run_id == "q025", ], 
            aes(x=year, y=dalys_total, group=run_id), 
            color = "black", size = 1) +
  geom_line(data = camp_clean[camp_clean$run_id == "q975", ], 
            aes(x=year, y=dalys_total, group=run_id), 
            color = "black", size = 1) +
  facet_wrap( ~ country, ncol = 3, scales = "free") +
  theme_minimal()
dev.off()





### Do it again for no vaccination scenario (still borrowing "camp")
camp_clean <- nova %>% 
  group_by(country, run_id, year) %>% 
  summarise(cases_total = sum(cases, na.rm = TRUE), 
            deaths_total = sum(deaths, na.rm = TRUE), 
            dalys_total = sum(dalys, na.rm = TRUE))
camp_clean$run_id <- as.character(camp_clean$run_id)

camp_clean_mean <- camp_clean %>% 
  group_by(country, year) %>% 
  summarise(run_id = "mean",
            cases_total_tmp = mean(cases_total, na.rm = TRUE), 
            deaths_total_tmp = mean(deaths_total, na.rm = TRUE), 
            dalys_total_tmp = mean(dalys_total, na.rm = TRUE)) %>% 
  rename(cases_total = cases_total_tmp,
         deaths_total = deaths_total_tmp,
         dalys_total = dalys_total_tmp) %>% 
  select(names(camp_clean))

camp_clean_min <- camp_clean %>% 
  group_by(country, year) %>% 
  summarise(run_id = "q025",
            cases_total_tmp = quantile(cases_total, 0.025, na.rm = TRUE), 
            deaths_total_tmp = quantile(deaths_total, 0.025, na.rm = TRUE), 
            dalys_total_tmp = quantile(dalys_total, 0.025, na.rm = TRUE)) %>% 
  rename(cases_total = cases_total_tmp,
         deaths_total = deaths_total_tmp,
         dalys_total = dalys_total_tmp) %>% 
  select(names(camp_clean))

camp_clean_max <- camp_clean %>% 
  group_by(country, year) %>% 
  summarise(run_id = "q975",
            cases_total_tmp = quantile(cases_total, 0.975, na.rm = TRUE), 
            deaths_total_tmp = quantile(deaths_total, 0.975, na.rm = TRUE), 
            dalys_total_tmp = quantile(dalys_total, 0.975, na.rm = TRUE)) %>% 
  rename(cases_total = cases_total_tmp,
         deaths_total = deaths_total_tmp,
         dalys_total = dalys_total_tmp) %>% 
  select(names(camp_clean))

camp_clean <- rbind(camp_clean, camp_clean_mean, camp_clean_min, camp_clean_max)

### Plot and save
#camp_clean <- camp_clean[camp_clean$country == "BEN", ]
#pdf(paste0(current_working_dir, "/novac_stochastic_setting", setting_num, ".pdf"), width = 9*3, height = 9*3)
pdf(paste0(current_working_dir, "/novac_stochastic_setting_", setting_num, ".pdf"), width = 6*3, height = 30*3)
camp_clean %>%
  ggplot( aes(x=year, y=(cases_total), group=run_id)) +
  geom_line(color = "lightblue") + 
  geom_line(data = camp_clean[camp_clean$run_id == "mean", ], 
            aes(x=year, y=(cases_total), group=run_id), 
            color = "red", size = 2) +
  geom_line(data = camp_clean[camp_clean$run_id == "q025", ], 
            aes(x=year, y=(cases_total), group=run_id), 
            color = "black", size = 1) +
  geom_line(data = camp_clean[camp_clean$run_id == "q975", ], 
            aes(x=year, y=(cases_total), group=run_id), 
            color = "black", size = 1) +
  facet_wrap( ~ country, ncol = 3, scales = "free") +
  theme_minimal()
camp_clean %>%
  ggplot( aes(x=year, y=deaths_total, group=run_id)) +
  geom_line(color = "lightblue") + 
  geom_line(data = camp_clean[camp_clean$run_id == "mean", ], 
            aes(x=year, y=deaths_total, group=run_id), 
            color = "red", size = 2) +
  geom_line(data = camp_clean[camp_clean$run_id == "q025", ], 
            aes(x=year, y=deaths_total, group=run_id), 
            color = "black", size = 1) +
  geom_line(data = camp_clean[camp_clean$run_id == "q975", ], 
            aes(x=year, y=deaths_total, group=run_id), 
            color = "black", size = 1) +
  facet_wrap( ~ country, ncol = 3, scales = "free") +
  theme_minimal()
camp_clean %>%
  ggplot( aes(x=year, y=dalys_total, group=run_id)) +
  geom_line(color = "lightblue") + 
  geom_line(data = camp_clean[camp_clean$run_id == "mean", ], 
            aes(x=year, y=dalys_total, group=run_id), 
            color = "red", size = 2) +
  geom_line(data = camp_clean[camp_clean$run_id == "q025", ], 
            aes(x=year, y=dalys_total, group=run_id), 
            color = "black", size = 1) +
  geom_line(data = camp_clean[camp_clean$run_id == "q975", ], 
            aes(x=year, y=dalys_total, group=run_id), 
            color = "black", size = 1) +
  facet_wrap( ~ country, ncol = 3, scales = "free") +
  theme_minimal()
dev.off()


