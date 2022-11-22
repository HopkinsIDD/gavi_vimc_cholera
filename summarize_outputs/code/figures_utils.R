# scripts of functions used in figures.R

## get_tp_eff_median -----------------------------------------------------
# Function to summarize median targeted pop and median efficiency of each country 
get_tp_eff_median <- function(rc){
  
  df_tp_cumu <- rc %>% 
    filter(is_target == 1) %>%
    group_by(ISO, year, run_id, threshold, confirmation_lens, admin_level) %>%
    summarize(target_pop = sum(actual_fvp)) %>%
    group_by(ISO, run_id, threshold, confirmation_lens, admin_level) %>%
    mutate(target_pop_cumu = cumsum(target_pop)) %>% # get sum tp of all years
    dplyr::select(-target_pop) %>%
    group_by(ISO, year, threshold, confirmation_lens, admin_level) %>%
    summarize(tp_cumu_median = median(target_pop_cumu, na.rm = TRUE)) %>% # only keep median for each ISO
    dplyr::select(year, threshold, confirmation_lens, tp_cumu_median, admin_level) %>%
    group_by(ISO, threshold, confirmation_lens, admin_level) %>%
    summarize(tp_cumu_median = max(tp_cumu_median))
  
  # averted true cases per 1000 fvp 
  df_eff <- rc %>% 
    group_by(ISO, year, run_id, threshold, confirmation_lens, admin_level) %>%
    summarize(true_ac = sum(true_averted_cases)) %>%
    group_by(ISO, run_id, threshold, confirmation_lens, admin_level) %>%
    mutate(true_ac_cumu = cumsum(true_ac)) %>%
    dplyr::select(-true_ac) %>%
    # filter(year == max_year)
    group_by(ISO, run_id, threshold, confirmation_lens, admin_level) %>%
    summarize(true_ac_cumu = max(true_ac_cumu)) %>%
    right_join(
      rc %>% 
        group_by(ISO, year, run_id, threshold, confirmation_lens, admin_level) %>%
        summarize(fvp = sum(actual_fvp)) %>%
        group_by(ISO, run_id, threshold, confirmation_lens, admin_level) %>%
        mutate(fvp_cumu = cumsum(fvp)) %>%
        dplyr::select(-fvp) %>%
        group_by(ISO, run_id, threshold, confirmation_lens, admin_level) %>%
        summarize(fvp_cumu = max(fvp_cumu)) ,
      by = c("ISO", "run_id", "threshold", "confirmation_lens", "admin_level")
    )%>%
    # filter out all those fvp == 0 rows to get rid of Inf
    filter(fvp_cumu != 0) %>%
    # calculate efficiency
    mutate(efficiency = true_ac_cumu / fvp_cumu * 1000) %>%
    group_by(ISO, threshold, confirmation_lens, admin_level) %>%
    summarize(efficiency_median = median(efficiency, na.rm = TRUE)) %>% # keep only median for each country
    dplyr::select(ISO, threshold, confirmation_lens,admin_level, efficiency_median)
  
  
  # join two tables 
  df_result <- df_tp_cumu %>%
    left_join(df_eff, by = c("ISO", "threshold", "confirmation_lens", "admin_level")) %>%
    # convert unit of tp to "million"
    mutate(doses = tp_cumu_median/1000000 * 2) %>%
    rename(efficiency = efficiency_median) 
  
  return(df_result)
}

## make_boxplot -------------------------------------------------
# function to make boxplots 
make_boxplot <- function(df_sum, 
                         item_plotted, # item_plotted = "annual_doses"/"all_years_doses"/"efficiency"
                         palette = "Set1"){ 
  
  # rename
  df_sum <- df_sum %>%
    mutate(confirmation_lens = case_when(confirmation_lens == "no-estimate" ~ "Clinical case defintion",
                                         confirmation_lens == "global-estimate" ~ "National labs",
                                         confirmation_lens == "district-estimate" ~ "District labs")) %>%
    mutate(admin_level = case_when(admin_level == "admin1" ~ "OCV targeting at admin 1 level",
                                   admin_level == "admin2" ~ "OCV targeting at admin 2 level")) %>%
    mutate(threshold = case_when(threshold == 0.001 ~ "Threshold: 10/10,000",
                                 threshold == 2e-04 ~ "Threshold: 2/10,000",
                                 threshold == 1e-04 ~ "Threshold: 1/10,000")) %>%
    mutate(threshold = factor(threshold, levels = c("Threshold: 10/10,000", "Threshold: 2/10,000", "Threshold: 1/10,000")))
  
  # set color options 
  num_colors <- length(unique(df_sum$confirmation_lens))
  if(item_plotted == "all_years_doses"){ylab_name = "Doses of OCV administered from 2022-2030 (million)"}
  if(item_plotted == "annual_doses"){ylab_name = "Annual doses administered (million)"}
  if(item_plotted == "efficiency"){ylab_name = "OCV efficiency (averted true cases / 1,000 fvp)"}
  
  # get the column to be plotted
  if(item_plotted == "all_years_doses"){df_sum <- df_sum %>% rename(all_year_doses = doses)}
  if(item_plotted == "annual_doses"){df_sum <- df_sum %>% mutate(annual_doses = doses / 14)}
  
  plt <- 
    ggplot(df_sum, aes(x = as.factor(confirmation_lens), y = get(item_plotted), fill = as.factor(confirmation_lens))) +
    geom_boxplot(alpha = 0.8, width = 0.5, outlier.shape = NA) +
    # without outliers 
    geom_jitter(aes(color = as.factor(confirmation_lens)),
                alpha = 0.8, width = 0.1, height = 0.01) +
    theme_bw() +
    facet_grid(admin_level ~ threshold, scales = "free", space = "free") +
    # scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
    scale_color_manual(values=brewer.pal(n=num_colors, name = palette)) +
    scale_fill_brewer(palette = palette) +
    theme(axis.text = element_text(size = 10),
          legend.position = 'none') +
    xlab('Cholera confirmation') +
    ylab(ylab_name)
    
  return(plt)
}


# plot_threshold ---------------------------------------------------
# which incidence threshold had the highest mean efficiency (for admin2 targeting + district-estimate only) for a given country
# problem is: nearly all countries have highest efficiency at 1/1000 threshold 
plot_threshold <- function(df){
  
  plt <- 
    df %>%
    mutate(threshold = as.factor(threshold)) %>%
    mutate(country_name = countrycode(ISO, origin = "iso3c", destination = "country.name")) %>%
    ggplot(aes(x = incid_cv, y = log(baseline_incidence), color = threshold))+
    geom_label(aes(label = country_name), size = 4, show.legend = F, label.size = NA, label.padding = unit(0, "lines")) +
    geom_point()+
    theme_bw()+
    #geom_text(size=0.3, nudge_x = -1, nudge_y = 0.00001, hjust = 0.000001) +
    ylab("log(National Baseline Incidence Rate)")+
    xlab("Coefficient of Variation of district level incidence rate")+
    labs(color = "Threshold")
  
  return(plt)
}


## plot_gained_eff_with_size  ----------------------------------------------------
plot_gained_eff_with_size <- function(df_sum,
                                      baseline_estimate,
                                      improved_estimate,
                                      threshold_chosen,
                                      shade = "no"){
  
  if(threshold_chosen == 0.001){title = "Threshold: 10/10,000"}
  if(threshold_chosen == 2e-04){title = "Threshold: 2/10,000"}
  if(threshold_chosen == 1e-04){title = "Threshold: 1/10,000"}
  
  # filter out countries if there is no data for either of the estimates
  incomplete_ISO <- (df_sum %>%
                       filter(admin_level == "admin2") %>%
                       filter(threshold == threshold_chosen) %>%
                       filter(confirmation_lens == baseline_estimate | confirmation_lens == improved_estimate) %>%
                       group_by(ISO) %>% count() %>% filter(n < 2))$ISO
  
  df_sum <- df_sum %>% filter(!ISO %in% incomplete_ISO)
  
  df_wide <- df_sum %>% 
    filter(admin_level == "admin2" & threshold == threshold_chosen) %>%
    dplyr::select(ISO, threshold, confirmation_lens, efficiency) %>%
    # mutate(threshold = case_when(threshold == 0.001 ~ "Threshold: 1/1,000",
    #                              threshold == 2e-04 ~ "Threshold: 1/5,000",
    #                              threshold == 1e-04 ~ "Threshold: 1/10,000")) %>%
    # mutate(threshold = factor(threshold, levels = c("Threshold: 1/1,000", "Threshold: 1/5,000", "Threshold: 1/10,000"))) %>%
    pivot_wider(names_from = "confirmation_lens", values_from = "efficiency") %>%
    mutate(eff_gained = get(improved_estimate) - get(baseline_estimate)) %>%
    arrange(eff_gained) 
  
  df_dose <- df_sum %>% 
    filter(admin_level == "admin2" & threshold == threshold_chosen) %>%
    dplyr::select(ISO, threshold, confirmation_lens, doses) %>%
    pivot_wider(names_from = "confirmation_lens", values_from = c("doses")) %>%
    rename(`dose_district-estimate` = `district-estimate`,
           `dose_global-estimate` = `global-estimate`,
           `dose_no-estimate` = `no-estimate`)
  
  df_plt <- df_wide %>%
    left_join(df_dose, by = c("ISO", "threshold")) %>%
    mutate(country_name = countrycode(ISO, origin = "iso3c", destination = "country.name"))
    
  # arrange the ISO by eff_gained
  levels <- df_plt$country_name
  
  # dose legned
  dose_break1 <- c(5, 10, 15, 20)
  dose_break2 <- c(25, 50, 75, 100)
  dose_lab1 <- c("5", "10", "15", "20")
  dose_lab2 <- c("25", "50", "75", "100")
  dose_range1 <- c(1/2, 8/2)
  dose_range2 <- c(1,8)
  if(threshold_chosen == 0.001){
    dose_break <- dose_break1
    dose_lab <- dose_lab1
    dose_range <- dose_range1
  }else{
    dose_break <- dose_break2
    dose_lab <- dose_lab2
    dose_range <- dose_range2
  }
  
  plt <-
    df_plt %>%
    filter(!is.na(eff_gained)) %>%
    mutate(country_name = factor(country_name, levels = levels)) %>%
    ggplot()
  
  
  # add shade
  if(shade == "yes" & threshold_chosen == 0.001){
    plt <- plt + geom_rect(xmin = -Inf, xmax = 6, ymin = -Inf, ymax = Inf, fill = '#efefef', alpha = 0.1)
  }
  if(shade == "yes" & threshold_chosen != 0.001){
    plt <- plt + geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = '#efefef', alpha = 0.1)
  }
  
  # legend labels
  if(baseline_estimate == "no-estimate"){leg_labs_baseline <- "Clinical case definition"}
  if(baseline_estimate == "global-estimate"){leg_labs_baseline <- "National labs"}
  if(improved_estimate == "global-estimate"){leg_labs_improved <- "National labs"}
  if(improved_estimate == "district-estimate"){leg_labs_improved <- "District labs"}
  
  # label moving distance
  if(threshold_chosen == 0.001){label_distance = -0.7}
  if(threshold_chosen != 0.001){label_distance = -0.5}
  
  plt <- plt +
    geom_segment(aes(x = get(baseline_estimate), xend = get(improved_estimate),
                     y = country_name, yend = country_name)) +
    #geom_label(aes(label =  country_name, x = (get(baseline_estimate) - ISO_strlength)/30, y = ISO), size=3, label.size = NA, label.padding = unit(0, "lines")) +
    #geom_label(aes(label =  country_name, x = label_distance, y = ISO), size=3, label.size = NA, label.padding = unit(0, "lines")) +
    geom_point(aes(y = country_name, x = get(improved_estimate), color = "#fdb462", size = get(paste0("dose_", improved_estimate)))) +
    geom_point(aes(y = country_name, x = get(baseline_estimate), color = "#80b1d3", size = get(paste0("dose_", baseline_estimate)))) +
    xlab("OCV Efficiency") + 
    theme_bw() +
    theme(axis.text.y = element_text(size = 10),
          #axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          legend.text=element_text(size = 13),
          legend.title=element_text(size = 13)) +
    ggtitle(title) +
    scale_color_manual(name="OCV Efficiency", values=c("#fdb462", "#80b1d3"), labels = c(leg_labs_baseline, leg_labs_improved)) +
    scale_size_continuous(breaks=dose_break, labels=dose_lab, range = dose_range) + 
    guides(color = guide_legend(override.aes = list(size = c(1.5, 1.5) ) )) +
    labs(size = "Total doses administered (million)")
    
  return(plt)
}


## plot_eff_over_inc -------------------------------------
plot_eff_over_inc <- function(df_sum, 
                              df_incidence){
  plt <- df_sum %>%
    left_join(df_incidence, by = "ISO") %>%
    mutate(confirmation_lens = case_when(confirmation_lens == "district-estimate" ~ "District labs",
                                         confirmation_lens == "global-estimate" ~ "National lab",
                                         confirmation_lens == "no-estimate" ~ "Clinical definition")) %>%
    mutate(threshold = case_when(threshold == 1e-03 ~ "10/10,000",
                                 threshold == 2e-04 ~ "2/10,000",
                                 threshold == 1e-04 ~ "1/10,000")) %>%
    mutate(threshold = factor(threshold, levels = c("10/10,000", "2/10,000", "1/10,000"))) %>%
    mutate(admin_level = case_when(admin_level == "admin1" ~ "Province-level campaigns",
                                   admin_level == "admin2" ~ "District-level campaigns ")) %>%
    ggplot(aes(x = baseline_incidence, y = efficiency, color = as.factor(threshold))) +
    geom_point() + 
    geom_smooth(alpha = 0.3) + 
    facet_grid(confirmation_lens ~ admin_level, scales = "free", space = "free") + 
    theme_bw() +
    xlab("Incidence rate of clinical cases") + 
    ylab("OCV Efficiency") +
    labs(color = "Threshold") +
    scale_color_brewer(palette="Set2")

  return(plt)
}



## plot_gained_eff_over_inc  ----------------------------------------------------
plot_gained_eff_over_inc <- function(df_sum,
                                     baseline_estimate, 
                                     improved_estimate,
                                     df_incidence = df_incid0){
  df_wide <- df_sum %>% 
    dplyr::select(ISO, threshold, confirmation_lens, efficiency, admin_level) %>%
    # mutate(threshold = case_when(threshold == 0.001 ~ "Threshold: 1/1,000",
    #                              threshold == 2e-04 ~ "Threshold: 1/5,000",
    #                              threshold == 1e-04 ~ "Threshold: 1/10,000")) %>%
    # mutate(threshold = factor(threshold, levels = c("Threshold: 1/1,000", "Threshold: 1/5,000", "Threshold: 1/10,000"))) %>%
    pivot_wider(names_from = "confirmation_lens", values_from = "efficiency") %>%
    mutate(eff_gained = get(improved_estimate) - get(baseline_estimate)) %>%
    arrange(eff_gained) %>% 
    # link with national baseline incidence
    left_join(df_incidence, by = "ISO")
  
  # make scatter plot
  plt <- df_wide %>%
    ggplot(aes(x = baseline_incidence, y = eff_gained, color = as.factor(threshold))) + 
    geom_point() + 
    geom_smooth() +
    facet_grid( ~ admin_level, scales = "free", space = "free") + 
    theme_bw() +
    labs(color = "Threshold") + 
    xlab("") + ylab("") +
    ggtitle(paste0("Gained efficiency from ", baseline_estimate, " to ", improved_estimate))

  return(plt)
  
}

## plot_eff_over_cv -------------------------------------
plot_eff_over_cv <- function(df_sum, 
                             df_incidence){
  
  df_incidence <- df_incidence %>%
    group_by(ISO) %>% summarize(incid_cv = cv(baseline_incidence, na.rm = T) / 100)
  
  plt <- df_sum %>%
    left_join(df_incidence, by = "ISO") %>%
    ggplot(aes(x = incid_cv, y = efficiency, color = as.factor(threshold))) +
    geom_point() +
    geom_smooth() +
    facet_grid(confirmation_lens ~ admin_level, scales = "free", space = "free") + 
    theme_bw() +
    xlab("CV of district-level incidence rate") + 
    ylab("Efficiency") +
    labs(color = "Threshold")
  
  return(plt)
}



## plot_gained_eff_over_cv  ----------------------------------------------------
plot_gained_eff_over_cv <- function(df_sum,
                                    baseline_estimate, 
                                    improved_estimate,
                                    df_incidence){
  df_incidence <- df_incidence %>%
    group_by(ISO) %>% summarize(incid_cv = cv(baseline_incidence, na.rm = T) / 100)
  
  df_wide <- df_sum %>% 
    dplyr::select(ISO, threshold, confirmation_lens, efficiency, admin_level) %>%
    # mutate(threshold = case_when(threshold == 0.001 ~ "Threshold: 1/1,000",
    #                              threshold == 2e-04 ~ "Threshold: 1/5,000",
    #                              threshold == 1e-04 ~ "Threshold: 1/10,000")) %>%
    # mutate(threshold = factor(threshold, levels = c("Threshold: 1/1,000", "Threshold: 1/5,000", "Threshold: 1/10,000"))) %>%
    pivot_wider(names_from = "confirmation_lens", values_from = "efficiency") %>%
    mutate(eff_gained = get(improved_estimate) - get(baseline_estimate)) %>%
    arrange(eff_gained) %>% 
    # link with national baseline incidence
    left_join(df_incidence, by = "ISO")
  
  # make scatter plot
  plt <- df_wide %>%
    ggplot(aes(x = incid_cv, y = eff_gained, color = as.factor(threshold))) + 
    geom_point() + 
    geom_smooth() +
    facet_grid( ~ admin_level, scales = "free", space = "free") + 
    theme_bw() +
    labs(color = "Threshold") + 
    xlab("") + ylab("") +
    ggtitle(paste0("Gained efficiency from ", baseline_estimate, " to ", improved_estimate))
  
  return(plt)
  
}




# ## plot_gained_eff ----------------------------------------------------
# plot_gained_eff <- function(df_sum,
#                             baseline_estimate,
#                             improved_estimate,
#                             threshold_chosen){
#   
#   if(threshold_chosen == 0.001){title = "Threshold: 1/1,000"}
#   if(threshold_chosen == 2e-04){title = "Threshold: 1/5,000"}
#   if(threshold_chosen == 1e-04){title = "Threshold: 1/10,000"}
#   
#   # filter out countries if there is no data for either of the estimates
#   incomplete_ISO <- (df_sum %>%
#                        filter(admin_level == "admin2") %>%
#                        filter(threshold == threshold_chosen) %>%
#                        filter(confirmation_lens == baseline_estimate | confirmation_lens == improved_estimate) %>%
#                        group_by(ISO) %>% count() %>% filter(n < 2))$ISO
#   
#   df_sum <- df_sum %>% filter(!ISO %in% incomplete_ISO)
#   
#   df_wide <- df_sum %>% 
#     filter(admin_level == "admin2" & threshold == threshold_chosen) %>%
#     dplyr::select(ISO, threshold, confirmation_lens, efficiency) %>%
#     # mutate(threshold = case_when(threshold == 0.001 ~ "Threshold: 1/1,000",
#     #                              threshold == 2e-04 ~ "Threshold: 1/5,000",
#     #                              threshold == 1e-04 ~ "Threshold: 1/10,000")) %>%
#     # mutate(threshold = factor(threshold, levels = c("Threshold: 1/1,000", "Threshold: 1/5,000", "Threshold: 1/10,000"))) %>%
#     pivot_wider(names_from = "confirmation_lens", values_from = "efficiency") %>%
#     mutate(eff_gained = get(improved_estimate) - get(baseline_estimate)) %>%
#     arrange(eff_gained)
#   
#   # arrange the ISO by eff_gained
#   levels <- df_wide$ISO
#   
#   plt <-
#     df_wide %>%
#     mutate(country_name = countrycode(ISO, origin = "iso3c", destination = "country.name")) %>%
#     mutate(ISO = factor(ISO, levels = levels)) %>%
#     ggplot() +
#     geom_segment(aes(x = get(baseline_estimate), xend = get(improved_estimate),
#                      y = ISO, yend = ISO)) +
#     geom_label(aes(label =  country_name, x=get(baseline_estimate) - label_distance, y = ISO), size=3, label.size = NA, label.padding = unit(0, "lines")) +
#     geom_point(aes(y = ISO, x = get(improved_estimate), color = "#ED0000FF"), size = 2) +
#     geom_point(aes(y = ISO, x = get(baseline_estimate), color = "#00468BFF"), size = 2) +
#     xlab("OCV Efficiency") +
#     theme_bw() +
#     theme(axis.text.y = element_blank(),
#           axis.title.y = element_blank()) +
#     ggtitle(title) +
#     scale_color_manual(name="OCV Efficiency", values=c("#ED0000FF", "#00468BFF"), labels = c(baseline_estimate, improved_estimate)) +
#     guides(color = guide_legend(override.aes = list(size = c(1.5, 1.5) ) ) )
#   
#   return(plt)
# }

