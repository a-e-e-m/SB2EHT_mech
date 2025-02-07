# generating publication tables 

# packages
library(tidyverse)
library(tidybayes)
library(rstan)
library(gt)
library(webshot2)

# load runs
models2run <- read.csv(file = './models2run.csv', sep = ';', stringsAsFactors = FALSE)

# index of model run to fetch sigma_v, mu_x and sigma_x estimates from
i <- 4
dir_i <- paste0("fitting/", models2run$run[i]) # folder

# index of corresponding model run with int dose BA data only
i_BAonly <- 3
dir_i_BAonly <- paste0("fitting/", models2run$run[i_BAonly]) # folder


## grouping

  # group_number lookup joint model
  B_all <- readRDS(file.path(dir_i, 'B_all.rds'))
  
  group_number_lookup <- B_all |>
    filter(treat == 1) |>
    group_by(group_number, country, site, year, int_dose_available, Trial_code) |>
    summarise(insecticide = paste(unique(insecticide), collapse = ', '),
              test_type = paste(unique(test_type), collapse = ', '),
              .groups = "keep"
    )
  
  # group_number lookup BAonly model
  B_all_BAonly <- readRDS(file.path(paste0("fitting/", models2run$run[i_BAonly]) , 'B_all.rds'))
  group_number_conversion <- B_all_BAonly |> 
    filter(treat == 1) |>
    group_by(group_number, country, site, year, int_dose_available, Trial_code) |>
    summarise(insecticide = paste(unique(insecticide), collapse = ', '),
              test_type = paste(unique(test_type), collapse = ', '),
              .groups = "keep"
    ) |>
    dplyr::rename(group_number_conversion = group_number)
  
  # group number conversion
  group_number_lookup <- group_number_lookup |>
    left_join(group_number_conversion)
  
  # reorder groups
  # first the INT dose groups, then other groups by country, year, insecticide (with this ordering)
  group_number_lookup <- group_number_lookup |> 
    arrange(group_number_conversion, country, site, year, insecticide) |>
    rowid_to_column(var = "group_number_pub") 

## load data
  # load data joint model (BA and EHT)
  fit_i <- readRDS(file.path(dir_i, 'fit.rds'))
  # load data BAonly model
  fit_i_BAonly <- readRDS(file.path(dir_i_BAonly, 'fit.rds'))

  
## data groups table as csv
data_summary <- read_csv(file.path(dir_i, 'data_summary.csv'))

data_groups <- as_tibble(ungroup(data_summary)) |>
  left_join(as_tibble(ungroup(group_number_lookup)), 
            by = c("Group number" = "group_number", "Country" = "country", "Trial_code" , "test_type", "EHT site" = "site", "Year" = "year")) |>
  select(!c(insecticide, group_number_conversion)) |>
  mutate(Country = case_when(Country == "BurkinaFaso" ~ "Burkina Faso",
                             TRUE ~ Country)) |>
  rename("Site" = "EHT site",
         "Group" = "group_number_pub",
         "SB protocol" = test_type,
         "Hut type" = hut_type,
         "Insecticides \nin SB" = "Insecticides tested in resistance bio assay",
         "Insecticide \non ITN" = "Insecticides in nets",
         "ITN products" = "Net products",
         "Total in treated SB" = "Total number of mosquitoes in intervention resistance bio assays",
         "Total in control SB" = "Total number of mosquitoes in control resistance bio assays",
         "Total in intervention hut" = "Total number of mosquitoes in intervention arm of EHT",
         "Total in control hut" = "Total number of mosquitoes in control arm of EHT",
         "Group number old" = "Group number") |>
  arrange(Group) |>
  relocate(Group, "Group number old", Reference, Country, Site, Year, "SB protocol", "Hut type", "Insecticides \nin SB", "Insecticide \non ITN", "ITN products", "Total in treated SB", "Total in control SB", "Total in intervention hut", "Total in control hut")

  #save as csv to then add citation keys in excel, then convert to latex by online tool and copy paste into overleaf
  write_csv2(data_groups, file = file.path(dir_i, "Data_groups_table.csv"))

## parameter draws and diagnostics tables as image (use gtsave("tab_1.png", expand = 10)) with tidybayes package
  # joint model
    draws_i <- gather_draws(fit_i, p_b_control, sigma_v, mu_d[group_number], sigma_d[group_number], p_h_control[group_number], mu_x, sigma_x)
    summary_i <- tidybayes::summarise_draws(draws_i) |>
      left_join(group_number_lookup) |>
      arrange(!is.na(group_number_pub), match(.variable, c("p_b_control", "sigma_v", "mu_x", "sigma_x", "mu_d", "sigma_d", "p_h_control")), group_number_pub) |>
      select(!variable) |>
      ungroup() |>
      rename("Parameter" = ".variable",
             "Group" = "group_number_pub",
             "Mean" = "mean",
             "Median" = "median",
             "SD" = "sd",
             "MAD" = "mad",
             "Rhat" = "rhat",
             "ESS_bulk" = "ess_bulk",
             "ESS_tail" = "ess_tail",
             "5% Quantile" = "q5",
             "95% Quantile" = "q95") |>
      select(Parameter, Group, Mean, Median, SD, MAD, "5% Quantile", "95% Quantile", Rhat, ESS_bulk, ESS_tail)
    
    # polish table: change row names, round, spanners
    T_summary_i <- gt(summary_i) |>
      tab_header(
        title = "Parameter estimates and sampling diagnostics",
        subtitle = "from joint SB and EHT model"
      ) |>
      cols_label(Group = "Assay pair") |> 
      text_case_match(
        "p_b_control" ~ "p<sup>(0)</sup><sub>SB</sub>",
        "p_h_control" ~ "p<sup>(0)</sup><sub>EHT</sub>",
        "sigma_v" ~ "&sigma;<sub>V</sub>",
        "mu_x" ~ "&mu;<sub>X</sub>",
        "sigma_x" ~ "&sigma;<sub>X</sub>",
        "mu_d" ~ "&mu;<sub>T</sub>",
        "sigma_d" ~ "&sigma;<sub>T</sub>",
        "NA" ~ "-"
        #.locations = cells_stub()
      ) |>
      fmt_number(
        columns = 3:11,
        #decimals = 2,
        n_sigfig = 4
      ) |>
      tab_footnote(
        footnote = "Standard deviation",
        locations = cells_column_labels(columns = SD)
      ) |>
      tab_footnote(
        footnote = "Median absolute deviation",
        locations = cells_column_labels(columns = MAD)
      ) |>
      tab_footnote(
        footnote = "Rhat MCMC convergence diagnostic",
        locations = cells_column_labels(columns = Rhat)
      ) |>
      tab_footnote(
        footnote = "Bulk effective sample size",
        locations = cells_column_labels(columns = ESS_bulk)
      ) |>
      tab_footnote(
        footnote = "Tail effective sample size",
        locations = cells_column_labels(columns = ESS_tail)
      ) |>
      tab_source_note(
        source_note = md("Reference for sampling diagnostics: Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and Paul-Christian Bürkner (2021). Rank-normalization, folding, and localization: An improved R-hat for assessing convergence of MCMC (with discussion). *Bayesian Data Analysis*. 16(2), 667-–718. doi:10.1214/20-BA1221")
      ) 
    
    #save
    gtsave(T_summary_i, filename = "table_parameter_estimates_SBandEHT.html", path = dir_i)
    gtsave(T_summary_i, filename = "table_parameter_estimates_SBandEHT.html", path = file.path("tables_mechmodel_article_SI") ) 
    
  # BA_only model
    draws_i_BAonly <- gather_draws(fit_i_BAonly, p_b_control[Group], sigma_v, mu_d[Group], sigma_d[Group])
    summary_i_BAonly2 <- tidybayes::summarise_draws(draws_i_BAonly) |>
      arrange(!is.na(Group), match(.variable, c("p_b_control", "sigma_v", "mu_d", "sigma_d")), Group) |>
      select(!variable) |>
      ungroup() |>
      rename("Parameter" = ".variable",
             "Mean" = "mean",
             "Median" = "median",
             "SD" = "sd",
             "MAD" = "mad",
             "Rhat" = "rhat",
             "ESS_bulk" = "ess_bulk",
             "ESS_tail" = "ess_tail",
             "ESS_tail" = "ess_tail",
             "5% Quantile" = "q5",
             "95% Quantile" = "q95") |>
      select(Parameter, Group, Mean, Median, SD, MAD, "5% Quantile", "95% Quantile", Rhat, ESS_bulk, ESS_tail)

    # polish table: change row names, round, spanners
    T_summary_i_BAonly2 <- gt(summary_i_BAonly2) |>
      tab_header(
        title = "Parameter estimates and sampling diagnostics",
        subtitle = "from SB only model"
      ) |>
      cols_label(Group = "Assay pair") |> 
      text_case_match(
        "p_b_control" ~ "p<sup>(0)</sup><sub>SB</sub>",
        "sigma_v" ~ "&sigma;<sub>V</sub>",
        "mu_d" ~ "&mu;<sub>T</sub>",
        "sigma_d" ~ "&sigma;<sub>T</sub>",
        "NA" ~ "-"
        #.locations = cells_stub()
      ) |>
      fmt_number(
        columns = 3:11,
        #decimals = 2,
        n_sigfig = 4
      )  |>
      tab_footnote(
        footnote = "Standard deviation",
        locations = cells_column_labels(columns = SD)
      ) |>
      tab_footnote(
        footnote = "Median absolute deviation",
        locations = cells_column_labels(columns = MAD)
      ) |>
      tab_footnote(
        footnote = "Rhat MCMC convergence diagnostic",
        locations = cells_column_labels(columns = Rhat)
      ) |>
      tab_footnote(
        footnote = "Bulk effective sample size",
        locations = cells_column_labels(columns = ESS_bulk)
      ) |>
      tab_footnote(
        footnote = "Tail effective sample size",
        locations = cells_column_labels(columns = ESS_tail)
      ) |>
      tab_source_note(
        source_note = md("Reference for sampling diagnostics: Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and Paul-Christian Bürkner (2021). Rank-normalization, folding, and localization: An improved R-hat for assessing convergence of MCMC (with discussion). *Bayesian Data Analysis*. 16(2), 667-–718. doi:10.1214/20-BA1221")
      )
    
    #save
    gtsave(T_summary_i_BAonly2, filename = "table_parameter_estimates_SBonly.html", path = dir_i)
    gtsave(T_summary_i_BAonly2, filename = "table_parameter_estimates_SBonly.html", path = file.path("tables_mechmodel_article_SI") ) 
    
# ITN killing effect predictions 
    # joint model
    draws_i2 <- spread_draws(fit_i, p_b_control, sigma_v, mu_d[group_number], sigma_d[group_number], p_h_control[group_number], mu_x, sigma_x) |>
      left_join(group_number_lookup) |> 
      rename("Group" = "group_number_pub") |>
      ungroup() |>
      group_by(Group) |>
      mutate(LD50 = exp(mu_d),
             ITN_killing_effect = 0 + ( 1 - 0) * pnorm( (mu_x - mu_d) / sqrt(sigma_x^2 + sigma_d^2 )),
             lethal_dose = exp(rnorm(n(), mu_d, sigma_d))
      ) 
    
    pred_i <- draws_i2 |>
      median_qi(ITN_killing_effect, LD50, mu_d, sigma_d, .width = 0.95) |>
      left_join(summarise(draws_i2, corr = cor(mu_d, sigma_d, method = "pearson"))) |>
      select(!c(".width", ".point", ".interval")) 
    
    pred_i_lethaldose <- draws_i2 |>
      median_qi(lethal_dose, .width = c(0.95, 0.8, 0.5))
    
    # html table formatting
    T_pred_i <- gt(pred_i)  |>
      tab_header(
        title = "Model predictions per group",
        subtitle = "from joint SB and EHT model"
      ) |>
      tab_spanner(
        label = "ITN killing effect",
        columns = c(ITN_killing_effect, ITN_killing_effect.lower, ITN_killing_effect.upper)
      ) |>
      tab_spanner(
        label = html("LD<sub>50</sub>"),
        columns = c(LD50, LD50.lower, LD50.upper)
      )  |>
      tab_spanner(
        label = html("&mu;<sub>T</sub>"),
        columns = c(mu_d, mu_d.lower, mu_d.upper)
      ) |>
      tab_spanner(
        label = html("&sigma;<sub>T</sub>"),
        columns = c(sigma_d, sigma_d.lower, sigma_d.upper)
      )  |>
      tab_spanner(
        label = html("Cor(&mu;<sub>T</sub>, &sigma;<sub>T</sub>)"),
        columns = c(corr)
      )  |>
      cols_label(
        Group = "Assay pair",
        ITN_killing_effect = "median",
        ITN_killing_effect.lower = "q025",
        ITN_killing_effect.upper = "q975",
        LD50 = "median",
        LD50.lower = "q025",
        LD50.upper = "q975",
        mu_d = "median",
        mu_d.lower = "q025",
        mu_d.upper = "q975",
        sigma_d = "median",
        sigma_d.lower = "q025",
        sigma_d.upper = "q975",
        corr = ""
      ) |>
      fmt_percent(
        columns = c("ITN_killing_effect", "ITN_killing_effect.lower", "ITN_killing_effect.upper"),
        decimals = 2
      ) |>
      fmt_scientific(
        columns = c("mu_d", "mu_d.lower", "mu_d.upper", "LD50", "LD50.lower", "LD50.upper", "sigma_d", "sigma_d.lower", "sigma_d.upper"),
        decimals = 2
      ) |>
      fmt_number(
        columns = c(corr),
        decimals = 2
      ) |>
      tab_footnote(
        footnote = html("Population median of lethal dose, equal to exp(&mu;<sub>T</sub>) (scale: multiples of discriminating dose)"),
        locations = cells_column_spanners(spanners = "LD<sub>50</sub>")
      ) |>
      tab_footnote(
        footnote = "Mean of log lethal dose (scale: log of multiples of discriminating dose)",
        locations = cells_column_spanners(spanners = "&mu;<sub>T</sub>")
      ) |>
      tab_footnote(
        footnote = "Standard deviation of log lethal dose (scale: log of multiples of discriminating dose)",
        locations = cells_column_spanners(spanners = "&sigma;<sub>T</sub>")
      ) |>
      tab_footnote(
        footnote = "Pearson correlation coefficient",
        locations = cells_column_spanners(spanners = "Cor(&mu;<sub>T</sub>, &sigma;<sub>T</sub>)")
      ) 
      
    #save
    gtsave(T_pred_i, filename = "table_pred_SBandEHT.html", path = dir_i) 
    gtsave(T_pred_i, filename = "table_pred_SBandEHT.html", path = file.path("tables_mechmodel_article_SI") ) 
    
  
    # BA_only model
    draws_i_BAonly2 <- spread_draws(fit_i_BAonly, p_b_control[Group], sigma_v, mu_d[Group], sigma_d[Group]) |>
      mutate(LD50 = exp(mu_d)
      ) 
    
    pred_i_BAonly <- draws_i_BAonly2 |>
      median_qi(LD50, sigma_d, .width = 0.95) |>
      left_join(summarise(draws_i_BAonly2, corr = cor(mu_d, sigma_d, method = "pearson"))) |>
      select(!c(".width", ".point", ".interval")) 
    
    # html table formatting
    T_pred_i_BAonly <- gt(pred_i_BAonly)  |>
      tab_header(
        title = "Model predictions per group",
        subtitle = "from SB only model"
      ) |>
      tab_spanner(
        label = html("LD<sub>50</sub>"),
        columns = c(LD50, LD50.lower, LD50.upper)
      )  |>
      tab_spanner(
        label = html("&sigma;<sub>T</sub>"),
        columns = c(sigma_d, sigma_d.lower, sigma_d.upper)
      )  |>
      tab_spanner(
        label = html("Cor(&mu;<sub>T</sub>, &sigma;<sub>T</sub>)"),
        columns = c(corr)
      )  |>
      cols_label(
        LD50 = "median",
        LD50.lower = "q025",
        LD50.upper = "q975",
        sigma_d = "median",
        sigma_d.lower = "q025",
        sigma_d.upper = "q975",
        corr = ""
      ) |>
      fmt_scientific(
        columns = c("LD50", "LD50.lower", "LD50.upper", "sigma_d", "sigma_d.lower", "sigma_d.upper"),
        decimals = 2
      ) |>
      fmt_number(
        columns = c(corr),
        decimals = 2
      )
    
    #save 
    gtsave(T_pred_i_BAonly, filename = "table_pred_SBonly.html", path = dir_i_BAonly)  
    
   

    
    
    
    
    
    