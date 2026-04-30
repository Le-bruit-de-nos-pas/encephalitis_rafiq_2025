

library(readxl)
library(tidyverse)
library(data.table)

# Baseline summary info -----------

Kappa_df <- read_xlsx(path="../data/Kappa_NFL20250721_clean.xlsx", trim_ws = TRUE)

age_df <- read_xlsx(path="../data/codeDDN_IndexKappaNfL.xlsx", trim_ws = TRUE)
age_df$ddn <- as.Date(age_df$ddn)

Kappa_df$IndeK_clean <- as.numeric(Kappa_df$IndeK_clean)
Kappa_df$KappaPlasma <- as.numeric(Kappa_df$KappaPlasma)
Kappa_df$AlbuminPlasma <- as.numeric(Kappa_df$AlbuminPlasma)
Kappa_df$AlbuminCSF <- as.numeric(Kappa_df$AlbuminCSF)
Kappa_df$KappaCSF <- as.numeric(Kappa_df$KappaCSF)
Kappa_df$NfL <- as.numeric(Kappa_df$NfL)

length(unique(Kappa_df$Code_P)) # 50

Kappa_df %>% select(Code_P, Anticorps) %>% distinct() # 50

Kappa_df %>% group_by(Code_P) %>% count() %>%
  ungroup() %>%
  summarise(mean=mean(n),
            sd=sd(n),
            median=median(n),
            q1=quantile(n, 0.25),
            q3=quantile(n,0.75))

#    mean    sd median    q1    q3
#   <dbl> <dbl>  <dbl> <dbl> <dbl>
# 1  2.72 0.858      3     2     3


Kappa_df <- Kappa_df %>% mutate(NLR=`Neutrophiles - nb absolu (sang) clean`/`Lymphocytes - nb absolu (sang) clean`)

firsts_df <- Kappa_df %>% arrange(Code_P, Date_IndexK) %>% group_by(Code_P) %>% slice(1)

firsts_df %>% group_by(Anticorps) %>% count()

#  1 CASPR2               7
#  2 CV2                  2
#  3 DPPX                 1
#  4 GAD                  4
#  5 GFAP                 2
#  6 Hu                   1
#  7 Hu post ICI          2
#  8 IgLON5               5
#  9 KELCH11              1
# 10 LGI1                 8
# 11 MA2 post ICI         1
# 12 MAP1B                1
# 13 NMDAr                7
# 14 SERONEG              3
# 15 SERONEG post ICI     1
# 16 Yo                   4


variables_to_track <- c("IndeK_clean", "KappaPlasma", "AlbuminPlasma" , "AlbuminCSF" , "KappaCSF", 
                        "NfL", "Case",  "mRS", "abnormal MRI onset", "néo associé" ,"ICI", 
                        "CSF leuco count" , "CSF prot", "BOC_count"  )

summary_table <- firsts_df %>%
  ungroup() %>%
  select(all_of(variables_to_track)) %>%
  summarise(across(everything(), 
                   ~paste0(
                     round(mean(.x, na.rm = TRUE), 2), 
                     " (±", round(sd(.x, na.rm = TRUE), 2), ") | ",
                     round(median(.x, na.rm = TRUE), 2), 
                     " [", round(quantile(.x, 0.25, na.rm = TRUE), 2), 
                     "-", round(quantile(.x, 0.75, na.rm = TRUE), 2), "]"
                   )
  ))


# Optional: convert to long format for a tidy table
summary_table_long <- summary_table %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Summary")


summary_table_long

firsts_df %>% ungroup() %>% group_by(`abnormal MRI onset`) %>% count()


Kappa_df$`date symptômes` <- as.Date(as.numeric(Kappa_df$`date symptômes`), origin = "1899-12-30")
Kappa_df$`date dg (T0)`  <- as.Date(Kappa_df$`date dg (T0)`  )
Kappa_df$`Date_sampling`  <- as.Date(Kappa_df$`Date_sampling`  )


library(lubridate)

Kappa_df <- Kappa_df %>%
  mutate(elapsed_sympt = time_length(interval(`date symptômes`, Date_sampling ), "days")) %>%
  mutate(elapsed_dx = time_length(interval(`date dg (T0)`, Date_sampling ), "days")) 

firsts_df <- Kappa_df %>% arrange(Code_P, Date_IndexK) %>% group_by(Code_P) %>% slice(1)


firsts_df %>% ungroup() %>%
  mutate(elapsed_sympt=elapsed_sympt/30.44) %>%
  ungroup() %>%
  summarise(mean=mean(elapsed_sympt),
            sd=sd(elapsed_sympt),
            median=median(elapsed_sympt),
            q1=quantile(elapsed_sympt, 0.25),
            q3=quantile(elapsed_sympt,0.75))




antibodies <- tibble::tribble(
  ~Antibodies,     ~n,
  "CASPR2",        7,
  "CV2",           2,
  "DPPX",          1,
  "GAD",           4,
  "GFAP",          2,
  "Hu (ANNA-1)",   3,
  "IgLON5",        5,
  "KELCH11",       1,
  "LGI1",          8,
  "MA2",           1,
  "MAP1B",         1,
  "NMDAr",         7,
  "SERONEG",       4,
  "Yo / PCA-1",    4
)


# Add percentages
antibodies <- antibodies %>%
  mutate(percentage = n / sum(n) * 100)

# Reorder antibody levels by percentage (descending)
antibodies <- antibodies %>%
  arrange(desc(percentage)) %>%
  mutate(Antibodies = factor(Antibodies, levels = Antibodies))

# Stacked bar (sorted)
to_plot <- ggplot(antibodies, aes(x = "Autoantibodies", y = percentage, fill = Antibodies)) +
  geom_bar(stat = "identity", width = 0.6) +
  labs(x = "", y = "Cohort percentage (%) \n") +
  theme_minimal() +
    scale_fill_viridis_d(option = "cividis") +   # "plasma", "cividis", "magma", etc
   theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt")) +
  theme(legend.position = "right")

ggsave(file="antibodies_break.svg", plot=to_plot, width=3, height=4)


variables_to_track <- c("IndeK_clean", "KappaPlasma", "AlbuminPlasma" , "AlbuminCSF" , "KappaCSF", 
                        "NfL", "Case",  "mRS", 
                        "CSF leuco count" , "CSF prot", "BOC_count" )


names(firsts_df)



vars <- c(
  "elapsed_dx",
  "elapsed_sympt",
  "BOC_count",
  "CSF prot",
  "CSF leuco count",
"IndeK_clean",
    "mRS",
  "Case",
  "NfL",
  "AlbuminCSF",
  "AlbuminPlasma",
  "KappaCSF",
  "KappaPlasma"
)


extreme_outliers <- function(x) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  
  lower <- q1 - 3 * iqr
  upper <- q3 + 3 * iqr
  
  is_outlier <- x < lower | x > upper
  
  list(
    lower = lower,
    upper = upper,
    n_outliers = sum(is_outlier, na.rm = TRUE),
    outlier_idx = which(is_outlier)
  )
}

outlier_summary <- map_dfr(vars, function(v) {
  x <- firsts_df[[v]]
  res <- extreme_outliers(x)
  
  tibble(
    variable = v,
    n = sum(!is.na(x)),
    n_outliers = res$n_outliers,
    pct_outliers = round(100 * res$n_outliers / sum(!is.na(x)), 1),
    lower_bound = res$lower,
    upper_bound = res$upper
  )
})

outlier_summary

firsts_df_no_outliers <- firsts_df

for (v in vars) {
  x <- firsts_df_no_outliers[[v]]
  res <- extreme_outliers(x)
  
  firsts_df_no_outliers[[v]][res$outlier_idx] <- NA
}




to_plot <- firsts_df_no_outliers %>% mutate(elapsed_dx=elapsed_dx/30.44) %>% ggplot(aes(y=`elapsed_dx`, x="")) +
  geom_boxplot( color="#747475", alpha=0, notch = TRUE, width=0.5, linewidth=1.5) +
  geom_jitter(width=0.2, alpha=0.6, size=2, shape=1, stroke=2, color="#00204D") +
    labs(x = "First evaluation", 
       y = "Diagnosis duration (months)\n", 
       title="Diagnosis duration") +
   theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt"))
to_plot
ggsave(file="Diagnosis_first.svg", plot=to_plot, width=4, height=5)



to_plot <- firsts_df_no_outliers %>% mutate(elapsed_sympt=elapsed_sympt/30.44) %>% ggplot(aes(y=`elapsed_sympt`, x="")) +
  geom_boxplot( color="#747475", alpha=0, notch = TRUE, width=0.5, linewidth=1.5) +
  geom_jitter(width=0.2, alpha=0.6, size=2, shape=1, stroke=2, color="#00204D") +
    labs(x = "First evaluation", 
       y = "Symptom/disease duration (months)\n", 
       title="Symptom/disease duration") +
   theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt"))
to_plot
ggsave(file="elapsed_sympt_first.svg", plot=to_plot, width=4, height=5)




to_plot <- ggplot(firsts_df_no_outliers, aes(y=BOC_count, x="")) +
  geom_boxplot( color="#747475", alpha=0, notch = TRUE, width=0.5, linewidth=1.5) +
  geom_jitter(width=0.2, alpha=0.6, size=2, shape=1, stroke=2, color="#00204D") +
    labs(x = "First evaluation", 
       y = "Oligoclonal bands count\n", 
       title="OBC #") +
   theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt"))
to_plot
ggsave(file="OBC_first.svg", plot=to_plot, width=4, height=5)



to_plot <- ggplot(firsts_df_no_outliers, aes(y=`CSF prot`, x="")) +
  geom_boxplot( color="#747475", alpha=0, notch = TRUE, width=0.5, linewidth=1.5) +
  geom_jitter(width=0.2, alpha=0.6, size=2, shape=1, stroke=2, color="#00204D") +
    labs(x = "First evaluation", 
       y = "Protein CSF (g/L)\n", 
       title="Protein CSF") +
   theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt"))
to_plot
ggsave(file="protein_CSF_first.svg", plot=to_plot, width=4, height=5)





to_plot <- ggplot(firsts_df_no_outliers, aes(y=`CSF leuco count`, x="")) +
  geom_boxplot( color="#747475", alpha=0, notch = TRUE, width=0.5, linewidth=1.5) +
  geom_jitter(width=0.2, alpha=0.6, size=2, shape=1, stroke=2, color="#00204D") +
    labs(x = "First evaluation", 
       y = "Leukocytes CSF (cells/µL)\n", 
       title="Leukocytes CSF") +
   theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt"))
to_plot
ggsave(file="Leukocytes_CSF_first.svg", plot=to_plot, width=4, height=5)



to_plot <- ggplot(firsts_df_no_outliers, aes(y=mRS, x="")) +
  geom_boxplot( color="#747475", alpha=0, notch = TRUE, width=0.5, linewidth=1.5) +
  geom_jitter(width=0.2, alpha=0.6, size=2, shape=1, stroke=2, color="#00204D") +
    labs(x = "First evaluation", 
       y = "Modified Rankin Scale (mRS)\n", 
       title="mRS") +
   theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt"))
to_plot
ggsave(file="mrs_first.svg", plot=to_plot, width=4, height=5)


to_plot <- ggplot(firsts_df_no_outliers, aes(y=Case, x="")) +
  geom_boxplot( color="#747475", alpha=0, notch = TRUE, width=0.5, linewidth=1.5) +
  geom_jitter(width=0.2, alpha=0.6, size=2, shape=1, stroke=2, color="#00204D") +
    labs(x = "First evaluation", 
       y = "Clinical Assessment Scale in Autoimmune Encephalitis (CASE)\n", 
       title="CASE") +
   theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt"))
to_plot
ggsave(file="case_first.svg", plot=to_plot, width=4, height=5)




to_plot <- ggplot(firsts_df_no_outliers, aes(y=NfL, x="")) +
  geom_boxplot( color="#747475", alpha=0, notch = TRUE, width=0.5, linewidth=1.5) +
  geom_jitter(width=0.2, alpha=0.6, size=2, shape=1, stroke=2, color="#00204D") +
    labs(x = "First evaluation", 
       y = "Neurofilament light chain (NfL, pg/mL)\n", 
       title="Nfl") +
   theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt"))
to_plot
ggsave(file="NfL_first.svg", plot=to_plot, width=4, height=5)


to_plot <- ggplot(firsts_df_no_outliers, aes(y=AlbuminCSF, x="")) +
  geom_boxplot( color="#747475", alpha=0, notch = TRUE, width=0.5, linewidth=1.5) +
  geom_jitter(width=0.2, alpha=0.6, size=2, shape=1, stroke=2, color="#00204D") +
    labs(x = "First evaluation", 
       y = "Albumin CSF (mg/L)\n", 
       title="Albumin CSF") +
   theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt"))
to_plot
ggsave(file="AlbuminCSF_first.svg", plot=to_plot, width=4, height=5)


to_plot <- ggplot(firsts_df_no_outliers, aes(y=AlbuminPlasma, x="")) +
  geom_boxplot( color="#747475", alpha=0, notch = TRUE, width=0.5, linewidth=1.5) +
  geom_jitter(width=0.2, alpha=0.6, size=2, shape=1, stroke=2, color="#00204D") +
    labs(x = "First evaluation", 
       y = "Albumin Plasma (g/L)\n", 
       title="Albumin Plasma") +
   theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt"))
to_plot
ggsave(file="AlbuminPlasma_first.svg", plot=to_plot, width=4, height=5)



to_plot <- ggplot(firsts_df_no_outliers, aes(y=KappaCSF, x="")) +
  geom_boxplot( color="#747475", alpha=0, notch = TRUE, width=0.5, linewidth=1.5) +
  geom_jitter(width=0.2, alpha=0.6, size=2, shape=1, stroke=2, color="#00204D") +
    labs(x = "First evaluation", 
       y = "Kappa FLC CSF (mg/L)\n", 
       title="Kappa FLC CSF") +
   theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt"))
to_plot
ggsave(file="KappaCSF_first.svg", plot=to_plot, width=4, height=5)


to_plot <- ggplot(firsts_df_no_outliers, aes(y=KappaPlasma, x="")) +
  geom_boxplot( color="#747475", alpha=0, notch = TRUE, width=0.5, linewidth=1.5) +
  geom_jitter(width=0.2, alpha=0.6, size=2, shape=1, stroke=2, color="#00204D") +
    labs(x = "First evaluation", 
       y = "Kappa FLC Plasma (mg/L) \n", 
       title="Kappa FLC Plasma") +
   theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt"))
to_plot
ggsave(file="KappaPlasma_first.svg", plot=to_plot, width=4, height=5)


to_plot <- ggplot(firsts_df_no_outliers, aes(y=IndeK_clean , x="")) +
  geom_boxplot( color="#747475", alpha=0, notch = TRUE, width=0.5, linewidth=1.5) +
  geom_jitter(width=0.2, alpha=0.6, size=2, shape=1, stroke=2, color="#00204D") +
    labs(x = "First evaluation", 
       y = "Kappa FLC Index \n", 
       title="Kappa FLC Index") +
   theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt"))
to_plot
ggsave(file="Kappaindex_first.svg", plot=to_plot, width=4, height=5)















library(lubridate)


Kappa_df %>% group_by(Code_P) %>% filter(Date_sampling==min(Date_sampling)) %>%
  select(Code_P, Date_sampling) %>% mutate(Date_sampling=as.Date(Date_sampling)) %>%
  inner_join(age_df %>% rename("Code_P"="code")) %>%
    mutate(age = time_length(interval(ddn, Date_sampling ), "years"))  %>%
    ungroup() %>%
  summarise(mean=mean(age),
            sd=sd(age),
            median=median(age),
            q1=quantile(age, 0.25),
            q3=quantile(age,0.75))




Kappa_df <- Kappa_df %>%
  mutate(elapsed_sympt = time_length(interval(`date symptômes`, Date_sampling ), "days")) %>%
  mutate(elapsed_dx = time_length(interval(`date dg (T0)`, Date_sampling ), "days")) 





to_plot <-  Kappa_df %>% group_by(Code_P) %>% filter(Date_sampling==min(Date_sampling)) %>%
  select(Code_P, Date_sampling) %>% mutate(Date_sampling=as.Date(Date_sampling)) %>%
  inner_join(age_df %>% rename("Code_P"="code")) %>%
    mutate(age = time_length(interval(ddn, Date_sampling ), "years"))  %>%
    ungroup() %>%
  ggplot(aes(y=age, x="")) +
  geom_boxplot( color="#747475", alpha=0, notch = TRUE, width=0.5, linewidth=1.5) +
  geom_jitter(width=0.2, alpha=0.6, size=2, shape=1, stroke=2, color="#00204D") +
    labs(x = "First evaluation", 
       y = "Age baseline (years)\n", 
       title="Age baseline") +
   theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt"))

ggsave(file="age_first.svg", plot=to_plot, width=4, height=5)



# --------------


# Correlations ---------------

Kappa_df <- read_xlsx(path="../data/Kappa_NFL20250721_clean.xlsx", trim_ws = TRUE)

Kappa_df$IndeK_clean <- as.numeric(Kappa_df$IndeK_clean)
Kappa_df$KappaPlasma <- as.numeric(Kappa_df$KappaPlasma)
Kappa_df$AlbuminPlasma <- as.numeric(Kappa_df$AlbuminPlasma)
Kappa_df$AlbuminCSF <- as.numeric(Kappa_df$AlbuminCSF)
Kappa_df$KappaCSF <- as.numeric(Kappa_df$KappaCSF)
Kappa_df$NfL <- as.numeric(Kappa_df$NfL)

length(unique(Kappa_df$Code_P)) # 50

Kappa_df <- Kappa_df %>% mutate(NLR=`Neutrophiles - nb absolu (sang) clean`/`Lymphocytes - nb absolu (sang) clean`)


Kappa_df$`date symptômes` <- as.Date(as.numeric(Kappa_df$`date symptômes`), origin = "1899-12-30")
Kappa_df$`date dg (T0)`  <- as.Date(Kappa_df$`date dg (T0)`  )
Kappa_df$`Date_sampling`  <- as.Date(Kappa_df$`Date_sampling`  )

library(lubridate)



age_df <- read_xlsx(path="../data/codeDDN_IndexKappaNfL.xlsx", trim_ws = TRUE)
age_df$ddn <- as.Date(age_df$ddn)

Kappa_df <- Kappa_df %>% inner_join(age_df %>% rename("Code_P"="code"))


Kappa_df <- Kappa_df %>%
  mutate(elapsed_sympt = time_length(interval(`date symptômes`, Date_sampling ), "months")) %>%
  mutate(elapsed_dx = time_length(interval(`date dg (T0)`, Date_sampling ), "months"))  %>%
  mutate(age = time_length(interval(ddn, Date_sampling ), "years"))


variables_to_track <- c("Code_P",  "IndeK_clean", "KappaPlasma", "AlbuminPlasma" , "AlbuminCSF" , "KappaCSF", 
                        "NfL", "Case",  "mRS", "abnormal MRI onset", "néo associé" ,"ICI", 
                        "CSF leuco count" , "CSF prot", "BOC_count" , "CRP" , "Neutrophiles - nb absolu (sang) clean" ,
                        "Lymphocytes - nb absolu (sang) clean" , "NLR", "elapsed_sympt", "elapsed_dx", "age")

Kappa_df <- Kappa_df %>% select(all_of(variables_to_track))




vars <- c(
  "elapsed_dx",
  "elapsed_sympt",
  "BOC_count",
  "CSF prot",
  "CSF leuco count",
  "IndeK_clean",
    "mRS",
  "Case",
  "NfL",
  "AlbuminCSF",
  "AlbuminPlasma",
  "KappaCSF",
  "KappaPlasma"
)


extreme_outliers <- function(x) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  
  lower <- q1 - 3 * iqr
  upper <- q3 + 3 * iqr
  
  is_outlier <- x < lower | x > upper
  
  list(
    lower = lower,
    upper = upper,
    n_outliers = sum(is_outlier, na.rm = TRUE),
    outlier_idx = which(is_outlier)
  )
}


firsts_df_no_outliers <- firsts_df

for (v in vars) {
  x <- firsts_df_no_outliers[[v]]
  res <- extreme_outliers(x)
  
  firsts_df_no_outliers[[v]][res$outlier_idx] <- NA
}


names(firsts_df_no_outliers)

target_var <- "elapsed_sympt"
label <- "Disease duration"

to_plot <- firsts_df_no_outliers %>%  mutate(elapsed_sympt=elapsed_sympt/30.44) %>%
  select(all_of(c(target_var, "Case", "mRS"))) %>%
  pivot_longer(cols = c(Case, mRS),
               names_to = "Outcome",
               values_to = "Value")

ploted <- ggplot(to_plot, aes_string(x = target_var , y = "Value")) +
  geom_jitter(width=0.2, alpha=0.6, size=2, shape=1, stroke=2, color="#00204D") +
  geom_smooth(method = "lm", se = TRUE, color = "#747475", fill="#747475") +
  facet_wrap(~Outcome, ncol = 1, scales = "free_y") +
  labs(
    x = paste0(label),
    y = "Case / mRS values",
    title = paste(label)
  ) +
   theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold", hjust = 0),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt")) 

ploted
ggsave(file=paste0(label, "_all.svg"), plot=ploted, width=5, height=7)


options(scipen = 999)


names(firsts_df_no_outliers)



vars <- firsts_df %>% ungroup() %>%
  select(-Code_P, -`Delay_T1T2(m)`, -`Delay_T2T3(m)`,
  -`Delay_T3T4(m)`, -`Delay_T4T5(m)`,  -`Delay_T5T6(m)`,
  -CRP , -`Neutrophiles - nb absolu (sang) clean`, 
 -`Lymphocytes - nb absolu (sang) clean` , -`Traitement 1ère ligne (O/N)`, 
 -`Traitement 2nd ligne (O/N)`, -NLR) %>%
  select(where(is.numeric)) %>%
  names()


# Select only your variables
df_corr <- firsts_df_no_outliers %>%
  ungroup() %>%
  select(all_of(vars)) 

# Compute Spearman correlation matrix
cor_matrix <- cor(df_corr, method = "spearman", use = "pairwise.complete.obs")

round(cor_matrix,2)

colnames(cor_matrix)


library(Hmisc)

res <- rcorr(as.matrix(df_corr), type = "pearson")

cor_matrix <- res$r
p_matrix   <- res$P

cor_long <- as.data.frame(cor_matrix) %>%
  rownames_to_column("Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "Correlation")

p_long <- as.data.frame(p_matrix) %>%
  rownames_to_column("Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "p_value")

# Merge
heatmap_long <- left_join(cor_long, p_long, by = c("Var1", "Var2"))


heatmap_long <- heatmap_long %>%
  mutate(
    stars = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    label = sprintf("%.2f%s", Correlation, stars)
  )


var_labels <- c(
  "IndeK_clean"        = "κ FLC Index",
  "KappaPlasma"        = "κ FLC Plasma",
  "AlbuminPlasma"      = "Albumin Plasma",
  "AlbuminCSF"         = "Albumin CSF",
  "KappaCSF"           = "κ FLC CSF",
  "NfL"                = "NfL",
  "Case"               = "CASE",
  "mRS"                = "mRS",
  "abnormal MRI onset" = "Abnormal MRI",
  "néo associé"        = "Neoplasia",
  "ICI"                = "Post-ICI",
  "CSF leuco count"    = "Leuko CSF",
  "CSF prot"           = "CSF Protein",
  "BOC_count"          = "OCB count",
  "elapsed_sympt"      = "Symptom duration",
  "elapsed_dx"         = "Diagnosis duration"
)

heatmap_long <- heatmap_long %>%
  mutate(
    Var1 = recode(Var1, !!!var_labels),
    Var2 = recode(Var2, !!!var_labels)
  )


var_order_clean <- c(
  "mRS",
  "CASE",
  "Abnormal MRI",
  "Neoplasia",
  "Post-ICI",
  "Symptom duration",
  "Diagnosis duration",
  "κ FLC Index",
  "κ FLC CSF",
  "NfL",
  "OCB count",
  "Leuko CSF",
  "CSF Protein",
  "κ FLC Plasma",
  "Albumin CSF",
  "Albumin Plasma"
)

heatmap_long <- heatmap_long %>%
  mutate(
    Var1 = factor(Var1, levels = (var_order_clean) ),
    Var2 = factor(Var2, levels = rev(var_order_clean) )
  )


plot <- ggplot(heatmap_long, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient2(
    low = "#00204D",
    mid = "white",
    high = "firebrick",
    midpoint = 0,
    name = "Pearson\nrho"
  ) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plot
ggsave(file="correlations.svg", plot=plot, width=9, height=8) 









# -----------


# XGBoost CASE --------------
Kappa_df <- read_xlsx(path="../data/Kappa_NFL20250721_clean.xlsx", trim_ws = TRUE)


age_df <- read_xlsx(path="../data/codeDDN_IndexKappaNfL.xlsx", trim_ws = TRUE)
age_df$ddn <- as.Date(age_df$ddn)


Kappa_df$IndeK_clean <- as.numeric(Kappa_df$IndeK_clean)
Kappa_df$KappaPlasma <- as.numeric(Kappa_df$KappaPlasma)
Kappa_df$AlbuminPlasma <- as.numeric(Kappa_df$AlbuminPlasma)
Kappa_df$AlbuminCSF <- as.numeric(Kappa_df$AlbuminCSF)
Kappa_df$KappaCSF <- as.numeric(Kappa_df$KappaCSF)
Kappa_df$NfL <- as.numeric(Kappa_df$NfL)

length(unique(Kappa_df$Code_P)) # 50

Kappa_df <- Kappa_df %>% mutate(NLR=`Neutrophiles - nb absolu (sang) clean`/`Lymphocytes - nb absolu (sang) clean`)


Kappa_df$`date symptômes` <- as.Date(as.numeric(Kappa_df$`date symptômes`), origin = "1899-12-30")
Kappa_df$`date dg (T0)`  <- as.Date(Kappa_df$`date dg (T0)`  )
Kappa_df$`Date_sampling`  <- as.Date(Kappa_df$`Date_sampling`  )

library(lubridate)


Kappa_df <- Kappa_df %>%
   inner_join(age_df %>% rename("Code_P"="code")) %>%
    mutate(age = time_length(interval(ddn, Date_sampling ), "years"))  %>%
  mutate(elapsed_sympt = time_length(interval(`date symptômes`, Date_sampling ), "months")) %>%
  mutate(elapsed_dx = time_length(interval(`date dg (T0)`, Date_sampling ), "months"))  



variables_to_track <- c("Code_P",  "IndeK_clean", "KappaPlasma", "AlbuminPlasma" , "AlbuminCSF" , "KappaCSF", 
                        "NfL", "Case",  "mRS", "abnormal MRI onset", "néo associé" ,"ICI", 
                        "CSF leuco count" , "CSF prot", "BOC_count" ,  "elapsed_sympt", "elapsed_dx", "age")

Kappa_df <- Kappa_df %>% select(all_of(variables_to_track))


Kappa_df <- Kappa_df %>% rename("Leukocytes CSF"="CSF leuco count") %>%
  rename("Protein CSF"="CSF prot") %>%
  rename("κ FLC Index"="IndeK_clean") %>%
  rename("Albumin Plasma"="AlbuminPlasma") %>%
  rename("κ Plasma"="KappaPlasma") %>%
  rename("Albumin CSF"="AlbuminCSF") %>%
  rename("Oligoclonal bands"="BOC_count") %>%
  rename("κ CSF"="KappaCSF") %>%
  rename("Abnormal MRI"="abnormal MRI onset") %>%
  rename("Neoplasia-assoc"="néo associé") 



Kappa_df <- Kappa_df %>% group_by(Code_P) %>% mutate(future_case=lead(Case)) %>% ungroup() %>%
  select(-elapsed_dx, -elapsed_sympt, -Code_P, -mRS, -Case) %>%
  filter(!is.na(future_case))




library(xgboost)
library(dplyr)

# Remove target and convert tibble to data.frame
df <- Kappa_df %>% as.data.frame()

# Prepare data: predictors and target
X <- df %>% select(-future_case) # all columns except target
y <- df$future_case



X_mat <- as.matrix(X)
n <- nrow(X_mat)

# Create 10 folds indices
folds <- sample(rep(1:10, length.out = n))

preds_full <- numeric(n)  # store predictions for all samples

for (fold in 1:10) {
  train_idx <- which(folds != fold)
  test_idx <- which(folds == fold)
  
  dtrain <- xgb.DMatrix(data = X_mat[train_idx, , drop=FALSE], 
                        label = y[train_idx], missing = NA)
  dtest <- xgb.DMatrix(data = X_mat[test_idx, , drop=FALSE], missing = NA)
  
  # Train model
  model <- xgb.train(params = list(objective = "reg:squarederror"),
                     data = dtrain,
                     nrounds = 100,
                     verbose = 0)
  
  # Predict test fold
  preds_full[test_idx] <- predict(model, dtest)
}


rmse_cv <- y %>% bind_cols(preds_full) %>% filter(`...1`<20& `...2`<20)

# Calculate RMSE on all samples
error <- sqrt(mean((rmse_cv$`...1` - rmse_cv$`...2`)^2))
print(paste("10-fold CV RMSE:", error))



# Plot observed vs predicted for all samples
plot_df <- data.frame(observed = rmse_cv$`...1`, predicted = rmse_cv$`...2`)




#plot_df_inc <- plot_df %>% mutate(group="Including Prev mRS|CASE")


plot_df_exc <- plot_df %>% mutate(group="Excluding Prev mRS|CASE")


to_plot <- plot_df_inc %>% bind_rows(plot_df_exc) %>%
  ggplot(aes(observed, predicted, colour=group, fill=group)) +
  geom_jitter(width=0.2, alpha=0.5, size=1.5, shape=1, stroke=2) +
  geom_smooth(method = "lm", se = FALSE) +
  coord_cartesian(xlim = c(0, 15), ylim = c(0, 15)) +
  scale_colour_manual(values=c("#00204D", "firebrick")) +
  scale_fill_manual(values=c("#00204D", "firebrick")) +
  labs(title = "",
       x = "\n Observed Future CASE Score", y = "Predicted Future CASE Score \n") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt"))

to_plot
ggsave(file="xgboost.svg", plot=to_plot, width=7, height=5)




# Train final model on full data for feature importance
dtrain_full <- xgb.DMatrix(data = X_mat, label = y, missing = NA)

final_model <- xgb.train(params = list(objective = "reg:squarederror"),
                         data = dtrain_full,
                         nrounds = 100,
                         verbose = 0)

importance <- xgb.importance(model = final_model)

xgb.plot.importance(importance)

to_plot <- importance %>%
  arrange(Importance) %>%
  mutate(Feature=factor(Feature, levels=Feature)) %>%
  ggplot(aes(Feature, Importance)) +
  geom_col(fill="#00204D", colour="#00204D", alpha=0.8) +
  coord_flip() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt")) +
  ylab("\n Feature")+ xlab("Relative Feature Importance (a.u.) \n")


to_plot

ggsave(file="xgboost_im_inc.svg", plot=to_plot, width=5, height=5)

# -----------

# XGBoost mRS --------------
Kappa_df <- read_xlsx(path="../data/Kappa_NFL20250721_clean.xlsx", trim_ws = TRUE)


age_df <- read_xlsx(path="../data/codeDDN_IndexKappaNfL.xlsx", trim_ws = TRUE)
age_df$ddn <- as.Date(age_df$ddn)

Kappa_df$IndeK_clean <- as.numeric(Kappa_df$IndeK_clean)
Kappa_df$KappaPlasma <- as.numeric(Kappa_df$KappaPlasma)
Kappa_df$AlbuminPlasma <- as.numeric(Kappa_df$AlbuminPlasma)
Kappa_df$AlbuminCSF <- as.numeric(Kappa_df$AlbuminCSF)
Kappa_df$KappaCSF <- as.numeric(Kappa_df$KappaCSF)
Kappa_df$NfL <- as.numeric(Kappa_df$NfL)

length(unique(Kappa_df$Code_P)) # 50

Kappa_df <- Kappa_df %>% mutate(NLR=`Neutrophiles - nb absolu (sang) clean`/`Lymphocytes - nb absolu (sang) clean`)


Kappa_df$`date symptômes` <- as.Date(as.numeric(Kappa_df$`date symptômes`), origin = "1899-12-30")
Kappa_df$`date dg (T0)`  <- as.Date(Kappa_df$`date dg (T0)`  )
Kappa_df$`Date_sampling`  <- as.Date(Kappa_df$`Date_sampling`  )

library(lubridate)

Kappa_df <- Kappa_df %>%
   inner_join(age_df %>% rename("Code_P"="code")) %>%
    mutate(age = time_length(interval(ddn, Date_sampling ), "years"))  %>%
  mutate(elapsed_sympt = time_length(interval(`date symptômes`, Date_sampling ), "months")) %>%
  mutate(elapsed_dx = time_length(interval(`date dg (T0)`, Date_sampling ), "months"))  



variables_to_track <- c("Code_P",  "IndeK_clean", "KappaPlasma", "AlbuminPlasma" , "AlbuminCSF" , "KappaCSF", 
                        "NfL", "Case",  "mRS", "abnormal MRI onset", "néo associé" ,"ICI", 
                        "CSF leuco count" , "CSF prot", "BOC_count" ,  "elapsed_sympt", "elapsed_dx", "age")

Kappa_df <- Kappa_df %>% select(all_of(variables_to_track))


Kappa_df <- Kappa_df %>% rename("Leukocytes CSF"="CSF leuco count") %>%
  rename("Protein CSF"="CSF prot") %>%
  rename("κ FLC Index"="IndeK_clean") %>%
  rename("Albumin Plasma"="AlbuminPlasma") %>%
  rename("κ Plasma"="KappaPlasma") %>%
  rename("Albumin CSF"="AlbuminCSF") %>%
  rename("Oligoclonal bands"="BOC_count") %>%
  rename("κ CSF"="KappaCSF") %>%
  rename("Abnormal MRI"="abnormal MRI onset") %>%
  rename("Neoplasia-assoc"="néo associé") 



Kappa_df <- Kappa_df %>% group_by(Code_P) %>% mutate(future_mRS=lead(mRS)) %>% ungroup() %>%
  select(-elapsed_dx, -elapsed_sympt, -Code_P, -mRS, -Case) %>%
  filter(!is.na(future_mRS))




library(xgboost)
library(dplyr)

# Remove target and convert tibble to data.frame
df <- Kappa_df %>% as.data.frame()

# Prepare data: predictors and target
X <- df %>% select(-future_mRS) # all columns except target
y <- df$future_mRS



X_mat <- as.matrix(X)
n <- nrow(X_mat)

# Create 10 folds indices
folds <- sample(rep(1:10, length.out = n))

preds_full <- numeric(n)  # store predictions for all samples

for (fold in 1:10) {
  train_idx <- which(folds != fold)
  test_idx <- which(folds == fold)
  
  dtrain <- xgb.DMatrix(data = X_mat[train_idx, , drop=FALSE], 
                        label = y[train_idx], missing = NA)
  dtest <- xgb.DMatrix(data = X_mat[test_idx, , drop=FALSE], missing = NA)
  
  # Train model
  model <- xgb.train(params = list(objective = "reg:squarederror"),
                     data = dtrain,
                     nrounds = 100,
                     verbose = 0)
  
  # Predict test fold
  preds_full[test_idx] <- predict(model, dtest)
}


rmse_cv <- y %>% bind_cols(preds_full) %>% filter(`...1`<20& `...2`<20)

# Calculate RMSE on all samples
error <- sqrt(mean((rmse_cv$`...1` - rmse_cv$`...2`)^2))
print(paste("10-fold CV RMSE:", error))


# 0.7  1.22

# Plot observed vs predicted for all samples
plot_df <- data.frame(observed = rmse_cv$`...1`, predicted = rmse_cv$`...2`)




#plot_df_inc <- plot_df %>% mutate(group="Including Prev mRS|CASE")

plot_df_exc <- plot_df %>% mutate(group="Excluding Prev mRS|CASE")


to_plot <- plot_df_inc %>% bind_rows(plot_df_exc) %>%
  ggplot(aes(observed, predicted, colour=group, fill=group)) +
  geom_jitter(width=0.2, alpha=0.5, size=1.5, shape=1, stroke=2) +
  geom_smooth(method = "lm", se = FALSE) +
  coord_cartesian(xlim = c(0, 6), ylim = c(0, 6)) +
  scale_colour_manual(values=c("#00204D", "firebrick")) +
  scale_fill_manual(values=c("#00204D", "firebrick")) +
  labs(title = "",
       x = "\n Observed Future mRS Score", y = "Predicted Future mRS Score \n") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt"))

ggsave(file="xgboost.svg", plot=to_plot, width=7, height=5)




# Train final model on full data for feature importance
dtrain_full <- xgb.DMatrix(data = X_mat, label = y, missing = NA)
final_model <- xgb.train(params = list(objective = "reg:squarederror"),
                         data = dtrain_full,
                         nrounds = 100,
                         verbose = 0)

importance <- xgb.importance(model = final_model)
xgb.plot.importance(importance)

to_plot <- importance %>%
  arrange(Importance) %>%
  mutate(Feature=factor(Feature, levels=Feature)) %>%
  ggplot(aes(Feature, Importance)) +
  geom_col(fill="#00204D", colour="#00204D", alpha=0.8) +
  coord_flip() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt")) +
  ylab("\n Feature")+ xlab("Relative Feature Importance (a.u.) \n")


ggsave(file="xgboost_im_inc.svg", plot=to_plot, width=5, height=5)

# -----------


# Observed clinical/biochemical deltas and physician assessment -----------

age_df <- read_xlsx(path="../data/codeDDN_IndexKappaNfL.xlsx", trim_ws = TRUE)
age_df$ddn <- as.Date(age_df$ddn)

Kappa_df <- read_xlsx(path="../data/Kappa_NFL20250721_clean.xlsx", trim_ws = TRUE)

Kappa_df$IndeK_clean <- as.numeric(Kappa_df$IndeK_clean)
Kappa_df$KappaPlasma <- as.numeric(Kappa_df$KappaPlasma)
Kappa_df$AlbuminPlasma <- as.numeric(Kappa_df$AlbuminPlasma)
Kappa_df$AlbuminCSF <- as.numeric(Kappa_df$AlbuminCSF)
Kappa_df$KappaCSF <- as.numeric(Kappa_df$KappaCSF)
Kappa_df$NfL <- as.numeric(Kappa_df$NfL)

length(unique(Kappa_df$Code_P)) # 50

Kappa_df <- Kappa_df %>% mutate(NLR=`Neutrophiles - nb absolu (sang) clean`/`Lymphocytes - nb absolu (sang) clean`)


Kappa_df$`date symptômes` <- as.Date(as.numeric(Kappa_df$`date symptômes`), origin = "1899-12-30")
Kappa_df$`date dg (T0)`  <- as.Date(Kappa_df$`date dg (T0)`  )
Kappa_df$`Date_sampling`  <- as.Date(Kappa_df$`Date_sampling`  )

library(lubridate)

Kappa_df <- Kappa_df %>%
   inner_join(age_df %>% rename("Code_P"="code")) %>%
    mutate(age = time_length(interval(ddn, Date_sampling ), "years"))  %>%
  mutate(elapsed_sympt = time_length(interval(`date symptômes`, Date_sampling ), "months")) %>%
  mutate(elapsed_dx = time_length(interval(`date dg (T0)`, Date_sampling ), "months"))  

names(Kappa_df)

variables_to_track <- c("Code_P",  "IndeK_clean", "KappaPlasma", "AlbuminPlasma" , "AlbuminCSF" , "KappaCSF", 
                        "NfL", "Case",  "mRS",  "Date_sampling",
                        "CSF leuco count" , "CSF prot", "BOC_count" ,  "Statut clinique", "age")



Kappa_df <- Kappa_df %>% select(all_of(variables_to_track))

names(Kappa_df)

Kappa_df <- Kappa_df %>% rename("CSF_leuco_count"="CSF leuco count") %>%
  rename("CSF_prot"="CSF prot") %>%
  rename("Statut_clinique"="Statut clinique")


#Kappa_df <- Kappa_df %>% select(-Date_sampling)

num_vars <- Kappa_df %>%
  select(where(is.numeric)) %>%
  names()


library(dplyr)

Kappa_delta <- Kappa_df %>%
  group_by(Code_P) %>%
  arrange(Code_P, Date_sampling) %>%   
  mutate(
    across(
      all_of(num_vars),
      ~ .x - lag(.x),
      .names = "delta_{.col}"
    )
  ) %>%
  ungroup()



Kappa_delta <- Kappa_delta %>% select(Code_P, Statut_clinique, delta_IndeK_clean:delta_age) %>% 
  ungroup() %>% drop_na()


delta_scaled <- Kappa_delta %>%
  mutate(across(where(is.numeric), scale))


unique(delta_scaled$Statut_clinique)

delta_scaled <- delta_scaled %>%
  mutate(
    Statut_clinique = str_to_lower(Statut_clinique),
    Statut_clinique = case_when(
      str_detect(Statut_clinique, "améli") ~ "amelioration",
      str_detect(Statut_clinique, "stable") ~ "stable",
      str_detect(Statut_clinique, "worst|aggrav|détério") ~ "worsening",
      TRUE ~ Statut_clinique
    ),
    Statut_clinique = factor(Statut_clinique, levels = c("worsening", "stable", "amelioration"),
      ordered = TRUE)
  )


delta_vars <- delta_scaled %>%
  select(starts_with("delta_")) %>%
  names()

kw_results <- lapply(delta_vars, function(v) {

  test <- kruskal.test(
    delta_scaled[[v]] ~ delta_scaled$Statut_clinique
  )

  data.frame(
    variable = v,
    p_value = test$p.value
  )
}) %>%
  bind_rows() %>%
  arrange(p_value)

kw_results

#                 variable      p_value
# 1             delta_Case 2.499387e-09
# 2              delta_mRS 1.068793e-05
# 3         delta_KappaCSF 3.769250e-02
# 4         delta_CSF_prot 7.161946e-02
# 5       delta_AlbuminCSF 1.886066e-01
# 6  delta_CSF_leuco_count 2.890404e-01
# 7      delta_KappaPlasma 3.188212e-01
# 8      delta_IndeK_clean 3.962789e-01
# 9              delta_age 4.281019e-01
# 10             delta_NfL 9.014621e-01
# 11   delta_AlbuminPlasma 9.131346e-01
# 12       delta_BOC_count 9.596542e-01

plot_df <- delta_scaled %>%
  group_by(Statut_clinique) %>%
  summarise(
    mean_delta = mean(delta_Case, na.rm = TRUE),
    se_delta = sd(delta_Case, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

plot <- plot_df %>%
  ggplot(aes(x = Statut_clinique, y = mean_delta, colour=Statut_clinique, fill=Statut_clinique)) +
  geom_col(alpha=0.8) +
  geom_segment(
    aes(
      x = Statut_clinique,
      xend = Statut_clinique,
      y = mean_delta - se_delta,
      yend = mean_delta + se_delta
    ),
    linewidth = 3.2,
    lineend = "round", alpha=0.6
  ) +
  labs(
    x = "\n  \nPhysician-based clinical trajectory",
    y = "Δ CASE \n [Standardized Mean ± SE] \n",
    title="Δ CASE & Δ Clinical Trajectory"
  ) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt")) +
  scale_fill_manual(values=c("#00204D", "darkgray", "firebrick")) +
  scale_colour_manual(values=c("#00204D", "darkgray", "firebrick")) 


ggsave(file="case_delta.svg", plot=plot, width=4, height=5)

plot_df <- delta_scaled %>%
  group_by(Statut_clinique) %>%
  summarise(
    mean_delta = mean(delta_mRS, na.rm = TRUE),
    se_delta = sd(delta_mRS, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


plot <- plot_df %>%
  ggplot(aes(x = Statut_clinique, y = mean_delta, colour=Statut_clinique, fill=Statut_clinique)) +
  geom_col(alpha=0.8) +
  geom_segment(
    aes(
      x = Statut_clinique,
      xend = Statut_clinique,
      y = mean_delta - se_delta,
      yend = mean_delta + se_delta
    ),
    linewidth = 3.2,
    lineend = "round", alpha=0.6
  ) +
  labs(
    x = "\n  \nPhysician-based clinical trajectory",
    y = "Δ mRS \n [Standardized Mean ± SE] \n",
    title="Δ mRS & Δ Clinical Trajectory"
  ) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt")) +
  scale_fill_manual(values=c("#00204D", "darkgray", "firebrick")) +
  scale_colour_manual(values=c("#00204D", "darkgray", "firebrick")) 


ggsave(file="case_mRS.svg", plot=plot, width=4, height=5)


plot_df <- delta_scaled %>%
  group_by(Statut_clinique) %>%
  summarise(
    mean_delta = mean(delta_CSF_prot, na.rm = TRUE),
    se_delta = sd(delta_CSF_prot, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


plot <- plot_df %>%
  ggplot(aes(x = Statut_clinique, y = mean_delta, colour=Statut_clinique, fill=Statut_clinique)) +
  geom_col(alpha=0.8) +
  geom_segment(
    aes(
      x = Statut_clinique,
      xend = Statut_clinique,
      y = mean_delta - se_delta,
      yend = mean_delta + se_delta
    ),
    linewidth = 3.2,
    lineend = "round", alpha=0.6
  ) +
  labs(
    x = "\n  \nPhysician-based clinical trajectory",
    y = "Δ Protein CSF \n [Standardized Mean ± SE] \n",
    title="Δ Protein CSF & Δ Clinical Trajectory"
  ) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt")) +
  scale_fill_manual(values=c("#00204D", "darkgray", "firebrick")) +
  scale_colour_manual(values=c("#00204D", "darkgray", "firebrick")) 

plot
ggsave(file="proteincsf.svg", plot=plot, width=4, height=5)




set.seed(123)


library(xgboost)
library(dplyr)

X <- delta_scaled %>%
  select(starts_with("delta_")) %>%
  as.matrix()

y <- delta_scaled$Statut_clinique
y_num <- as.numeric(y) - 1   # worsening=0, stable=1, amelioration=2


patients <- unique(delta_scaled$Code_P)

pred_all <- c()
true_all <- c()

prob_all <- NULL

for (p in patients) {

  train <- delta_scaled %>% filter(Code_P != p)
  test  <- delta_scaled %>% filter(Code_P == p)

  X_train <- train %>%
    select(starts_with("delta_")) %>%
    as.matrix()

  X_test <- test %>%
    select(starts_with("delta_")) %>%
    as.matrix()

  y_train <- as.numeric(train$Statut_clinique) - 1
  y_test  <- as.numeric(test$Statut_clinique) - 1

  dtrain <- xgb.DMatrix(X_train, label = y_train)
  dtest  <- xgb.DMatrix(X_test)

  model <- xgb.train(
    params = list(
      objective = "multi:softprob",
      num_class = 3,
      eval_metric = "mlogloss"
    ),
    data = dtrain,
    nrounds = 150,
    verbose = 0
  )

  pred <- predict(model, dtest)
  pred_mat <- matrix(pred, ncol = 3, byrow = TRUE)
  prob_all <- rbind(prob_all, pred_mat)
  pred_class <- max.col(pred_mat) - 1

  pred_all <- c(pred_all, pred_class)
  true_all <- c(true_all, y_test)
}


mean(pred_all == true_all)
mean(abs(pred_all - true_all))

table(Predicted = pred_all, True = true_all)

cor(pred_all, true_all, method = "spearman")


cm <- table(Predicted = pred_all, True = true_all)

cm_df <- as.data.frame(cm)

ggplot(cm_df, aes(x = True, y = Predicted, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 5) +
  scale_fill_gradient(low = "white", high = "firebrick") +
  labs(
    x = "True class",
    y = "Predicted class",
    fill = "Count"
  ) +
  theme_minimal()


df_pred <- data.frame(
  true = true_all,
  pred = pred_all
)

plot <- ggplot(df_pred, aes(x = true, y = pred)) +
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.5, size=2, shape=1, stroke=2) +
  geom_smooth(method = "lm", se = T, color = "firebrick", fill="firebrick", alpha=0.2) +
  scale_x_continuous(breaks = 0:2, labels = c("worsening", "stable", "improve")) +
  scale_y_continuous(breaks = 0:2, labels = c("worsening", "stable", "improve")) +
  labs(
    x = "\n True class",
    y = "Predicted class \n"
  ) +
   theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt"))

ggsave(file="agree.svg", plot=plot, width=5, height=5)




dtrain_full <- xgb.DMatrix(X, label = y_num)

model_full <- xgb.train(
  params = list(
    objective = "multi:softprob",
    num_class = 3
  ),
  data = dtrain_full,
  nrounds = 150
)

importance <- xgb.importance(model = model_full)

xgb.plot.importance(importance)


to_plot <- importance %>%
  arrange(Importance) %>%
  mutate(Feature=factor(Feature, levels=Feature)) %>%
  ggplot(aes(Feature, Importance)) +
  geom_col(fill="#00204D", colour="#00204D", alpha=0.8) +
  coord_flip() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt")) +
  ylab("\n Feature")+ xlab("Relative Feature Importance (a.u.) \n")


to_plot

ggsave(file="xgboost_im_inc.svg", plot=to_plot, width=5, height=5)




prob_df <- as.data.frame(prob_all)
colnames(prob_df) <- c("worsening", "stable", "improvement")

prob_df$true <- factor(true_all, labels = c("worsening", "stable", "improvement"))

prob_long <- prob_df %>%
  pivot_longer(cols = worsening:improvement,
               names_to = "predicted_class",
               values_to = "probability")

plot <- ggplot(prob_long, aes(x = predicted_class, y = probability, fill = true, colour=true)) +
  geom_boxplot(outliers = F, alpha=0.6, notch=F) +
  labs(
    x = "Predicted class\n",
    y = "\n Predicted probability",
    fill = "True class", color= "True class"
  ) +
  coord_flip() +
   theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt"))  +
  scale_fill_manual(values=c("#00204D", "darkgray", "firebrick")) +
  scale_colour_manual(values=c("#00204D", "darkgray", "firebrick")) 

ggsave(file="xgboost_im_inc.svg", plot=plot, width=5, height=5)


# -------
# Correlations forestplots ---------------

Kappa_df <- read_xlsx(path="../data/Kappa_NFL20250721_clean.xlsx", trim_ws = TRUE)

Kappa_df$IndeK_clean <- as.numeric(Kappa_df$IndeK_clean)
Kappa_df$KappaPlasma <- as.numeric(Kappa_df$KappaPlasma)
Kappa_df$AlbuminPlasma <- as.numeric(Kappa_df$AlbuminPlasma)
Kappa_df$AlbuminCSF <- as.numeric(Kappa_df$AlbuminCSF)
Kappa_df$KappaCSF <- as.numeric(Kappa_df$KappaCSF)
Kappa_df$NfL <- as.numeric(Kappa_df$NfL)

length(unique(Kappa_df$Code_P)) # 50

Kappa_df <- Kappa_df %>% mutate(NLR=`Neutrophiles - nb absolu (sang) clean`/`Lymphocytes - nb absolu (sang) clean`)


Kappa_df$`date symptômes` <- as.Date(as.numeric(Kappa_df$`date symptômes`), origin = "1899-12-30")
Kappa_df$`date dg (T0)`  <- as.Date(Kappa_df$`date dg (T0)`  )
Kappa_df$`Date_sampling`  <- as.Date(Kappa_df$`Date_sampling`  )

library(lubridate)



age_df <- read_xlsx(path="../data/codeDDN_IndexKappaNfL.xlsx", trim_ws = TRUE)
age_df$ddn <- as.Date(age_df$ddn)

Kappa_df <- Kappa_df %>% inner_join(age_df %>% rename("Code_P"="code"))


Kappa_df <- Kappa_df %>%
  mutate(elapsed_sympt = time_length(interval(`date symptômes`, Date_sampling ), "months")) %>%
  mutate(elapsed_dx = time_length(interval(`date dg (T0)`, Date_sampling ), "months"))  %>%
  mutate(age = time_length(interval(ddn, Date_sampling ), "years"))


variables_to_track <- c("Code_P",  "IndeK_clean", "KappaPlasma", "AlbuminPlasma" , "AlbuminCSF" , "KappaCSF", 
                        "NfL", "Case",  "mRS", "abnormal MRI onset", "néo associé" ,"ICI", 
                        "CSF leuco count" , "CSF prot", "BOC_count" , "CRP" , "Neutrophiles - nb absolu (sang) clean" ,
                        "Lymphocytes - nb absolu (sang) clean" , "NLR", "elapsed_sympt", "elapsed_dx", "age")

Kappa_df <- Kappa_df %>% select(all_of(variables_to_track))


# Variables to test against Case
vars_to_test <- c("IndeK_clean", "KappaPlasma", "AlbuminPlasma", "AlbuminCSF", "KappaCSF",
                  "NfL", "mRS", "CSF leuco count", "CSF prot", "BOC_count", "CRP",
                  "Neutrophiles - nb absolu (sang) clean",
                  "Lymphocytes - nb absolu (sang) clean", "NLR", "elapsed_sympt", "elapsed_dx")

results_list <- list()

library(lmerTest)
library(rmcorr)

results_list <- list()

for (v in vars_to_test) {
  
  df_sub <- Kappa_df %>% 
    select(Code_P, Case, age, all_of(v)) %>% 
    drop_na()
  
  if(nrow(df_sub) < 5) next
  
  df_sub$age <- scale(df_sub$age)
  df_sub$Case_scaled <- scale(df_sub$Case)
  df_sub$Var_scaled  <- scale(df_sub[[v]])
  
  df_sub$Case_rank_scaled <- scale(rank(df_sub$Case))
  df_sub$Var_rank_scaled  <- scale(rank(df_sub[[v]]))
  
  # Age-adjusted models
  raw_model <- lmer(Case_scaled ~ Var_scaled + age + (1 | Code_P), data = df_sub)
  raw_coef <- summary(raw_model)$coefficients[2, c(1,2,5)]
  
  rank_model <- lmer(Case_rank_scaled ~ Var_rank_scaled + age + (1 | Code_P), data = df_sub)
  rank_coef <- summary(rank_model)$coefficients[2, c(1,2,5)]
  
  # Age-adjusted residual rmcorr
  df_sub <- df_sub %>%
  group_by(Code_P) %>%
  mutate(
    Case_resid = resid(lm(reformulate("age", response = "Case"), data = cur_data())),
    Var_resid  = resid(lm(reformulate("age", response = v), data = cur_data()))
  ) %>%
  ungroup()
  
  rm <- tryCatch({
    rmcorr(participant = "Code_P",
           measure1 = "Case_resid",
           measure2 = "Var_resid",
           dataset = df_sub)
  }, error = function(e) NULL)
  
  if (!is.null(rm)) {
    rm_vals <- c(r = rm$r, CI_low = rm$CI[1], CI_high = rm$CI[2], p = rm$p)
  } else {
    rm_vals <- c(r = NA, CI_low = NA, CI_high = NA, p = NA)
  }
  
  results_list[[v]] <- data.frame(
    Variable = v,
    Raw_Beta = raw_coef[1],
    Raw_SE   = raw_coef[2],
    Raw_p    = raw_coef[3],
    Rank_Beta = rank_coef[1],
    Rank_SE   = rank_coef[2],
    Rank_p    = rank_coef[3],
    r = rm_vals["r"],
    CI_low = rm_vals["CI_low"],
    CI_high = rm_vals["CI_high"],
    r_p = rm_vals["p"],
    n = nrow(df_sub),
    stringsAsFactors = FALSE
  )
}


# Collapse into one table
results_df <- do.call(rbind, results_list)

# Look at results
print(results_df)


results_df <- results_df %>% filter(Variable != "NLR" & 
                                      Variable !="Lymphocytes - nb absolu (sang) clean" & 
                                      Variable!="Neutrophiles - nb absolu (sang) clean" & 
                                      Variable!="CRP")

plot_df <- results_df %>%
  mutate(
    # 95% CI for raw beta
    Raw_low  = Raw_Beta - 1.96 * Raw_SE,
    Raw_high = Raw_Beta + 1.96 * Raw_SE,
    
    # 95% CI for rank beta
    Rank_low  = Rank_Beta - 1.96 * Rank_SE,
    Rank_high = Rank_Beta + 1.96 * Rank_SE
  ) %>%
  
  # reshape to long format
  select(
    Variable,
    Raw_Beta, Raw_low, Raw_high,
    Rank_Beta, Rank_low, Rank_high,
    r, CI_low, CI_high
  ) %>%
  
  pivot_longer(
    cols = -Variable,
    names_to = "metric",
    values_to = "value"
  )


plot_df <- plot_df %>%
  mutate(
    type = case_when(
      grepl("Raw_", metric) ~ "[1] Mixed Effects Std. Raw beta",
      grepl("Rank_", metric) ~ "[2] Mixed Effects Rank beta",
      metric %in% c("r", "CI_low", "CI_high") ~ "[3] Repeated Measures Correlation"
    ),
    
    stat = case_when(
      grepl("Beta$", metric) ~ "estimate",
      grepl("_low$", metric) | metric == "CI_low" ~ "low",
      grepl("_high$", metric) | metric == "CI_high" ~ "high",
      metric == "r" ~ "estimate"
    )
  ) %>%
  
  select(-metric) %>%
  
  pivot_wider(
    names_from = stat,
    values_from = value
  )


plot <- ggplot(plot_df, aes(x = estimate, y = reorder(Variable, estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2) +
  geom_linerange(aes(xmin = low, xmax = high, y = reorder(Variable, estimate)),
                 size = 2, alpha = 0.7, lineend = "round", color = "#32435d") +
  facet_wrap(~ type, ncol=1) +
  labs(
    x = "\n Effect size \n [Age-adjusted]",
    y = "Feature \n", title="CASE Score"
  ) +
theme_minimal() +
  coord_cartesian(xlim=c(-0.5,1)) +
geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dotted", alpha = 0.3) +
  geom_vline(xintercept = c(-0.4, 0.4), linetype = "dotted", alpha = 0.3) +
    geom_vline(xintercept = c(-0.6, 0.6), linetype = "dotted", alpha = 0.3) +
    geom_vline(xintercept = c(-0.8, 0.8), linetype = "dotted", alpha = 0.3) +
      geom_vline(xintercept = c(-1.0, 1.0), linetype = "dotted", alpha = 0.3) +
  theme(axis.text.y = element_text(size = 10, face = "bold", hjust = 1),
        axis.ticks.y = element_blank(),
        legend.position = "bottom") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 12, angle = 0, vjust = -0.1),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0, size = 14),
        plot.subtitle = element_text(hjust = 0, size = 11),
        plot.margin = margin(10, 10, 10, 10, "pt")) 

ggsave(file="plot1.svg", plot=plot, width=7, height=9)







# Variables to test against Case
vars_to_test <- c("IndeK_clean", "KappaPlasma", "AlbuminPlasma", "AlbuminCSF", "KappaCSF",
                  "NfL", "Case", "CSF leuco count", "CSF prot", "BOC_count", "CRP",
                  "Neutrophiles - nb absolu (sang) clean",
                  "Lymphocytes - nb absolu (sang) clean", "NLR", "elapsed_sympt", "elapsed_dx")

results_list <- list()

library(lmerTest)
library(rmcorr)

results_list <- list()

for (v in vars_to_test) {
  
  df_sub <- Kappa_df %>% 
    select(Code_P, mRS, age, all_of(v)) %>% 
    drop_na()
  
  if (nrow(df_sub) < 5) next
  
  # Scale variables
  df_sub$age <- scale(df_sub$age)
  df_sub$mRS_scaled <- scale(df_sub$mRS)
  df_sub$Var_scaled  <- scale(df_sub[[v]])
  
  df_sub$mRS_rank_scaled <- scale(rank(df_sub$mRS))
  df_sub$Var_rank_scaled  <- scale(rank(df_sub[[v]]))
  
  # --- Age-adjusted mixed models ---
  raw_model <- lmer(mRS_scaled ~ Var_scaled + age + (1 | Code_P), data = df_sub)
  raw_coef <- summary(raw_model)$coefficients[2, c(1,2,5)]
  
  rank_model <- lmer(mRS_rank_scaled ~ Var_rank_scaled + age + (1 | Code_P), data = df_sub)
  rank_coef <- summary(rank_model)$coefficients[2, c(1,2,5)]
  
  # --- Age-adjusted residuals for repeated-measures correlation ---
  df_sub <- df_sub %>%
    group_by(Code_P) %>%
    mutate(
      mRS_resid = resid(lm(reformulate("age", response = "mRS"), data = cur_data())),
      Var_resid = resid(lm(reformulate("age", response = v), data = cur_data()))
    ) %>%
    ungroup()
  
  rm <- tryCatch({
    rmcorr(participant = "Code_P",
           measure1 = "mRS_resid",
           measure2 = "Var_resid",
           dataset = df_sub)
  }, error = function(e) NULL)
  
  if (!is.null(rm)) {
    rm_vals <- c(r = rm$r, CI_low = rm$CI[1], CI_high = rm$CI[2], p = rm$p)
  } else {
    rm_vals <- c(r = NA, CI_low = NA, CI_high = NA, p = NA)
  }
  
  # --- Store results ---
  results_list[[v]] <- data.frame(
    Variable = v,
    Raw_Beta = raw_coef[1],
    Raw_SE   = raw_coef[2],
    Raw_p    = raw_coef[3],
    Rank_Beta = rank_coef[1],
    Rank_SE   = rank_coef[2],
    Rank_p    = rank_coef[3],
    r = rm_vals["r"],
    CI_low = rm_vals["CI_low"],
    CI_high = rm_vals["CI_high"],
    r_p = rm_vals["p"],
    n = nrow(df_sub),
    stringsAsFactors = FALSE
  )
}



# Collapse into one table
results_df <- do.call(rbind, results_list)

# Look at results
print(results_df)






results_df <- results_df %>% filter(Variable != "NLR" & 
                                      Variable !="Lymphocytes - nb absolu (sang) clean" & 
                                      Variable!="Neutrophiles - nb absolu (sang) clean" & 
                                      Variable!="CRP")

plot_df <- results_df %>%
  mutate(
    # 95% CI for raw beta
    Raw_low  = Raw_Beta - 1.96 * Raw_SE,
    Raw_high = Raw_Beta + 1.96 * Raw_SE,
    
    # 95% CI for rank beta
    Rank_low  = Rank_Beta - 1.96 * Rank_SE,
    Rank_high = Rank_Beta + 1.96 * Rank_SE
  ) %>%
  
  # reshape to long format
  select(
    Variable,
    Raw_Beta, Raw_low, Raw_high,
    Rank_Beta, Rank_low, Rank_high,
    r, CI_low, CI_high
  ) %>%
  
  pivot_longer(
    cols = -Variable,
    names_to = "metric",
    values_to = "value"
  )


plot_df <- plot_df %>%
  mutate(
    type = case_when(
      grepl("Raw_", metric) ~ "[1] Mixed Effects Std. Raw beta",
      grepl("Rank_", metric) ~ "[2] Mixed Effects Rank beta",
      metric %in% c("r", "CI_low", "CI_high") ~ "[3] Repeated Measures Correlation"
    ),
    
    stat = case_when(
      grepl("Beta$", metric) ~ "estimate",
      grepl("_low$", metric) | metric == "CI_low" ~ "low",
      grepl("_high$", metric) | metric == "CI_high" ~ "high",
      metric == "r" ~ "estimate"
    )
  ) %>%
  
  select(-metric) %>%
  
  pivot_wider(
    names_from = stat,
    values_from = value
  )


plot <- ggplot(plot_df, aes(x = estimate, y = reorder(Variable, estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2) +
  geom_linerange(aes(xmin = low, xmax = high, y = reorder(Variable, estimate)),
                 size = 2, alpha = 0.7, lineend = "round", color = "firebrick") +
  facet_wrap(~ type, ncol=1) +
  labs(
    x = "\n Effect size \n [Age-adjusted]",
    y = "Feature \n", title="mRS Score"
  ) +
theme_minimal() +
 coord_cartesian(xlim=c(-0.6,1)) +
geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dotted", alpha = 0.3) +
  geom_vline(xintercept = c(-0.4, 0.4), linetype = "dotted", alpha = 0.3) +
    geom_vline(xintercept = c(-0.6, 0.6), linetype = "dotted", alpha = 0.3) +
    geom_vline(xintercept = c(-0.8, 0.8), linetype = "dotted", alpha = 0.3) +
      geom_vline(xintercept = c(-1.0, 1.0), linetype = "dotted", alpha = 0.3) +
  theme(axis.text.y = element_text(size = 10, face = "bold", hjust = 1),
        axis.ticks.y = element_blank(),
        legend.position = "bottom") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 12, angle = 0, vjust = -0.1),
        axis.title.x = element_text(size = 12, vjust = -0.5),
        axis.title.y = element_text(size = 12, vjust = -0.5),
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0, size = 14),
        plot.subtitle = element_text(hjust = 0, size = 11),
        plot.margin = margin(10, 10, 10, 10, "pt")) 

plot
ggsave(file="plot2.svg", plot=plot, width=7, height=9)


# --------------
