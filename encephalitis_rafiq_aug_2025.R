

library(readxl)
library(tidyverse)
library(data.table)

# Baseline summary info -----------

Kappa_df <- read_xlsx(path="../data/Kappa_NFL20250721_clean.xlsx", trim_ws = TRUE)

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
                        "CSF leuco count" , "CSF prot", "BOC_count" , "CRP" , "Neutrophiles - nb absolu (sang) clean" ,
                        "Lymphocytes - nb absolu (sang) clean" , "NLR")

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
                        "CSF leuco count" , "CSF prot", "BOC_count" , "CRP" , "Neutrophiles - nb absolu (sang) clean" ,
                        "Lymphocytes - nb absolu (sang) clean" , "NLR")


names(firsts_df)


to_plot <- firsts_df %>% mutate(elapsed_dx=elapsed_dx/30.44) %>% ggplot(aes(y=`elapsed_dx`, x="")) +
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

ggsave(file="Diagnosis_first.svg", plot=to_plot, width=4, height=5)



to_plot <- firsts_df %>% mutate(elapsed_sympt=elapsed_sympt/30.44) %>% ggplot(aes(y=`elapsed_sympt`, x="")) +
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

ggsave(file="elapsed_sympt_first.svg", plot=to_plot, width=4, height=5)




to_plot <- ggplot(firsts_df, aes(y=`NLR`, x="")) +
  geom_boxplot( color="#747475", alpha=0, notch = TRUE, width=0.5, linewidth=1.5) +
  geom_jitter(width=0.2, alpha=0.6, size=2, shape=1, stroke=2, color="#00204D") +
    labs(x = "First evaluation", 
       y = "Neutrophiles-to-lymphocytes ratio (NLR)\n", 
       title="NLR") +
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

ggsave(file="nlr_Serum_first.svg", plot=to_plot, width=4, height=5)

to_plot <- ggplot(firsts_df, aes(y=`Lymphocytes - nb absolu (sang) clean`, x="")) +
  geom_boxplot( color="#747475", alpha=0, notch = TRUE, width=0.5, linewidth=1.5) +
  geom_jitter(width=0.2, alpha=0.6, size=2, shape=1, stroke=2, color="#00204D") +
    labs(x = "First evaluation", 
       y = "Lymphocytes Serum (10⁹ cells/L)\n", 
       title="Lymphocytes Serum") +
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

ggsave(file="Lymphocytes_Serum_first.svg", plot=to_plot, width=4, height=5)


to_plot <- ggplot(firsts_df, aes(y=`Neutrophiles - nb absolu (sang) clean`, x="")) +
  geom_boxplot( color="#747475", alpha=0, notch = TRUE, width=0.5, linewidth=1.5) +
  geom_jitter(width=0.2, alpha=0.6, size=2, shape=1, stroke=2, color="#00204D") +
    labs(x = "First evaluation", 
       y = "Neutrophiles Serum (10⁹ cells/L)\n", 
       title="Neutrophiles Serum") +
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

ggsave(file="Neutrophiles_Serum_first.svg", plot=to_plot, width=4, height=5)

to_plot <- ggplot(firsts_df, aes(y=CRP, x="")) +
  geom_boxplot( color="#747475", alpha=0, notch = TRUE, width=0.5, linewidth=1.5) +
  geom_jitter(width=0.2, alpha=0.6, size=2, shape=1, stroke=2, color="#00204D") +
    labs(x = "First evaluation", 
       y = "C-reactive protein (CRP, mg/L)\n", 
       title="CRP") +
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

ggsave(file="crp_first.svg", plot=to_plot, width=4, height=5)


to_plot <- ggplot(firsts_df, aes(y=BOC_count, x="")) +
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

ggsave(file="OBC_first.svg", plot=to_plot, width=4, height=5)



to_plot <- ggplot(firsts_df, aes(y=`CSF prot`, x="")) +
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

ggsave(file="protein_CSF_first.svg", plot=to_plot, width=4, height=5)



to_plot <- ggplot(firsts_df, aes(y=`CSF leuco count`, x="")) +
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

ggsave(file="Leukocytes_CSF_first.svg", plot=to_plot, width=4, height=5)



to_plot <- ggplot(firsts_df, aes(y=mRS, x="")) +
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

ggsave(file="mrs_first.svg", plot=to_plot, width=4, height=5)


to_plot <- ggplot(firsts_df, aes(y=Case, x="")) +
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

ggsave(file="case_first.svg", plot=to_plot, width=4, height=5)




to_plot <- ggplot(firsts_df, aes(y=NfL, x="")) +
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

ggsave(file="NfL_first.svg", plot=to_plot, width=4, height=5)


to_plot <- ggplot(firsts_df, aes(y=AlbuminCSF, x="")) +
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

ggsave(file="AlbuminCSF_first.svg", plot=to_plot, width=4, height=5)


to_plot <- ggplot(firsts_df, aes(y=AlbuminPlasma, x="")) +
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

ggsave(file="AlbuminPlasma_first.svg", plot=to_plot, width=4, height=5)



to_plot <- ggplot(firsts_df, aes(y=KappaCSF, x="")) +
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

ggsave(file="KappaCSF_first.svg", plot=to_plot, width=4, height=5)


to_plot <- ggplot(firsts_df, aes(y=KappaPlasma, x="")) +
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

ggsave(file="KappaPlasma_first.svg", plot=to_plot, width=4, height=5)


to_plot <- ggplot(firsts_df, aes(y=KappaPlasma, x="")) +
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

ggsave(file="KappaPlasma_first.svg", plot=to_plot, width=4, height=5)



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

Kappa_df <- Kappa_df %>%
  mutate(elapsed_sympt = time_length(interval(`date symptômes`, Date_sampling ), "months")) %>%
  mutate(elapsed_dx = time_length(interval(`date dg (T0)`, Date_sampling ), "months")) 



variables_to_track <- c("Code_P",  "IndeK_clean", "KappaPlasma", "AlbuminPlasma" , "AlbuminCSF" , "KappaCSF", 
                        "NfL", "Case",  "mRS", "abnormal MRI onset", "néo associé" ,"ICI", 
                        "CSF leuco count" , "CSF prot", "BOC_count" , "CRP" , "Neutrophiles - nb absolu (sang) clean" ,
                        "Lymphocytes - nb absolu (sang) clean" , "NLR", "elapsed_sympt", "elapsed_dx")

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
    select(Code_P, Case, all_of(v)) %>% 
    drop_na()
  
  if(nrow(df_sub) < 5) next
  
  df_sub$Case_scaled <- scale(df_sub$Case)
  df_sub$Var_scaled  <- scale(df_sub[[v]])
  
  df_sub$Case_rank_scaled <- scale(rank(df_sub$Case))
  df_sub$Var_rank_scaled  <- scale(rank(df_sub[[v]]))
  
  # Raw model
  raw_model <- lmer(Case_scaled ~ Var_scaled + (1 | Code_P), data = df_sub)
  raw_coef <- summary(raw_model)$coefficients[2, c(1,2,5)]
  
  # Rank model
  rank_model <- lmer(Case_rank_scaled ~ Var_rank_scaled + (1 | Code_P), data = df_sub)
  rank_coef <- summary(rank_model)$coefficients[2, c(1,2,5)]
  
  # Rmcorr
  rm <- tryCatch({
    rmcorr(participant = "Code_P",
           measure1 = "Case",
           measure2 = v,
           dataset = df_sub)
  }, error = function(e) NULL)
  
  if (!is.null(rm)) {
    rm_vals <- c(r = rm$r, CI_low = rm$CI[1], CI_high = rm$CI[2], p = rm$p)
  } else {
    rm_vals <- c(r = NA, CI_low = NA, CI_high = NA, p = NA)
  }
  
  # Store results
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


rownames(results_df) <- NULL

results_df$Raw_Beta <- round(results_df$Raw_Beta, 2)
results_df$Raw_SE <- round(results_df$Raw_SE, 2)
results_df$Rank_Beta <- round(results_df$Rank_Beta, 2)
results_df$Rank_SE <- round(results_df$Rank_SE, 2)
results_df$r <- round(results_df$r, 2)
results_df$CI_low <- round(results_df$CI_low, 2)
results_df$CI_high <- round(results_df$CI_high, 2)

results_df %>%
  mutate(Raw_Beta=paste0(Raw_Beta, paste0(" ±", Raw_SE))) %>% select(-Raw_SE) %>%
  mutate(Rank_Beta=paste0(Rank_Beta, paste0(" ±", Rank_SE))) %>% select(-Rank_SE) %>%
  mutate(r=paste0(r, paste0(" [", paste0(CI_low , paste0(" , ", paste0(CI_high, "]")))))) %>% select(-CI_low , -CI_high)
  






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
    select(Code_P, mRS, all_of(v)) %>% 
    drop_na()
  
  if(nrow(df_sub) < 5) next
  
  df_sub$mRS_scaled <- scale(df_sub$mRS)
  df_sub$Var_scaled  <- scale(df_sub[[v]])
  
  df_sub$mRS_rank_scaled <- scale(rank(df_sub$mRS))
  df_sub$Var_rank_scaled  <- scale(rank(df_sub[[v]]))
  
  # Raw model
  raw_model <- lmer(mRS_scaled ~ Var_scaled + (1 | Code_P), data = df_sub)
  raw_coef <- summary(raw_model)$coefficients[2, c(1,2,5)]
  
  # Rank model
  rank_model <- lmer(mRS_rank_scaled ~ Var_rank_scaled + (1 | Code_P), data = df_sub)
  rank_coef <- summary(rank_model)$coefficients[2, c(1,2,5)]
  
  # Rmcorr
  rm <- tryCatch({
    rmcorr(participant = "Code_P",
           measure1 = "mRS",
           measure2 = v,
           dataset = df_sub)
  }, error = function(e) NULL)
  
  if (!is.null(rm)) {
    rm_vals <- c(r = rm$r, CI_low = rm$CI[1], CI_high = rm$CI[2], p = rm$p)
  } else {
    rm_vals <- c(r = NA, CI_low = NA, CI_high = NA, p = NA)
  }
  
  # Store results
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


rownames(results_df) <- NULL

results_df$Raw_Beta <- round(results_df$Raw_Beta, 2)
results_df$Raw_SE <- round(results_df$Raw_SE, 2)
results_df$Rank_Beta <- round(results_df$Rank_Beta, 2)
results_df$Rank_SE <- round(results_df$Rank_SE, 2)
results_df$r <- round(results_df$r, 2)
results_df$CI_low <- round(results_df$CI_low, 2)
results_df$CI_high <- round(results_df$CI_high, 2)

results_df %>%
  mutate(Raw_Beta=paste0(Raw_Beta, paste0(" ±", Raw_SE))) %>% select(-Raw_SE) %>%
  mutate(Rank_Beta=paste0(Rank_Beta, paste0(" ±", Rank_SE))) %>% select(-Rank_SE) %>%
  mutate(r=paste0(r, paste0(" [", paste0(CI_low , paste0(" , ", paste0(CI_high, "]")))))) %>% select(-CI_low , -CI_high)
  


Kappa_df 





names(Kappa_df)

target_var <- "elapsed_sympt"
label <- "Disease duration"

to_plot <- Kappa_df %>% # mutate(elapsed_sympt=elapsed_sympt/30.44) %>%
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

ggsave(file=paste0(label, "_all.svg"), plot=ploted, width=5, height=7)


options(scipen = 999)


# List of candidate variables (numeric only, excluding Case)
vars <- Kappa_df %>%
  select(-Code_P) %>%
  select(where(is.numeric)) %>%
  select(-mRS) %>%
  names()

# Run Spearman correlation with Case for each variable
# cor_results <- map_dfr(vars, function(var) {
#   df <- Kappa_df %>%
#     select(Case, !!sym(var)) %>%
#     drop_na()
#   
#   test <- cor.test(df$Case, df[[var]], method = "spearman")
#   
#   tibble(
#     Variable = var,
#     Spearman_rho = unname(test$estimate),
#     p_value = test$p.value
#   )
# })


cor_results <- map_dfr(vars, function(var) {
  df <- Kappa_df %>%
    select(mRS, !!sym(var)) %>%
    drop_na()
  
  # Compute 5% and 95% quantiles for the variable
  q_low  <- quantile(df[[var]], 0.05, na.rm = TRUE)
  q_high <- quantile(df[[var]], 0.95, na.rm = TRUE)
  
  # Keep only values within quantile range
  df_trimmed <- df %>%
    filter(df[[var]] >= q_low & df[[var]] <= q_high)
  
  # Run Spearman correlation on trimmed data
  test <- cor.test(df_trimmed$mRS, df_trimmed[[var]], method = "spearman")
  
  tibble(
    Variable = var,
    Spearman_rho = unname(test$estimate),
    p_value = test$p.value,
    n_used = nrow(df_trimmed)   # track sample size after trimming
  )
})



# Show results sorted by strongest correlation
cor_results <- cor_results %>%
  arrange(desc(abs(Spearman_rho)))

data.frame(cor_results)



# Enter the data in your specified order
variables <- c(
  "κ FLC Index", "κ FLC Plasma", "Albumin Plasma", "Albumin CSF", "κ FLC CSF",
  "NfL", "CASE/mRS", "Abnormal MRI", "Neoplasia", "ICI", "Leuko CSF", "CSF Protein",
  "OCB count", "CRP", "Neutrophiles", "Lymphocytes", "NLR", "Symptom duration", "Diagnosis duration"
)

case_corr <- c(
  0.03833395 , 0.08433665, -0.19190672 , 0.18753894 , 0.24460075 ,
  0.22547683 , 0.71201251 , 0.08834478 , 0.18806050 , 0.21956403 , 0.02921881 , 0.19970999 ,
  0.42053507 , 0.02058906 , 0.11301147 , -0.15598102 , 0.14754256 , -0.20242644 , -0.17180055 
)

mrs_corr <- c(
  0.13544338 , 0.03127982 , -0.28511294 , 0.03461965 , 0.27831722 ,
  0.34853861 , 0.65876671 , 0.16194733 , 0.26927762 , 0.27889101 , 0.03847854 , 0.10370737 ,
  0.39751244 , -0.12801217 , 0.10217663 , -0.03492053 , 0.01254586 , -0.30453936 , -0.23976206 
)

# Adjust variable names for display in plot (short, R-friendly)
plot_vars <- c(
  "κ FLC Index", "κ FLC Plasma", "Albumin Plasma", "Albumin CSF",
  "κ FLC CSF", "NfL", "CASE/mRS", "Abnormal MRI", "Neoplasia", "ICI",
  "Leuko CSF", "CSF Protein", "OCB count", "CRP", "Neutrophiles", "Lymphocytes",
  "NLR", "Symptom duration", "Diagnosis duration"
)

# Organize into a data frame for melting
heatmap_mat <- rbind(CASE = case_corr, mRS = mrs_corr)
colnames(heatmap_mat) <- plot_vars
rownames(heatmap_mat) <- c("CASE", "mRS")

heatmap_long <- melt(heatmap_mat)
colnames(heatmap_long) <- c("Scale", "Variable", "Correlation")


# Plot heatmap
to_plot <- ggplot(heatmap_long, aes(x = Variable, y = Scale, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Correlation)), size = 3) +
  scale_fill_gradient2(low = "#00204D", mid = "white", high = "firebrick", midpoint = 0, 
                       name = "Spearman\nrho") +
  labs(x = "Covariate", y = "CASE vs.mRS") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave(file="correlations.svg", plot=to_plot, width=10, height=2) 


# -----------

# Predict future/next CASE/mRS/Global Clinician appreciation -----------



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
  mutate(elapsed_sympt = time_length(interval(`date symptômes`, Date_sampling ), "months")) %>%
  mutate(elapsed_dx = time_length(interval(`date dg (T0)`, Date_sampling ), "months"))  %>%
  mutate(elapsed_dx = time_length(interval(`date dg (T0)`, Date_sampling ), "months"))  


names(Kappa_df)

variables_to_track <- c("Code_P",  "IndeK_clean", "KappaPlasma", "AlbuminPlasma" , "AlbuminCSF" , "KappaCSF", 
                        "NfL", "Case",  "mRS",  "Date_sampling",
                        "CSF leuco count" , "CSF prot", "BOC_count" , "CRP" , "Neutrophiles - nb absolu (sang) clean" ,
                        "Lymphocytes - nb absolu (sang) clean" , "NLR", "Statut clinique")

Kappa_df <- Kappa_df %>% select(all_of(variables_to_track))

names(Kappa_df)

Kappa_df <- Kappa_df %>% rename("CSF_leuco_count"="CSF leuco count") %>%
  rename("CSF_prot"="CSF prot") %>%
  rename("Neutrophiles"="Neutrophiles - nb absolu (sang) clean") %>%
  rename("Lymphocytes"="Lymphocytes - nb absolu (sang) clean") %>%
  rename("Statut_clinique"="Statut clinique")

library(dplyr)
library(lme4)
library(lmerTest)
library(purrr)
library(tidyr)

# Outcome can be "Case" or "mRS"
outcome_var <- "Case"

# Function to fit lagged models
fit_lagged_models <- function(df, outcome_var, covariates, scale_data = TRUE) {
  
  df <- df %>% arrange(Code_P, Date_sampling)
  
  # create lagged outcome
  df <- df %>%
    group_by(Code_P) %>%
    mutate(lag_outcome = lag(!!sym(outcome_var))) %>%
    ungroup()
  
  results <- map_dfr(covariates, function(var) {
    
    df_lag <- df %>%
      group_by(Code_P) %>%
      mutate(!!paste0("lag_", var) := lag(!!sym(var))) %>%
      ungroup() %>%
      drop_na(lag_outcome, !!sym(paste0("lag_", var))) # remove NAs
    
    # Optionally scale the lagged covariate (and outcome if needed)
    if (scale_data) {
      df_lag <- df_lag %>%
        mutate(
          !!paste0("lag_", var) := as.numeric(scale(!!sym(paste0("lag_", var)))),
          lag_outcome = as.numeric(scale(lag_outcome))
        )
    }
    
    # Unadjusted model
    model_unadj <- lmer(
      formula(paste0(outcome_var, " ~ lag_", var, " + (1 | Code_P)")),
      data = df_lag
    )
    
    # Adjusted model (lagged outcome)
    model_adj <- lmer(
      formula(paste0(outcome_var, " ~ lag_outcome + lag_", var, " + (1 | Code_P)")),
      data = df_lag
    )
    
    tibble(
      Covariate = var,
      Unadj_Estimate = summary(model_unadj)$coefficients[2,1],
      Unadj_SE       = summary(model_unadj)$coefficients[2,2],
      Unadj_p        = summary(model_unadj)$coefficients[2,5],
      Adj_Estimate   = summary(model_adj)$coefficients[3,1],
      Adj_SE         = summary(model_adj)$coefficients[3,2],
      Adj_p          = summary(model_adj)$coefficients[3,5]
    )
    
  })
  
  return(results)
}

# Run for selected variables
covariates_to_test <- c("IndeK_clean", "KappaPlasma", "AlbuminPlasma" , "AlbuminCSF" , 
                        "KappaCSF", "NfL", "BOC_count" , "CRP", "CSF_leuco_count" , "CSF_prot", 
                        "Neutrophiles" ,
                        "Lymphocytes" , "NLR")

lagged_results <- fit_lagged_models(Kappa_df, outcome_var, covariates_to_test)

lagged_results


#    Covariate       Unadj_Estimate Unadj_SE Unadj_p Adj_Estimate Adj_SE  Adj_p
#    <chr>                    <dbl>    <dbl>   <dbl>        <dbl>  <dbl>  <dbl>
#  1 IndeK_clean            -0.0529    0.236  0.823        0.291   0.124 0.0217
#  2 KappaPlasma             0.0866    0.396  0.828        0.128   0.221 0.563 
#  3 AlbuminPlasma           0.0371    0.229  0.872        0.0172  0.132 0.897 
#  4 AlbuminCSF             -0.208     0.327  0.528       -0.185   0.128 0.154 
#  5 KappaCSF                0.507     0.347  0.148        0.462   0.222 0.0404
#  6 NfL                     0.654     0.339  0.0578       0.266   0.233 0.258 
#  7 BOC_count               0.239     0.201  0.239        0.273   0.136 0.0485
#  8 CRP                     0.327     0.280  0.247        0.0812  0.130 0.535 
#  9 CSF_leuco_count         0.156     0.197  0.431        0.198   0.172 0.253 
# 10 CSF_prot                0.224     0.270  0.410       -0.176   0.144 0.226 
# 11 Neutrophiles            0.673     0.426  0.119        0.0501  0.271 0.854 
# 12 Lymphocytes             0.317     0.453  0.487        0.329   0.261 0.211 
# 13 NLR                     0.163     0.419  0.698       -0.235   0.271 0.390






# Outcome can be "Case" or "mRS"
outcome_var <- "mRS"

# Function to fit lagged models
fit_lagged_models <- function(df, outcome_var, covariates, scale_data=TRUE) {
  
  df <- df %>% arrange(Code_P, Date_sampling)
  
  # create lagged outcome
  df <- df %>%
    group_by(Code_P) %>%
    mutate(lag_outcome = lag(!!sym(outcome_var))) %>%
    ungroup()
  
  results <- map_dfr(covariates, function(var) {
    
    df_lag <- df %>%
      group_by(Code_P) %>%
      mutate(!!paste0("lag_", var) := lag(!!sym(var))) %>%
      ungroup() %>%
      drop_na(lag_outcome, !!sym(paste0("lag_", var)))  # remove NAs
    
    # Optionally scale the lagged covariate (and outcome if needed)
    if (scale_data) {
      df_lag <- df_lag %>%
        mutate(
          !!paste0("lag_", var) := as.numeric(scale(!!sym(paste0("lag_", var)))),
          lag_outcome = as.numeric(scale(lag_outcome))
        )
    }
    

    
    # Unadjusted model
    model_unadj <- lmer(
      formula(paste0(outcome_var, " ~ lag_", var, " + (1 | Code_P)")),
      data = df_lag
    )
    
    # Adjusted model (lagged outcome)
    model_adj <- lmer(
      formula(paste0(outcome_var, " ~ lag_outcome + lag_", var, " + (1 | Code_P)")),
      data = df_lag
    )
    
    tibble(
      Covariate = var,
      Unadj_Estimate = summary(model_unadj)$coefficients[2,1],
      Unadj_SE       = summary(model_unadj)$coefficients[2,2],
      Unadj_p        = summary(model_unadj)$coefficients[2,5],
      Adj_Estimate   = summary(model_adj)$coefficients[3,1],
      Adj_SE         = summary(model_adj)$coefficients[3,2],
      Adj_p          = summary(model_adj)$coefficients[3,5]
    )
    
  })
  
  return(results)
}

# Run for selected variables
covariates_to_test <- c("IndeK_clean", "KappaPlasma", "AlbuminPlasma" , "AlbuminCSF" , 
                        "KappaCSF", "NfL", "BOC_count" , "CRP", "CSF_leuco_count" , "CSF_prot", 
                        "Neutrophiles" ,
                        "Lymphocytes" , "NLR")

lagged_results <- fit_lagged_models(Kappa_df, outcome_var, covariates_to_test)


lagged_results


#    Covariate       Unadj_Estimate Unadj_SE Unadj_p Adj_Estimate Adj_SE     Adj_p
#    <chr>                    <dbl>    <dbl>   <dbl>        <dbl>  <dbl>     <dbl>
#  1 IndeK_clean            0.0261    0.0653   0.692      0.234   0.0525 0.0000268
#  2 KappaPlasma            0.00293   0.0811   0.971      0.156   0.0644 0.0179   
#  3 AlbuminPlasma          0.0182    0.0607   0.766     -0.0358  0.0631 0.572    
#  4 AlbuminCSF             0.00691   0.0894   0.939     -0.0517  0.0713 0.471    
#  5 KappaCSF               0.0332    0.0600   0.583      0.179   0.0574 0.00262  
#  6 NfL                    0.0249    0.0548   0.652      0.0477  0.0615 0.441    
#  7 BOC_count              0.0101    0.0566   0.860      0.191   0.0599 0.00204  
#  8 CRP                    0.136     0.106    0.204      0.00878 0.0579 0.880    
#  9 CSF_leuco_count        0.0846    0.0559   0.139      0.0992  0.0580 0.0927   
# 10 CSF_prot              -0.0166    0.0831   0.843     -0.0342  0.0738 0.645    
# 11 Neutrophiles           0.0824    0.0800   0.310     -0.0194  0.0712 0.786    
# 12 Lymphocytes            0.0852    0.0959   0.379      0.0760  0.0685 0.271    
# 13 NLR                   -0.0279    0.0706   0.696     -0.0714  0.0689 0.304  






library(ordinal)

# Assume OutcomeFlag is factor ordered: amélioration < Stable < aggravation
Kappa_df$Statut_clinique <- factor(Kappa_df$Statut_clinique, levels = c("amélioration", "stable", "aggravation"), ordered = TRUE)


# Lagged covariates
test_df <- Kappa_df %>%
  arrange(Code_P, Date_sampling) %>%
  group_by(Code_P) %>%
  mutate(
    lag_Statut_clinique = lag(Statut_clinique),
    lag_BOC_count  = lag(BOC_count )
    # add other covariates as needed
  ) %>% 
  mutate(lag_BOC_count=scale(lag_BOC_count)) %>%
  ungroup() %>%
  drop_na(lag_Statut_clinique, lag_BOC_count )

# Ordinal mixed model with random intercept for patient
model_flag <- clmm(Statut_clinique ~ lag_BOC_count  + (1 | Code_P), data = test_df)
summary(model_flag)


# Adjusted for prior outcome
model_flag_adj <- clmm(Statut_clinique ~ lag_Statut_clinique + lag_BOC_count  + (1 | Code_P), data = test_df)
summary(model_flag_adj)



# Lagged covariates
test_df <- Kappa_df %>%
  arrange(Code_P, Date_sampling) %>%
  group_by(Code_P) %>%
  mutate(
    lag_Statut_clinique = lag(Statut_clinique),
    lag_AlbuminCSF   = lag(AlbuminCSF  )
    # add other covariates as needed
  ) %>% 
  mutate(lag_AlbuminCSF =scale(lag_AlbuminCSF )) %>%
  ungroup() %>%
  drop_na(lag_Statut_clinique, lag_AlbuminCSF  )

# Ordinal mixed model with random intercept for patient
model_flag <- clmm(Statut_clinique ~ lag_AlbuminCSF   + (1 | Code_P), data = test_df)
summary(model_flag)


# Adjusted for prior outcome
model_flag_adj <- clmm(Statut_clinique ~ lag_Statut_clinique + lag_AlbuminCSF   + (1 | Code_P), data = test_df)
summary(model_flag_adj)




# ---------


# Predict categorical ---------------------------

# library(dplyr)
# library(lme4)
# library(irr)
# library(purrr)
# 
# #---------- Parameters
# tol <- 2   # tolerance for "stable" changes
# 
# # Covariates to test
# covariates_to_test <- c("IndeK_clean")
# 
# outcome_var <- "Case"   # could also loop for mRS
# 
#  # ---------- Helper Functions
# # Categorize change as improved/stable/worsened
# case_change_flag <- function(next_case, lag_case, tol = 0.25) {
#   diff <- next_case - lag_case
#   dplyr::case_when(
#     diff > tol  ~ "worsened",
#     diff < -tol ~ "improved",
#     TRUE        ~ "stable"
#   )
# }
# 
# # Fit lagged models and compute predicted flags
# fit_predict_flags <- function(df, outcome_var, covariate, tol = 0.25) {
#   
#   df <- df %>% arrange(Code_P, Date_sampling)
#   
#   # Create lagged outcome
#   df <- df %>%
#     group_by(Code_P) %>%
#     mutate(
#       lag_outcome = lag(!!sym(outcome_var)),
#       lag_cov = lag(!!sym(covariate))
#     ) %>%
#     ungroup() %>%
#     drop_na(lag_outcome, lag_cov)
#   
#   # Scale lagged covariate (optional)
#   df <- df %>%
#     mutate(
#       lag_cov_scaled = as.numeric(scale(lag_cov)),
#       lag_outcome_scaled = as.numeric(scale(lag_outcome))
#     )
#   
#   # --- Fit models
#   # Unadjusted: lagged covariate only
#   model_unadj <- lmer(
#     formula(paste0(outcome_var, " ~ lag_outcome_scaled + (1 | Code_P)")),
#     data = df
#   )
#   
#   # Adjusted: lagged outcome + lagged covariate
#   model_adj <- lmer(
#     formula(paste0(outcome_var, " ~ lag_outcome_scaled + lag_cov_scaled + (1 | Code_P)")),
#     data = df
#   )
#   
#   # ---- Predicted next CASE
#   df <- df %>%
#     mutate(
#       pred_unadj = predict(model_unadj, re.form = NULL),
#       pred_adj   = predict(model_adj, re.form = NULL),
#       observed_flag = case_change_flag(!!sym(outcome_var), lag_outcome, tol),
#       pred_flag_unadj = case_change_flag(pred_unadj, lag_outcome, tol),
#       pred_flag_adj   = case_change_flag(pred_adj, lag_outcome, tol)
#     )
#   
#   # ---- Agreement (Cohen's kappa
#   kappa_unadj <- irr::kappa2(data.frame(df$observed_flag, df$pred_flag_unadj))$value
#   kappa_adj   <- irr::kappa2(data.frame(df$observed_flag, df$pred_flag_adj))$value
#   
#   # ---- Collect coefficients
#   tibble(
#     Covariate = covariate,
#     Unadj_Estimate = summary(model_unadj)$coefficients[2,1],
#     Unadj_SE       = summary(model_unadj)$coefficients[2,2],
#     Unadj_p        = summary(model_unadj)$coefficients[2,5],
#     Adj_Estimate   = summary(model_adj)$coefficients[3,1],
#     Adj_SE         = summary(model_adj)$coefficients[3,2],
#     Adj_p          = summary(model_adj)$coefficients[3,5],
#     Kappa_Unadj    = kappa_unadj,
#     Kappa_Adj      = kappa_adj
#   )
#   
#   return(df)
# }
# 
# 
# fit_predict_flags(Kappa_df, "Case", "IndeK_clean", tol = tol) %>%
#   select(Code_P, IndeK_clean , Case, Statut_clinique, lag_outcome, 
#          lag_cov, lag_cov_scaled,lag_outcome_scaled, pred_unadj, 
#          pred_adj, observed_flag,pred_flag_unadj, pred_flag_adj)
 


ignore_df <- Kappa_df %>% select(Code_P, BOC_count, Case) %>% 
  group_by(Code_P) %>% mutate(lag_Case=lag(Case)) %>% mutate(lag_BOC_count=lag(BOC_count)) %>% 
  drop_na() %>% ungroup()



model_unadj <- lmer(
    formula(" Case ~ lag_Case + (1 | Code_P)"),
    data = ignore_df
  )


model_adj <- lmer(
    formula(" Case ~ lag_Case + lag_BOC_count + (1 | Code_P)"),
    data = ignore_df
  )

model_ind <- lmer(
    formula(" Case ~  lag_BOC_count + (1 | Code_P)"),
    data = ignore_df
  )


ignore_df <- ignore_df %>% 
  bind_cols(data.frame(predict(model_unadj, re.form = NULL)) )  %>% rename("preds_unadj"="predict.model_unadj..re.form...NULL.") %>%
  bind_cols(data.frame(predict(model_adj, re.form = NULL)))  %>% rename("preds_adj"="predict.model_adj..re.form...NULL.")  %>%
  bind_cols(data.frame(predict(model_ind, re.form = NULL)))  %>% rename("preds_ind"="predict.model_ind..re.form...NULL.") 

ignore_df <- ignore_df %>%
  mutate(preds_unadj=round(preds_unadj)) %>%
  mutate(preds_adj=round(preds_adj)) %>%
  mutate(preds_ind=round(preds_ind)) %>%
  mutate(OBSERVED=ifelse( (Case>lag_Case),"Inc",
                                      ifelse(Case<lag_Case, "Dec", "Stab"))) %>%
  mutate(UNADJ=ifelse( (preds_unadj >lag_Case),"Inc",
                                      ifelse(preds_unadj <lag_Case, "Dec", "Stab"))) %>% 
  mutate(ADJ=ifelse( (preds_adj  >lag_Case),"Inc",
                                      ifelse(preds_adj  <lag_Case, "Dec", "Stab"))) %>% 
  mutate(IND=ifelse( (preds_ind >lag_Case),"Inc",
                                      ifelse(preds_unadj <lag_Case, "Dec", "Stab"))) 


kappa_UNADJ <- irr::kappa2(data.frame(ignore_df$OBSERVED, ignore_df$UNADJ))
kappa_ADJ <- irr::kappa2(data.frame(ignore_df$OBSERVED, ignore_df$ADJ))
kappa_IND <- irr::kappa2(data.frame(ignore_df$OBSERVED, ignore_df$IND))



# -------------
 
# Get R2 for adjusted vs non adjusted -----------
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
  mutate(elapsed_sympt = time_length(interval(`date symptômes`, Date_sampling ), "months")) %>%
  mutate(elapsed_dx = time_length(interval(`date dg (T0)`, Date_sampling ), "months"))  %>%
  mutate(elapsed_dx = time_length(interval(`date dg (T0)`, Date_sampling ), "months"))  


names(Kappa_df)

variables_to_track <- c("Code_P",  "IndeK_clean", "KappaPlasma", "AlbuminPlasma" , "AlbuminCSF" , "KappaCSF", 
                        "NfL", "Case",  "mRS",  "Date_sampling",
                        "CSF leuco count" , "CSF prot", "BOC_count" , "CRP" , "Neutrophiles - nb absolu (sang) clean" ,
                        "Lymphocytes - nb absolu (sang) clean" , "NLR", "Statut clinique")

Kappa_df <- Kappa_df %>% select(all_of(variables_to_track))

names(Kappa_df)

Kappa_df <- Kappa_df %>% rename("CSF_leuco_count"="CSF leuco count") %>%
  rename("CSF_prot"="CSF prot") %>%
  rename("Neutrophiles"="Neutrophiles - nb absolu (sang) clean") %>%
  rename("Lymphocytes"="Lymphocytes - nb absolu (sang) clean") %>%
  rename("Statut_clinique"="Statut clinique")

library(dplyr)
library(lme4)
library(lmerTest)
library(purrr)
library(tidyr)
library(Metrics)



# Covariates to test
covariates_to_test <- c("IndeK_clean", "KappaPlasma", "AlbuminPlasma" , "AlbuminCSF" , 
                        "KappaCSF", "NfL", "BOC_count" , "CRP", "CSF_leuco_count" , 
                        "CSF_prot", "Neutrophiles", "Lymphocytes" , "NLR")

# Outcome variable
outcome_var <- "mRS"


# Function to fit lagged models and compute R² + MAE
fit_lagged_models_metrics <- function(df, outcome_var, covariates, scale_data = TRUE) {
  
  df <- df %>% arrange(Code_P, Date_sampling)
  
  # create lagged outcome
  df <- df %>%
    group_by(Code_P) %>%
    mutate(lag_outcome = lag(!!sym(outcome_var))) %>%
    ungroup()
  
  results <- map_dfr(covariates, function(var) {
    
    df_lag <- df %>%
      group_by(Code_P) %>%
      mutate(!!paste0("lag_", var) := lag(!!sym(var))) %>%
      ungroup() %>%
      drop_na(lag_outcome, !!sym(paste0("lag_", var))) # remove NAs
    
    # Optionally scale the lagged covariate (and outcome)
    if(scale_data){
      df_lag <- df_lag %>%
        mutate(
          !!paste0("lag_", var) := as.numeric(scale(!!sym(paste0("lag_", var)))),
          lag_outcome = as.numeric(scale(lag_outcome))
        )
    }
    
    # Model 1: Covariate-only
    model_cov <- lmer(
      formula(paste0(outcome_var, " ~ lag_", var, " + (1 | Code_P)")),
      data = df_lag
    )
    
    pred_cov <- predict(model_cov)
    var_fixed_cov <- var(predict(model_cov, re.form = NA))
    var_random_cov <- sum(as.numeric(VarCorr(model_cov)))
    var_resid_cov <- sigma(model_cov)^2
    R2_marginal_cov <- var_fixed_cov / (var_fixed_cov + var_random_cov + var_resid_cov)
    R2_conditional_cov <- (var_fixed_cov + var_random_cov) / (var_fixed_cov + var_random_cov + var_resid_cov)
    MAE_cov <- mae(df_lag[[outcome_var]], pred_cov)
    
    # Model 2: Outcome-only (lagged outcome)
    model_outcome <- lmer(
      formula(paste0(outcome_var, " ~ lag_outcome + (1 | Code_P)")),
      data = df_lag
    )
    pred_outcome <- predict(model_outcome)
    var_fixed_out <- var(predict(model_outcome, re.form = NA))
    var_random_out <- sum(as.numeric(VarCorr(model_outcome)))
    var_resid_out <- sigma(model_outcome)^2
    R2_marginal_out <- var_fixed_out / (var_fixed_out + var_random_out + var_resid_out)
    R2_conditional_out <- (var_fixed_out + var_random_out) / (var_fixed_out + var_random_out + var_resid_out)
    MAE_out <- mae(df_lag[[outcome_var]], pred_outcome)
    
    # Model 3: Adjusted (lagged outcome + covariate)
    model_adj <- lmer(
      formula(paste0(outcome_var, " ~ lag_outcome + lag_", var, " + (1 | Code_P)")),
      data = df_lag
    )
    pred_adj <- predict(model_adj)
    var_fixed_adj <- var(predict(model_adj, re.form = NA))
    var_random_adj <- sum(as.numeric(VarCorr(model_adj)))
    var_resid_adj <- sigma(model_adj)^2
    R2_marginal_adj <- var_fixed_adj / (var_fixed_adj + var_random_adj + var_resid_adj)
    R2_conditional_adj <- (var_fixed_adj + var_random_adj) / (var_fixed_adj + var_random_adj + var_resid_adj)
    MAE_adj <- mae(df_lag[[outcome_var]], pred_adj)
    
    tibble(
      Covariate = var,
      
      # Covariate-only
      Cov_only_Est = summary(model_cov)$coefficients[2,1],
      Cov_only_SE  = summary(model_cov)$coefficients[2,2],
      Cov_only_p   = summary(model_cov)$coefficients[2,5],
      R2_marginal_cov,
      R2_conditional_cov,
      MAE_cov,
      
      # Outcome-only
      Outcome_only_Est = summary(model_outcome)$coefficients[2,1],
      Outcome_only_SE  = summary(model_outcome)$coefficients[2,2],
      Outcome_only_p   = summary(model_outcome)$coefficients[2,5],
      R2_marginal_out,
      R2_conditional_out,
      MAE_out,
      
      # Adjusted
      Adj_Estimate = summary(model_adj)$coefficients[3,1],
      Adj_SE = summary(model_adj)$coefficients[3,2],
      Adj_p = summary(model_adj)$coefficients[3,5],
      R2_marginal_adj,
      R2_conditional_adj,
      MAE_adj
    )
    
  })
  
  return(results)
}

# Run the pipeline
covariates_to_test <- c("IndeK_clean", "KappaPlasma", "AlbuminPlasma", "AlbuminCSF", 
                        "KappaCSF", "NfL", "BOC_count", "CRP", "CSF_leuco_count", 
                        "CSF_prot", "Neutrophiles", "Lymphocytes", "NLR")

lagged_metrics <- fit_lagged_models_metrics(Kappa_df, "mRS", covariates_to_test)

lagged_metrics


# ----------
# Plot CASE and mRS over time ------------

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
  mutate(elapsed_sympt = time_length(interval(`date symptômes`, Date_sampling ), "months")) %>%
  mutate(elapsed_dx = time_length(interval(`date dg (T0)`, Date_sampling ), "months"))  %>%
  mutate(elapsed_dx = time_length(interval(`date dg (T0)`, Date_sampling ), "months"))  


names(Kappa_df)

variables_to_track <- c("Code_P",  "IndeK_clean", "KappaPlasma", "AlbuminPlasma" , "AlbuminCSF" , "KappaCSF", 
                        "NfL", "Case",  "mRS",  "Date_sampling",
                        "CSF leuco count" , "CSF prot", "BOC_count" , "CRP" , "Neutrophiles - nb absolu (sang) clean" ,
                        "Lymphocytes - nb absolu (sang) clean" , "NLR", "Statut clinique")

Kappa_df <- Kappa_df %>% select(all_of(variables_to_track))

names(Kappa_df)

Kappa_df <- Kappa_df %>% rename("CSF_leuco_count"="CSF leuco count") %>%
  rename("CSF_prot"="CSF prot") %>%
  rename("Neutrophiles"="Neutrophiles - nb absolu (sang) clean") %>%
  rename("Lymphocytes"="Lymphocytes - nb absolu (sang) clean") %>%
  rename("Statut_clinique"="Statut clinique")


to_plot <- Kappa_df %>% group_by(Code_P) %>% mutate(First=min(Date_sampling)) %>%
  select(Code_P,Case, First, Date_sampling) %>%
  mutate(elapsed=as.numeric(Date_sampling-First)/30.44 ) %>%
  ggplot(aes(elapsed, Case)) +
  geom_line(aes(group=Code_P), linewidth=1, alpha=0.2, color="#00204D") +
  geom_point( alpha=0.6, shape=1, stroke=2,  size=2, color="#00204D") +
  geom_smooth(color="#00204D", fill="#747475", alpha=0.2, size=2, method="loess") +
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
  ylab("CASE Score \n")+ xlab("\n Number of elapsed months from 1st evaluation ")

ggsave(file="CASE_over_time.svg", plot=to_plot, width=5, height=5)


to_plot <- Kappa_df %>% group_by(Code_P) %>% mutate(First=min(Date_sampling)) %>%
  select(Code_P,mRS, First, Date_sampling) %>%
  mutate(elapsed=as.numeric(Date_sampling-First)/30.44 ) %>%
  ggplot(aes(elapsed, mRS)) +
  geom_line(aes(group=Code_P), linewidth=1, alpha=0.2, color="#00204D") +
  geom_point( alpha=0.6, shape=1, stroke=2,  size=2, color="#00204D") +
  geom_smooth(color="#00204D", fill="#747475", alpha=0.2, size=2, method="loess") +
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
  ylab("mRS Score \n")+ xlab("\n Number of elapsed months from 1st evaluation ")

ggsave(file="mRS_over_time.svg", plot=to_plot, width=5, height=5)





to_plot <- Kappa_df %>% group_by(Code_P) %>% slice(1) %>%
  ungroup() %>% select(Code_P, Case) %>% rename("first_case"="Case") %>%
  left_join(Kappa_df, by="Code_P") %>%
  group_by(Code_P) %>% mutate(First=min(Date_sampling)) %>%
  select(Code_P,Case, first_case, First, Date_sampling) %>%
  mutate(elapsed=as.numeric(Date_sampling-First)/30.44 ) %>%
  mutate(Case_delta=Case-first_case) %>%
  ggplot(aes(elapsed, Case_delta)) +
  geom_line(aes(group=Code_P), linewidth=1, alpha=0.2, color="#00204D") +
  geom_point( alpha=0.6, shape=1, stroke=2,  size=2, color="#00204D") +
  geom_smooth(color="#00204D", fill="#747475", alpha=0.2, size=2, method="loess") +
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
  ylab("CASE Score Delta from 1st evaluation\n")+ xlab("\n Number of elapsed months from 1st evaluation ")

ggsave(file="CASE_deltas_over_time.svg", plot=to_plot, width=5, height=5)





to_plot <- Kappa_df %>% group_by(Code_P) %>% slice(1) %>%
  ungroup() %>% select(Code_P, mRS) %>% rename("first_mRS"="mRS") %>%
  left_join(Kappa_df, by="Code_P") %>%
  group_by(Code_P) %>% mutate(First=min(Date_sampling)) %>%
  select(Code_P,mRS, first_mRS, First, Date_sampling) %>%
  mutate(elapsed=as.numeric(Date_sampling-First)/30.44 ) %>%
  mutate(mRS_delta=mRS-first_mRS) %>%
  ggplot(aes(elapsed, mRS_delta)) +
  geom_line(aes(group=Code_P), linewidth=1, alpha=0.2, color="#00204D") +
  geom_point( alpha=0.6, shape=1, stroke=2,  size=2, color="#00204D") +
  geom_smooth(color="#00204D", fill="#747475", alpha=0.2, size=2, method="loess") +
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
  ylab("mRS Score Delta from 1st evaluation\n")+ xlab("\n Number of elapsed months from 1st evaluation ")

ggsave(file="mRS_deltas_over_time.svg", plot=to_plot, width=5, height=5)

# -----------

# XGBoost CASE --------------
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
  mutate(elapsed_sympt = time_length(interval(`date symptômes`, Date_sampling ), "months")) %>%
  mutate(elapsed_dx = time_length(interval(`date dg (T0)`, Date_sampling ), "months")) 



variables_to_track <- c("Code_P",  "IndeK_clean", "KappaPlasma", "AlbuminPlasma" , "AlbuminCSF" , "KappaCSF", 
                        "NfL", "Case",  "mRS", "abnormal MRI onset", "néo associé" ,"ICI", 
                        "CSF leuco count" , "CSF prot", "BOC_count" , "CRP" , "Neutrophiles - nb absolu (sang) clean" ,
                        "Lymphocytes - nb absolu (sang) clean" , "NLR", "elapsed_sympt", "elapsed_dx")

Kappa_df <- Kappa_df %>% select(all_of(variables_to_track))


Kappa_df <- Kappa_df %>% rename("Leukocytes CSF"="CSF leuco count") %>%
  rename("Protein CSF"="CSF prot") %>%
  rename("Neutrophiles Serum"="Neutrophiles - nb absolu (sang) clean") %>%
  rename("Lymphocytes Serum"="Lymphocytes - nb absolu (sang) clean")  %>%
  rename("κ FLC Index"="IndeK_clean") %>%
  rename("Albumin Plasma"="AlbuminPlasma") %>%
  rename("κ Plasma"="KappaPlasma") %>%
  rename("Albumin CSF"="AlbuminCSF") %>%
  rename("Oligoclonal bands"="BOC_count") %>%
  rename("κ CSF"="KappaCSF") %>%
  rename("Abnormal MRI"="abnormal MRI onset") %>%
  rename("Neoplasia-assoc"="néo associé") 



Kappa_df <- Kappa_df %>% group_by(Code_P) %>% mutate(future_case=lead(Case)) %>% ungroup() %>%
  select(-elapsed_dx, -elapsed_sympt, -Code_P, mRS, Case) %>%
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


#rmse_cv <- y %>% bind_cols(preds_full) %>% filter(`...1`<20& `...2`<20)

# Calculate RMSE on all samples
#error <- sqrt(mean((rmse_cv$`...1` - rmse_cv$`...2`)^2))
#print(paste("10-fold CV RMSE:", error))



# Plot observed vs predicted for all samples
#plot_df <- data.frame(observed = rmse_cv$`...1`, predicted = rmse_cv$`...2`)




#plot_df_inc <- plot_df %>% mutate(group="Including Prev mRS|CASE")

#plot_df_exc <- plot_df %>% mutate(group="Excluding Prev mRS|CASE")


# to_plot <- plot_df_inc %>% bind_rows(plot_df_exc) %>%
#   ggplot(aes(observed, predicted, colour=group, fill=group)) +
#   geom_jitter(width=0.2, alpha=0.5, size=1.5, shape=1, stroke=2) +
#   geom_smooth(method = "lm", se = FALSE) +
#   coord_cartesian(xlim = c(0, 15), ylim = c(0, 15)) +
#   scale_colour_manual(values=c("#00204D", "firebrick")) +
#   scale_fill_manual(values=c("#00204D", "firebrick")) +
#   labs(title = "",
#        x = "\n Observed Future CASE Score", y = "Predicted Future CASE Score \n") +
#   theme(axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         legend.position = "right") +
#   theme(panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.background = element_blank(),
#         strip.text = element_blank(),
#         axis.line = element_blank(),
#         axis.text.x = element_text(size = 10),
#         axis.text.y = element_text(size = 10),
#         axis.title.x = element_text(size = 12, vjust = -0.5),
#         axis.title.y = element_text(size = 12, vjust = -0.5),
#         plot.margin = margin(5, 5, 5, 5, "pt"))
# 
# ggsave(file="xgboost.svg", plot=to_plot, width=7, height=5)




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
  geom_col(fill="firebrick", colour="firebrick", alpha=0.8) +
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

# XGBoost mRS --------------
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
  mutate(elapsed_sympt = time_length(interval(`date symptômes`, Date_sampling ), "months")) %>%
  mutate(elapsed_dx = time_length(interval(`date dg (T0)`, Date_sampling ), "months")) 



variables_to_track <- c("Code_P",  "IndeK_clean", "KappaPlasma", "AlbuminPlasma" , "AlbuminCSF" , "KappaCSF", 
                        "NfL", "Case",  "mRS", "abnormal MRI onset", "néo associé" ,"ICI", 
                        "CSF leuco count" , "CSF prot", "BOC_count" , "CRP" , "Neutrophiles - nb absolu (sang) clean" ,
                        "Lymphocytes - nb absolu (sang) clean" , "NLR", "elapsed_sympt", "elapsed_dx")

Kappa_df <- Kappa_df %>% select(all_of(variables_to_track))


Kappa_df <- Kappa_df %>% rename("Leukocytes CSF"="CSF leuco count") %>%
  rename("Protein CSF"="CSF prot") %>%
  rename("Neutrophiles Serum"="Neutrophiles - nb absolu (sang) clean") %>%
  rename("Lymphocytes Serum"="Lymphocytes - nb absolu (sang) clean")  %>%
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

# Correlations first vs last visits --------------

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
  mutate(elapsed_sympt = time_length(interval(`date symptômes`, Date_sampling ), "months")) %>%
  mutate(elapsed_dx = time_length(interval(`date dg (T0)`, Date_sampling ), "months")) 



variables_to_track <- c("Code_P", "Date_sampling", "IndeK_clean", "KappaPlasma", "AlbuminPlasma" , "AlbuminCSF" , "KappaCSF", 
                        "NfL", "Case",  "mRS", "abnormal MRI onset", "néo associé" ,"ICI", 
                        "CSF leuco count" , "CSF prot", "BOC_count" , "CRP" , "Neutrophiles - nb absolu (sang) clean" ,
                        "Lymphocytes - nb absolu (sang) clean" , "NLR", "elapsed_sympt", "elapsed_dx", "Statut clinique")

Kappa_df <- Kappa_df %>% select(all_of(variables_to_track))


Kappa_df <- Kappa_df %>% rename("Leukocytes CSF"="CSF leuco count") %>%
  rename("Protein CSF"="CSF prot") %>%
  rename("Neutrophiles Serum"="Neutrophiles - nb absolu (sang) clean") %>%
  rename("Lymphocytes Serum"="Lymphocytes - nb absolu (sang) clean")  %>%
  rename("κ FLC Index"="IndeK_clean") %>%
  rename("Albumin Plasma"="AlbuminPlasma") %>%
  rename("κ Plasma"="KappaPlasma") %>%
  rename("Albumin CSF"="AlbuminCSF") %>%
  rename("Oligoclonal bands"="BOC_count") %>%
  rename("κ CSF"="KappaCSF") %>%
  rename("Abnormal MRI"="abnormal MRI onset") %>%
  rename("Neoplasia-assoc"="néo associé")  %>%
  rename("Clinical_status"="Statut clinique")  


ignore <- Kappa_df %>% select(`κ FLC Index`, Case, mRS)




# # subset and drop missing values
# df_kappaindex <- Kappa_df %>%
#   select(Code_P, Case, `κ FLC Index`) %>%
#   drop_na()
# 
# # plot
# df_kappaindex %>% group_by(Code_P) %>% slice(1)  %>% ungroup() %>%
#   ggplot(aes(x = `κ FLC Index`, y = Case, color = Code_P)) +
#   geom_point(size = 2, alpha = 0.6) +
#   ylim(0,15) +
#   xlim(0,60) +
#   geom_smooth(method = "lm", se = FALSE, aes(group = 1), color = "black", linetype = "dashed", size = 1.2) +
#   geom_smooth(method = "lm", se = FALSE, aes(group = Code_P), linetype = "solid", size = 0.7, alpha = 0.5) +
#   labs(
#     title = "`κ FLC Index` vs Case",
#     subtitle = "Dashed = overall regression (between subjects)\nColored lines = within-subject regressions",
#     x = "`κ FLC Index`",
#     y = "Case"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(legend.position = "none")


cor.test(ignore$`κ FLC Index`, ignore$Case, method = "spearman")

cor.test(ignore$`κ FLC Index`, ignore$mRS, method = "spearman")

names(Kappa_df)

test <- Kappa_df %>% group_by(Code_P) %>% filter(Date_sampling==min(Date_sampling)) %>%
  select(Code_P, `κ FLC Index`) %>%
  inner_join(Kappa_df %>% group_by(Code_P) %>% filter(Date_sampling==min(Date_sampling)) %>%
               select(Code_P, `Case`) %>% rename("First"="Case")) %>%
  inner_join(Kappa_df %>% group_by(Code_P) %>% filter(Date_sampling==max(Date_sampling)) %>%
               select(Code_P, `Case`) %>% rename("Last"="Case")) %>%
  drop_na() %>% ungroup() 

cor.test(test$`κ FLC Index`, test$First)
cor.test(test$`κ FLC Index`, test$Last)


test %>%
  gather("Time", "Case", First, Last) %>%
  filter(Case<=20) %>%
  ggplot(aes(x=`κ FLC Index`, y=Case, color=Time)) +
  facet_wrap(~Time) +
  geom_jitter() +
  geom_smooth(method="lm", se=T) 











# List of variables you want to keep
vars <- c("κ FLC Index", "κ Plasma", "Albumin Plasma", "Albumin CSF", 
          "κ CSF", "NfL", "Case", "mRS", "Abnormal MRI", "Neoplasia-assoc", 
          "ICI", "Leukocytes CSF", "Protein CSF", "Oligoclonal bands", 
          "CRP", "Neutrophiles Serum", "Lymphocytes Serum", "NLR")

# First visit per patient
first_df <- Kappa_df %>%
  group_by(Code_P) %>%
  filter(Date_sampling == min(Date_sampling)) %>%
  select(Code_P, all_of(vars)) %>%
  rename_with(~ paste0("First_", .), -Code_P)

# Last visit per patient
last_df <- Kappa_df %>%
  group_by(Code_P) %>%
  filter(Date_sampling == max(Date_sampling)) %>%
  select(Code_P, Case, mRS) %>%
  rename(Last_Case = Case, Last_mRS = mRS)

# Merge into one per-patient table
patient_df <- first_df %>%
  inner_join(last_df, by = "Code_P") %>%
  ungroup()



# Biomarkers to correlate
bio_vars <- setdiff(names(patient_df), c("Code_P", "First_Case", "First_mRS", "Last_Case", "Last_mRS"))

# Function to compute correlations
cor_fun <- function(var, outcome) {
  x <- patient_df[[var]]
  y <- patient_df[[outcome]]
  if (is.numeric(x) & is.numeric(y)) {
    res <- cor.test(x, y, method="spearman")
    tibble(Variable = var,
           Outcome = outcome,
           Correlation = res$estimate,
           p_value = res$p.value)
  } else {
    tibble(Variable = var, Outcome = outcome,
           Correlation = NA, p_value = NA)
  }
}

# Apply to all combinations (first biomarker vs outcomes)
results <- purrr::map_dfr(bio_vars, ~ cor_fun(.x, "First_Case"))
results <- bind_rows(results,
                     purrr::map_dfr(bio_vars, ~ cor_fun(.x, "Last_Case")))
results <- bind_rows(results,
                     purrr::map_dfr(bio_vars, ~ cor_fun(.x, "First_mRS")))
results <- bind_rows(results,
                     purrr::map_dfr(bio_vars, ~ cor_fun(.x, "Last_mRS")))


results %>% arrange(Variable)


# -------
# Baseline summary info per Ab location type -----------

Kappa_df <- read_xlsx(path="../data/Kappa_NFL20250721_clean.xlsx", trim_ws = TRUE)

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


Kappa_df %>% group_by(Code_P, `Ab type`) %>% count() %>%
  ungroup() %>% group_by(`Ab type`) %>%
  summarise(mean=mean(n),
            sd=sd(n),
            median=median(n),
            q1=quantile(n, 0.25),
            q3=quantile(n,0.75))


test <- Kappa_df %>% group_by(Code_P, `Ab type`) %>% count() %>%
  ungroup() %>% rename("Ab_type"=`Ab type`)

kruskal.test(n ~ Ab_type, data = test)




Kappa_df <- Kappa_df %>% mutate(NLR=`Neutrophiles - nb absolu (sang) clean`/`Lymphocytes - nb absolu (sang) clean`)

firsts_df <- Kappa_df %>% arrange(Code_P, Date_IndexK) %>% group_by(Code_P) %>% slice(1)

firsts_df %>% group_by(`Ab type`) %>% count()

firsts_df %>% group_by(`Traitement 1ère ligne (O/N)`, `Traitement 2nd ligne (O/N)`) %>% count()


firsts_df <- firsts_df %>% rename("Ab_type"=`Ab type`)

variables_to_track <- c("IndeK_clean", "KappaPlasma", "AlbuminPlasma" , "AlbuminCSF" , "KappaCSF", 
                        "NfL", "Case",  "mRS", "abnormal MRI onset", "néo associé" ,"ICI", 
                        "CSF leuco count" , "CSF prot", "BOC_count" , "CRP" , "Neutrophiles - nb absolu (sang) clean" ,
                        "Lymphocytes - nb absolu (sang) clean" , "NLR", "Ab_type")



summary_table <- firsts_df %>%
  ungroup() %>%
  select(all_of(variables_to_track)) %>%
  group_by(Ab_type) %>%
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

data.frame(summary_table_long)

firsts_df %>% ungroup() %>% group_by(Ab_type, `abnormal MRI onset`) %>% count()
firsts_df %>% ungroup() %>% group_by(Ab_type, `néo associé` ) %>% count()
firsts_df %>% ungroup() %>% group_by(Ab_type, ICI) %>% count()

firsts_df %>% ungroup() %>% group_by(Ab_type, `Traitement 1ère ligne (O/N)`) %>% count()
firsts_df %>% ungroup() %>% group_by(Ab_type, `Traitement 2nd ligne (O/N)`) %>% count()


Kappa_df$`date symptômes` <- as.Date(as.numeric(Kappa_df$`date symptômes`), origin = "1899-12-30")
Kappa_df$`date dg (T0)`  <- as.Date(Kappa_df$`date dg (T0)`  )
Kappa_df$`Date_sampling`  <- as.Date(Kappa_df$`Date_sampling`  )


library(lubridate)

Kappa_df <- Kappa_df %>%
  mutate(elapsed_sympt = time_length(interval(`date symptômes`, Date_sampling ), "days")) %>%
  mutate(elapsed_dx = time_length(interval(`date dg (T0)`, Date_sampling ), "days")) 

firsts_df <- Kappa_df %>% arrange(Code_P, Date_IndexK) %>% group_by(Code_P) %>% slice(1)

firsts_df <- firsts_df %>% rename("Ab_type"=`Ab type`)

firsts_df %>% ungroup() %>%
  mutate(elapsed_sympt=elapsed_sympt/30.44) %>%
  ungroup() %>%
  group_by(Ab_type) %>%
  summarise(mean=mean(elapsed_sympt),
            sd=sd(elapsed_sympt),
            median=median(elapsed_sympt),
            q1=quantile(elapsed_sympt, 0.25),
            q3=quantile(elapsed_sympt,0.75))



firsts_df <- firsts_df %>% ungroup() %>%
  mutate(elapsed_sympt=elapsed_sympt/30.44)

names(firsts_df)


firsts_df <- firsts_df %>% rename("CSF_leuco_count"="CSF leuco count") %>%
  rename("CSF_prot"="CSF prot") %>%
  rename("Neutrophiles"="Neutrophiles - nb absolu (sang) clean") %>%
  rename("Lymphocytes"="Lymphocytes - nb absolu (sang) clean") %>%
  rename("Statut_clinique"="Statut clinique")


# List of variables to test
vars_to_test <- c(
  "IndeK_clean", "KappaPlasma", "AlbuminPlasma", "AlbuminCSF",
  "KappaCSF", "NfL", "Case", "mRS", "CSF_leuco_count", "CSF_prot",
  "BOC", "BOC_count", "CRP", 
  "Neutrophiles", 
  "Lymphocytes", 
  "NLR", "elapsed_sympt", "elapsed_dx"
)



# Initialize a data frame to store results
kruskal_results <- data.frame(
  Variable = character(),
  Chi_Sq = numeric(),
  df = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop over variables
for (var in vars_to_test) {
  # Prepare formula
  f <- as.formula(paste(var, "~ Ab_type"))
  
  # Run Kruskal-Wallis test, ignoring NAs
  test <- kruskal.test(f, data = firsts_df)
  
  # Store results
  kruskal_results <- rbind(
    kruskal_results,
    data.frame(
      Variable = var,
      Chi_Sq = test$statistic,
      df = test$parameter,
      p_value = test$p.value
    )
  )
}

# Adjust p-values for multiple comparisons across variables
kruskal_results$p_adj <- p.adjust(kruskal_results$p_value, method = "holm")  # or "fdr"

# Optional: order by adjusted p-value
kruskal_results <- kruskal_results[order(kruskal_results$p_adj), ]

# View results
kruskal_results








per Ab type




# List of binary variables
binary_vars <- c(
  "abnormal MRI onset", 
  "néo associé", 
  "ICI", 
  "Traitement 1ère ligne (O/N)", 
  "Traitement 2nd ligne (O/N)"
)

# Initialize results data frame
fisher_results <- data.frame(
  Variable = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop over variables
for (var in binary_vars) {
  # Create contingency table: variable vs Ab_type
  tbl <- table(firsts_df$Ab_type, firsts_df[[var]], useNA = "no")
  
  # Run Fisher's Exact Test
  test <- fisher.test(tbl)
  
  # Store results
  fisher_results <- rbind(
    fisher_results,
    data.frame(
      Variable = var,
      p_value = test$p.value
    )
  )
}

# Adjust p-values for multiple variables (e.g., Holm, Bonferroni, or FDR)
fisher_results$p_adj <- p.adjust(fisher_results$p_value, method = "holm")  

# Optional: order by adjusted p-value
fisher_results <- fisher_results[order(fisher_results$p_adj), ]

# View results
fisher_results


# --------------

