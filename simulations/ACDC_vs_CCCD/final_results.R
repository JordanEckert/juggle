library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(latticeExtra)
library(forcats)
library(gridExtra)
library(grid)

path <- "~/Documents/GitHub/juggle/simulations/ACDC_vs_CCCD/results"

#### Overlapped Count Data ####
# List all files containing "counts" in the filename
count_files <- list.files(path = path,
                          pattern = "counts_.*\\.txt$", full.names = TRUE)

counts_data <- NULL

for (file in count_files){
  # Read the file
  counts <- read.csv(file, sep="")

  counts$V15 <- as.factor(1:5)
  rownames(counts) <- NULL

  # Create factors for dim, delta, ir
  counts$V1 <- as.factor(counts$V1)
  counts$V2 <- as.factor(counts$V2)
  counts$V3 <- as.factor(counts$V3)

  counts_data <- rbind(counts_data, counts)
}

colnames(counts_data) <- c("dim", "delta", "ir", "tau0", "tau.1", "tau.2", "tau.3", "tau.4", "tau.5", "tau.6", "tau.7", "tau.8", "tau.9", "tau1", "k")

### Count
df_long <- counts_data %>%
  pivot_longer(
    cols = starts_with("tau"),
    names_to = "tau_label",
    values_to = "count"
  )


df_long <- df_long %>%
  mutate(tau_label = gsub("tau", "", tau_label))

df_long$tau_label <- as.factor(df_long$tau_label)
# write.csv(df_long, "../long_count_data.csv")

# Only overlap
overlap_only <- df_long %>%
  filter(ir == 1) %>%
  mutate(proportion = count / 200)

overlap_only$tau_label <- fct_relevel(overlap_only$tau_label, "0")
overlap_only$dim <- fct_relevel(overlap_only$dim, "2", "3", "5", "10")
overlap_only$tau_label <- fct_recode(overlap_only$tau_label,
                                     "0.0" = "0",
                                     "0.1" = ".1",
                                     "0.2" = ".2",
                                     "0.3" = ".3",
                                     "0.4" = ".4",
                                     "0.5" = ".5",
                                     "0.6" = ".6",
                                     "0.7" = ".7",
                                     "0.8" = ".8",
                                     "0.9" = ".9",
                                     "1.0" = "1")

# Bar plot of counts for tau_label
ggplot(overlap_only,
       aes(x = tau_label, y = proportion, fill = as.factor(k))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(dim ~ delta, scales = "free") +
  labs(title = "Proportion of Optimal Tau by k",
       x = expression("Tau" ~ (tau)), y = "Proportion", fill = "k") +
  theme_minimal()

#### Overlapped Results Data ####
result_files <- list.files(path = path,
                           pattern = "results_.*\\.txt$",
                           full.names = TRUE)

results_data <- NULL

for(file in result_files){
  results <- read.csv(file, sep="")
  rownames(results) <- NULL
  colnames(results) <- c("dim", "delta", "ir", "PCCCD", "RWCCCD",
                         "ACDC1", "ACDC2", "ACDC3", "ACDC4", "ACDC5")

  # AUC Boundary
  results[,4:10] <- lapply(results[,4:10],
                           function(x) ifelse(x < .5, .5, x))

  results_data <- rbind(results_data, results)

}

results_long <- results_data %>%
  pivot_longer(cols = PCCCD:ACDC5,  # Select all model performance columns
               names_to = "model",  # New column for model names
               values_to = "performance")

overlap_results <- results_long %>%
  filter(ir == 1)

ggplot(overlap_results, aes(x = delta, y = performance, color = model)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ dim) +  # Separate plots for each dimension
  theme_minimal() +
  labs(title = "Model Performance Across Delta and Dimension",
       x = expression("Delta" ~ (Delta)),
       y = "AUC",
       color = "Model")

results_long_2 <- results_data %>%
  filter(ir == 1) %>%
  mutate("P-ACDC1" = PCCCD - ACDC1,
         "P-ACDC2" = PCCCD - ACDC2,
         "P-ACDC3" = PCCCD - ACDC3,
         "P-ACDC4" = PCCCD - ACDC4,
         "P-ACDC5" = PCCCD - ACDC5,
         "RW-ACDC1" = RWCCCD - ACDC1,
         "RW-ACDC2" = RWCCCD - ACDC2,
         "RW-ACDC3" = RWCCCD - ACDC3,
         "RW-ACDC4" = RWCCCD - ACDC4,
         "RW-ACDC5" = RWCCCD - ACDC5) %>%
  pivot_longer(cols = 'P-ACDC1':'RW-ACDC5',
               names_to = "model",  # New column for model names
               values_to = "performance")

ggplot(results_long_2, aes(x = delta, y = performance, color = model)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ dim) +  # Separate plots for each dimension
  theme_minimal() +
  labs(title = "Model Performance Differences Across Delta and Dim",
       x = expression("Delta" ~ (Delta)),
       y = "Difference in AUC",
       color = "Model")

#### Overlapped and Imbalanced Count Data ####
# List all files containing "counts" in the filename
count_files <- list.files(path = path,
                          pattern = "counts_.*\\.txt$", full.names = TRUE)

counts_data <- NULL

for (file in count_files){
  # Read the file
  counts <- read.csv(file, sep="")

  counts$V15 <- as.factor(1:5)
  rownames(counts) <- NULL

  # Create factors for dim, delta, ir
  counts$V1 <- as.factor(counts$V1)
  counts$V2 <- as.factor(counts$V2)
  counts$V3 <- as.factor(counts$V3)

  counts_data <- rbind(counts_data, counts)
}

colnames(counts_data) <- c("dim", "delta", "ir", "tau0", "tau.1", "tau.2", "tau.3", "tau.4", "tau.5", "tau.6", "tau.7", "tau.8", "tau.9", "tau1", "k")

### Count
df_long <- counts_data %>%
  pivot_longer(
    cols = starts_with("tau"),
    names_to = "tau_label",
    values_to = "count"
  )


df_long <- df_long %>%
  mutate(tau_label = gsub("tau", "", tau_label))

df_long$tau_label <- as.factor(df_long$tau_label)
# write.csv(df_long, "../long_count_data.csv")

df_long$tau_label <- fct_relevel(df_long$tau_label, "0")
df_long$dim <- fct_relevel(df_long$dim, "2", "3", "5", "10")
df_long$tau_label <- fct_recode(df_long$tau_label,
                                     "0.0" = "0",
                                     "0.1" = ".1",
                                     "0.2" = ".2",
                                     "0.3" = ".3",
                                     "0.4" = ".4",
                                     "0.5" = ".5",
                                     "0.6" = ".6",
                                     "0.7" = ".7",
                                     "0.8" = ".8",
                                     "0.9" = ".9",
                                     "1.0" = "1")

df_long_imbalance_1 <- df_long %>%
  mutate(proportion = count / 200) %>%
  filter (ir != 1) %>%
  filter (delta == .25)

df_long_imbalance_2 <- df_long %>%
  mutate(proportion = count / 200) %>%
  filter (ir != 1) %>%
filter (delta == .50)

df_long_imbalance_3 <- df_long %>%
  mutate(proportion = count / 200) %>%
  filter (ir != 1) %>%
  filter (delta == .75)

ggplot(df_long_imbalance_1,
       aes(x = tau_label, y = proportion, fill = as.factor(k))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(dim ~ ir, scales = "free") +
  labs(title = expression("Proportion of Optimal Tau" ~ (Delta ~ '= 0.25')),
       x = expression("Tau" ~ (tau)), y = "Proportion", fill = "k") +
  theme_minimal()

ggplot(df_long_imbalance_2,
       aes(x = tau_label, y = proportion, fill = as.factor(k))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(dim ~ ir, scales = "free") +
  labs(title = expression("Proportion of Optimal Tau" ~ (Delta ~ '= 0.50')),
       x = expression("Tau" ~ (tau)), y = "Proportion", fill = "k") +
  theme_minimal()

ggplot(df_long_imbalance_3,
       aes(x = tau_label, y = proportion, fill = as.factor(k))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(dim ~ ir, scales = "free") +
  labs(title = expression("Proportion of Optimal Tau" ~ (Delta ~ '= 0.75')),
       x = expression("Tau" ~ (tau)), y = "Proportion", fill = "k") +
  theme_minimal()


#### Overlapped and Imbalanced Results Data ####
result_files <- list.files(path = path,
                           pattern = "results_.*\\.txt$",
                           full.names = TRUE)

results_data <- NULL

for(file in result_files){
  results <- read.csv(file, sep="")
  rownames(results) <- NULL
  colnames(results) <- c("dim", "delta", "ir", "PCCCD", "RWCCCD",
                         "ACDC1", "ACDC2", "ACDC3", "ACDC4", "ACDC5")

  # AUC Boundary
  results[,4:10] <- lapply(results[,4:10],
                           function(x) ifelse(x < .5, .5, x))

  results_data <- rbind(results_data, results)

}

acdc_cols <- grep("ACDC", names(results_data), value = TRUE)

df_new <- results_data %>%
  mutate(across(all_of(acdc_cols), ~ PCCCD - .x, .names = "P-{.col}")) %>%
  mutate(across(all_of(acdc_cols), ~ RWCCCD - .x, .names = "RW-{.col}")) %>%
  select(dim, delta, ir, starts_with("P-"), starts_with("RW-"))

subset_data_2 <- df_new %>%
  filter(dim == 2)

subset_data_3 <- df_new %>%
  filter(dim == 3)

subset_data_5 <- df_new %>%
  filter(dim == 5)

subset_data_10 <- df_new %>%
  filter(dim == 10)

p_acdc1_2 <- levelplot(`P-ACDC1` ~ delta * ir, df_new,
                       panel = panel.levelplot.points, cex = 1.2,
                       xlab = expression("Shifting Parameter " ~ (Delta)),
                       ylab = "Global Imbalance Ratio (q)") +
  layer_(panel.2dsmoother(..., n = 200))

p_acdc2_2 <- levelplot(`P-ACDC2` ~ delta * ir, subset_data_2,
                       panel = panel.levelplot.points, cex = 1.2,
                       xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                       ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

p_acdc3_2 <- levelplot(`P-ACDC3` ~ delta * ir, subset_data_2,
                       panel = panel.levelplot.points, cex = 1.2,
                       xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                       ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

p_acdc4_2 <- levelplot(`P-ACDC4` ~ delta * ir, subset_data_2,
                       panel = panel.levelplot.points, cex = 1.2,
                       xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                       ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

p_acdc5_2 <- levelplot(`P-ACDC5` ~ delta * ir, subset_data_2,
                       panel = panel.levelplot.points, cex = 1.2,
                       xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                       ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

p_acdc1_3 <- levelplot(`P-ACDC1` ~ delta * ir, subset_data_3,
                       panel = panel.levelplot.points, cex = 1.2,
                       xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                       ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

p_acdc2_3 <- levelplot(`P-ACDC2` ~ delta * ir, subset_data_3,
                       panel = panel.levelplot.points, cex = 1.2,
                       xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                       ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

p_acdc3_3 <- levelplot(`P-ACDC3` ~ delta * ir, subset_data_3,
                       panel = panel.levelplot.points, cex = 1.2,
                       xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                       ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

p_acdc4_3 <- levelplot(`P-ACDC4` ~ delta * ir, subset_data_3,
                       panel = panel.levelplot.points, cex = 1.2,
                       xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                       ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

p_acdc5_3 <- levelplot(`P-ACDC5` ~ delta * ir, subset_data_3,
                       panel = panel.levelplot.points, cex = 1.2,
                       xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                       ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

p_acdc1_5 <- levelplot(`P-ACDC1` ~ delta * ir, subset_data_5,
                       panel = panel.levelplot.points, cex = 1.2,
                       xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                       ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

p_acdc2_5 <- levelplot(`P-ACDC2` ~ delta * ir, subset_data_5,
                       panel = panel.levelplot.points, cex = 1.2,
                       xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                       ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

p_acdc3_5 <- levelplot(`P-ACDC3` ~ delta * ir, subset_data_5,
                       panel = panel.levelplot.points, cex = 1.2,
                       xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                       ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

p_acdc4_5 <- levelplot(`P-ACDC4` ~ delta * ir, subset_data_5,
                       panel = panel.levelplot.points, cex = 1.2,
                       xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                       ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

p_acdc5_5 <- levelplot(`P-ACDC5` ~ delta * ir, subset_data_5,
                       panel = panel.levelplot.points, cex = 1.2,
                       xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                       ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

p_acdc1_10 <- levelplot(`P-ACDC1` ~ delta * ir, subset_data_10,
                        panel = panel.levelplot.points, cex = 1.2,
                        xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                        ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

p_acdc2_10 <- levelplot(`P-ACDC2` ~ delta * ir, subset_data_10,
                        panel = panel.levelplot.points, cex = 1.2,
                        xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                        ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

p_acdc3_10 <- levelplot(`P-ACDC3` ~ delta * ir, subset_data_10,
                        panel = panel.levelplot.points, cex = 1.2,
                        xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                        ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

p_acdc4_10 <- levelplot(`P-ACDC4` ~ delta * ir, subset_data_10,
                        panel = panel.levelplot.points, cex = 1.2,
                        xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                        ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

p_acdc5_10 <- levelplot(`P-ACDC5` ~ delta * ir, subset_data_10,
                        panel = panel.levelplot.points, cex = 1.2,
                        xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                        ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

rw_acdc1_10 <- levelplot(`RW-ACDC1` ~ delta * ir, subset_data_10,
                         panel = panel.levelplot.points, cex = 1.2,
                         xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                         ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

rw_acdc2_10 <- levelplot(`RW-ACDC2` ~ delta * ir, subset_data_10,
                         panel = panel.levelplot.points, cex = 1.2,
                         xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                         ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

rw_acdc3_10 <- levelplot(`RW-ACDC3` ~ delta * ir, subset_data_10,
                         panel = panel.levelplot.points, cex = 1.2,
                         xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                         ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

rw_acdc4_10 <- levelplot(`RW-ACDC4` ~ delta * ir, subset_data_10,
                         panel = panel.levelplot.points, cex = 1.2,
                         xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                         ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

rw_acdc5_10 <- levelplot(`RW-ACDC5` ~ delta * ir, subset_data_10,
                         panel = panel.levelplot.points, cex = 1.2,
                         xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                         ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

rw_acdc1_5 <- levelplot(`RW-ACDC1` ~ delta * ir, subset_data_5,
                        panel = panel.levelplot.points, cex = 1.2,
                        xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                        ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

rw_acdc2_5 <- levelplot(`RW-ACDC2` ~ delta * ir, subset_data_5,
                        panel = panel.levelplot.points, cex = 1.2,
                        xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                        ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

rw_acdc3_5 <- levelplot(`RW-ACDC3` ~ delta * ir, subset_data_5,
                        panel = panel.levelplot.points, cex = 1.2,
                        xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                        ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

rw_acdc4_5 <- levelplot(`RW-ACDC4` ~ delta * ir, subset_data_5,
                        panel = panel.levelplot.points, cex = 1.2,
                        xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                        ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

rw_acdc5_5 <- levelplot(`RW-ACDC5` ~ delta * ir, subset_data_5,
                        panel = panel.levelplot.points, cex = 1.2,
                        xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                        ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

rw_acdc1_3 <- levelplot(`RW-ACDC1` ~ delta * ir, subset_data_3,
                        panel = panel.levelplot.points, cex = 1.2,
                        xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                        ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

rw_acdc2_3 <- levelplot(`RW-ACDC2` ~ delta * ir, subset_data_3,
                        panel = panel.levelplot.points, cex = 1.2,
                        xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                        ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

rw_acdc3_3 <- levelplot(`RW-ACDC3` ~ delta * ir, subset_data_3,
                        panel = panel.levelplot.points, cex = 1.2,
                        xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                        ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

rw_acdc4_3 <- levelplot(`RW-ACDC4` ~ delta * ir, subset_data_3,
                        panel = panel.levelplot.points, cex = 1.2,
                        xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                        ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

rw_acdc5_3 <- levelplot(`RW-ACDC5` ~ delta * ir, subset_data_3,
                        panel = panel.levelplot.points, cex = 1.2,
                        xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                        ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

rw_acdc1_2 <- levelplot(`RW-ACDC1` ~ delta * ir, subset_data_2,
                        panel = panel.levelplot.points, cex = 1.2,
                        xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                        ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

rw_acdc2_2 <- levelplot(`RW-ACDC2` ~ delta * ir, subset_data_2,
                        panel = panel.levelplot.points, cex = 1.2,
                        xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                        ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

rw_acdc3_2 <- levelplot(`RW-ACDC3` ~ delta * ir, subset_data_2,
                        panel = panel.levelplot.points, cex = 1.2,
                        xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                        ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

rw_acdc4_2 <- levelplot(`RW-ACDC4` ~ delta * ir, subset_data_2,
                        panel = panel.levelplot.points, cex = 1.2,
                        xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                        ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))

rw_acdc5_2 <- levelplot(`RW-ACDC5` ~ delta * ir, subset_data_2,
                        panel = panel.levelplot.points, cex = 1.2,
                        xlab = expression("Shifting Parameter " ~ (Delta)),  # Adds the delta symbol
                        ylab = "Global Imbalance Ratio (q)") + layer_(panel.2dsmoother(..., n = 200))


grid_plots1 <- c(p_acdc1_10, p_acdc2_10, p_acdc3_10, p_acdc4_10, p_acdc5_10,
                 p_acdc1_5, p_acdc2_5, p_acdc3_5, p_acdc4_5, p_acdc5_5,
                 p_acdc1_3, p_acdc2_3, p_acdc3_3, p_acdc4_3, p_acdc5_3,
                 p_acdc1_2, p_acdc2_2, p_acdc3_2, p_acdc4_2, p_acdc5_2,
                 layout = c(5,4))

grid_plots1


grid.text("Heatmaps of P-CCCD - ACDC",
          x = 0.5, y = .98,  # Centered at the top
          gp = gpar(fontsize = 14, fontface = "bold"))


grid_plots2 <- c(rw_acdc1_10, rw_acdc2_10, rw_acdc3_10, rw_acdc4_10, rw_acdc5_10,
                 rw_acdc1_5, rw_acdc2_5, rw_acdc3_5, rw_acdc4_5, rw_acdc5_5,
                 rw_acdc1_3, rw_acdc2_3, rw_acdc3_3, rw_acdc4_3, rw_acdc5_3,
                 rw_acdc1_2, rw_acdc2_2, rw_acdc3_2, rw_acdc4_2, rw_acdc5_2,
                 layout = c(5,4))


grid_plots2

grid.text("Heatmaps of RW-CCCD - ACDC",
          x = 0.5, y = .98,  # Centered at the top
          gp = gpar(fontsize = 14, fontface = "bold"))

