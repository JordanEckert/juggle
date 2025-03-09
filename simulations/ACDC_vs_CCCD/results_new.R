library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)

# Import count_data and results_data first

### Count
df_long <- count_data %>%
  pivot_longer(
    cols = starts_with("tau"),
    names_to = "tau_label",
    values_to = "tau_value"
  )

df_long <- df_long %>%
  mutate(tau_label = gsub("tau", "", tau_label))

head(df_long)
df_long$tau_label <- as.factor(df_long$tau_label)

write.csv(df_long, "long_count_data.csv")

summary_table <- df_long %>%
  group_by(dim, delta, ir, k, tau_label) %>%
  summarize(
    mean_count = mean(tau_value),
    sd_count   = sd(tau_value),
    n          = n(),
    .groups = "drop"
  ) %>%
  arrange(dim, delta, ir, k, tau_label)

print(summary_table)

data <- df_long %>%
  mutate(proportion = tau_value / 200)

custom_labeller <- labeller(
  dim = function(x) paste("Dimension =", x),
  ir = function(x) paste("IR =", x)
)

ggplot(data, aes(x = as.numeric(as.character(tau_label)), y = proportion, color = as.factor(k))) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun = mean, geom = "point") +
  facet_grid(ir~dim, labeller = custom_labeller) +
  labs(
    title = "Performance of ACDC Across Different Settings",
    x = "Tau",
    y = "Proportion of Pilot Studies with Highest AUC",
    color = "k"
  ) +
  theme_minimal()


### Results
acdc_cols <- grep("ACDC", names(results_data), value = TRUE)

df_new <- results_data %>%
  mutate(across(all_of(acdc_cols), ~ PCCCD - .x, .names = "PCCCD-{.col}")) %>%
  mutate(across(all_of(acdc_cols), ~ RWCCCD - .x, .names = "RWCCCD-{.col}")) %>%
  select(dim, delta, ir, starts_with("PCCCD-"), starts_with("RWCCCD-"))

# View the transformed dataset
head(df_new)



# Reshape data to long format for easier plotting
df_long <- melt(df_new, id.vars = c("dim", "delta", "ir"),
                variable.name = "Comparison", value.name = "Difference")

df_long_PCCCD <- subset(df_long, grepl("PCCCD", df_long$Comparison))
df_long_RWCCCD <- subset(df_long, grepl("RWCCCD", df_long$Comparison))

ggplot(df_long_PCCCD, aes(x = delta, y = Difference, color = as.factor(ir), group = ir)) +
  geom_line() +
  facet_wrap(~ Comparison + dim, scales = "free") +
  labs(title = "Trends in PCCCD and RWCCCD Differences Across Delta",
       x = "Delta", y = "Difference", color = "IR") +
  theme_minimal()
