library(ggplot2)
library(tidyr)
library(dplyr)

# Define N and compute the parameters for each case
N <- 1:10
data <- data.frame(
  N = N,
  Parameters_3 = N^2 + 2*N,
  Parameters_2 = N^2 + 2*N - (N - 1),
  Parameters_1 = N^2 + 2 - 2*(N - 1)
)

# Reshape data to long format for ggplot
long_data <- pivot_longer(data, cols = starts_with("Parameters"),
                          names_to = "Case", values_to = "Parameters")

# Relabel the cases
long_data$Case <- recode(long_data$Case,
                         Parameters_3 = "3 State-dependent",
                         Parameters_2 = "2 State-dependent",
                         Parameters_1 = "1 State-dependent")

# Define shapes
shape_map <- c("3 State-dependent" = 15,  # square
               "2 State-dependent" = 17,  # triangle
               "1 State-dependent" = 16)  # circle

# Plot
ggplot(long_data, aes(x = N, y = Parameters, color = Case, shape = Case)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = c("3 State-dependent" = "#CC5500",
                                "2 State-dependent" = "#404080",
                                "1 State-dependent" = "#87A96B")) +
  scale_shape_manual(values = shape_map) +
  labs(x = "Number of States", y = "Number of Parameters", color = "", shape = "") +
  scale_x_continuous(breaks = 1:10) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "transparent", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.margin = margin(t = 2, r = 4, b = 2, l = 4),
    legend.spacing.y = unit(1, "pt"),
    legend.text = element_text(margin = margin(t = -2))
  )
