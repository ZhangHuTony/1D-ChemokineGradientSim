# Load required libraries
library(ggplot2)
library(dplyr)
library(gridExtra)

# ========== 1) Load Data ==========
# Adjust file paths if necessary
yfp <- read.csv("2D_Data/10mm_tumor/YFP_1D_projected.csv", header = TRUE)
rfp <- read.csv("2D_Data/10mm_tumor/RFP_1D_projected.csv", header = TRUE)

# Label and assign heat values: YFP = 1, RFP = 2
yfp$Type <- "YFP"
rfp$Type <- "RFP"
yfp$heat <- 1
rfp$heat <- 4

# Combine the cancer cell data
cancer_df <- rbind(yfp, rfp)
cancer_positions <- cancer_df$x

# ========== 2) Set Parameters ==========
D <- 10     # Diffusion-related parameter
delta <- 1  # Decay-related parameter

# ========== 3) Define x-range ==========
b_min <- min(cancer_positions)
b_max <- max(cancer_positions)
x_seq <- seq(b_min, b_max, length.out = 200)

# ========== 4) Define Concentration & Gradient Functions ==========

# c(x) = 1 / (2 * sqrt(delta * D)) * sum_i [ M_i * exp(-sqrt(delta/D)*|x - x_i|) ]
compute_concentration <- function(x, positions, heats, D, delta) {
  prefactor <- 1 / (2 * sqrt(delta * D))
  inside_sum <- heats * exp(-sqrt(delta / D) * abs(x - positions))
  prefactor * sum(inside_sum)
}

# ∇c(x) = -1/(2D) * sum_i [ M_i * exp(-sqrt(delta/D)*|x - x_i|) * sgn(x - x_i) ]
compute_gradient <- function(x, positions, heats, D, delta) {
  prefactor <- -1 / (2 * D)
  inside_sum <- heats * exp(-sqrt(delta / D) * abs(x - positions)) * sign(x - positions)
  prefactor * sum(inside_sum)
}

# ========== 5) Compute Values Over x_seq ==========
concentration_values <- sapply(
  x_seq,
  compute_concentration,
  positions = cancer_positions,
  heats = cancer_df$heat,
  D = D,
  delta = delta
)

gradient_values <- sapply(
  x_seq,
  compute_gradient,
  positions = cancer_positions,
  heats = cancer_df$heat,
  D = D,
  delta = delta
)

# Convert to data frames for plotting
conc_df <- data.frame(x = x_seq, concentration = concentration_values)
grad_df <- data.frame(x = x_seq, gradient = gradient_values)

# ========== 6) Determine Scaling Factors ==========
# We'll match each curve's maximum to the histogram's maximum density.
dens <- density(cancer_positions)
max_density <- max(dens$y)

max_conc <- max(abs(concentration_values))
max_grad <- max(abs(gradient_values))

scale_factor_conc <- max_density / max_conc
scale_factor_grad <- max_density / max_grad

# ========== 7) Build Plots with ggplot2 ==========

# --- (A) Plot for Chemokine Concentration ---
p_conc <- ggplot() +
  # Histograms for YFP (blue) and RFP (red)
  geom_histogram(data = filter(cancer_df, Type == "YFP"),
                 aes(x = x, y = ..density..),
                 fill = "green", alpha = 0.5, bins = 30) +
  geom_histogram(data = filter(cancer_df, Type == "RFP"),
                 aes(x = x, y = ..density..),
                 fill = "red", alpha = 0.5, bins = 30) +
  # Overlay the scaled concentration
  geom_line(data = conc_df,
            aes(x = x, y = concentration * scale_factor_conc),
            color = "black", size = 1) +
  scale_y_continuous(name = "Density",
                     sec.axis = sec_axis(~ . / scale_factor_conc,
                                         name = "Chemokine Concentration")) +
  ggtitle("Chemokine Concentration vs. Cancer Cell Distribution") +
  xlab("Position") +
  theme_minimal()

# --- (B) Plot for Chemokine Gradient ---
p_grad <- ggplot() +
  # Histograms for YFP (blue) and RFP (red)
  geom_histogram(data = filter(cancer_df, Type == "YFP"),
                 aes(x = x, y = ..density..),
                 fill = "green", alpha = 0.5, bins = 30) +
  geom_histogram(data = filter(cancer_df, Type == "RFP"),
                 aes(x = x, y = ..density..),
                 fill = "red", alpha = 0.5, bins = 30) +
  # Overlay the scaled gradient
  geom_line(data = grad_df,
            aes(x = x, y = gradient * scale_factor_grad),
            color = "black", size = 1) +
  scale_y_continuous(name = "Density",
                     sec.axis = sec_axis(~ . / scale_factor_grad,
                                         name = "Chemokine Gradient")) +
  ggtitle("Chemokine Gradient vs. Cancer Cell Distribution") +
  xlab("Position") +
  theme_minimal()

# ========== 8) Arrange and Display Both Plots ==========
grid.arrange(p_conc, p_grad, nrow = 2)


