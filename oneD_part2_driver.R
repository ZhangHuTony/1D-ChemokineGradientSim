# First, ensure you have these packages installed
library(Rcpp)
library(ggplot2)
library(gridExtra)

# Compile the simulation code from previous answer
Rcpp::sourceCpp('oneD_part2_sim.cpp') 

# ---- Data Loading ----
# Read all position files (adjust path as needed)
yfp <- read.csv("2D_Data/10mm_tumor/YFP_1D_projected.csv", header = TRUE)
rfp <- read.csv("2D_Data/10mm_tumor/RFP_1D_projected.csv", header = TRUE)
cd4 <- read.csv("2D_Data/10mm_tumor/CD4_1D_projected.csv", header = TRUE)

# ---- Parameter Setup ----
# Combine cancer cell positions with heats
cancer_x <- c(yfp$x, rfp$x)
cancer_M <- c(rep(1, nrow(yfp)), rep(2, nrow(rfp))) # YFP=1, RFP=2

# Determine simulation parameters
num_tcells <- nrow(cd4)
all_positions <- c(yfp$x, rfp$x, cd4$x)
b_min <- min(all_positions)
b_max <- max(all_positions)

# ---- Run Simulation ----
set.seed(123)  # For reproducibility
sim_results <- run_simulation(
  num_tcells = num_tcells,
  b_min = b_min,
  b_max = b_max,
  cancer_x = cancer_x,
  cancer_M = cancer_M,
  D = 1.5,      # May need to tune these parameters
  delta = 0.7,  # based on your system
  alpha = 2.5,
  num_steps = 100
)

# ---- Comparison Analysis ----
# Create comparison data frame
real_data <- data.frame(
  x = cd4$x,
  type = "Experimental"
)

sim_data <- data.frame(
  x = sim_results$final_positions,
  type = "Simulated"
)

combined_data <- rbind(real_data, sim_data)

# Visual comparison using density plots
ggplot(combined_data, aes(x = x, fill = type)) +
  geom_density(alpha = 0.5) +
  ggtitle("T-cell Position Distribution Comparison") +
  xlab("Position") + ylab("Density") +
  theme_minimal()


# Chemokine visualization
library(ggplot2)

# Create a grid of x-values
x_grid <- seq(b_min, b_max, length.out = 1000)

# Calculate concentration profile
conc_profile <- calculate_chemokine(
  x_grid = x_grid,
  cancer_x = cancer_x,
  cancer_M = cancer_M,
  D = 1.5,   # Match your simulation parameters
  delta = 0.7
)

# Create plot data
plot_data <- data.frame(
  x = x_grid,
  concentration = conc_profile
)

# Create cancer cell data
cancer_data <- data.frame(
  x = cancer_x,
  type = factor(rep(c("YFP", "RFP"), 
                    c(nrow(yfp), nrow(rfp)))),
  M = cancer_M
)
  

  