library(Rcpp)
library(ggplot2)

# Source the C++ simulation code
Rcpp::sourceCpp('oneD_part2_sim.cpp')

# ---- 1. Read in parameters for the chemokine gradient ----

# Load the datasets (adjust the file paths as necessary)
yfp <- read.csv("2D_Data/10mm_tumor/YFP_1D_projected.csv", header = TRUE)
rfp <- read.csv("2D_Data/10mm_tumor/RFP_1D_projected.csv", header = TRUE)
cd4 <- read.csv("2D_Data/10mm_tumor/CD4_1D_projected.csv", header = TRUE)

# Define cancer cell positions and corresponding heats:
# For example, assign YFP cells a heat of 1 and RFP cells a heat of 4.
cancer_x <- c(yfp$x, rfp$x)
cancer_M <- c(rep(1, nrow(yfp)), rep(4, nrow(rfp)))

# Gradient-specific parameters
D <- 100         # Diffusion coefficient
delta <- 0.1     # Decay rate

# Determine the domain from all positions (cancer cells and T-cells)
all_positions <- c(yfp$x, rfp$x, cd4$x)
b_min <- min(all_positions)
b_max <- max(all_positions)

# Apply padding to the domain
pad <- 50
grid_start <- b_min - pad
grid_end <- b_max + pad
grid_step <- 0.1

# ---- 2. Look for the CSV file of the precomputed gradient ----

# Create the directory if it doesn't exist
if(!dir.exists("precomputedGradients")){
  dir.create("precomputedGradients")
}

# Create a unique filename based on key gradient parameters
grad_filename <- sprintf("precomputedGradients/gradient_D%.1f_delta%.1f_bmin%.2f_bmax%.2f_resolution%.5f.csv", 
                         D, delta, b_min, b_max, grid_step)

# ---- 3. Compute and save the gradient table if the file does not exist ----

if(file.exists(grad_filename)){
  grad_table <- read.csv(grad_filename, header = TRUE)$gradient
  cat("Loaded precomputed gradient from file.\n")
} else {
  cat("Precomputed gradient not found. Computing gradient table...\n")
  grad_table <- precompute_gradient_table(grid_start, grid_end, grid_step,
                                          cancer_x, cancer_M, D, delta)
  write.csv(data.frame(gradient = grad_table), grad_filename, row.names = FALSE)
  cat("Gradient table computed and saved to file.\n")
}

# ---- 5. Run the T-cell simulation ----

# Simulation-specific parameters (only pertaining to the T-cell movement)
num_tcells <- nrow(cd4)  # number of T-cells from the CD4 dataset
alpha <- 0.01           # sensitivity parameter for T-cell movement
always_update <- FALSE   # flag to always update the T-cell velocity
num_steps <- 500        # number of simulation steps

set.seed(123)  # For reproducibility
sim_time <- system.time({
  sim_results <- run_simulation(
    num_tcells = num_tcells,
    b_min = b_min,
    b_max = b_max,
    cancer_x = cancer_x,
    cancer_M = cancer_M,
    D = D,
    delta = delta,
    alpha = alpha,
    num_steps = num_steps,
    always_update = always_update,
    gradient_table = grad_table,
    grid_start = grid_start,
    grid_step = grid_step,
    grid_end = grid_end  # pass grid_end for boundary checking
  )
})
cat("Simulation time: ", sim_time["elapsed"], " seconds.\n")

# ---- 6. Retrieve and display results ----
real_data <- data.frame(x = cd4$x, type = "Experimental")
sim_data <- data.frame(x = sim_results$final_positions, type = "Simulated")
combined_data <- rbind(real_data, sim_data)

ggplot(combined_data, aes(x = x, fill = type)) +
  geom_density(alpha = 0.5) +
  ggtitle("T-cell Position Distribution Comparison") +
  xlab("Position") + ylab("Density") +
  theme_minimal()


  

  