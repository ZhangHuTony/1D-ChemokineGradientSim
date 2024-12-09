







# Actual Data
CD4_cells = read.csv('Data/3mm_tumor/CrossSection1/CD4_XPos.csv')
YFP_cells = read.csv('Data/3mm_tumor/CrossSection1/YFP_XPos.csv')
RFP_cells = read.csv('Data/3mm_tumor/CrossSection1/RFP_XPos.csv')

#constants:
left_x = min(CD4_cells$pos_x, YFP_cells$pos_x, RFP_cells$pos_x, na.rm = TRUE)
right_x = max(CD4_cells$pos_x, YFP_cells$pos_x, RFP_cells$pos_x, na.rm = TRUE)


# Chemokine Gradient Calculation ----

#k: chemokine decay constant
#d: chemokine diffusion constant
#M: heat of cancer cell
#x: distance from cancer cell to current x position
calculate_diffusion_1d <- function(k, d, m, x) {
  
  ifelse(x < 0,
         (m / (2 * sqrt(d * k))) * exp(x * sqrt(k / d)),
         (m / (2 * sqrt(d * k))) * exp(-x * sqrt(k / d)))
}

step_size = 1 #minimum step size of agents

#k: chemokine decay constant
#d: chemokine diffusion constant
#r: rfp (cold) heat
#y: yfp (hot) heat 

#outputs: 
#'gradient': a dataframe with columns 'x_pos' and 'concentration' 

calculate_chemokine_gradient <- function(k = 0.2, d = 100, y = 2, r = 1){

  
  chemokine_x_positions = seq(from=left_x, to=right_x, by = step_size)
  
  chemokine_concentration = numeric(length(chemokine_x_positions))
  
  #calculate the chemokine for just the RFP cells
  for (row in 1:nrow(RFP_cells)) {
    diffusion_value <- calculate_diffusion_1d(k = k, d = d, m = r, x = (chemokine_x_positions - RFP_cells[row, "pos_x"]))
    chemokine_concentration =  chemokine_concentration + diffusion_value
  }
  
  # calculate the chemokine concentration for the YFP cells (adds it ontop of RFP)
  for (row in 1:nrow(YFP_cells)) {
    diffusion_value <- calculate_diffusion_1d(k = k, d = d, m = y, x = (chemokine_x_positions - YFP_cells[row, "pos_x"]))
    chemokine_concentration =  chemokine_concentration + diffusion_value
  }
  
  chemokine_gradient = data.frame(x_pos = chemokine_x_positions, concentration = chemokine_concentration)
  
  return(chemokine_gradient)
}

# T-cell Migration ----

## functions


#function which returns the next movement in which the cell will move
get_next_step <- function(left_conc, right_conc, stoch){

  
  higher_conc_direction = ifelse(right_conc>left_conc, 
                                 1, 
                                 -1)

  #calculates the probablity that a cell will move in the direction with a higher concentration
  prob_move_higher_conc = 1 / (1 + exp(-1 * stoch * abs(left_conc - right_conc)))
  
  rand_num = runif(1)
  
  #if randomly generated number is less than probability move in direction of higher concentration else move in opposite direction
  if(rand_num < prob_move_higher_conc){
    return(step_size * higher_conc_direction) 
  }else{
    return(step_size * (-1 * higher_conc_direction))
  }
  
}

#if the new x position goes past the boundaries, the cell will reappear randomly back in the tumor.
get_new_pos <- function(x_pos, delta_x){
  new_pos = x_pos + delta_x
  
  if((new_pos < left_x) | (new_pos > right_x)){
    new_pos = sample(left_x:right_x, 1)
  }
  
  return(new_pos)
}

#function which will return the x position after moving some delta x such that it will be treated as a torus. 
#If the new x position exceeds the right position will be placed on the left most position and vice versa.
get_new_torus_pos <- function(x_pos, delta_x){
  new_pos = x_pos + delta_x
  
  if(new_pos < left_x){
    new_pos = right_x
  }
  if(new_pos > right_x){
    new_pos = left_x
  }
  
  return(new_pos)
}

## Generate test T cell data
#'n': number of t cells
#'gradient': chemokine gradient that the t cells will live in.
generate_agents <- function(n, gradient){
  x_positions = sample(left_x:right_x, n, replace = TRUE) #randomly place tcells on tumor
  
  agents <- data.frame(
    id = 1:n,
    frame = rep(1, n),
    x_pos = x_positions
  )
  
  agents = merge(agents, gradient[,c("x_pos","concentration")], by = "x_pos", all.x = TRUE)
  
}




#function which runs the t_cell migration and keeps all x positions and concentrations for every frame.
run_migration_by_frame <- function(n, gradient, iterations, stoch){
  
  agents = generate_agents(n, gradient)
  
  column_names = c("frame", "id", "x_pos","concentration")
  
  
  for (i in 1:iterations){
    
    previous_frame <- max(agents$frame)
    next_frame <- previous_frame + 1
    
    
    if (next_frame <= iterations){
      
      previous_frame_agents <- agents[agents$frame == previous_frame,]
      
      
      next_frame_agents <- data.frame(matrix(ncol = length(column_names), nrow = 0))
      colnames(next_frame_agents) <- column_names
      
      for (id in 1:n){
        
        previous_frame_agent <- previous_frame_agents[previous_frame_agents$id == id,]
        
        previous_pos_x <- previous_frame_agent$x_pos

        
        # Get concentration to right
        previous_pos_x_right <- get_new_torus_pos(previous_pos_x, step_size)
        previous_sum_conc_right <- gradient$concentration[gradient$x_pos == previous_pos_x_right]
        
        # Get concentration to left
        previous_pos_x_left <- get_new_torus_pos(previous_pos_x, -1*step_size)
        previous_sum_conc_left <- gradient$concentration[gradient$x_pos == previous_pos_x_left]
        
       
        delta_x <- get_next_step(previous_sum_conc_left, previous_sum_conc_right, stoch)

        
        new_pos_x = get_new_pos(previous_pos_x, delta_x)
          
        next_agent <- data.frame(
          id = id,
          frame = next_frame,
          x_pos = new_pos_x,
          concentration = gradient$concentration[gradient$x_pos == new_pos_x]
        )
        next_frame_agents <- rbind(next_frame_agents, next_agent)
      }
      
      agents <- rbind(agents, next_frame_agents)
      
    }
  }
  
  return(agents)
  
}

# Function: run_migration_by_frame
# Purpose: Simulates the movement of agents over multiple frames based on a chemokine gradient
# Inputs:
#   - n: Number of agents
#   - gradient: Data representing the chemokine gradient
#   - iterations: Number of frames to simulate
#   - stoch: Level of stochasticity in agent movement
# Output: Returns the x positions of all agents after the final frame

run_migration <- function(n, gradient, iterations, stoch){
  
  # Generate initial agents with their properties and associate with the gradient
  agents <- generate_agents(n, gradient)
  agents$frame <- 1
  
  # Define the column names for consistency
  column_names = c("frame", "id", "x_pos","concentration")
  
  # Iterate through the specified number of iterations
  for (i in 2:iterations){
    
    # Get the agents from the previous frame
    previous_frame_agents <- agents
    # Calculate the next frame number
    next_frame <- max(previous_frame_agents$frame) + 1
    
    # Update each agent's position and concentration based on gradient and stochasticity
    for (id in 1:n){
      
      # Retrieve the agent's data from the previous frame
      previous_frame_agent <- previous_frame_agents[previous_frame_agents$id == id,]
      
      # Get the agent's current position
      previous_pos_x <- previous_frame_agent$x_pos
      
      
      # Get concentration to the right of the current position
      previous_pos_x_right <- get_new_torus_pos(previous_pos_x, step_size)
      previous_sum_conc_right <- gradient$concentration[gradient$x_pos == previous_pos_x_right]
      
      # Get concentration to the left of the current position
      previous_pos_x_left <- get_new_torus_pos(previous_pos_x, -1*step_size)
      previous_sum_conc_left <- gradient$concentration[gradient$x_pos == previous_pos_x_left]
      
      # Calculate the change in position based on gradient and stochasticity
      delta_x <- get_next_step(previous_sum_conc_left, previous_sum_conc_right, stoch)
      
      # Update the agent's position
      new_pos_x <- get_new_pos(previous_pos_x, delta_x)
      
      # Update the agent's data in the current frame
      previous_frame_agents[previous_frame_agents$id == id, "x_pos"] <- new_pos_x
      previous_frame_agents[previous_frame_agents$id == id, "concentration"] <- gradient$concentration[gradient$x_pos == new_pos_x]
      previous_frame_agents[previous_frame_agents$id == id, "frame"] <- next_frame
    }
    
    # Replace the agents with updated positions for the next iteration
    agents <- previous_frame_agents
  }
  
  # Return only the x positions of all agents
  return(agents$x_pos)
}




# KL-divergence ----



library(philentropy) #library which has KL calculation
#defines function which takes vectors of x positions and returns kl.
get_kl <- function(pos_p, pos_q){
  hist_breaks = seq(min(c(pos_p, pos_q)), max(c(pos_p, pos_q)), by = 1)
  
  hist_p = hist(pos_p, breaks= hist_breaks, plot = FALSE)
  hist_q = hist(pos_q, breaks = hist_breaks, plot = FALSE)
  
  dist_p = hist_p$counts / sum(hist_p$counts)
  dist_q = hist_q$counts / sum(hist_q$counts)
  
  x = rbind(dist_p, dist_q)
  
  kl = as.numeric(KL(x, unit = "log"))
  return(kl)
}


# Function to calculate KL divergence for each frame in agents data frame
calculate_kl_for_frames <- function(observed, simulated) {
  
  
  frames <- unique(simulated$frame)
  
  
  kl_divergences <- numeric(length(frames))
  
  
  for (i in seq_along(frames)) {
    frame <- frames[i]
    
    # Get x positions of agents in the current frame
    pos_q <- simulated$x_pos[simulated$frame == frame]
    
    
    kl_divergences[i] <- get_kl(observed$pos_x, pos_q)
  }
  
  
  kl_data <- data.frame(frame = frames, kl_divergence = kl_divergences)
  
  return(kl_data)
}



#Visualizing the data ----
library(ggplot2)

#See concentration and t_cells at specific frame
makePlotAtFrame <- function(agents, concentrations, Frame){
  agents_filtered = agents[agents$frame == Frame, ]
  
  ggplot()+
    geom_line(data = concentrations, aes(x = x_pos, y = concentration), color = "gray", linewidth = 1) +
    
    geom_point(data = agents_filtered, aes(x = x_pos, y = concentration), color = agents_filtered$id, size = 3) +
    
    labs(title = paste("Frame:", Frame), x = "X Position", y = "Chemokine Concentration")
}


library(gganimate)
library(av)

#Make animation of tcells over time
makeAnimation <- function(agents, concentrations) {
  # Ensure the 'frame' column is numeric
  agents$frame <- as.numeric(agents$frame)
  
  # Create the plot
  p <- ggplot() +
    # Line plot in the background from concentrations dataframe
    geom_line(data = concentrations, aes(x = x_pos, y = concentration), color = "gray", size = 1) +
    
    # Scatter plot for agents that will be animated
    geom_point(data = agents, aes(x = x_pos, y = concentration, color = factor(id)), size = 3) +
    
    # Add labels
    labs(title = "Frame: {frame}", x = "X Position", y = "Concentration") +
    
    # Animate over the 'frame' column
    transition_time(frame) +               
    ease_aes('linear')
  
  # Save the animation as a GIF
  animate(p, renderer = av_renderer(), nframes = length(unique(agents$frame)), fps = 5)
  anim_save("animation.mp4")
}

plot_distributions <- function(observed, simulated){
  # Find the final frame for 'df' and subset the data
  final_frame <- max(simulated$frame)
  print(final_frame)
  final_frame_data <- simulated[simulated$frame == final_frame, ]
  
  # Plot the density of 'x_pos' for the final frame from 'df'
  plot(density(final_frame_data$x_pos), main = "Density of x_pos and pos_x for Final Frame", 
       xlab = "Position", ylab = "Density", col = "blue")
  
  # Add the density of 'pos_x' from 'df2' to the plot
  lines(density(observed$pos_x), col = "red")
  

}



# Driver ---- 

# Function: calculate_migration_and_kl
# Purpose: Combines chemokine gradient calculation, agent migration simulation, and KL divergence calculation
# Inputs:
#   - k: Chemokine decay constant
#   - d: Chemokine diffusion constant
#   - m: Heat (intensity) of the cancer cell
#   - stoch: Level of stochasticity in agent movement
#   - observed_positions: X positions of the observed data (can be a vector or a single-column data frame)
#   - iterations: Number of frames to simulate (default = 100)
#   - num_runs: Number of times to repeat the migration simulation (default = 10)
# Output: Returns the average KL divergence between the simulated and observed distributions

calculate_migration_and_kl <- function(k, d, m, stoch, observed_positions, iterations = 50, num_runs = 10) {
  
  # If observed_positions is a data frame, extract the first column as a vector
  if (is.data.frame(observed_positions)) {
    observed_positions <- observed_positions[[1]]
  }
  
  # Derive the number of agents from the input positions
  n <- length(observed_positions)
  
  
  # Step 1: Calculate the chemokine gradient
  gradient <- calculate_chemokine_gradient(k, d, m)
  
  # Initialize a variable to store the total KL divergence
  total_kl_divergence <- 0
  
  # Step 2: Run the migration simulation multiple times
  for (i in 1:num_runs) {
    # Simulate the migration of agents
    simulated_positions <- run_migration(n, gradient, iterations, stoch)
    
    # Calculate the KL divergence for this run
    kl_divergence <- get_kl(observed_positions, simulated_positions)
    
    # Accumulate the KL divergence
    total_kl_divergence <- total_kl_divergence + kl_divergence
  }
  
  # Calculate the average KL divergence
  average_kl_divergence <- total_kl_divergence / num_runs
  
  return(average_kl_divergence)
}


# Function: optimize_parameters
# Purpose: Optimize k, d, m, and stoch to minimize KL divergence
# Inputs:
#   - input_positions: X positions of the observed data (can be a vector or a single-column data frame)
#   - iterations: Number of frames to simulate (default = 50)
#   - num_runs: Number of times to repeat the migration simulation (default = 10)
#   - bounds: A list specifying bounds for k, d, m, and stoch
# Output: Returns the optimized values of k, d, m, and stoch

optimize_parameters <- function(input_positions, iterations = 50, num_runs = 10, bounds) {
  
  # Objective function for optimization
  objective_function <- function(params) {
    k <- params[1]
    d <- params[2]
    m <- params[3]
    stoch <- params[4]
    
    calculate_migration_and_kl(k, d, m, stoch, input_positions, iterations, num_runs)
  }
  
  # Set up optimization using L-BFGS-B method with bounds
  result <- optim(
    par = c(mean(bounds$k), mean(bounds$d), mean(bounds$m), mean(bounds$stoch)),
    fn = objective_function,
    method = "L-BFGS-B",
    lower = c(bounds$k[1], bounds$d[1], bounds$m[1], bounds$stoch[1]),
    upper = c(bounds$k[2], bounds$d[2], bounds$m[2], bounds$stoch[2])
  )
  
  # Return optimized parameters
  return(list(k = result$par[1], d = result$par[2], m = result$par[3], stoch = result$par[4]))
}


# Define bounds for parameters (adjust these based on your context)
bounds <- list(
  k = c(0.01, 1.0),      # Bounds for chemokine decay constant
  d = c(50, 150),      # Bounds for chemokine diffusion constant
  m = c(1.0, 10.0),      # Bounds for cancer cell heat
  stoch = c(0.0, 2.0)    # Bounds for stochasticity
)

# Optimize the parameters
optimized_params <- optimize_parameters(
  input_positions = CD4_cells,
  iterations = 50,
  num_runs = 10,
  bounds = bounds
)

# Print the optimized parameters
cat("Optimized Parameters:\n")
cat("k:", optimized_params$k, "\n")
cat("d:", optimized_params$d, "\n")
cat("m:", optimized_params$m, "\n")
cat("stoch:", optimized_params$stoch, "\n")


print(
  calculate_migration_and_kl(
  k = optimized_params$k,
  d = optimized_params$d,
  m = optimized_params$m,
  stoch = optimized_params$stoch,
  observed_positions = CD4_cells,
  num_runs = 1000
  )
)



total_kl = 0
for(i in 1:1000){
  total_kl = total_kl + get_kl(CD4_cells$pos_x,sample(left_x:right_x, 36, replace = TRUE) )
}

print(total_kl/1000)

# Rename 'pos_x' to 'x_pos' to facilitate merging
CD4_cells <- CD4_cells %>%
  rename(x_pos = pos_x)

# Merge df2 with df1 to get the corresponding concentration values for df2 x_pos
df2 <- merge(CD4_cells, grad, by = "x_pos")

# Create the line plot with points
ggplot() +
  geom_line(data = grad, aes(x = x_pos, y = concentration), color = 'black', size = 0.5) +  # Line from df1
  geom_point(data = df2, aes(x = x_pos, y = concentration), color = 'blue', size = 2) +  # Points from df2 on the line
  labs(title = "Concentration vs X Position with Points from Second DataFrame",
       x = "X Position",
       y = "Concentration") +
  theme_minimal()  # Clean minimal theme

