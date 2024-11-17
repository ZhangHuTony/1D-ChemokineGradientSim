







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
chem_grad = calculate_chemokine_gradient()

simulated_cells = run_migration_by_frame(n = nrow(CD4_cells), gradient = chem_grad, iterations = 250, stoch = 1)

kl_by_frame = calculate_kl_for_frames(CD4_cells, simulated_cells)

plot(kl_by_frame$frame, kl_by_frame$kl_divergence, type = 'l')


