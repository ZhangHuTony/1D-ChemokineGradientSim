#parameters
hot_heat = 2
cold_heat = 1
step_size = 1






# # Test data
# 
# ## parameters:
# 
# hot_size = 25
# cold_size = 0
# left_x = 0
# right_x = 1000
# 
# ## very hot data
# 
# 
# hot_data = data.frame(
#   pos_x = seq(0,1000, length.out = hot_size),
#   heat = rep(hot_heat,hot_size)
# )
# 
# 
# ## cold data
# 
# 
# cold_data = data.frame(
#   pos_x = seq(0,1000, length.out = cold_size),
#   heat = rep(cold_heat,cold_size)
# )
# 
# ## combine data
# cancer_pos = rbind(hot_data, cold_data)

# Actual Data
library(readxl) #package for reading excel files
tumor_3_raw = read_excel('Data/6235_3mm X coordinates for cross sections.xlsx', sheet = 1)

##function to split dataset into cancer and CD4 dataframes (removes any rows with multiple markers)
library(dplyr)
clean_raw_data <- function(df) {
  # Clean the data: Trim spaces and convert to lowercase
  df <- df %>%
    mutate(
      `AF647_clean` = trimws(tolower(`AF647?`)),
      `YFP_clean` = trimws(tolower(`YFP?`)),
      `RFP_clean` = trimws(tolower(`RFP?`))
    )
  
  # Filter CD4 cells (AF647)
  CD4_Cells <- df %>%
    filter( AF647_clean == "af647" & is.na(YFP_clean) ) #keep rows with only AF647
  
  # Filter Cancer cells (YFP or RFP)
  Cancer_Cells <- df %>%
    filter(
      (YFP_clean == "yfp" & is.na(RFP_clean) & is.na(AF647_clean) ) #keep rows with only YFP
      | # or
      (is.na(YFP_clean) & RFP_clean == "rfp" & is.na(AF647_clean) ) #keep rows with only RFP
    )
  
  return(list(CD4_df = CD4_Cells, Cancer_df = Cancer_Cells))
}

test = clean_raw_data(tumor_3_raw)
cancer_test = test$Cancer_df
cd4_test = test$CD4_df

# Chemokine Gradient Calculation

  #k: chemokine decay constant
k = 0.2
  #d: chemokine diffusion constant
d = 100  
  #M: heat of cancer cell
  #x: distance from cancer cell to current x position
calculate_diffusion_1d <- function(k, d, m, x) {
  
  ifelse(x < 0,
         (m / (2 * sqrt(d * k))) * exp(x * sqrt(k / d)),
         (m / (2 * sqrt(d * k))) * exp(-x * sqrt(k / d)))
}

chemokine_x_positions = seq(from=left_x, to=right_x, by = step_size)

chemokine_concentration = numeric(length(chemokine_x_positions))


for (row in 1:nrow(cancer_pos)) {
  y <- calculate_diffusion_1d(k = k, d = d, m = cancer_pos[row, "heat"], x = (chemokine_x_positions - cancer_pos[row, "pos_x"]))
  chemokine_concentration =  chemokine_concentration + y 
}

chemokine_gradient = data.frame(x_pos = chemokine_x_positions, concentration = chemokine_concentration)



# Simulation

## parameters
stoch = 1
iterations = 500

## functions


#function which returns the next movement in which the cell will move
get_next_step <- function(left_conc, right_conc){
  
  
  higher_conc_direction = ifelse(right_conc>left_conc, 
                                 1, 
                                 -1)

  #calculates the probablity that a cell will move in the direction with a higher concentration
  prob_move_higher_conc = 1 / (1 + exp(-1 * stoch * abs(left_conc - right_conc)))
  
  rand_num = runif(1)
  
  #if randomly generated number is less than probablity move in direction of higher concentration else move in opposite direction
  if(rand_num < prob_move_higher_conc){
    return(step_size * higher_conc_direction) 
  }else{
    return(step_size * (-1 * higher_conc_direction))
  }
  
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
num_tcells = 11

x_positions <- seq(left_x, right_x, length.out = num_tcells)

tcell_df <- data.frame(
  frame = rep(1, num_tcells),
  id = 1:num_tcells,
  x_pos = x_positions
)

tcell_df = merge(tcell_df, chemokine_gradient[,c("x_pos","concentration")], by = "x_pos", all.x = TRUE)

column_names = c("frame", "id", "x_pos","concentration")


run_sim <- function(agents){
  for (i in 1:iterations){
    
    previous_frame <- max(agents$frame)
    next_frame <- previous_frame + 1
    
    
    if (next_frame <= iterations){
      
      previous_frame_agents <- agents[agents$frame == previous_frame,]
      
      
      next_frame_agents <- data.frame(matrix(ncol = length(column_names), nrow = 0))
      colnames(next_frame_agents) <- column_names
      
      for (id in 1:num_tcells){
        
        previous_frame_agent <- previous_frame_agents[previous_frame_agents$id == id,]
        
        previous_pos_x <- previous_frame_agent$x_pos

        
        # Get concentration to right
        previous_pos_x_right <- get_new_torus_pos(previous_pos_x, step_size)
        previous_sum_conc_right <- chemokine_gradient$concentration[chemokine_gradient$x_pos == previous_pos_x_right]
        
        # Get concentration to left
        previous_pos_x_left <- get_new_torus_pos(previous_pos_x, -1*step_size)
        previous_sum_conc_left <- chemokine_gradient$concentration[chemokine_gradient$x_pos == previous_pos_x_left]
        
       
        delta_x <- get_next_step(previous_sum_conc_left, previous_sum_conc_right)

        
        new_pos_x = get_new_torus_pos(previous_pos_x, delta_x)
          
        next_agent <- data.frame(
          id = id,
          frame = next_frame,
          x_pos = new_pos_x,
          concentration = chemokine_gradient$concentration[chemokine_gradient$x_pos == new_pos_x]
        )
        next_frame_agents <- rbind(next_frame_agents, next_agent)
      }
      
      agents <- rbind(agents, next_frame_agents)
      
    }
  }
  
  return(agents)
  
}

#Visualizing the data
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
  animate(p, renderer = av_renderer(), nframes = iterations, fps = 5)
  anim_save("animation.mp4")
}



