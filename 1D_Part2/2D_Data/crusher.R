library(ggplot2)
library(gridExtra)

# Read the CSV files
CD4_df <- read.csv("10mm_tumor/CD4_Pos.csv")
RFP_df <- read.csv("10mm_tumor/RFP_Pos.csv")
YFP_df <- read.csv("10mm_tumor/YFP_Pos.csv")

# Function to project a point onto a given line y = mx + b
project_onto_line <- function(x, y, m, b) {
  # Compute projection scalar
  t <- (x + m * (y - b)) / (1 + m^2)
  return(t)  # Return vector of 1D projection coordinates
}

# Define the new line equation
m <- 0  # Slope
b <- 2400  # Intercept

# Plot each dataset separately
ggplot() +
  geom_point(data = YFP_df, aes(x = pos_x, y = pos_y), color = "green", alpha = 0.8) +
  geom_point(data = RFP_df, aes(x = pos_x, y = pos_y), color = "red", alpha = 0.8) +
  geom_point(data = CD4_df, aes(x = pos_x, y = pos_y), color = "blue", alpha = 0.8) +
  #geom_abline(slope = m, intercept = b, color = "black") +
  theme_minimal() +
  labs(title = "Scatterplot of CD4, RFP, and YFP Positions",
       x = "Position X",
       y = "Position Y") +
  theme(legend.position = "none")

# # Apply projection to each dataset
# CD4_proj <- project_onto_line(CD4_df$pos_x, CD4_df$pos_y, m, b)
# RFP_proj <- project_onto_line(RFP_df$pos_x, RFP_df$pos_y, m, b)
# YFP_proj <- project_onto_line(YFP_df$pos_x, YFP_df$pos_y, m, b)
# 
# # write.csv(data.frame(CD4_proj), "5mm_tumor/CD4_1D_projected.csv", row.names = FALSE, col.names = FALSE)
# # write.csv(data.frame(RFP_proj), "5mm_tumor/RFP_1D_projected.csv", row.names = FALSE, col.names = FALSE)
# # write.csv(data.frame(YFP_proj), "5mm_tumor/YFP_1D_projected.csv", row.names = FALSE, col.names = FALSE)
# 
# 
# # Convert to dataframes
# CD4_1D <- data.frame(proj_1D = CD4_proj, Type = "CD4")
# RFP_1D <- data.frame(proj_1D = RFP_proj, Type = "RFP")
# YFP_1D <- data.frame(proj_1D = YFP_proj, Type = "YFP")
# 
# # Combine projected points into a single dataframe for visualization
# proj_df <- rbind(CD4_1D, RFP_1D, YFP_1D)
# 
# # Plot the 1D probability density function
# p2 <- ggplot(proj_df, aes(x = proj_1D, fill = Type)) +
#   geom_density(alpha = 0.5, position = "identity") + 
#   theme_minimal() +
#   labs(title = "Histogram of 1D Projections onto Custom Line",
#        x = "Projected Coordinate",
#        y = "Frequency") +
#   scale_fill_manual(values = c("CD4" = "blue", "RFP" = "red", "YFP" = "green"))
# # Arrange plots side by side
# grid.arrange(p1, p2, ncol = 2)


