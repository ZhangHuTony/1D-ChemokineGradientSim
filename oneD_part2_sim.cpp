#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List run_simulation(int num_tcells, double b_min, double b_max,
                    NumericVector cancer_x, NumericVector cancer_M,
                    double D, double delta, double alpha, int num_steps) {
  
  // Input validation
  if (cancer_x.size() != cancer_M.size()) {
    stop("Cancer cell positions and heats must have equal length");
  }
  
  const int n_cancer = cancer_x.size();
  const double sqrt_d_over_D = sqrt(delta/D);
  const double gradient_factor = -1.0 / (2.0 * D);
  
  // Initialize T-cells
  NumericVector positions(num_tcells);
  NumericVector velocities(num_tcells);
  RNGScope rngScope;
  
  for(int i = 0; i < num_tcells; ++i) {
    positions[i] = R::runif(b_min, b_max);
    velocities[i] = R::runif(-1.0, 1.0);
  }
  
  // Main simulation loop
  for(int step = 0; step < num_steps; ++step) {
    NumericVector new_positions(num_tcells);
    NumericVector new_velocities(num_tcells);
    
    for(int i = 0; i < num_tcells; ++i) {
      const double current_x = positions[i];
      double gradient = 0.0;
      
      // Calculate chemokine gradient
      for(int j = 0; j < n_cancer; ++j) {
        const double dx = current_x - cancer_x[j];
        const double abs_dx = std::abs(dx);
        const double sgn_dx = (dx > 0) ? 1.0 : ((dx < 0) ? -1.0 : 0.0);
        gradient += cancer_M[j] * exp(-sqrt_d_over_D * abs_dx) * sgn_dx;
      }
      gradient *= gradient_factor;
      
      // Calculate velocity update
      const double abs_grad = std::abs(gradient);
      const double sigma = exp(alpha * abs_grad) / (1.0 + exp(alpha * abs_grad));
      const double p = R::runif(0.0, 1.0);
      double new_vel = velocities[i];  // Default to previous velocity
      
      if(p < sigma) {
        new_vel = (gradient > 0) ? sigma : -sigma;
      }
      
      // Update position and velocity
      new_positions[i] = current_x + new_vel;
      new_velocities[i] = new_vel;
    }
    
    positions = new_positions;
    velocities = new_velocities;
  }
  
  return List::create(
    _["final_positions"] = positions,
    _["final_velocities"] = velocities
  );
}
// [[Rcpp::export]]
NumericVector calculate_chemokine(const NumericVector& x_grid,
                                  const NumericVector& cancer_x,
                                  const NumericVector& cancer_M,
                                  double D, double delta) {
  const int n = x_grid.size();
  const int n_cancer = cancer_x.size();
  const double coeff = 1.0 / (2.0 * std::sqrt(D * delta));
  const double exponent_term = std::sqrt(delta/D);
  
  NumericVector concentration(n);
  
#pragma omp parallel for
  for(int i = 0; i < n; ++i) {
    double sum = 0.0;
    const double x = x_grid[i];
    
    for(int j = 0; j < n_cancer; ++j) {
      const double dx = std::abs(x - cancer_x[j]);
      sum += cancer_M[j] * std::exp(-exponent_term * dx);
    }
    
    concentration[i] = coeff * sum;
  }
  
  return concentration;
}