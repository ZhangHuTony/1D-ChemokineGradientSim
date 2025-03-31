// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Tumor class: precomputes a lookup table and provides gradient lookups
class Tumor {
public:
  NumericVector cancer_x;
  NumericVector cancer_M;
  int n_cancer;
  double sqrt_d_over_D;
  double gradient_factor;
  
  // Lookup table parameters
  double lookup_start;
  double lookup_step;
  std::vector<double> lookup_table;
  
  Tumor(NumericVector cancer_x_, NumericVector cancer_M_, double D, double delta) {
    if(cancer_x_.size() != cancer_M_.size()){
      stop("Cancer cell positions and heats must have equal length");
    }
    cancer_x = cancer_x_;
    cancer_M = cancer_M_;
    n_cancer = cancer_x.size();
    sqrt_d_over_D = sqrt(delta / D);
    gradient_factor = -1.0 / (2.0 * D);
  }
  
  // Initialize the lookup table using precomputed gradients
  void initialize_lookup(double grid_start, double grid_step, NumericVector gradient_table) {
    this->lookup_start = grid_start;
    this->lookup_step = grid_step;
    int n = gradient_table.size();
    lookup_table.resize(n);
    for (int i = 0; i < n; i++) {
      lookup_table[i] = gradient_table[i];
    }
  }
  
  // Return the chemokine gradient at position x via lookup
  double get_gradient(double x) {
    int index = static_cast<int>(std::round((x - lookup_start) / lookup_step));
    if(index < 0) index = 0;
    if(index >= static_cast<int>(lookup_table.size()))
      index = lookup_table.size() - 1;
    return lookup_table[index];
  }
  
  // Compute sigma given a gradient and the sensitivity parameter alpha
  double get_sigma(double gradient, double alpha) {
    double abs_grad = std::abs(gradient);
    double expo = alpha * abs_grad;
    double cap = 700.0; // cap to avoid overflow in exp
    if (expo > cap) {
      expo = cap;
    }
    double exp_val = exp(expo);
    double sigma = exp_val / (1.0 + exp_val);
    return sigma;
  }
};

// TCell class: represents a T cell that updates its state based on the local gradient
class TCell {
public:
  double position;
  double velocity;
  
  TCell(double pos, double vel) : position(pos), velocity(vel) {}
  
  // Update T cell using the gradient lookup from Tumor
  void update(Tumor &tumor, double alpha, bool always_update, double p = -1.0) {
    double gradient = tumor.get_gradient(position);
    double sigma = tumor.get_sigma(gradient, alpha);
    double new_vel = velocity;
    
    if (always_update) {
      new_vel = (gradient > 0) ? sigma : -sigma;
    } else {
      if (p < 0) {
        p = R::runif(0.0, 1.0);
      }
      if (p < sigma) {
        new_vel = (gradient > 0) ? sigma : -sigma;
      }
    }
    velocity = new_vel;
    position += new_vel;
  }
};

// Precompute the chemokine gradient lookup table over a grid
// [[Rcpp::export]]
NumericVector precompute_gradient_table(double grid_start, double grid_end, double grid_step,
                                        NumericVector cancer_x, NumericVector cancer_M,
                                        double D, double delta) {
  Tumor tumor(cancer_x, cancer_M, D, delta);
  int n_points = std::ceil((grid_end - grid_start) / grid_step) + 1;
  NumericVector table(n_points);
  
  for (int i = 0; i < n_points; i++) {
    double x = grid_start + i * grid_step;
    double gradient = 0.0;
    for (int j = 0; j < tumor.n_cancer; j++) {
      double dx = x - tumor.cancer_x[j];
      double abs_dx = std::abs(dx);
      double sgn_dx = (dx > 0) ? 1.0 : ((dx < 0) ? -1.0 : 0.0);
      gradient += tumor.cancer_M[j] * exp(-tumor.sqrt_d_over_D * abs_dx) * sgn_dx;
    }
    gradient *= tumor.gradient_factor;
    table[i] = gradient;
  }
  return table;
}

// Run the T-cell simulation using the precomputed gradient lookup table.
// A new parameter grid_end is added so we can enforce the boundary condition.
// [[Rcpp::export]]
List run_simulation(int num_tcells, double b_min, double b_max,
                    NumericVector cancer_x, NumericVector cancer_M,
                    double D, double delta, double alpha, int num_steps,
                    bool always_update,
                    NumericVector gradient_table,
                    double grid_start, double grid_step, double grid_end) {
  RNGScope rngScope;
  
  Tumor tumor(cancer_x, cancer_M, D, delta);
  tumor.initialize_lookup(grid_start, grid_step, gradient_table);
  
  // Initialize T-cells with random positions (within b_min to b_max) and random velocities.
  std::vector<TCell> tcells;
  tcells.reserve(num_tcells);
  for (int i = 0; i < num_tcells; i++){
    double pos = R::runif(b_min, b_max);
    double vel = R::runif(-1.0, 1.0);
    tcells.push_back(TCell(pos, vel));
  }
  
  // Simulation loop
  for (int step = 0; step < num_steps; step++){
    std::vector<double> randoms;
    if (!always_update){
      randoms.resize(num_tcells);
      for (int i = 0; i < num_tcells; i++){
        randoms[i] = R::runif(0.0, 1.0);
      }
    }
    
#pragma omp parallel for schedule(static)
    for (int i = 0; i < num_tcells; i++){
      if (always_update)
        tcells[i].update(tumor, alpha, always_update);
      else
        tcells[i].update(tumor, alpha, always_update, randoms[i]);
    }
    
    // Boundary condition: if a T-cell moves outside grid_start or grid_end,
    // assign it a new random position between b_min and b_max.
    for (int i = 0; i < num_tcells; i++){
      if (tcells[i].position < grid_start || tcells[i].position > grid_end){
        tcells[i].position = R::runif(b_min, b_max);
      }
    }
  }
  
  // Collect final positions and velocities
  NumericVector final_positions(num_tcells);
  NumericVector final_velocities(num_tcells);
  for (int i = 0; i < num_tcells; i++){
    final_positions[i] = tcells[i].position;
    final_velocities[i] = tcells[i].velocity;
  }
  
  return List::create(_["final_positions"] = final_positions,
                      _["final_velocities"] = final_velocities);
}




