// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

class Tumor {
public:
  NumericVector cancer_x;
  NumericVector cancer_M;
  int n_cancer;
  double sqrt_d_over_D;
  double gradient_factor;
  double lookup_start, lookup_step;
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
  
  void initialize_lookup(double grid_start, double grid_step, NumericVector gradient_table) {
    this->lookup_start = grid_start;
    this->lookup_step = grid_step;
    int n = gradient_table.size();
    lookup_table.resize(n);
    for (int i = 0; i < n; i++) {
      lookup_table[i] = gradient_table[i];
    }
  }
  
  double get_gradient(double x) {
    int index = static_cast<int>(std::round((x - lookup_start) / lookup_step));
    if(index < 0) index = 0;
    if(index >= static_cast<int>(lookup_table.size())) index = lookup_table.size() - 1;
    return lookup_table[index];
  }
  
  double get_sigma(double gradient, double alpha) {
    double abs_grad = std::abs(gradient);
    double expo = alpha * abs_grad;
    double cap = 700.0;
    if (expo > cap) expo = cap;
    double exp_val = exp(expo);
    double sigma = exp_val / (1.0 + exp_val);
    return sigma;
  }
};

class TCell {
public:
  double position;
  double velocity;
  double initial_position;
  double path_length;
  
  TCell(double pos, double vel) : position(pos), velocity(vel), initial_position(pos), path_length(0.0) {}
  
  void update(Tumor &tumor, double alpha, double speed, bool always_update, double p = -1.0) {
    double gradient = tumor.get_gradient(position);
    double sigma = tumor.get_sigma(gradient, alpha);
    double new_vel = velocity;
    
    if (always_update) {
      new_vel = (gradient > 0) ? sigma*speed : -sigma*speed;
    } else {
      if (p < 0) p = R::runif(0.0, 1.0);
      if (p < sigma) new_vel = (gradient > 0) ? sigma*speed : -sigma*speed;
    }
    velocity = new_vel;
    path_length += std::abs(new_vel);
    position += new_vel;
  }
  
  double get_sinuosity() {
    double displacement = std::abs(position - initial_position);
    if (displacement == 0.0) return NA_REAL;
    return path_length / displacement;
  }
};

// [[Rcpp::export]]
NumericVector precompute_gradient_table(double grid_start, double grid_end, double grid_step,
                                        NumericVector cancer_x, NumericVector cancer_M,
                                        double D, double delta) {
  auto t1 = std::chrono::high_resolution_clock::now();
  
  Tumor tumor(cancer_x, cancer_M, D, delta);
  int n_points = std::ceil((grid_end - grid_start) / grid_step) + 1;
  NumericVector table(n_points);
  
#pragma omp parallel for schedule(static)
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
  
  auto t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = t2 - t1;
  Rcpp::Rcout << "Gradient computation time: " << duration.count() << " seconds.\n";
  
  return table;
}

// [[Rcpp::export]]
List run_simulation(int num_tcells, double b_min, double b_max,
                    NumericVector cancer_x, NumericVector cancer_M,
                    double D, double delta, double alpha, double speed, int num_steps,
                    bool always_update,
                    NumericVector gradient_table,
                    double grid_start, double grid_step, double grid_end) {
  RNGScope rngScope;
  
  Tumor tumor(cancer_x, cancer_M, D, delta);
  tumor.initialize_lookup(grid_start, grid_step, gradient_table);
  
  std::vector<TCell> tcells;
  tcells.reserve(num_tcells);
  for (int i = 0; i < num_tcells; i++){
    double pos = R::runif(b_min, b_max);
    double vel = R::runif(-1.0, 1.0);
    tcells.push_back(TCell(pos, vel));
  }
  
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
        tcells[i].update(tumor, alpha, speed, always_update);
      else
        tcells[i].update(tumor, alpha, speed, always_update, randoms[i]);
    }
    
    for (int i = 0; i < num_tcells; i++){
      if (tcells[i].position < grid_start || tcells[i].position > grid_end){
        tcells[i].position = R::runif(b_min, b_max);
      }
    }
  }
  
  NumericVector final_positions(num_tcells);
  NumericVector final_velocities(num_tcells);
  NumericVector sinuosities(num_tcells);
  
  for (int i = 0; i < num_tcells; i++){
    final_positions[i] = tcells[i].position;
    final_velocities[i] = tcells[i].velocity;
    sinuosities[i] = tcells[i].get_sinuosity();
  }
  
  return List::create(_["final_positions"] = final_positions,
                      _["final_velocities"] = final_velocities,
                      _["sinuosities"] = sinuosities);
}



