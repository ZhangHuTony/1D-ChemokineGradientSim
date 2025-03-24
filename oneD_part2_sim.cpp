// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Tumor class handles chemokine gradient and sigma calculation
class Tumor {
public:
  NumericVector cancer_x;
  NumericVector cancer_M;
  int n_cancer;
  double sqrt_d_over_D;
  double gradient_factor;
  
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
  
  // Calculate chemokine gradient at a given position x
  double get_gradient(double x) {
    double gradient = 0.0;
    for (int j = 0; j < n_cancer; ++j) {
      double dx = x - cancer_x[j];
      double abs_dx = std::abs(dx);
      double sgn_dx = (dx > 0) ? 1.0 : ((dx < 0) ? -1.0 : 0.0);
      gradient += cancer_M[j] * exp(-sqrt_d_over_D * abs_dx) * sgn_dx;
    }
    gradient *= gradient_factor;
    return gradient;
  }
  
  // Compute sigma based on the gradient and sensitivity parameter alpha
  double get_sigma(double gradient, double alpha) {
    double abs_grad = std::abs(gradient);
    double sigma = exp(alpha * abs_grad) / (1.0 + exp(alpha * abs_grad));
    return sigma;
  }
};

// TCell class handles a T cell's position and velocity
class TCell {
public:
  double position;
  double velocity;
  
  TCell(double pos, double vel) : position(pos), velocity(vel) {}
  
  // Update the T cell's state.
  // If always_update is true, the T cell always updates its velocity.
  // If false, it uses the provided random number (p) to decide whether to update.
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

// [[Rcpp::export]]
List run_simulation(int num_tcells, double b_min, double b_max,
                    NumericVector cancer_x, NumericVector cancer_M,
                    double D, double delta, double alpha, int num_steps,
                    bool always_update) {
  RNGScope rngScope;
  
  // Create a Tumor instance to manage chemokine gradient calculations
  Tumor tumor(cancer_x, cancer_M, D, delta);
  
  // Initialize T cells with random positions and velocities
  std::vector<TCell> tcells;
  tcells.reserve(num_tcells);
  for (int i = 0; i < num_tcells; ++i) {
    double pos = R::runif(b_min, b_max);
    double vel = R::runif(-1.0, 1.0);
    tcells.push_back(TCell(pos, vel));
  }
  
  // Main simulation loop: update each T cell for num_steps iterations
  for (int step = 0; step < num_steps; ++step) {
    std::vector<double> randoms;
    if (!always_update) {
      randoms.resize(num_tcells);
      // Pre-generate random numbers on the main thread to avoid thread safety issues
      for (int i = 0; i < num_tcells; ++i) {
        randoms[i] = R::runif(0.0, 1.0);
      }
    }
    
    // Parallelize the T cell updates using OpenMP
#pragma omp parallel for schedule(static)
    for (int i = 0; i < num_tcells; ++i) {
      if (always_update) {
        tcells[i].update(tumor, alpha, always_update);
      } else {
        tcells[i].update(tumor, alpha, always_update, randoms[i]);
      }
    }
  }
  
  // Collect final positions and velocities
  NumericVector final_positions(num_tcells);
  NumericVector final_velocities(num_tcells);
  for (int i = 0; i < num_tcells; ++i) {
    final_positions[i] = tcells[i].position;
    final_velocities[i] = tcells[i].velocity;
  }
  
  return List::create(
    _["final_positions"] = final_positions,
    _["final_velocities"] = final_velocities
  );
}


