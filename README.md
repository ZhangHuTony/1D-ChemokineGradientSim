

# 1D and 1D\_Part2 Chemokine Gradient Simulation

This repository contains two subprojects focused on simulating T-cell behavior in response to chemokine gradients in heterogeneous tumor environments.

## Project Motivation

Understanding the spatial organization of T-cells in tumors provides insights into immunotherapy response. These simulations are part of a broader effort to computationally model immune-tumor interactions.

## Folders

* `1D/`: The original R-based simulation of chemokine gradients and T-cell positioning.
* `1D_Part2/`: A second-generation implementation using C++ for faster computation, gradient precomputation, and derivative-based updates.

## New in 1D\_Part2

`1D_Part2/` addresses performance and scalability limitations from the original implementation. Major improvements include:

* **C++ Core**: Rewritten simulation logic in C++ for significantly faster runtimes.
* **Precomputed Gradient**: The chemokine field gradient is calculated once and stored to avoid redundant calculations.
* **Derivative-Based T-cell Movement**: T-cell movement is now influenced by spatial chemokine derivatives instead of discrete sampling.



## Usage Instructions

Follow the individual README files in each subdirectory for detailed setup and usage information. In short:

1. Navigate to either `1D/` or `1D_Part2/`
2. Run the simulation using R or compile with `make` for the C++ version
3. Analyze results and compare T-cell distributions

## Summary

This project aims to simulate a chemokine gradient within heterogeneous tumors to better understand and replicate the spatial distribution of T-cells. The `1D` implementation uses R to model chemokine diffusion and compare simulated distributions against real tumor data using KL divergence as a metric. The follow-up implementation in `1D_Part2` enhances performance by leveraging C++ and introduces gradient precomputation and derivative-based cell movement.



## Acknowledgments

* Prof. Fred Adler
* Montana Ferita
* Dr. Melissa Reeves & Robert Letchworth
* Noah Moffat
* UROP, University of Utah ([link](https://our.utah.edu/research-scholarship-opportunities/urop/))

## References

> Tanaka, Miho et al. "Tumor cell heterogeneity drives spatial organization of the intratumoral immune response in squamous cell skin carcinoma." *bioRxiv*, doi:10.1101/2023.04.25.538140.

## Contact

For issues or questions, contact Tony Zhang at [zhanghutony@gmail.com](mailto:zhanghutony@gmail.com)
