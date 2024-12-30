# 1D Chemokine Gradient Simulation in Heterogeneous Tumors

This repository contains the R code and data files for simulating a 1D chemokine gradient in heterogeneous tumors, aiming to replicate the observed T-cell distribution within tumors. This simulation is based on the research presented in ["Tumor cell heterogeneity drives spatial organization of the intratumoral immune response in squamous cell skin carcinoma"](https://pmc.ncbi.nlm.nih.gov/articles/PMC10168251/).

## Repository Structure

### 1. **R Scripts**
- **`TumorDataModel.R`**: The primary R script containing the simulation code. This script models the chemokine gradient, compares, and optimizes the model to the observed data.


### 2. **Data Folder**
- **`Data/`**: Contains datasets of the positions of tumor cells and their classified cell type as well as images of the tumor cross sections. These datasets are used to validate the simulation model.
- **`Data/DataTransformation.R`**: A utility script for transforming raw tumor data into the required format for `TumorDataModel.R`.

### 3. **Project Summary**
- **`ProjectSummary.pdf`**: A crude summary of the project, detailing the objectives, methods, results, and conclusions.

## Dependencies


- **R**: Ensure you have R installed on your system.
- **Required Packages**: The following R packages are needed for the simulation installable via `install.packages()`:
  - `ggplot2`
  - `philentropy`
  - `gganimate`
  - `av`




## Results
- The simulation aims to replicate the observed T-cell distributions within heterogeneous tumors. Outputs include plots and metrics to compare the simulated and actual distributions using the KL divergence.

## Acknowledgments
Special thanks to the following contributors:
- **Professor Fred Adler**: For guidance and support throughout this project.
- **Graduate Student Montana Ferita**: For assistance in developing the R code and analytical methods.
- **Dr. Melissa Reeves and Robert Letchworth**: For providing the tumor data used in this study.
- **Noah Moffat**: For laying the groundwork.
- **Undergraduate Research Opportunity Program (UROP), University of Utah**: For funding this project. Learn more about UROP [here](https://our.utah.edu/research-scholarship-opportunities/urop/).



## Citation
> Tanaka, Miho et al. “Tumor cell heterogeneity drives spatial organization of the intratumoral immune response in squamous cell skin carcinoma.” bioRxiv : the preprint server for biology 2023.04.25.538140. 21 Jun. 2023, doi:10.1101/2023.04.25.538140. Preprint.

If you use this code please give credit by citing this repository. If you use the data please give credit to the [Reeve's Lab](https://www.reeveslab.com/).

---
For any questions or suggestions, please contact Tony Zhang at zhanghutony@gmail.com.

