# HDO_Data_Analysis

This repository contains the data and analysis files used to analyze the hydrodeoxygenation (HDO) of guaiacol on activated carbon supports. Files include the following:

- **HDO_Data.csv**
  - A CSV file containing the experimental HDO data.
- **Data_Import.jl**
  - Import the experimental data into a convenient data structure for analysis. 
- **Functions.jl**
  - Create ODE systems for the experimental reactions, as well as objective functions for optimization purposes.
- **BlackBox_Optimization.jl**
  - Loop through experimental data and perform black-box optimization using the resampling inheritence memetic search method.
- **Sensitivity_Analysis.jl**
  - Loop through sets of black-box optimization results to determine the system sensitivity to each parameter.
- **Plot_Data.xlsx**
  - An XLSX file containing organized data that will be used to generate plots.
- **Plotting.jl**
  - Create plots showing model-predicted state trajectories versus the experimental data, using the rate constants obtained through black-box optimization.