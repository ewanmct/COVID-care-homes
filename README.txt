# The zip file (COVID_CH) contains R code and Rdata files to generate the results from the manuscript 
# "COVID-19 in Scottish care homes: A metapopulation model of spread among residents and staff".
#
#
# -------------------------------------------------------------------------------------------------------
#  
# -------------------------------------------------------------------------------------------------------
#
# IN SHORT...
#
# 	The following R scripts generate the data for the figures in the manuscript
#
# 	- 3_CALIBRATION_FIXED_SEEDING.R
# 	- 5_CREATE_BARPLOT_DATA.R
# 	- 7_CREATE_ALL_HEATMAP_DATA.R
#
#       However, note that scripts 3 and 7 will not calculate the correct aggregated sum of squares as the care home 
#       case and death data has been randomised. The real care home data are available on request from Prof Bruce
#       Guthrie https://doi.org/10.1016/S2666-7568(20)30012-X. You must use that in place of the randomised data in
#       RANDOMISED_Care_Home_7DayAverageCases_NHS_Lothian.csv and RANDOMISED_Care_Home_WeeklyDeaths_NHS_Lothian 
#       to create the RData files used to make our figures in the manuscript. 
#       
# 	These correct RData files are included in the zip file if you wish to plot the figures from the 
# 	actual manuscript. Have all files from the zip file in the same directory, then run the following
# 	scripts to output the figures
#
# 	- 4_COMPARE_FIXED_SEEDING_FITS.R
# 	- 6_PLOT_BARPLOT.R
# 	- 8_PLOT_EPSILON_GAMMA_DELTA_HEATMAPS.R
# 	- 9_PLOT_OTHER_HEATMAPS.R.
#
# -------------------------------------------------------------------------------------------------------
#
# -------------------------------------------------------------------------------------------------------
#
# BRIEF DESCRIPTION OF EACH R SCRIPT
#
#
# READ_DATA.R :
#
#       Reads in surveillance data and manipulates it into the format for our use
#      - NHS Lothian weekly cases
#      - NHS Lothian weekly deaths
#      - NHS Lothian weekly cases in care home residents  (this data has been randomised. Real data available on request from Prof Bruce Guthrie https://doi.org/10.1016/S2666-7568(20)30012-X )
#      - NHS Lothian weekly cases in care home residents  (this data has been randomised. Real data available on request from Prof Bruce Guthrie https://doi.org/10.1016/S2666-7568(20)30012-X )
#
#
# 1_SEIRD_FUNCTIONS.R :
#
#      Contains functions used in 2_SEIRD.R
#
#
# 2_SEIRD.R :
#
#      Model parameters are fixed here. Running this script generates SEIARD time series data for 
#      each subpopulation
#
#
# 3_CALIBRATION_FIXED_SEEDING.R :
#
#      Running this script performs our calibration process and fit the model to data. 
#
#      The number of homes seeded is fixed here. The script then simulates (solves the ode system) over a list of 
#      permutations of the unfixed parameters (see manuscript for fixed parameters and assumptions). We use a 
#      HPC (ARCHIE-WeSt www.archie-west.ac.uk) and parallelisation. These parameter sets/scenarios are ranked 
#      by measuring the aggregated sum of squared error of the model output against the data. 
#
#      Running this script saves 
#     - all the parameter sets simulated over (scenario_matrix.RData)
#     - the aggregated least squares for every parameter set (Least_Squares_Errors.RData)
#     - model solution for the top ten best fitting scenarios (ones with lowest aggregated sum of 
#          squares error)
#     - model solution, time-dependent parameters, fitted parameter values for top scenario + all relevant data from
#         calibration process (e.g., seeded_4home_E_best_fit.RData) plus some extra for later use
#
# 
# 4_COMPARE_FIXED_SEEDING_FITS.R :
#
#      This script is to compare fitted models and their parameters (for each home seeded). It plots data made using 
#      3_CALIBRATION_FIXED_SEEDING.R.
#
#      Running the script outputs
#      - a plot of surveillance data and best-fitting model, "seeded_4home_plot.png" (note the randomised care home data)
#      - plot of corresponding Rt and visitation curve for the best fit (fitted time-dependent parameters), "Rt_and_visitation.png"
#      - a plot of the distribution of fitted parameters as a function of homes seeded "all_pars.png"
#      - a plot of the quality of fit as a function of homes seeded "Least_Square_Error_plot.png"
#
#
# 5_CREATE_BARPLOT_DATA.R :
#
#      Running this script makes data for plotting barplots (tornado plots) showing the model's sensitivity to its parameters.
#
#      All parameters are first set to their fitted/base case values for four homes seeded. 
#      We perturb individual parameters and measure the change in each population 
#      deaths. These results are then visualised as a barplot.
#
#      The script outputs
#      - final_deaths_change_storage_C_seeding_varied.RData", containing the barplot data
#
#
# 6_PLOT_BARPLOT.R :
#
#      Running this script outputs the plotted barplots ("single_fit_CSR_barplot_all_vertical.png") using data made in
#      5_CREATE_BARPLOT_DATA.R.
# 
#
# 7_CREATE_ALL_HEATMAP_DATA.R :
# 
#      Running this script creates data for the heatmap figures in the manuscript
#      - "results_epsilon_delta.RData", "results_epsilon_gamma.RData", "results_gamma_delta.RData", used for the heatmap
#             looking at sensitivity of the final resident deaths to the time-share/mixing parameters (ω^γ_high,δ,ε) 
#      - "results_beta_rr_end_beta_c_end.RData", used for the heatmap looking at natural log of aggregated sum of squared
#             error in a ω^C_end−ω^G_end parameter space.       
#
# 
# 8_PLOT_EPSILON_GAMMA_DELTA_HEATMAPS.R :
#
#      Running this script outputs the plotted heatmap looking at the sensitivity of the final resident deaths to the 
#      time-share/mixing parameters ("epsilon_gamma_delta_plots_tiles_contours_tag.png")
#
#
# 9_PLOT_OTHER_HEATMAPS.R :
#
#      Running this script outputs the plotted heatmap looking at the natural log of the aggregated sum of squared error in 
#      a ω^C_end−ω^G_end parameter space ("beta_rr_end_beta_c_end_LS_plot.png")
#
#
#
#
# END -----------------------------------------------------------------------------------------------------------------
