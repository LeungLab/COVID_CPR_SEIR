# COVID_CPR_SEIR

Description of code and .csv files accompanying “Modeling reductions in SARS-CoV-2 transmission and hospital burden achieved by prioritizing testing using a 
clinical prediction rule” by Jody R. Reimer, Sharia M. Ahmed, Benjamin Brintz, Rashmee U. Shah, Lindsay T. Keegan, Matthew J. Ferrari, Daniel T. Leung. 

Contact either Jody Reimer (reimer@math.utah.edu) or Sharia Ahmed (Sharia.M.Ahmed@utah.edu) with questions. 

PatientData.csv 
This .csv file contains (unrealistic) data for 100 theoretical patients, and is require to run CPR.R 

CPR.R 
R script used to fit the clinical prediction rules (main, ModelA, ModelB, ModelC). Runs with PatientData.csv. 

MainModel.csv, ModelA.csv, ModelB.csv, ModelC.csv 
Each of these .csv files is the output of the cross-validated dataset from a clinical prediction rule for the main model 
used in the manuscript as well as the alternate models considered (models A, B, and C). One of these is required to run both 
Clinical_impacts.R and SEIR_with_prioritized_testing.R 

Clinical_impacts.R R script used to produce Fig. 2. To run this script, either MainModel.csv, ModelA.csv, ModelB.csv, or ModelC.csv is required to be in the same folder. 

SEIR_with_prioritized_testing.R 
R script used to produce Fig. 3, Fig. 4 and Table 1. To run this script, either MainModel.csv, ModelA.csv, ModelB.csv, or ModelC.csv is required to be in the same folder. 
Parameters in this script were changed to also produce Table S2.
