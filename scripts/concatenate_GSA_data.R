# Concatenate all the data from the Global Sensitivity Analysis found in GSA_results_final/data

# Create dataframe ------------------------------------------
m = 2000

RF_SBeq_total_data = matrix(NA, nrow = m, ncol = 9)
RF_CVeq_total_data = matrix(NA, nrow = m, ncol = 9)
RF_age_total_data = matrix(NA, nrow = m, ncol = 9)

#load in all data from slurm job array and concatenate results from 
# different slurm jobs
for(i in 0:19) {
  load(paste0("GSA_data_0523/RF_SBeq_data_", i ,".Rdata"))
  RF_SBeq_total_data[(1:100) + (i*100), ] = as.matrix(RF_SBeq_data)
}
colnames(RF_SBeq_total_data) = colnames(RF_SBeq_data)

for(i in 0:19) {
  load(paste0("GSA_data_0523/RF_cv_data_", i ,".Rdata"))
  RF_CVeq_total_data[(1:100) + (i*100), ] = as.matrix(RF_cv_data)
}
colnames(RF_CVeq_total_data) = colnames(RF_cv_data)

for(i in 0:19) {
  load(paste0("GSA_data_0523/RF_age_data_", i ,".Rdata"))
  RF_age_total_data[(1:100) + (i*100), ] = as.matrix(RF_age_data)
}
colnames(RF_age_total_data) = colnames(RF_age_data)

# merge spawning biomass stability, age-structure, and spawning biomass outputs into one dataframe
GSA_total_df = as.data.frame(RF_CVeq_total_data)
GSA_total_df$rockfish_avg = RF_SBeq_total_data[,1]
GSA_total_df$rockfish_age = RF_age_total_data[,1]

#save(GSA_total_df, file = "GSA_results_final/GSA_total_df.Rdata")
#save(RF_SBeq_total_data, file = "GSA_results_final/RF_SBeq_total_data.Rdata")
#save(RF_CVeq_total_data, file = "GSA_results_final/RF_CVeq_total_data.Rdata")
#save(RF_age_total_data, file = "GSA_results_final/RF_age_total_data.Rdata")
