#!/usr/bin/env bash
#SBATCH -N 1 -n 20
#SBATCH -p hebbe
#SBATCH -t 12:00:00
#SBATCH -A C3SE2017-1-12

module load matlab/R2016a
matlab -nodesktop -nosplash -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -r run_sim_study_compare_px_rc &

wait