#!/usr/bin/env bash
#SBATCH -N 1 -n 20
#SBATCH -p hebbe
#SBATCH -t 12:05:00
#SBATCH -A C3SE2017-1-12

module load Matlab/2016a
module add Matlab/2016a

matlab -nodesktop -nosplash -nojvm -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -nojvm -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -nojvm -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -nojvm -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -nojvm -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -nojvm -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -nojvm -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -nojvm -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -nojvm -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -nojvm -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -nojvm -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -nojvm -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -nojvm -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -nojvm -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -nojvm -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -nojvm -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -nojvm -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -nojvm -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -nojvm -r run_sim_study_compare_px_rc &
matlab -nodesktop -nosplash -nojvm -r run_sim_study_compare_px_rc &

wait