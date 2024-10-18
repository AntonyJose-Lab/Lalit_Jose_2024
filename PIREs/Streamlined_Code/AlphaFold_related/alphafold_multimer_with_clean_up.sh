#!/bin/bash
#SBATCH -N 1
#SBATCH -t 2-00:00:00
#SBATCH -c 8
#SBATCH --mem-per-cpu=6G
#SBATCH -p gpu
#SBATCH --gpus=a100:1
#SBATCH -o "slurm-%j.out"
#SBATCH --oversubscribe # This is the default on Zaratan
#SBATCH --mail-type=BEGIN,END,FAIL # Send email notifications
#SBATCH --mail-user=amjose@umd.edu # email TO   

. ~/.bashrc

WORKDIR="/home/amjose/scratch.amjose-prj"
cd $WORKDIR

## variable
FASTA_SEQ="name"

FASTA_FILE="${WORKDIR}/data/${FASTA_SEQ}.fasta"

## make directories
RESULTSROOT="${WORKDIR}/results"
if [ ! -d $RESULTSROOT ]; then
        mkdir $RESULTSROOT
fi
RESULTSDIR="results/${FASTA_SEQ}"
mkdir $RESULTSDIR

ALPHAFOLD_DIR=/software/alphafold/v2.3.2/alphafold
ALPHAFOLD_MINICONDA=/software/alphafold/v2.3.2/miniconda3
ALPHAFOLD_DATABASE_DIR=/scratch/zt1/software/alphafold/v4/alphafold_db

echo "Starting job $SLURM_JOB_ID"
date
hostname
pwd

# Activate the 'alphafold' conda environment:
source ${ALPHAFOLD_MINICONDA}/etc/profile.d/conda.sh
conda activate alphafold


module purge
module load hpcc/zaratan
module load gcc/9.4.0
module load cuda
module load openmm

# Generate features:
python ${ALPHAFOLD_DIR}/run_alphafold.py \
--fasta_paths=${FASTA_FILE} \
--output_dir=${RESULTSROOT} \
--data_dir=${ALPHAFOLD_DATABASE_DIR} \
--uniref90_database_path=${ALPHAFOLD_DATABASE_DIR}/uniref90/uniref90.fasta \
--uniref30_database_path=${ALPHAFOLD_DATABASE_DIR}/uniref30/UniRef30_2021_03 \
--mgnify_database_path=${ALPHAFOLD_DATABASE_DIR}/mgnify/mgy_clusters_2022_05.fa \
--template_mmcif_dir=${ALPHAFOLD_DATABASE_DIR}/pdb_mmcif/mmcif_files  \
--max_template_date=2100-01-01 \
--obsolete_pdbs_path=${ALPHAFOLD_DATABASE_DIR}/pdb_mmcif/obsolete.dat \
--use_gpu_relax=True \
--model_preset=multimer \
--uniprot_database_path= \
--use_precomputed_msas=True \
--models_to_relax=best \
--pdb_seqres_database_path= \
--bfd_database_path=${ALPHAFOLD_DATABASE_DIR}/bfd/bfd_metaclust_clu_complete_id30_$
--pdb_seqres_database_path=${ALPHAFOLD_DATABASE_DIR}/pdb_seqres/pdb_seqres.txt \
--uniprot_database_path=${ALPHAFOLD_DATABASE_DIR}/uniprot/uniprot.fasta

# Clean up to remove models with lower ranks.
cd results/${FASTA_SEQ}
rm ranked_*
rm unrelaxed_*
relaxed_root_name=$(find . -name "relaxed_*" | tail -c +11 | head -c -5)
mv result_${relaxed_root_name}.pkl _result_${relaxed_root_name}.pkl
rm result_*
mv _result_${relaxed_root_name}.pkl result_${relaxed_root_name}.pkl
cd $WORKDIR
