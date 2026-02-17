# SLURM Cluster Configuration

## Account Information
- **Standard account:** `tdlong_lab`
- **GPU account:** `tdlong_lab_gpu`

## Module Loading
```bash
module load python/3.10.2
module load R/4.2.2
```

## Conda / Mamba Activation (BBDuk)

BBDuk (part of BBMap/BBTools) is installed via conda. Activate before use:

```bash
source ~/miniforge3/etc/profile.d/conda.sh
source ~/miniforge3/etc/profile.d/mamba.sh
mamba activate bbmap
```

This is already handled in `scripts/run_te_kmer_count.sh`.

## CPU Partitions

| Partition | Max memory per core | Default / Max runtime |
|-----------|--------------------|-----------------------|
| standard  | 6 GB               | 2 day / 14 day        |
| highmem   | 10 GB              | 2 day / 14 day        |
| hugemem   | 18 GB              | 2 day / 14 day        |

### Key Rules
1. **Memory:** Request explicitly with `--mem-per-cpu=XG`
2. **Scaling memory:** Add more cores (`--cpus-per-task=N`) or move to a higher partition
3. **Runtime:** Default is 2 days; extend with `--time=` (max 14 days)

### Resource Request Examples
```bash
# Standard partition: 2 cores × 6 GB = 12 GB total
--mem-per-cpu=6G --cpus-per-task=2

# Highmem partition: 4 cores × 10 GB = 40 GB total
--mem-per-cpu=10G --cpus-per-task=4

# Hugemem partition: 3 cores × 18 GB = 54 GB total
--mem-per-cpu=18G --cpus-per-task=3
```

### CPU Job Template
```bash
#!/bin/bash
#SBATCH --job-name=cpu_analysis
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=4:00:00

module load python/3.10.2

# Your commands here
```

### CPU Array Job Template
```bash
#!/bin/bash
#SBATCH --job-name=array_job
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=4:00:00
#SBATCH --array=1-20

module load python/3.10.2

# Use SLURM_ARRAY_TASK_ID in your script
python scripts/process.py --task_id ${SLURM_ARRAY_TASK_ID}
```

## GPU Partitions

| Partition  | Account          | Interruptible? | Notes                         |
|------------|------------------|----------------|-------------------------------|
| gpu        | tdlong_lab_gpu   | No             | ~$0.70/hr per GPU             |
| free-gpu   | tdlong_lab       | Yes            | Good for short runs (<5 min)  |
| gpu-debug  | tdlong_lab_gpu   | No             | ~15 min limit, often instant  |

### GPU Resource Notes
- GPU jobs get 1 CPU (3 GB) by default -- request more if preprocessing on CPU
- Request `--cpus-per-task=4-8` and `--mem-per-cpu=6G` for CPU-heavy preprocessing
- Use `module load python/3.10.2` for PyTorch (pip-installed includes CUDA)

### GPU Job Template
```bash
#!/bin/bash
#SBATCH --job-name=train_gpu
#SBATCH -A tdlong_lab_gpu
#SBATCH -p gpu
#SBATCH --gres=gpu:V100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=6G
#SBATCH --time=12:00:00

module load python/3.10.2

python scripts/train_model.py --epochs 100
```

## Quick Reference

```bash
# Standard CPU job
sbatch -A tdlong_lab -p standard --cpus-per-task=1 --mem-per-cpu=5G slurm_jobs/job.sh

# Higher memory CPU job (highmem)
sbatch -A tdlong_lab -p highmem --cpus-per-task=4 --mem-per-cpu=10G slurm_jobs/job.sh

# Production GPU job
sbatch -A tdlong_lab_gpu -p gpu --gres=gpu:V100:1 --cpus-per-task=8 --mem-per-cpu=6G -t 12:00:00 slurm_jobs/train.sh
```
