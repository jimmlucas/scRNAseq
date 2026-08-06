# Server execution

The production server profile uses SLURM plus Singularity. It keeps scheduler
settings outside the scientific analysis plan so the same approved policies
can be reproduced on another installation.

## Required software

- Nextflow 26.04.4 or newer
- Java supported by that Nextflow release
- SLURM client commands available on the submission node
- Singularity or Apptainer available on compute nodes
- Read access to the input data and write access to results, work and cache
  directories from every compute node
- Access to GHCR, or pre-pulled SIF images on shared storage

## Environment

Configure the installation before launching the workflow:

```bash
export BBUDY_SLURM_QUEUE=compute
export BBUDY_SLURM_ACCOUNT=my_account
export BBUDY_WORK_DIR=/shared/scratch/my_project/nextflow-work
export BBUDY_QUEUE_SIZE=20
export BBUDY_DOUBLET_CPUS=1
export NXF_SINGULARITY_CACHEDIR=/shared/containers/nextflow
```

`BBUDY_SLURM_QUEUE` and `BBUDY_SLURM_ACCOUNT` are optional when the cluster
provides suitable defaults. `BBUDY_WORK_DIR` should normally point to scratch
or another high-throughput shared filesystem. Do not put `work` on a small
home filesystem.

`BBUDY_DOUBLET_CPUS` controls the CPU allocation inside each DoubletFinder
scoring task. Its default is 1. There is no `maxForks` cap on these tasks:
Nextflow can run several samples concurrently and the local executor or SLURM
decides how many fit according to requested CPUs, memory, queue limits and
current capacity. Increase `BBUDY_DOUBLET_CPUS` only when internal
DoubletFinder sweep parallelism is preferred over sample-level parallelism.

## Preflight

Resolve the final configuration without submitting jobs:

```bash
nextflow config -profile slurm .
tests/server_profile/run_tests.sh
```

Confirm that the compute environment can pull both images. If compute nodes
cannot reach GHCR, pre-pull the images into `NXF_SINGULARITY_CACHEDIR` from an
allowed node.

## Run

The `slurm` profile enables Singularity itself, so a production invocation is:

```bash
nextflow run main.nf \
  --mode analysis_plan_qc_process \
  --analysis_plan configs/analysis_plan_test.yaml \
  --project_id PROJECT_ID \
  --project_name PROJECT_NAME \
  --biological_run_id RUN_ID \
  --execution_id EXECUTION_ID \
  --outdir /shared/results/PROJECT_ID \
  -profile slurm \
  -resume
```

For a non-SLURM server, use `-profile singularity` or `-profile apptainer` and
provide the executor in a site-specific configuration file:

```bash
nextflow run main.nf -profile apptainer -c conf/site.config ...
```

## Resource behavior

The default resource labels now start at:

| Label | Initial memory | Initial time |
|---|---:|---:|
| `lightweight` | 512 MB | 30 min |
| `medium_memory` | 8 GB | 8 h |
| `high_memory` | 16 GB | 12 h |
| `very_high_memory` | 32 GB | 72 h |

Memory failures (exit codes 137-140) and termination/time-limit signals (143)
are retried up to two times. Memory grows with every attempt and time doubles:
for example, a 72-hour task retries with 144 and then 288 hours. A scheduler
cannot continue a killed R process from the exact instruction where it
stopped; the task restarts, while completed Nextflow tasks remain reusable via
`-resume`. Cluster limits still take precedence.

## Operational checks

Before accepting a real run:

1. Verify every task is submitted to SLURM rather than executed on the login
   node.
2. Confirm `work`, results and the container cache are visible to compute
   nodes.
3. Interrupt a small run and verify that `-resume` reuses completed tasks.
4. Inspect the trace for peak RSS, CPU efficiency and runtimes.
5. Tune label resources using observed trace values before the full dataset.
