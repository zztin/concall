jobs: 50
cluster: "sbatch -t {resources.time_min} --mem={resources.mem_mb}M -c {resources.cpus} -o logs_slurm/{rule}_{wildcards}.stout.log -e logs_slurm/{rule}_{wildcards}.error.log"
use-conda: true
default-resources: [cpus=1, mem_mb=4000, time_min=30]
latency-wait: 240
retry-times: 3
keep-going: true
rerun-incomplete: true
