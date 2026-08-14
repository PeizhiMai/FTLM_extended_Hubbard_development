# 4x4 U=7 compressibility campaign

This directory fixes the production parameters used for the 16-twist CADES
campaign:

- `lx=ly=4`, `U=7`, `tx=ty=1`, `tp=V=0`
- `m=80`, exact momentum-block threshold `256`
- initial target `R=16` (with retained `R=8` reductions)
- `beta=2.857142857142857,12.5` and 281 points on `mu=[-3,4]`
- the 4x4 twist grid in `twists.csv`
- one exclusive `high_mem` node per twist, 16 cores, 220 GiB, four hours

`twist_005` is `(phix,phiy)=(0.25,0.25)` and has permanent seed
`5012360`. The largest block is `(Nup,Ndown,mx,my)=(8,8,0,0)`, with
dimension `10,353,252`.

## Why OpenMP, MKL, and no MPI

The production binary uses OpenMP over random vectors so all workers share one
sector basis and one sparse momentum-block Hamiltonian. MPI ranks on the same
node would independently rebuild those dominant data structures and multiply
memory use, so MPI is intentionally not used. Intel MKL is linked in sequential
mode (`MKL_NUM_THREADS=1`) for LAPACK work; the large-block cost is the sparse
Lanczos matvec, not dense LAPACK. The CADES build gate can compare Intel and GCC,
but production defaults to Intel oneAPI 2024.1 plus MKL 2024.1.

## Deployment and gate

Deployment refuses a dirty git tree and installs a commit-addressed source and
binary snapshot:

```bash
campaigns/lx4_ly4_u7/cades/deploy_to_cades.sh intel
```

On CADES, submit the largest-block probe:

```bash
CAMPAIGN_ROOT=/lustre/or-scratch24/scratch/9pm/ftlm_codex/lx4_ly4_u7_campaign \
  source/campaigns/lx4_ly4_u7/cades/submit_resource_probe.sh
```

After the job finishes, apply the gate:

```bash
python3 source/campaigns/lx4_ly4_u7/cades/validate_resource_probe.py
```

The validator requires all sample IDs `0-15`, dimension `10,353,252`, every
sample below 180 minutes, and peak RSS at most 210 GiB. RSS below 180 GiB keeps
16 concurrent samples; 180-210 GiB selects batches of eight. It writes
`RESOURCE_GATE_PASSED.json`, without which production submission is refused.

## Production waves and extension

Submit every incomplete twist concurrently:

```bash
python3 source/campaigns/lx4_ly4_u7/cades/submit_wave.py --samples 16
python3 source/campaigns/lx4_ly4_u7/cades/campaign_status.py --samples 16
```

Each four-hour job stops internally after 210 minutes, appends every completed
sample durably, and exits cleanly. Run the same command again for another wave.
To extend later, change only the target:

```bash
python3 source/campaigns/lx4_ly4_u7/cades/submit_wave.py --samples 32
```

Sample IDs `0-15` are reused; only `16-31` are generated. The same checkpoint
can still be reduced exactly at `R=8` or `R=16`.

## Reduction, averaging, and the x axis

Twists are combined at fixed `(beta,mu)`: average `n` and `kappa` over the 16
twists, then set `x=1-mean(n)`. This handles the fact that individual twists
have different densities at the same chemical potential without interpolating
or averaging at mismatched thermodynamic control parameters.

Once all twist CSVs exist:

```bash
python3 source/campaigns/lx4_ly4_u7/finalize_campaign.py \
  --campaign-root "$CAMPAIGN_ROOT" --samples 16
```

The finalizer requires two betas, 281 rows per beta, monotone density,
nonnegative compressibility, and complete `x=[0,0.35]` coverage. It writes the
equal-observable twist average and a publication-style PNG/PDF. Because beta and
mu are reducer inputs rather than checkpoint metadata, retained checkpoints can
produce new temperature or chemical-potential grids without more Lanczos work.
