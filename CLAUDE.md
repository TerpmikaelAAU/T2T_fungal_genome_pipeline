# T2T fungal genome pipeline

Snakemake 9 workflow, run on AAU BioCloud through the SLURM executor plugin.
README.md is the user documentation; read it before changing a rule.

## Config
- `config/config.yaml` holds settings only, no comments. Every key is
  explained in README.md, "Configuring a run". A new key gets a row in the
  README's tables, not a comment in the config.
- The committed config must point at `example_data/`. Real samples (cluster
  paths) stay uncommitted: CI dry-runs the committed config and fails on a
  path that only exists on the cluster.

## Checking a change
Nothing here can run the tools; the cluster does that. What can be checked
locally:
- CI (`.github/workflows/dry-run.yml`): `py_compile` on `workflow/scripts/*.py`,
  then `snakemake -n` on the default config, with every optional feature on,
  and with every feature off. Run the same steps locally with Snakemake
  9.26.1 in a `uv venv` before pushing.
- CI never executes a script. To test one, e.g. `summarize.py`, `exec` it
  with a stub `snakemake` object (`input`, `params`, `output` as dicts with
  attribute access) against real output files the user provides.

## Gotchas
- Most intermediates are `temp()`. Adding a new consumer of one to a
  finished run makes Snakemake rebuild it, and that can cascade into
  rerunning hifiasm for every grid cell. Say so in the PR, and have the user
  `snakemake -n` before running.
- One job per grid cell, not one job over all cells, when a rule reads the
  per-cell FASTQs: otherwise every cell's reads must sit on disk at once.
