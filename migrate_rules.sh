#!/usr/bin/env bash
# Migrate workflow/rules/*.smk from Snakemake 7 to Snakemake 9 conventions.
# Run ONCE from the repo root:   bash migrate_rules.sh
#
# What it changes:
#   1. resources[...]["time"]        -> resources[...]["runtime"]   (KeyError fix)
#   2. partition="gpu"/"general"/"shared" lines  -> removed (dead code on BioCloud)
#   3. gpus=1                        -> gres=GPU_GRES               (--gres=gpu:a10:1)
#   4. dorado_polish borrows porechop budget -> its own entry
#   5. envs/seqkit                   -> envs/seqkit.yml             (missing extension)

set -euo pipefail

if [[ ! -d workflow/rules ]]; then
  echo "ERROR: run this from the repo root (workflow/rules not found)." >&2
  exit 1
fi

# Keep a backup so this is reversible
backup="workflow/rules.backup.$(date +%Y%m%d-%H%M%S)"
cp -r workflow/rules "$backup"
echo "Backup written to: $backup"

python3 - <<'PY'
import re, os

d = "workflow/rules"
changed = []

for fn in sorted(os.listdir(d)):
    if not fn.endswith(".smk"):
        continue
    p = os.path.join(d, fn)
    orig = open(p).read()
    t = orig

    # 1. ["time"] -> ["runtime"]
    t = re.sub(r'resources\["(\w+)"\]\["time"\]', r'resources["\1"]["runtime"]', t)

    # 2. remove dead partition lines
    t = re.sub(r'\n\s*partition="(?:gpu|general|shared)",?', '', t)

    # 3. gpus=1 -> gres=GPU_GRES
    t = re.sub(r'\n(\s*)gpus=1,', r'\n\1gres=GPU_GRES,', t)

    # 5. seqkit env missing .yml
    t = t.replace('"../envs/seqkit"', '"../envs/seqkit.yml"')

    if t != orig:
        open(p, "w").write(t)
        changed.append(fn)

# 4. dorado_polish gets its own resource budget
p = os.path.join(d, "6_3_dorado_polish.smk")
if os.path.exists(p):
    t = open(p).read()
    t = t.replace('resources["porechop_api"]["mem_mb"]',
                  'resources["dorado_polish"]["mem_mb"]')
    t = t.replace('resources["porechop_api"]["runtime"]',
                  'resources["dorado_polish"]["runtime"]')
    open(p, "w").write(t)

print("Updated %d rule files:" % len(changed))
for c in changed:
    print("  " + c)
PY

echo
echo "Now verify with a dry run:"
echo "  snakemake -n --profile profile"
