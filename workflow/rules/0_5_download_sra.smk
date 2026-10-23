# ============================================================================
#  sra samples: download the reads from ENA
# ============================================================================
# A sample with `type: sra` names one or more run accessions instead of a
# path. This rule fetches their FASTQs, md5-checked, as one .fastq.gz at the
# path the Snakefile gave the sample, and from there on it is a fastq sample:
# prepare_fastq (0_3_decompress.smk) reads this file like any user's FASTQ.
#
# The FASTQs come from ENA, which mirrors every public SRA (and DDBJ) run
# within a few days and serves them as ready-made .fastq.gz -- no
# fasterq-dump conversion, and nothing to install beyond curl. The script
# refuses runs that aren't Oxford Nanopore (hifiasm runs with --ont), and
# says so when ENA has no FASTQ for a run yet. ENA's file server is slow
# and drops long transfers, so the script resumes each file where it
# stopped, from data/sra/<sample>.fastq.gz.partial.<file>; those files survive
# a failed job, so a retry or a later run resumes too. See
# workflow/scripts/download_sra.sh.
#
# Kept, not temp(): a dataset is several GB, and deleting it after
# prepare_fastq would mean downloading it again whenever a later run needs
# the reads. Delete data/sra/ by hand once a sample's results/ are final.
rule download_sra:
    output:
        a = "data/sra/{input}.fastq.gz"
    params:
        script = os.path.join(workflow.basedir, "workflow/scripts/download_sra.sh"),
        accessions = lambda w: " ".join(SAMPLES[w.input]["accession"]),
    threads:
        1
    resources:
        mem_mb = 2000,
        runtime = 1440,
    log:
        "logs/download_sra/{input}.log"
    shell:
        """
        bash "{params.script}" "{output.a}" {params.accessions} > {log} 2>&1
        """
