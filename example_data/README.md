# Example data

Tiny, synthetic files for testing that the pipeline runs end-to-end --
random sequence, not a real organism. Useful for checking your setup
(containers pull, jobs schedule, rules connect) before pointing the
pipeline at real data. Don't expect a real assembly from these.

| File                      | Entry point | Notes                          |
|---------------------------|-------------|---------------------------------|
| `sample_a.fastq.gz`       | `fastq`     | 20 reads, gzipped              |
| `sample_b.fastq`          | `fastq`     | 20 reads, plain text            |
| `sample_c.bam`            | `bam`       | 20 unaligned reads, no move table |
| `sample_d_pod5/reads.pod5`| `pod5`      | 5 reads, random signal          |

All reads are random bases/signal, so they're far too small and too low
coverage for hifiasm to assemble (it will correctly skip them -- see
`hifiasm_min_input_mb` in `config/config.yaml`). `sample_c.bam` has no
real move table, so polishing on it will not produce a meaningful result.
`sample_d_pod5` has no real signal, so basecalling it will not produce
meaningful reads either.
