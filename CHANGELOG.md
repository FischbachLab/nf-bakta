# Change Log

## v0.0.3

- updated Nextflow workflow to DSL2
- updated bakta version to the v1.9.3
- updated the database to v5.1.
- added a hmmsearch step against the PFamA v36 database.
- by default, all bakta outputs are produced.

## v0.0.4

- updated hmmsearch output is filtered to only report hits with a minimum domain coverage of 0.4. Any overlapping hits are compared and the one with the lowest e-value is reported.
- added a stats file in hmmsearch to track if any features were completely removed due to the filtering step.
- updated the bakta version to v1.9.4
- updated the pyhmmer version from v0.10.11 to v0.10.15