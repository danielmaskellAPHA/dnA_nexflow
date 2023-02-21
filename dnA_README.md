======V I 6 - N F - A S S E M B L Y - P I P E L I N E - V 1.0 ======
	
	** v1.0 - 08/02/2022 ** -- working build
	** v1.0.1 - 21/02/2022 ** -- fixed final consensus ordering and output issues

Adapted from original denovoAssembly.sh by Dan Maskell VI6 (daniel.maskell@APHA.gov.uk)

For use within VI6, or as proof-of-concept for other workgroups.

Note: as this is intended for use within VI6 there are several architecture-specific paths that are specific to our institution. These will be added as config options at a later date. Furthermore, there are a few dependency scripts that are not included with this repo.

====================================================================

Using paired-end illumina short-reads, this pipeline will:

- Map to specified host
- Extract predicted host-reads
- Run SPAdes de novo assembly pipeline
- BLAST contigs against local influenza database
- Construct hypothetical reference (segment mix & match from previous, real strains)
- Map reads to hypothetical reference
- Call tentative consensus
- Map reads to consensus
- Repeat for 4 iterations

Will publish:

- Contigs
- BLAST results
- 3rd iteration consensus
- Final map with stats
- Final consensus

Poor/non-Influenza reads will halt when attempting to construct a hypothetical reference.

Nextflow produces a "work" directory as it runs. This pipeline is written to remove sizeable mapping files from this directory as soon as they are obsolete. 
This means runs cannot be resumed. 
Unfortunately due to current system architecture this is an unavoidable workaround (running many reads in parrallel requires enormous disk-space otherwise).

While the work directory is likely to remain fairly small, it is good practise to delete it after completed runs. Do not delete it during a run.


-- Dependencies --

VI6 "denovoAssembly-v3" conda environment
Nextflow (contained in above env)

If nextflow is not root installed, this pipeline must be run with the conda env active.
Regardless, env must installed as the pipeline will try to activate the env as it runs.

Access to VI6Storage FSx drive
Access to VI6Bioinformatics FSx drive

-- Usage --

nextflow run dnA.nf --reads [s3 location of reads] --submission [VI6-NGS-XX-XX] 

Failing to input reads or submission will result in an immediate error

Options:
--user - $USER path override -- likely to break everything
--host - host database file (.fna.gz) -- default chicken
--ref  - virus database directory -- default influenza

-- Future Builds --

Host/Ref English input from list of local databases (e.g. --host Human --ref Coronaviridae)
