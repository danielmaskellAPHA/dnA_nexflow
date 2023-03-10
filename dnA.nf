#!usr/bin/env nextflow

params.outdir = null
params.host = "Chicken"
params.ref = "Influenza"
user = "$USER"
params.reads = null

include { MAP as MAP1 } from "/home/${user}/mnt/VI6Bioinformatics/Nextflow/denovoAssembly/modules/map/main.nf"
include { MAP as MAP2 } from "/home/${user}/mnt/VI6Bioinformatics/Nextflow/denovoAssembly/modules/map/main.nf"
include { MAP as MAP3 } from "/home/${user}/mnt/VI6Bioinformatics/Nextflow/denovoAssembly/modules/map/main.nf"
include { MAP as MAP4 } from "/home/${user}/mnt/VI6Bioinformatics/Nextflow/denovoAssembly/modules/map/main.nf"
include { VCF as VCF1 } from "/home/${user}/mnt/VI6Bioinformatics/Nextflow/denovoAssembly/modules/map/main.nf"
include { VCF as VCF2 } from "/home/${user}/mnt/VI6Bioinformatics/Nextflow/denovoAssembly/modules/map/main.nf"
include { VCF as VCF3 } from "/home/${user}/mnt/VI6Bioinformatics/Nextflow/denovoAssembly/modules/map/main.nf"

check_params()

if( !params.reads || params.reads instanceof Boolean ) error "ERROR: Missing --reads parameter, check --help for usage"
if( !params.outdir || params.outdir instanceof Boolean ) error "ERROR: Missing --outdir parameter, check --help for usage"
if( params.host == true ) error "ERROR: Host parameter raised but not specified. Input a host, or omit this for default (Chicken)"
if( params.ref == true ) error "ERROR: Reference parameter raised but not specified. Input a reference, or omit this for default (Influenza)"

pattern = "${params.reads}/*_{S*_R1,S*_R2}*.fastq.gz"

if( params.host == "Chicken" || params.host == "chicken" )
    host = "/home/${user}/mnt/VI6Bioinformatics/Central_Pipelines/denovoAssembly/Host_Genomes/GallusGallusGenome/GCA_000002315.5_GRCg6a_genomic.fna.gz"
else
    host = params.host
    
if( params.ref == "Influenza" || params.ref == "influenza" || params.ref == "flu" || params.ref == "Flu" )
    ref = "/home/${user}/mnt/VI6Bioinformatics/Central_Pipelines/denovoAssembly/Viral_Reference_Databases/Influenza_Database/"
else
    ref = params.ref
    


log.info """\
    V I 6 - N F - P I P E L I N E
    =============================
    outdir    : ${params.outdir}
    host      : ${params.host}
    hostdir   : ${host}
    virus     : ${params.ref}
    virusdir  : ${ref}
    reads     : ${params.reads}
    user      : ${user}
    """
    .stripIndent()
    
Channel
    .fromFilePairs( pattern )
        .set{ read_pairs_ch }

    
process HOSTMAP {

    tag "Mapping ${sampleid} to host"
     
    errorStrategy 'ignore'
    
    conda "/home/${user}/miniconda3/envs/denovoAssembly-v2"
    
    input:
    tuple val(sampleid), path(reads)
    
    output:
    tuple val(sampleid), path("${sampleid}_NotHostReads.bam"), path(reads)
    
    script:
    """
    bwa mem -t ${task.cpus} ${host} ${reads[0]} ${reads[1]} | samtools view -@ ${task.cpus} -b -f 4 -o "$sampleid"_NotHostReads.bam -;
    """
}
    
process SUBSAMPLE {

    tag "Subsampling ${sampleid} non-host reads"
    
    errorStrategy 'ignore'
    
    conda "/home/${user}/miniconda3/envs/denovoAssembly-v2"
    
    input:
    tuple val(sampleid), path(nonhost), path(reads)
    
    output:
    tuple val(sampleid), path("*NotHostReads*"), path(reads)    
    
    shell:
    '''
    NonHostReads=$(samtools view -c !{nonhost})
    MaxReads=1000000
    if (( NonHostReads > MaxReads ));
        then	
            samtools view -@ 4 -b -s 0.30 !{nonhost} > !{sampleid}_nonHost_subsample.bam
            BAMFile=!{sampleid}_nonHost_subsample.bam
            rm "$(readlink -f !{nonhost})"
        else
            BAMFile=!{nonhost}
    fi
    
    samtools sort -@ !{task.cpus} -n -O BAM $BAMFile | samtools fastq -@ !{task.cpus} -1 !{sampleid}_NotHostReads1.fastq -2 !{sampleid}_NotHostReads2.fastq -s !{sampleid}_NotHostSingletons.fastq -;
    rm "$(readlink -f $BAMFile)"
    '''
}

process SPADES {

    publishDir "${params.outdir}/${sampleid}", mode: 'copy', pattern: 'SPAdes/*.fasta'
    
    tag "Running SPAdes on ${sampleid}"
    
    errorStrategy 'ignore'
    
    conda "/home/${user}/miniconda3/envs/denovoAssembly-v2"
    
    input:
    tuple val(sampleid), path(nh_reads), path(reads)
    
    output:
    tuple val(sampleid), path("SPAdes/*contigs.fasta"), path(reads)
    
    script:
    """
    spades.py --threads ${task.cpus} --only-assembler --pe1-1 ${nh_reads[0]} --pe1-2 ${nh_reads[1]} -o SPAdes
    
    mv SPAdes/contigs.fasta SPAdes/${sampleid}_SPAdes_contigs.fasta
    """

}

process BLAST {
    
    conda "/home/${user}/miniconda3/envs/denovoAssembly-v2"
    
    errorStrategy 'ignore'
    
    publishDir "${params.outdir}/${sampleid}/BLAST_Hits", mode: 'copy', pattern: '*_crunch.txt'
    
    input:
    tuple val(sampleid), path(contigs), path(reads)
    
    output:
    tuple val(sampleid), path("*_crunch.txt"), path(reads)
    
    shell:
    '''
    for file in !{ref}*.fasta
    do
        seg=$(basename $file)
        segname=${seg%%.*}
        
        blastn -db $file -query !{contigs} -out !{sampleid}_"$segname"_crunch.txt -evalue 0.0001 -max_target_seqs 5 -outfmt 6 -num_threads 4 &
    done
   '''
}

process TOPMATCH {
    
    conda "/home/${user}/miniconda3/envs/denovoAssembly-v2"
    
    errorStrategy 'ignore'
    
    input:
    tuple val(sampleid), path(crunches), path(reads) 
    
    output:
    tuple val(sampleid), path("*_match.fas"), path(reads)  
    
    shell:
    '''
    for crunch in !{crunches}
    do
	topseg=${crunch%%.*}
	part=${topseg#!{sampleid}_}
	segname=${part%_*}

	sort -k12,12 -rn "$crunch" | head -1 - | awk '{print $2}' - > "$topseg"_match.txt

	    blastdbcmd -db !{ref}*"$segname".fasta -entry_batch "$topseg"_match.txt > "$segname"_match.fas &
    done
    '''    
}

process CAT {

    conda "/home/${user}/miniconda3/envs/denovoAssembly-v2"
    
    publishDir "${params.outdir}/${sampleid}", mode: 'copy', pattern: '*.fa'
    
    errorStrategy 'ignore'
    
    input:
    tuple val(sampleid), path(matches), path(reads)
    
    output:
    tuple val(sampleid), path("*.fa"), path(reads)
    
    shell:
    '''
    cat !{matches} > top_matches.fa
    
    if [ ! -s top_matches.fa ] ; then
    echo "*****Error - top_matches.fa is empty. No contigs had blast matches*****"
    exit
    fi   
    '''
}

process CONSENSUS {

    conda "/home/${user}/miniconda3/envs/denovoAssembly-v2"
    
    publishDir "${params.outdir}/${sampleid}", mode: 'copy'
    
    input:
    tuple val(sampleid), path(realign), path(rfile), path(reads)
    val count
    
    output:
    tuple val(sampleid), path("*.fasta")
    tuple val(sampleid), path(rfile)
    path("*.txt")
    path("*mapOnly*")
    
    
    shell:
    '''
    rfile=!{rfile}

	ref=$(basename "$rfile")
	refname=${ref%%_*}
	reffile=${ref%%.*}
	
    samplename=!{sampleid}
    
    samtools sort -n -@ "$threads" !{realign} -o - | \
    samtools fixmate -r -m 	-@ "$threads" - - | \
    samtools sort -@ "$threads" -o - | \
    samtools markdup -r -s -O BAM -@ "$threads" - "$samplename"_iter!{count}_clean_mapOnly.bam ;
		
    samtools index -@ "$threads" "$samplename"_iter!{count}_clean_mapOnly.bam

    mkdir genconsensus_Results
    
    cd genconsensus_Results
    
    python /home/!{user}/mnt/VI6Bioinformatics/Central_Pipelines/Utils/genconsensus/genconsensus.py  -t "0" -m "1" -n "-" -r ../"$rfile" -b ../"$samplename"_iter!{count}_clean_mapOnly.bam
    mv final_consensus.fasta ../"$samplename"_iter!{count}_consensus.fasta
    
    cd ..		
    
    rm -r genconsensus_Results
    
    samtools flagstat "$samplename"_iter!{count}_realign.bam > "$samplename"_iter!{count}_MappingStats.txt
    
    rm "$(readlink -f !{realign})"
    rm "$(readlink -f !{reads[0]})"
    rm "$(readlink -f !{reads[1]})"
    '''
    
}

workflow{
    HOSTMAP(read_pairs_ch)
    (bam, reads) = HOSTMAP.out
    ext_reads = SUBSAMPLE(bam)
    dn_contigs = SPADES(ext_reads)
    matches = BLAST(dn_contigs)
    top_matches = TOPMATCH(matches)
    dn_ref = CAT(top_matches)
    map1 = MAP1(dn_ref, 1)
    con1 = VCF1(map1, 1)
    map2 = MAP2(con1, 2)
    con2 = VCF2(map2, 2)
    map3 = MAP3(con2, 3)
    con3 = VCF3(map3, 3)
    map_final = MAP4(con3, 4)
    con_final = CONSENSUS(map_final, 4)            
}
   
def check_params() {

    if( params.remove('help') ) {
        println """\
        	--outdir   : specify output directory
        	--host     : specify host genome location or preset (default: Chicken, --hosts to view)
        	--ref      : specify reference genome location or preset (default: Influenza, --refs to view)
        	--reads    : specify location of reads
        	"""
        	.stripIndent()
        exit 0
    }
    
    if( params.remove('hosts') ) {
        println """\
        	     H O S T     L I S T     
        	=============================
        	'Chicken'    : Gallus Gallus
        	
        	that's it so far, sorry :) - 27/02/23
        	"""
        	.stripIndent()
        exit 0
    }
    
    if( params.remove('refs') ) {
        println """\
        	R E F E R E N C E     L I S T
        	=============================
        	'Influenza'    : Influenza A
        	
        	that's it so far, sorry :) - 27/02/23 DM
        	"""
        	.stripIndent()
        exit 0
    }

}
