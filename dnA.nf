#!usr/bin/env nextflow

params.submission = null
params.host = "Chicken"
params.ref = "Influenza"
params.user = "$USER"
params.reads = null

include { MAP as MAP1 } from '/home/danmaskell/modules/map'
include { MAP as MAP2 } from '/home/danmaskell/modules/map'
include { MAP as MAP3 } from '/home/danmaskell/modules/map'
include { MAP as MAP4 } from '/home/danmaskell/modules/map'
include { VCF as VCF1 } from '/home/danmaskell/modules/map'
include { VCF as VCF2 } from '/home/danmaskell/modules/map'
include { VCF as VCF3 } from '/home/danmaskell/modules/map'

if( !params.reads || params.reads instanceof Boolean ) error "ERROR: Missing --reads parameter"
if( !params.submission || params.submission instanceof Boolean ) error "ERROR: Missing --submission parameter"

pattern = "${params.reads}*_{S*_R1,S*_R2}*.fastq.gz"

if( params.host == "Chicken" )
    host = "/home/${params.user}/mnt/VI6Bioinformatics/Central_Pipelines/denovoAssembly/Host_Genomes/GallusGallusGenome/GCA_000002315.5_GRCg6a_genomic.fna.gz"
else
    host = params.host
    
if( params.ref == "Influenza" )
    ref = "/home/${params.user}/mnt/VI6Bioinformatics/Central_Pipelines/denovoAssembly/Viral_Reference_Databases/Influenza_Database/"
else
    ref = params.ref

log.info """\
    V I 6 - N F - P I P E L I N E
    =============================
    submission: ${params.submission}
    host      : ${params.host}
    hostdir   : ${host}
    virus     : ${params.ref}
    virusdir  : ${ref}
    reads     : ${params.reads}
    user      : ${params.user}
    """
    .stripIndent()
    
Channel
    .fromFilePairs( pattern )
        .set{ read_pairs_ch }

    
process HOSTMAP {

    tag "Mapping $sampleid to host"
     
    conda "/home/${params.user}/miniconda3/envs/denovoAssembly-v2"
    
    input:
    tuple val(sampleid), path(reads)
    
    output:
    tuple val(sampleid), path("${sampleid}_NotHostReads.bam")
    
    script:
    """
    bwa mem -t 4 ${host} ${reads[0]} ${reads[1]} | samtools view -@ 4 -b -f 4 -o "$sampleid"_NotHostReads.bam -;
    """
}
    
process SUBSAMPLE {

    tag "Subsampling ${sampleid} non-host reads"
    
    conda "/home/${params.user}/miniconda3/envs/denovoAssembly-v2"
    
    input:
    tuple val(sampleid), path(nonhost)
    
    output:
    tuple val(sampleid), path("*NotHostReads*")    
    
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
    
    samtools sort -@ 4 -n -O BAM $BAMFile | samtools fastq -@ 4 -1 !{sampleid}_NotHostReads1.fastq -2 !{sampleid}_NotHostReads2.fastq -s !{sampleid}_NotHostSingletons.fastq -;
    rm "$(readlink -f $BAMFile)"
    '''
}

process SPADES {

    publishDir "/home/${params.user}/mnt/VI6Storage/Raw_NGS_Data/${params.submission}/${sampleid}", mode: 'copy'
    
    conda "/home/${params.user}/miniconda3/envs/denovoAssembly-v2"
    
    input:
    tuple val(sampleid), path(reads)
    
    output:
    tuple val(sampleid), path("SPAdes/*contigs.fasta")
    
    script:
    """
    spades.py --threads 4 --only-assembler --pe1-1 ${reads[0]} --pe1-2 ${reads[1]} -o SPAdes
    
    mv SPAdes/contigs.fasta SPAdes/${sampleid}_SPAdes_contigs.fasta
    """

}

process BLAST {
    
    conda "/home/${params.user}/miniconda3/envs/denovoAssembly-v2"
    
    publishDir "/home/${params.user}/mnt/VI6Storage/Raw_NGS_Data/${params.submission}/${sampleid}/BLAST_Hits", mode: 'copy'
    
    input:
    tuple val(sampleid), path(contigs)
    
    output:
    tuple val(sampleid), path("*_crunch.txt")
    
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
    
    conda "/home/${params.user}/miniconda3/envs/denovoAssembly-v2"
    
    input:
    tuple val(sampleid), path(crunches) 
    
    output:
    tuple val(sampleid), path("*_match.fas")  
    
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

    conda "/home/${params.user}/miniconda3/envs/denovoAssembly-v2"
    
    publishDir "/home/${params.user}/mnt/VI6Storage/Raw_NGS_Data/${params.submission}/${sampleid}", mode: 'copy'
    
    errorStrategy 'finish'
    
    input:
    tuple val(sampleid), path(matches)
    
    output:
    tuple val(sampleid), path("*.fa")
    
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

    conda "/home/${params.user}/miniconda3/envs/denovoAssembly-v2"
    
    publishDir "/home/${params.user}/mnt/VI6Storage/Raw_NGS_Data/${params.submission}/${sampleid}", mode: 'copy'
    
    input:
    tuple val(sampleid), path(realign)
    tuple val(sampleid), path(rfile)
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
    threads=$(grep -c ^processor /proc/cpuinfo)
    
    samtools sort -n -@ "$threads" !{realign} -o - | \
    samtools fixmate -r -m 	-@ "$threads" - - | \
    samtools sort -@ "$threads" -o - | \
    samtools markdup -r -s -O BAM -@ "$threads" - "$samplename"_iter!{count}_clean_mapOnly.bam ;
		
    samtools index -@ "$threads" "$samplename"_iter!{count}_clean_mapOnly.bam

    mkdir genconsensus_Results
    
    cd genconsensus_Results
    
    python /home/!{params.user}/mnt/VI6Bioinformatics/Central_Pipelines/Utils/genconsensus/genconsensus.py  -t "0" -m "1" -n "-" -r ../"$rfile" -b ../"$samplename"_iter!{count}_clean_mapOnly.bam
    mv final_consensus.fasta ../"$samplename"_iter!{count}_consensus.fasta
    
    cd ..		
    
    rm -r genconsensus_Results
    
    samtools flagstat "$samplename"_iter!{count}_realign.bam > "$samplename"_iter!{count}_MappingStats.txt
    
    rm "$(readlink -f !{realign})"
    '''
    
}

workflow{
    read_pairs_ch.view()
    bam = HOSTMAP(read_pairs_ch)
    bam.view{ it }
    ext_reads = SUBSAMPLE(bam)
    ext_reads.view{ it }
    dn_contigs = SPADES(ext_reads)
    matches = BLAST(dn_contigs)
    top_matches = TOPMATCH(matches)
    dn_ref = CAT(top_matches)
    map1 = MAP1(dn_ref, read_pairs_ch, 1)
    con1 = VCF1(map1, dn_ref, 1)
    map2 = MAP2(con1, read_pairs_ch, 2)
    con2 = VCF2(map2, con1, 2)
    map3 = MAP3(con2, read_pairs_ch, 3)
    con3 = VCF3(map3, con2, 3)
    map_final = MAP4(con3, read_pairs_ch, 4)
    con_final = CONSENSUS(map_final, con3, 4)            
}

   
