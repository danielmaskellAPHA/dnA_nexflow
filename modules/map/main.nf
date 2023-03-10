user = "$USER"

process MAP {

    conda "/home/${user}/miniconda3/envs/denovoAssembly-v2"
    
    maxForks 5
    
    input:
    tuple val(sampleid), path(rfile), path(reads)
    val count
    
    output:
    tuple val(sampleid), path("*_realign.bam"), path(rfile), path(reads)

    shell:
    def avail_mem = task.memory ? ((task.memory.toBytes() - 6000000000) / task.cpus) : false
    def sort_mem = avail_mem && avail_mem > 2000000000 ? "-m $avail_mem" : ''
    '''
    rfile=!{rfile}

	ref=$(basename "$rfile")
	refname=${ref%%_*}
	reffile=${ref%%.*}
	
    samplename=!{sampleid}
    
    bwa index !{rfile}
    
    samtools faidx !{rfile}
    
    picard -XX:ParallelGCThreads=!{task.cpus} CreateSequenceDictionary R="$rfile" O=${rfile%%.*}.dict

    bwa mem -T10 -t !{task.cpus} -k 19 -B 4 -O 6 -R "$(echo "@RG\\tID:!{sampleid}\\tSM:!{sampleid}\\tLB:!{sampleid}")" "$rfile" !{reads[0]} !{reads[1]} | samtools view -@ !{task.cpus} -u - | samtools sort -@ !{task.cpus} !{sort_mem} -o "$samplename"_iter!{count}_map_sorted.bam
    
    samtools index -@ !{task.cpus} "$samplename"_iter!{count}_map_sorted.bam

    python /home/!{user}/mnt/VI6Bioinformatics/Central_Pipelines/denovoAssembly/gatk.py -T RealignerTargetCreator -nt !{task.cpus} -R "$rfile" -I "$samplename"_iter!{count}_map_sorted.bam -o indel!{count}.list
    
    python /home/!{user}/mnt/VI6Bioinformatics/Central_Pipelines/denovoAssembly/gatk.py -T IndelRealigner -R "$rfile" -I "$samplename"_iter!{count}_map_sorted.bam -targetIntervals indel!{count}.list -maxReads 50000 -o "$samplename"_iter!{count}_realign.bam
		
    rm "$samplename"_iter!{count}_map_sorted.bam
    rm "$samplename"_iter!{count}_map_sorted.bam.bai
    '''
    
}

process VCF {
    conda "/home/${user}/miniconda3/envs/denovoAssembly-v2"
    
    input:
    tuple val(sampleid), path(realign), path(rfile), path(reads)
    val count
    
    output:
    tuple val(sampleid), path("*.fasta"), path(reads)
    
    shell:
    '''
    rfile=!{rfile}

	ref=$(basename "$rfile")
	refname=${ref%%_*}
	reffile=${ref%%.*}
	
    samplename=!{sampleid}
    
    bcftools mpileup -L 10000 -Q1 -AEpf "$rfile" !{realign} | bcftools call --threads !{task.cpus} -c - > "$samplename"_iter!{count}.vcf
    perl /home/!{user}/mnt/VI6Bioinformatics/Central_Pipelines/Utils/vcf2consensus.pl  consensus -f "$rfile" "$samplename"_iter!{count}.vcf | sed '/^>/ s/-iter[0-9]//;/^>/ s/$/'-iter!{count}'/' - > "$samplename"_iter!{count}_consensus.fasta

    samtools flagstat !{realign} > "$samplename"_iter!{count}_MappingStats.txt
    
    rm "$(readlink -f !{realign})"
   
    rm "$(readlink -f !{rfile})"
    '''
}
