process MAP {

    conda "/home/${params.user}/miniconda3/envs/denovoAssembly-v2"
    
    input:
    tuple val(sampleid), path(top_match)
    tuple val(sampleid_og), path(reads)
    val count
    
    output:
    tuple val(sampleid), path("*_realign.bam")

    shell:
    '''
    rfile=!{top_match}

	ref=$(basename "$rfile")
	refname=${ref%%_*}
	reffile=${ref%%.*}
	
    samplename=!{sampleid}
    threads=$(grep -c ^processor /proc/cpuinfo)
    
    
    bwa index !{top_match}
    
    samtools faidx !{top_match}
    
    picard -XX:ParallelGCThreads="$threads" CreateSequenceDictionary R="$rfile" O=${rfile%%.*}.dict

    bwa mem -T10 -t "$threads" -k 19 -B 4 -O 6 -R "$(echo "@RG\\tID:!{sampleid}\\tSM:!{sampleid}\\tLB:!{sampleid}")" "$rfile" !{reads[0]} !{reads[1]} | samtools view -@ "$threads" -u - | samtools sort -@ "$threads" -o "$samplename"_iter!{count}_map_sorted.bam
    
    samtools index -@ "$threads" "$samplename"_iter!{count}_map_sorted.bam

    python '/home/!{params.user}/mnt/VI6Bioinformatics/Central_Pipelines/denovoAssembly/gatk.py' -T RealignerTargetCreator -nt "$threads" -R "$rfile" -I "$samplename"_iter!{count}_map_sorted.bam -o indel!{count}.list
    
    python '/home/!{params.user}/mnt/VI6Bioinformatics/Central_Pipelines/denovoAssembly/gatk.py' -T IndelRealigner -R "$rfile" -I "$samplename"_iter!{count}_map_sorted.bam -targetIntervals indel!{count}.list -maxReads 50000 -o "$samplename"_iter!{count}_realign.bam
		
    rm "$samplename"_iter!{count}_map_sorted.bam
    rm "$samplename"_iter!{count}_map_sorted.bam.bai
    
    if [ !{count} == 4 ]
    then
    	rm "$(readlink -f !{reads[0]})"
    	rm "$(readlink -f !{reads[1]})"
   	fi
    
    '''
    
}

process VCF {
    conda "/home/${params.user}/miniconda3/envs/denovoAssembly-v2"
    
    input:
    tuple val(sampleid), path(realign)
    tuple val(sampleid), path(rfile)
    val count
    
    output:
    tuple val(sampleid), path("*.fasta")
    
    shell:
    '''
    rfile=!{rfile}

	ref=$(basename "$rfile")
	refname=${ref%%_*}
	reffile=${ref%%.*}
	
    samplename=!{sampleid}
    threads=$(grep -c ^processor /proc/cpuinfo)
    
    bcftools mpileup -L 10000 -Q1 -AEpf "$rfile" !{realign} | bcftools call --threads "$threads" -c - > "$samplename"_iter!{count}.vcf
    perl /home/!{params.user}/mnt/VI6Bioinformatics/Central_Pipelines/Utils/vcf2consensus.pl  consensus -f "$rfile" "$samplename"_iter!{count}.vcf | sed '/^>/ s/-iter[0-9]//;/^>/ s/$/'-iter!{count}'/' - > "$samplename"_iter!{count}_consensus.fasta

    samtools flagstat !{realign} > "$samplename"_iter!{count}_MappingStats.txt
    
    rm "$(readlink -f !{realign})"
   
    rm "$(readlink -f !{rfile})"
    '''
}
