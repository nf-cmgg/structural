process BUILD_PLOIDY_MATRIX {
    tag "build_ploidy_matrix"
    label 'process_low'

    conda 'bioconda::tabix=1.11'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'quay.io/biocontainers/tabix:1.11--hdfd78af_0' }"

    input:
    path bincov_matrix

    output:
    path "*_ploidy_matrix.bed.gz" , emit: ploidy_matrix
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "batch"

    def binsize = task.ext.binsize ?: 1000000

    """
    zcat ${bincov_matrix} \\
        | awk ' \\
        function printRow() \\
            {printf "%s\t%d\t%d",chr,start,stop; \\
            for(i=4;i<=nf;++i) {printf "\t%d",vals[i]; vals[i]=0}; \\
            print ""} \\
        BEGIN {binSize=${binsize}} \\
        NR==1 {print substr(\$0,1)} \\
        NR==2 {chr=\$1; start=\$2; stop=start+binSize; nf=NF; for(i=4;i<=nf;++i) {vals[i]=\$i}} \\
        NR>2  {if(\$1!=chr){printRow(); chr=\$1; start=\$2; stop=start+binSize} \\
                else if(\$2>=stop) {printRow(); while(\$2>=stop) {start=stop; stop=start+binSize}} \\
                for(i=4;i<=nf;++i) {vals[i]+=\$i}} \\
        END   {if(nf!=0)printRow()}' \\
        | bgzip > ${prefix}_ploidy_matrix.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
