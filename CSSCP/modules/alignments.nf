process runMafftAdjust_1 {
    publishDir params.outdir, mode: 'copy'

    input:
    path pairfile

    output:
    path "adjusted_1/*.fa"

    script:
    """
    mkdir -p adjusted_1
    mafft --thread 1 --adjustdirectionaccurately --localpair --maxiterate 1000 ${pairfile} | \
        seqkit seq -g - | seqtk seq -l 10000000 - > adjusted_1/\$(basename ${pairfile})
    """
}

process runMafftAdjust_2 {
    publishDir params.outdir, mode: 'copy'

    input:
    path pairfile

    output:
    path "adjusted_2/*.fa"

    script:
    """
    mkdir -p adjusted_2
    mafft --thread 1 --adjustdirectionaccurately --localpair --maxiterate 1000 ${pairfile} | \
        seqkit seq -g - | seqtk seq -l 10000000 - > adjusted_2/\$(basename ${pairfile})
    """
}


process runMafftAdjust_3 {
    publishDir params.outdir, mode: 'copy'

    input:
    path pairfile

    output:
    path "adjusted_3/*.fa"

    script:
    """
    mkdir -p adjusted_3
    mafft --thread 1 --adjustdirectionaccurately --localpair --maxiterate 1000 ${pairfile} | \
        seqkit seq -g - | seqtk seq -l 10000000 - > adjusted_3/\$(basename ${pairfile})
    """
}
process runMarsRotate_1 {
    publishDir params.outdir, mode: 'copy'

    input:
    path pairfile

    output:
    path "rotated_1/*.fa"

    script:
    """
    mkdir -p rotated_1
    mars -i ${pairfile} -a DNA -m 0 -q 10 -T 1 -o rotated_1/\$(basename ${pairfile})
    """
}

process runMarsRotate_2 {
    publishDir params.outdir, mode: 'copy'

    input:
    path pairfile

    output:
    path "rotated_2/*.fa"

    script:
    """
    mkdir -p rotated_2
    mars -i ${pairfile} -a DNA -m 0 -q 10 -T 1 -o rotated_2/\$(basename ${pairfile})
    """
}

process runMarsRotate_3 {
    publishDir params.outdir, mode: 'copy'

    input:
    path pairfile

    output:
    path "rotated_3/*.fa"

    script:
    """
    mkdir -p rotated_3
    mars -i ${pairfile} -a DNA -m 0 -q 10 -T 1 -o rotated_3/\$(basename ${pairfile})
    """
}

process runMafftAlignment {
    publishDir params.outdir, mode: 'copy'

    input:
    path pairfile

    output:
    path "alignments/*.fa"


    script:
    """
    mkdir -p alignments
    mafft --thread 1 --localpair --maxiterate 1000 ${pairfile} > alignments/\$(basename ${pairfile})
    """
}