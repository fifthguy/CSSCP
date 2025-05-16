

// Include processes from modules
include { runMafftAdjust_1; runMarsRotate_1; runMafftAdjust_2; runMarsRotate_2; runMarsRotate_3; runMafftAdjust_3; runMafftAlignment } from './modules/alignments.nf'
include { prepareMultifasta; splitSequences; computeAndMergeSimilarity; cleanupAlignments; cleanupIntermediaries} from './modules/utilities.nf' 


workflow dataprep {

    take:
    fastaFile

    main:   
    preppedFasta = prepareMultifasta(fastaFile)
    splitFiles = splitSequences(preppedFasta).flatten()

    emit:
    splitFiles
}

workflow adjustMafft {

    take:
    files

    main:
    adjusted = runMafftAdjust_1(files)

    emit:
    adjusted

}

workflow adjustAndRotate {

    take:
    files

    main:
    adjusted = runMafftAdjust_1(files)
    rotated = runMarsRotate_1(adjusted)

    emit:
    rotated

}

workflow adjustAndRotate_round2 {

    take:
    files

    main:
    adjusted = runMafftAdjust_2(files)
    rotated = runMarsRotate_2(adjusted)

    emit:
    rotated

}


workflow adjustAndRotate_round3 {

    take:
    files

    main:
    adjusted = runMafftAdjust_3(files)
    rotated = runMarsRotate_3(adjusted)

    emit:
    rotated

}


workflow diagnostics {

    take:
    input

    main:
    aligned = runMafftAlignment(input)
    similarity_tsvs = computeAndMergeSimilarity(aligned.collect())

    cleanupIntermediaries("dummy")

    if (!params.keep) {
        cleanupAlignments("dummy")
    }
}


workflow simple {

    take:
    data
    
    emit:
    data
}

workflow adjust {

    take:
    data
    
    main:
    adjusted = adjustMafft(data)

    emit:
    adjusted

}



workflow adjust_rotate {

    take:
    data
    
    main:
    adjusted = adjustAndRotate(data)

    emit:
    adjusted
}


workflow adjust_rotate_2x {

    take:
    data
    
    main:
    adjusted = adjustAndRotate(data)
    adjusted2 = adjustAndRotate_round2(adjusted)

    emit:
    adjusted2
}


workflow adjust_rotate_3x {

    take:
    data
    
    main:
    adjusted = adjustAndRotate(data)
    adjusted2 = adjustAndRotate_round2(adjusted)
    adjusted3 = adjustAndRotate_round3(adjusted2)

    emit:
    adjusted3
}