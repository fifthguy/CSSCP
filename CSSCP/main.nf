#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.in = "$baseDir/data/sample.fa"
params.outdir = "results"

params.workflow = "simple"

params.keep = true


// Include workflows
include { dataprep; diagnostics; simple; adjust; adjust_rotate; adjust_rotate_2x; adjust_rotate_3x } from './workflows.nf'

/*
 * Main workflow logic
 */
workflow {
    filein = Channel.fromPath(params.in)
    data = dataprep(filein)
    
    // Conditional execution based on the workflow parameter
    if (params.workflow == 'simple') {
        out = simple(data)
    } else if (params.workflow == 'adjust') {
        out = adjust(data)
    } else if (params.workflow == 'adjust_rotate') {
        out = adjust_rotate(data)
    } else if (params.workflow == 'adjust_rotate_2x') {
        out = adjust_rotate_2x(data)    
    } else if (params.workflow == 'adjust_rotate_3x') {
        out = adjust_rotate_3x(data)
    } else {
        error "Unknown workflow: ${params.workflow}"
    }

    diagnostics(out)
}
