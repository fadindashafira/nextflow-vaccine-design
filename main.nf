#!/usr/bin/env nextflow

/*
========================================================================================
    In Silico Vaccine Design Pipeline
========================================================================================
    Github : https://github.com/username/in-silico-vaccine-design
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

// Load modules
include { retrieveSequence } from './modules/retrieveSequence'
include { predictBCellEpitopes } from './modules/predictEpitopes'
include { predictTCellEpitopesI } from './modules/predictEpitopes'
include { predictTCellEpitopesII } from './modules/predictEpitopes'
include { combineEpitopes } from './modules/combineEpitopes'
include { designVaccineConstruct } from './modules/designVaccineConstruct'
include { evaluateVaccineConstruct } from './modules/evaluateVaccine'
include { molecularDynamics } from './modules/molecularDynamics'

// Load params from config
params.accession = "BAL61222.1"
params.outdir = "results"
params.run_md = false
params.help = false

// Show help message
if (params.help) {
    log.info"""
    =============================================
    IN SILICO VACCINE DESIGN PIPELINE
    =============================================
    Usage:
    nextflow run main.nf --accession <accession> --outdir <output_directory>
    
    Parameters:
      --accession     Protein accession number from NCBI (default: ${params.accession})
      --outdir        Output directory (default: ${params.outdir})
      --run_md        Run molecular dynamics simulation (default: ${params.run_md})
      --help          Show this help message
    """.stripIndent()
    exit 0
}

log.info"""
=============================================
 IN SILICO VACCINE DESIGN PIPELINE
=============================================
 Accession: ${params.accession}
 Output directory: ${params.outdir}
 Run molecular dynamics: ${params.run_md}
"""

// Main workflow
workflow {
    // Step 1: Retrieve sequence
    retrieveSequence(params.accession)
    
    // Step 2: Predict epitopes
    predictBCellEpitopes(retrieveSequence.out.fasta)
    predictTCellEpitopesI(retrieveSequence.out.fasta)
    predictTCellEpitopesII(retrieveSequence.out.fasta)
    
    // Step 3: Combine epitopes and design vaccine
    combineEpitopes(
        predictBCellEpitopes.out.bcell_epitopes,
        predictTCellEpitopesI.out.tcell_i_epitopes,
        predictTCellEpitopesII.out.tcell_ii_epitopes
    )
    designVaccineConstruct(combineEpitopes.out.combined_epitopes)
    
    // Step 4: Evaluate vaccine
    evaluateVaccineConstruct(designVaccineConstruct.out.vaccine_construct)
    
    // Step 5: Molecular dynamics (optional)
    if (params.run_md) {
        molecularDynamics(designVaccineConstruct.out.vaccine)
    }
}

// On completion
workflow.onComplete {
    log.info"""
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    """
}