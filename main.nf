#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/associations
========================================================================================
    Github : https://github.com/nf-core/associations
    Website: https://nf-co.re/associations
    Slack  : https://nfcore.slack.com/channels/associations
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

//params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

//WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { ASSOCIATIONS } from './workflows/associations'

//
// WORKFLOW: Run main nf-core/associations analysis pipeline
//
workflow NFCORE_ASSOCIATIONS {
    ASSOCIATIONS ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_ASSOCIATIONS ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
