nextflow.enable.dsl=2

// STAAR docker image 
//                  docker pull zilinli/staarpipeline:0.9.6

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/../assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

    process vcf2gds {                     //vcf2gds.R            R -> seqArray
        input:
        file vcf_chr
        output:
        file gds_out
        script:
        """
        seqVCF2GDS(vcf_chr, gds_out)
        """
    }
    process annotationVariantList {     //Varinfo_gds.R

    }
    process annotationFAVOR {           //Annotate.R

    }
    process gds2agds {                  //gds2agds.R

    } 


/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/



/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow {

    // first create a channel for each of the chromosomes.
    input_channel = Channel.fromPath(params.input_tsv, followLinks: true, checkIfExists: true)
    channel_input_data_table
        .splitCsv(header: true, sep: params.input_tables_column_delimiter).map{row->tuple(row.vcf,row.chr )}
	    .set{vcf_chr}

    vcf_chr.view()

    //STEP 0: VCF to GDS (need only a R function seqVCF2GDS(vcf, gds))
    vcf2gds(vcf_chr)

    //STEP 1: Generate the variants list to be annotated.
        //Script: Varinfo_gds.R
    annotationVariantList()

    //STEP 2: Annotate the variants using the FAVOR database through xsv software.
        //Script: Annotate.R
    annotationFAVOR()

    //STEP 3: Generate the annotated GDS file.
        //Script: gds2agds.R
    gds2agds()


    // extract the tar chromosome notation files
    // tar -xzvf /lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/data/FAVOR/chr20.tar.gz -C ./

    // Run the 3 codes to create aGDS files - these are the outputs.
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

/*
workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}
*/

/*
========================================================================================
    THE END
========================================================================================
*/
