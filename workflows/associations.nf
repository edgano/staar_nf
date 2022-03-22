// STAAR docker image 
//                  docker pull zilinli/staarpipeline:0.9.6
//                  #docker run -v /media/edgano/24F8A36CF8A33AC6/projects/staar_nf-main/assets:/source -ti zilinli/staarpipeline:0.9.6
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

//ch_multiqc_config        = file("$projectDir/../assets/multiqc_config.yaml", checkIfExists: true)
//ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/
    process vcf2gds {                   //vcf2gds.R
        input:
        file vcf_chr
        output:
        path '*.gds', emit: gds_out
        script:
        """
#!/usr/bin/env Rscript

library(SeqArray)
library(Rcpp)

seqVCF2GDS("${vcf_chr}", "${vcf_chr}_out.gds")
        """
    }
    // Input: GDS files of each chromosome and the FAVOR database information FAVORdatabase_chrsplit.csv. 
    //        For more details, please see the R script.
    // Output: CSV files of the variants list. For each chromosome, the number of CSV files is listed in FAVORdatabase_chrsplit.csv.
/*
    process annotationVariantList {     //Varinfo_gds.R
        input:
        file ddbb_split
        file gds_in                     //this is the output of vcf2gds
        output:
        file "*.csv", emit: variantListCSV_out
        script:
        """

        """
    }
    // Input: CSV files of the variants list to be annotated, the FAVOR database information FAVORdatabase_chrsplit.csv,
    //        the FAVOR database, and the directory xsv software. For more details, please see the R script.
    // Output: CSV files of the annotated variants list.
    process annotationFAVOR {           //Annotate.R
        input:
        file ddbb_split
        path ddbb_path
        output:
        file "*.csv", emit: annotatedVariantList_out
        script:
        """

        """
    }
    // Input: GDS files and the csv files of annotated variants list (Anno_chrXX.csv or Anno_chrXX_STAARpipeline.csv). 
    //        For more details, please see the R script.
    // Output: aGDS files including both the genotype and annotation information
    process gds2agds {                  //gds2agds.R
        input:
        gds_in
        annootatedVL_in
        output:
        file "*.agds" , emit: aGDS_out
        script:
        """

        """
    } 
*/
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
// Info -> https://github.com/xihaoli/STAARpipeline-Tutorial

// Info required for completion email and summary
def multiqc_report = []

// TODO -> this var move at the "validate inputs"
inputVCF_ch = Channel.fromPath("$projectDir/assets/Interval_WGS_chr20_TF_binding_site_test.vcf", checkIfExists: true)
favorDDBB_split_ch = Channel.fromPath('$projectDir/assets/FAVORdatabase_chrsplit.csv.vcf')
ddbb_path = "/lustre/scratch123/hgi/teams/hgi/mo11/associations/FAVOR/n/holyscratch01/xlin/xihao_zilin/FAVORAnnotatorDB"
//

workflow ASSOCIATIONS {

    // first create a channel for each of the chromosomes.
    //input_channel = Channel.fromPath(params.input_tsv, followLinks: true, checkIfExists: true)
    /*channel_input_data_table
        .splitCsv(header: true, sep: params.inputVCF).map{row->tuple(row.vcf,row.chr )}
	    .set{vcf_chr}
        vcf_chr.view()

    */

    //STEP 0: VCF to GDS (need only a R function seqVCF2GDS(vcf, gds))
    vcf2gds(inputVCF_ch)

    //STEP 1: Generate the variants list to be annotated.
    //annotationVariantList(favorDDBB_split_ch, vcf2gds.out.gds_out)

    //STEP 2: Annotate the variants using the FAVOR database through xsv software.
    //annotationFAVOR(favorDDBB_split_ch,ddbb_path)

    //STEP 3: Generate the annotated GDS file.
    //gds2agds(gds_out,annotatedVariantList_out)

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
