
// STAAR docker image 
//                  docker pull zilinli/staarpipeline:0.9.6
//  https://github.com/xihaoli/STAARpipeline-Tutorial#association-analysis-using-staarpipeline
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
    // Step 0: Preparation for association analysis of whole-genome/whole-exome sequencing studies
        // Input: aGDS files of all 22 chromosomes. For more details, please see the R script.
        // Output: agds_dir.Rdata, Annotation_name_catalog.Rdata, jobs_num.Rdata.
        // Script: Association_Analysis_PreStep.r

    process analysisPreStep {     
        input:
        file aGDS
        output:
        file "*_dir.Rdata", emit: agds_dir
        file "*_catalog.Rdata", emit: annotation_name
        file "*_num.Rdata", emit: jobs_num
        script:
        """
#!/usr/bin/env Rscript

## load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)


        """
    }

    // Step 1: Fit STAAR null model
        // Input: Phenotype data and (sparse) genetic relatedness matrix. For more details, please see the R scripts.
        // Output: a Rdata file of the STAAR null model.
        // Script: STAARpipeline_Null_Model.r or STAARpipeline_Null_Model_GENESIS.r

    process fitNullModel {     
        input:
        file phenotypeCsv
        file sGRM
        output:
        file "*_nullmodel.Rdata", emit: objNullModel
        script:
        """
#!/usr/bin/env Rscript

library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
        """
    }

    // Step 2: Individual analysis
        // Input: aGDS files and the STAAR null model. For more details, please see the R script.
        // Output: Rdata files with the user-defined names.
        // Script: STAARpipeline_Individual_Analysis.r
/*
    process individualAnalysis {     
        input:

        output:

        script:
        """

        """
    }
*/
    // Step 3.1: Gene-centric coding analysis
        // Input: aGDS files and the STAAR null model. For more details, please see the R scripts.
        // Output: 381 Rdata files with the user-defined names. For more details, please see the R scripts.
        // Script: STAARpipeline_Gene_Centric_Coding.r and STAARpipeline_Gene_Centric_Coding_Long_Masks.r
/*
    process geneCentricCoding {     
        input:

        output:

        script:
        """

        """
    }
*/
    // Step 3.2: Gene-centric noncoding analysis
        // Input: aGDS files and the STAAR null model. For more details, please see the R scripts.
        // Output: 387 Rdata files with the user-defined names for protein-coding genes and 223 Rdata files with the user-defined names for ncRNA genes. For more details, please see the R scripts.
        // Script: STAARpipeline_Gene_Centric_Noncoding.r, STAARpipeline_Gene_Centric_Noncoding_Long_Masks.r, 
        //          STAARpipeline_Gene_Centric_ncRNA.r and STAARpipeline_Gene_Centric_ncRNA_Long_Masks.r
/*
    process geneCentricNoCoding {     
        input:

        output:

        script:
        """

        """
    }
*/
    // Step 4: Sliding window analysis
        // Input: aGDS files and the STAAR null model. For more details, please see the R script.
        // Output: Rdata files with the user-defined names.
        // Script: STAARpipeline_Sliding_Window.r
/*
    process slidingWindow {     
        input:

        output:

        script:
        """

        """
    }
*/
    // Step 5.0: Obtain SCANG-STAAR null model
        // Input: STAAR null model. For more details, please see the R script.
        // Output: a Rdata file of the SCANG-STAAR null model.
        // Script: STAARpipeline_STAAR2SCANG.r
/*
    process staar2scang {     
        input:

        output:

        script:
        """

        """
    }
*/
    // Step 5: Dynamic window analysis using SCANG-STAAR
        // Input: SCANG-STAAR null model. For more details, please see the R script.
        // Output: Rdata files with the user-defined names.
        // Script: STAARpipeline_Dynamic_Window.r
/*
    process dynamicWindowSCANG {     
        input:

        output:

        script:
        """

        """
    }
*/