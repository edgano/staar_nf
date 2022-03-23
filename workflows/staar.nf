
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
        ###############################
        #           Input
        ###############################
        ## Phenotype file
        # phenotype <- read.csv("/path_to_the_file/pheno.csv")   <- nf input phenotypeCsv   
        ## (sparse) GRM file
        # sgrm <- get(load("/path_to_the_file/sGRM.Rdata"))     <- nf input sGRM 

        ## output file name
        output_name <- "obj_nullmodel.Rdata"    ## not sure if its needed

        ###############################
        #        Main Function
        ###############################
        ### fit null model
        obj_nullmodel <- fit_nullmodel(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+as.factor(study_race),data=$phenotypeCsv,
                                    kins=$sGRM,use_sparse=TRUE,kins_cutoff=0.022,id="sample.id",family=gaussian(link="identity"),verbose=TRUE)

        # TODO -> check if safe obj is correct
        save(obj_nullmodel,file=paste0(".",output_name))

        """
    }
    process fitNullModelGenesis {     
        input:
        file phenotypeCsv
        file sGRM
        output:
        file "*_GENESIS.Rdata", emit: objNullModel
        script:
        """
        #!/usr/bin/env Rscript

        library(GENESIS)
        library(STAAR)
        library(STAARpipeline)
        ###############################
        #           Input
        ###############################
        ## Phenotype file
        # phenotype <- read.csv("/path_to_the_file/pheno.csv")   <- nf input phenotypeCsv   
        ## (sparse) GRM file
        # sgrm <- get(load("/path_to_the_file/sGRM.Rdata"))     <- nf input sGRM 

        ## output file name
        output_name <- "obj_nullmodel_GENESIS.Rdata"    ## not sure if its needed

        ###############################
        #        Main Function
        ###############################
        ### fit null model using GENESIS 
        data_GENESIS <- as($phenotypeCsv,"AnnotatedDataFrame")          # Make AnnotatedDataFrame (specifically required by GENESIS)
        obj_nullmodel_GENESIS <- fitNullModel(data_GENESIS,outcome="LDLadj.norm",
                                            covars=c("age","age2","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","study_race"),
                                            cov.mat=$sGRM,group.var="study_race",AIREML.tol=1e-4,verbose=TRUE)

        ### convert GENESIS null model to STAAR null model
        obj_nullmodel <- genesis2staar_nullmodel(obj_nullmodel_GENESIS)
        # TODO -> check if safe obj is correct
        save(obj_nullmodel,file=paste0(".",output_name))

        """
    }
    // Step 2: Individual analysis
        // Input: aGDS files and the STAAR null model. For more details, please see the R script.
        // Output: Rdata files with the user-defined names.
        // Script: STAARpipeline_Individual_Analysis.r

        //TODO -> needed the jobs_num.Rdata TOO
        //        arrayid <- as.numeric(commandArgs(TRUE)[1])        

    process individualAnalysis {     
        input:
        file aGDS           // file or path??
        file nullModel
        output:

        script:
        """
        #!/usr/bin/env Rscript

        ## load required packages
        library(gdsfmt)
        library(SeqArray)
        library(SeqVarTools)
        library(STAAR)
        library(STAARpipeline)
        ###############################
        #           Input
        ###############################
        ## Number of jobs for each chromosome
        jobs_num <- get(load("/path_to_the_file/jobs_num.Rdata"))
        ## aGDS directory
        agds_dir <- get(load(\"$aGDS\"))
        ## Null model
        obj_nullmodel <- get(load(\"$nullModel\"))

        ## QC_label
        QC_label <- "annotation/filter"
        ## variant_type
        variant_type <- "variant"
        ## geno_missing_imputation
        geno_missing_imputation <- "mean"

        ## output path
            ## output_path <- "/path_to_the_output_file/"
        ## output file name
        output_file_name <- "TOPMed_F5_LDL_results_individual_analysis"
        ## input array id from batch file (Harvard FAS RC cluster)
        arrayid <- as.numeric(commandArgs(TRUE)[1])

        ###############################
        #        Main Function
        ###############################
        chr <- which.max(arrayid <= cumsum(jobs_num\$individual_analysis_num))
        group.num <- jobs_num\$individual_analysis_num[chr]

        if (chr == 1){
        groupid <- arrayid
        }else{
        groupid <- arrayid - cumsum(jobs_num\$individual_analysis_num)[chr-1]
        }

        ### gds file
        gds.path <- agds_dir[chr]
        genofile <- seqOpen(gds.path)

        start_loc <- (groupid-1)*10e6 + jobs_num\$start_loc[chr]
        end_loc <- start_loc + 10e6 - 1
        end_loc <- min(end_loc,jobs_num\$end_loc[chr])

        a <- Sys.time()
        results_individual_analysis <- c()
        if(start_loc < end_loc)
        {
            results_individual_analysis <- Individual_Analysis(chr=chr,start_loc=start_loc,end_loc=end_loc,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                            QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation)
        }
        b <- Sys.time()
        b - a

        save(results_individual_analysis,file=paste0(".",output_file_name,"_",arrayid,".Rdata"))

        seqClose(genofile)
        """
    }

    // Step 3.1: Gene-centric coding analysis
        // Input: aGDS files and the STAAR null model. For more details, please see the R scripts.
        // Output: 381 Rdata files with the user-defined names. For more details, please see the R scripts.
        // Script: STAARpipeline_Gene_Centric_Coding.r and STAARpipeline_Gene_Centric_Coding_Long_Masks.r

        // TODO -> gene_num_in_array <- 50  
        //         table(genes_info[,2])
        //  ### exclude large genes     <-- magic numbers
    process geneCentricCoding {     
        input:
        file aGDS
        file nullModel
        output:

        script:
        """
         #!/usr/bin/env Rscript

        ###############################
        #           Input
        ###############################
        
        ###############################
        #        Main Function
        ###############################

        """
    }
/*
    process geneCentricCodingLongMask {     
        input:
        file aGDS
        file nullModel
        output:

        script:
        """
        #!/usr/bin/env Rscript

        ###############################
        #           Input
        ###############################
        
        ###############################
        #        Main Function
        ###############################
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
        file aGDS
        file nullModel
        output:

        script:
        """
        #!/usr/bin/env Rscript

        ###############################
        #           Input
        ###############################
        
        ###############################
        #        Main Function
        ###############################
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
        file aGDS
        file nullModel
        output:

        script:
        """
        #!/usr/bin/env Rscript

        ###############################
        #           Input
        ###############################
        
        ###############################
        #        Main Function
        ###############################
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
        file staarNullModel
        output:
        file "*.Rdata", emit: scangStaarNullModel
        script:
        """
        #!/usr/bin/env Rscript

        ###############################
        #           Input
        ###############################
        
        ###############################
        #        Main Function
        ###############################
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
        file staarNullModel
        output:
        
        script:
        """
        #!/usr/bin/env Rscript

        ###############################
        #           Input
        ###############################
        
        ###############################
        #        Main Function
        ###############################
        """
    }
*/
/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
// Info -> https://github.com/xihaoli/STAARpipeline-Tutorial

// Info required for completion email and summary
//def multiqc_report = []


workflow STAAR {

    // first create a channel

    //Step 0: Preparation for association analysis of whole-genome/whole-exome sequencing studies
    analysisPreStep(inputAgds)

    //Step 1: Fit STAAR null model
    fitNullModel(phenoCSV, sGRM)

    //Step 2: Individual analysis
    individualAnalysis(aGDS, fitNullModel.out.objNullModel)

    //Step 3.1: Gene-centric coding analysis
    geneCentricCoding(aGDS, fitNullModel.out.objNullModel)

    //Step 3.2: Gene-centric noncoding analysis
    geneCentricNoCoding(aGDS, fitNullModel.out.objNullModel)

    //Step 4: Sliding window analysis
    slidingWindow(aGDS, fitNullModel.out.objNullModel)

    // Step 5.0: Obtain SCANG-STAAR null model
    staar2scang(fitNullModel.out.objNullModel)

    //Step 5: Dynamic window analysis using SCANG-STAAR
    dynamicWindowSCANG(fitNullModel.out.objNullModel)

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
