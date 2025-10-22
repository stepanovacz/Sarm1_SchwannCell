# Sarm1_SchwannCell
This code supports the analysis of single-nucleus RNA sequencing (snRNA-seq) data from wild-type (WT) and Sarm1 knockout (KO) sciatic nerves following transection injury. Samples include sham, 2 hours post-injury, and 1 day post-injury time points.

The code was developed for the manuscript currently under revision:
Sarm1 Gates the Transition from Protective to Repair Schwann Cell States Following Nerve Injury
**Stepanova, E.¹⁶, Hunter-Chang, S.¹², Lee, J.¹, Pavelec, C.⁷⁸, Tripathi, A.¹, Cho, C.¹, Vegiraju, T.¹, Kim-Aun, C.¹, Kucenas, S.¹²⁴⁶, Leitinger, N.⁷⁸, Coutinho-Budd, J.²⁴⁶, Campbell, J.¹²⁶, Deppmann, C.¹–⁶

Session info for snRNAseq
> sessionInfo()
R version 4.4.1 (2024-06-14)
Platform: x86_64-pc-linux-gnu
Running under: Rocky Linux 8.10 (Green Obsidian)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggrepel_0.9.6               org.Mm.eg.db_3.20.0         AnnotationDbi_1.68.0       
 [4] clusterProfiler_4.14.6      Libra_1.0.0                 edgeR_4.4.2                
 [7] limma_3.62.2                ggridges_0.5.6              SeuratWrappers_0.4.0       
[10] monocle3_1.4.26             SingleCellExperiment_1.28.1 SummarizedExperiment_1.36.0
[13] GenomicRanges_1.58.0        GenomeInfoDb_1.42.3         IRanges_2.40.1             
[16] S4Vectors_0.44.0            MatrixGenerics_1.18.1       matrixStats_1.5.0          
[19] Biobase_2.66.0              BiocGenerics_0.52.0         gridExtra_2.3              
[22] future_1.67.0               RColorBrewer_1.1-3          pheatmap_1.0.13            
[25] patchwork_1.3.2             BUSpaRse_1.20.0             BiocManager_1.30.26        
[28] lubridate_1.9.3             forcats_1.0.0               stringr_1.5.1              
[31] purrr_1.1.0                 readr_2.1.5                 tidyr_1.3.1                
[34] tibble_3.3.0                ggplot2_3.5.2               tidyverse_2.0.0            
[37] dplyr_1.1.4                 Seurat_5.3.0                SeuratObject_5.1.0         
[40] sp_2.2-0                   

loaded via a namespace (and not attached):
  [1] R.methodsS3_1.8.2        progress_1.2.3           goftest_1.2-3            Biostrings_2.74.1       
  [5] TH.data_1.1-4            vctrs_0.6.5              ggtangle_0.0.7           spatstat.random_3.4-1   
  [9] digest_0.6.37            png_0.1-8                proxy_0.4-27             plyranges_1.26.0        
 [13] deldir_2.0-4             parallelly_1.45.1        MASS_7.3-61              reshape2_1.4.4          
 [17] httpuv_1.6.16            qvalue_2.38.0            withr_3.0.2              ggrastr_1.0.2           
 [21] ggfun_0.2.0              survival_3.7-0           memoise_2.0.1            ggbeeswarm_0.7.2        
 [25] emmeans_1.11.2-8         gson_0.1.0               tidytree_0.4.6           zoo_1.8-14              
 [29] pbapply_1.7-4            R.oo_1.27.1              prettyunits_1.2.0        KEGGREST_1.46.0         
 [33] promises_1.3.3           httr_1.4.7               restfulr_0.0.16          globals_0.18.0          
 [37] fitdistrplus_1.2-4       rstudioapi_0.16.0        UCSC.utils_1.2.0         miniUI_0.1.1.1          
 [41] generics_0.1.4           DOSE_4.0.1               curl_7.0.0               zlibbioc_1.52.0         
 [45] polyclip_1.10-7          GenomeInfoDbData_1.2.13  SparseArray_1.6.2        xtable_1.8-4            
 [49] S4Arrays_1.6.0           BiocFileCache_2.14.0     hms_1.1.3                irlba_2.3.5.1           
 [53] colorspace_2.1-1         filelock_1.0.3           hdf5r_1.3.11             ROCR_1.0-11             
 [57] reticulate_1.43.0        spatstat.data_3.1-8      magrittr_2.0.3           lmtest_0.9-40           
 [61] later_1.4.3              viridis_0.6.5            ggtree_3.14.0            lattice_0.22-6          
 [65] spatstat.geom_3.5-0      future.apply_1.11.2      scattermore_1.2          XML_3.99-0.19           
 [69] cowplot_1.2.0            RcppAnnoy_0.0.22         pillar_1.11.0            nlme_3.1-166            
 [73] compiler_4.4.1           RSpectra_0.16-2          stringi_1.8.7            tensor_1.5.1            
 [77] minqa_1.2.8              GenomicAlignments_1.42.0 plyr_1.8.9               crayon_1.5.3            
 [81] abind_1.4-8              BiocIO_1.16.0            blme_1.0-6               gridGraphics_0.5-1      
 [85] locfit_1.5-9.12          bit_4.6.0                sandwich_3.1-1           fastmatch_1.1-6         
 [89] codetools_0.2-20         multcomp_1.4-28          plotly_4.11.0            mime_0.13               
 [93] splines_4.4.1            Rcpp_1.1.0               fastDummies_1.7.5        dbplyr_2.5.0            
 [97] blob_1.2.4               AnnotationFilter_1.30.0  lme4_1.1-37              fs_1.6.6                
[101] listenv_0.9.1            Rdpack_2.6.4             ggplotify_0.1.2          estimability_1.5.1      
[105] Matrix_1.7-0             statmod_1.5.0            tzdb_0.5.0               pkgconfig_2.0.3         
[109] tools_4.4.1              cachem_1.1.0             rbibutils_2.3            RSQLite_2.4.3           
[113] viridisLite_0.4.2        DBI_1.2.3                numDeriv_2016.8-1.1      fastmap_1.2.0           
[117] scales_1.4.0             grid_4.4.1               pbmcapply_1.5.1          ica_1.0-3               
[121] Rsamtools_2.22.0         dotCall64_1.2            RANN_2.6.2               farver_2.1.2            
[125] reformulas_0.4.1         mgcv_1.9-1               yaml_2.3.10              rtracklayer_1.66.0      
[129] cli_3.6.5                tester_0.2.0             lifecycle_1.0.4          uwot_0.2.3              
[133] glmmTMB_1.1.10           mvtnorm_1.3-3            BiocParallel_1.40.2      timechange_0.3.0        
[137] gtable_0.3.6             rjson_0.2.23             progressr_0.14.0         parallel_4.4.1          
[141] ape_5.8-1                jsonlite_2.0.0           RcppHNSW_0.6.0           bitops_1.0-9            
[145] bit64_4.6.0-1            assertthat_0.2.1         Rtsne_0.17               yulab.utils_0.2.1       
[149] spatstat.utils_3.1-5     GOSemSim_2.32.0          zeallot_0.2.0            spatstat.univar_3.1-4   
[153] R.utils_2.13.0           lazyeval_0.2.2           shiny_1.11.1             htmltools_0.5.8.1       
[157] enrichplot_1.26.6        GO.db_3.20.0             sctransform_0.4.2        rappdirs_0.3.3          
[161] ensembldb_2.30.0         glue_1.8.0               spam_2.11-1              httr2_1.0.5             
[165] XVector_0.46.0           RCurl_1.98-1.17          treeio_1.30.0            BSgenome_1.74.0         
[169] boot_1.3-31              igraph_2.1.4             TMB_1.9.16               R6_2.6.1                
[173] DESeq2_1.46.0            labeling_0.4.3           GenomicFeatures_1.58.0   cluster_2.1.8.1         
[177] aplot_0.2.8              nloptr_2.2.1             DelayedArray_0.32.0      tidyselect_1.2.1        
[181] vipor_0.4.7              ProtGenerics_1.38.0      xml2_1.3.6               rsvd_1.0.5              
[185] KernSmooth_2.23-24       data.table_1.17.8        htmlwidgets_1.6.4        fgsea_1.32.4            
[189] biomaRt_2.62.1           rlang_1.1.6              spatstat.sparse_3.1-0    spatstat.explore_3.5-2  
[193] lmerTest_3.1-3           remotes_2.5.0            beeswarm_0.4.0


Session info for bulk RNAseq and bulk vs pseudobulk comparisons
R version 4.4.2 (2024-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26100)

Matrix products: default


locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
[4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] gridExtra_2.3               dplyr_1.1.4                 ggplot2_3.5.2               reshape2_1.4.4             
 [5] viridis_0.6.5               viridisLite_0.4.2           pheatmap_1.0.13             org.Mm.eg.db_3.20.0        
 [9] AnnotationDbi_1.68.0        DESeq2_1.46.0               SummarizedExperiment_1.36.0 Biobase_2.66.0             
[13] MatrixGenerics_1.18.1       matrixStats_1.5.0           GenomicRanges_1.58.0        GenomeInfoDb_1.42.3        
[17] IRanges_2.40.1              S4Vectors_0.44.0            BiocGenerics_0.52.0         tximport_1.34.0            

loaded via a namespace (and not attached):
 [1] KEGGREST_1.46.0         gtable_0.3.6            lattice_0.22-6          tzdb_0.5.0              vctrs_0.6.5            
 [6] tools_4.4.2             generics_0.1.4          parallel_4.4.2          tibble_3.2.1            RSQLite_2.3.9          
[11] blob_1.2.4              pkgconfig_2.0.3         Matrix_1.7-1            RColorBrewer_1.1-3      lifecycle_1.0.4        
[16] GenomeInfoDbData_1.2.13 stringr_1.5.1           compiler_4.4.2          farver_2.1.2            Biostrings_2.74.1      
[21] codetools_0.2-20        pillar_1.11.0           crayon_1.5.3            BiocParallel_1.40.0     DelayedArray_0.32.0    
[26] cachem_1.1.0            abind_1.4-8             tidyselect_1.2.1        locfit_1.5-9.11         stringi_1.8.4          
[31] labeling_0.4.3          fastmap_1.2.0           grid_4.4.2              colorspace_2.1-1        cli_3.6.3              
[36] SparseArray_1.6.1       magrittr_2.0.3          S4Arrays_1.6.0          withr_3.0.2             readr_2.1.5            
[41] scales_1.4.0            UCSC.utils_1.2.0        bit64_4.6.0-1           XVector_0.46.0          httr_1.4.7             
[46] bit_4.5.0.1             hms_1.1.3               png_0.1-8               memoise_2.0.1           rlang_1.1.5            
[51] Rcpp_1.0.14             glue_1.8.0              DBI_1.2.3               vroom_1.6.5             rstudioapi_0.17.1      
[56] jsonlite_1.8.9          plyr_1.8.9              R6_2.6.1                zlibbioc_1.52.0        
