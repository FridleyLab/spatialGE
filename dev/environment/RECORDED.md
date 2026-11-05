# Recorded build environment

## sessionInfo()

```
R version 4.5.3 (2026-03-11)
Platform: aarch64-apple-darwin20.0.0
Running under: macOS Sequoia 15.7.9

Matrix products: default
BLAS:   /opt/anaconda3/envs/spatialge-dev/lib/libblas.3.9.0.dylib 
LAPACK: /opt/anaconda3/envs/spatialge-dev/lib/liblapack.3.9.0.dylib  LAPACK version 3.9.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: system (macOS)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
  [1] SpatialExperiment_1.20.0    IRanges_2.44.0             
  [3] scSpatialSIM_0.1.4          R.methodsS3_1.8.2          
  [5] GSEABase_1.72.0             progress_1.2.3             
  [7] tiff_0.1-12                 goftest_1.2-3              
  [9] Biostrings_2.78.0           HDF5Array_1.38.0           
 [11] vctrs_0.7.3                 spatstat.random_3.5-1      
 [13] digest_0.6.39               png_0.1-9                  
 [15] shape_1.4.6.1               proxy_0.4-29               
 [17] registry_0.5-1              ggrepel_0.9.8              
 [19] deldir_2.0-4                parallelly_1.48.0          
 [21] magick_2.9.1                MASS_7.3-66                
 [23] reshape2_1.4.5              gstat_2.1-6                
 [25] foreach_1.5.2               BiocGenerics_0.56.0        
 [27] withr_3.0.3                 xfun_0.60                  
 [29] ggpubr_1.0.0                ellipsis_0.3.3             
 [31] memoise_2.0.1               s2_1.1.11                  
 [33] bench_1.1.4                 intervals_0.15.5           
 [35] janitor_2.2.1               Seqinfo_1.0.0              
 [37] spatstat_3.6-2              profvis_0.4.0              
 [39] ggpolypath_0.4.0            systemfonts_1.3.2          
 [41] zoo_1.9-0                   GlobalOptions_0.1.4        
 [43] pbapply_1.7-5               R.oo_1.27.1                
 [45] spData_2.3.5                Formula_1.2-6              
 [47] prettyunits_1.2.0           KEGGREST_1.50.0            
 [49] otel_0.2.0                  httr_1.4.9                 
 [51] rstatix_1.1.0               rhdf5filters_1.22.0        
 [53] globals_0.19.1              spaMM_4.7.0                
 [55] rhdf5_2.54.1                rstudioapi_0.19.0          
 [57] units_1.0-1                 generics_0.1.4             
 [59] sparsesvd_0.2-3             concaveman_1.2.0           
 [61] curl_8.0.0                  spacetime_1.3-3            
 [63] S4Vectors_0.48.1            sfsmisc_1.1-25             
 [65] h5mread_1.2.1               ScaledMatrix_1.18.0        
 [67] polyclip_1.10-7             SparseArray_1.10.10        
 [69] fftwtools_0.9-11            RcppEigen_0.3.4.0.2        
 [71] xtable_1.8-8                stringr_1.6.0              
 [73] doParallel_1.0.17           evaluate_1.0.5             
 [75] S4Arrays_1.10.1             GenomicRanges_1.62.1       
 [77] hms_1.1.4                   irlba_2.3.7                
 [79] colorspace_2.1-3            hdf5r_1.3.16               
 [81] readxl_1.5.0                spatstat.data_3.1-9        
 [83] magrittr_2.0.5              snakecase_0.11.1           
 [85] spdep_1.4-2                 readr_2.2.0                
 [87] lattice_0.23-1              spatstat.geom_3.8-2        
 [89] future.apply_1.20.2         iotools_0.4-0.1            
 [91] XML_3.99-0.24               matrixStats_1.5.0          
 [93] xts_0.14.3                  class_7.3-23               
 [95] pillar_1.11.1               nlme_3.1-171               
 [97] iterators_1.0.14            EBImage_4.52.0             
 [99] compiler_4.5.3              beachmat_2.26.0            
[101] stringi_1.8.9               sf_1.1-2                   
[103] SummarizedExperiment_1.40.0 tensor_1.5.1               
[105] minqa_1.2.8                 lubridate_1.9.5            
[107] devtools_2.5.2              plyr_1.8.9                 
[109] msigdbr_26.1.1              crayon_1.5.3               
[111] abind_1.4-8                 locfit_1.5-9.12            
[113] sp_2.2-3                    bit_4.6.0                  
[115] waldo_0.6.2                 dplyr_1.2.1                
[117] codetools_0.2-20            textshaping_1.0.5          
[119] BiocSingular_1.26.1         slam_0.1-56                
[121] e1071_1.7-17                GetoptLong_1.1.1           
[123] splines_4.5.3               circlize_0.4.18            
[125] Rcpp_1.1.2                  sparseMatrixStats_1.22.0   
[127] cellranger_1.1.0            knitr_1.52                 
[129] blob_1.3.0                  clue_0.3-68                
[131] fs_2.1.0                    listenv_1.0.0              
[133] checkmate_2.3.4             DelayedMatrixStats_1.32.0  
[135] Rdpack_2.6.6                pkgbuild_1.4.8             
[137] GSVA_2.4.9                  ggsignif_0.6.4             
[139] tibble_3.3.1                Matrix_1.7-6               
[141] tzdb_0.5.0                  svglite_2.2.2              
[143] tweenr_2.0.3                pkgconfig_2.0.3            
[145] tools_4.5.3                 cachem_1.1.0               
[147] rbibutils_2.4.1             RSQLite_3.53.3             
[149] viridisLite_0.4.3           DBI_1.3.0                  
[151] numDeriv_2016.8-1.1         wordspace_0.2-9            
[153] fastmap_1.2.0               rmarkdown_2.32             
[155] scales_1.4.0                grid_4.5.3                 
[157] usethis_3.2.2               broom_1.0.13               
[159] FNN_1.1.4.1                 BiocManager_1.30.27        
[161] dotCall64_1.2               graph_1.88.1               
[163] carData_3.0-6               rpart_4.1.27               
[165] farver_2.1.2                geoR_1.9-6                 
[167] reformulas_0.4.4            mgcv_1.9-4                 
[169] wk_0.9.5                    roxygen2_8.1.0             
[171] MatrixGenerics_1.22.0       cli_3.6.6                  
[173] purrr_1.2.2                 stats4_4.5.3               
[175] lifecycle_1.0.5             uwot_0.2.5                 
[177] Biobase_2.70.0              arrow_25.0.0               
[179] sessioninfo_1.2.4           backports_1.5.1            
[181] BiocParallel_1.44.0         annotate_1.88.0            
[183] timechange_0.4.0            gtable_0.3.6               
[185] rjson_0.2.23                progressr_1.0.0            
[187] parallel_4.5.3              spatstat.linnet_3.5-3      
[189] testthat_3.3.2              jsonlite_2.0.0             
[191] rex_1.2.2                   bitops_1.1-0               
[193] ggplot2_4.0.3               kableExtra_1.4.1           
[195] bit64_4.8.6                 assertthat_0.2.1           
[197] brio_1.1.5                  spatstat.utils_3.2-5       
[199] SeuratObject_5.4.0.9001     spatstat.univar_3.2-0      
[201] R.utils_2.13.0              dynamicTreeCut_1.63-1      
[203] htmltools_0.5.9             sctransform_0.4.3          
[205] glue_1.8.1                  tcltk_4.5.3                
[207] spam_2.11-4                 XVector_0.50.0             
[209] RCurl_1.98-1.20             ROI_1.0-2                  
[211] classInt_0.4-11             jpeg_0.1-11                
[213] gridExtra_2.3.1             boot_1.3-32                
[215] R6_2.6.1                    SingleCellExperiment_1.32.0
[217] tidyr_1.3.2                 RcppProgress_0.4.2         
[219] cluster_2.1.8.3             Rhdf5lib_1.32.0            
[221] pkgload_1.5.3               memuse_4.2-3               
[223] DelayedArray_0.36.1         tidyselect_1.2.1           
[225] covr_3.6.5                  ggforce_0.5.0              
[227] xml2_1.6.0                  car_3.1-5                  
[229] AnnotationDbi_1.72.0        future_1.75.0              
[231] khroma_1.17.0               rsvd_1.0.5                 
[233] KernSmooth_2.23-26          S7_0.2.2                   
[235] spatstat.model_3.7-2        splancs_2.01-45            
[237] data.table_1.18.6.1         htmlwidgets_1.6.4          
[239] ComplexHeatmap_2.26.1       RColorBrewer_1.1-3         
[241] rlang_1.3.0                 spatstat.sparse_3.2-0      
[243] tidyverse_2.0.0             spatstat.explore_3.8-2     
[245] remotes_2.5.0              
```

## extSoftVersion()

```
                                                       zlib 
                                                    "1.3.2" 
                                                      bzlib 
                                       "1.0.8, 13-Jul-2019" 
                                                         xz 
                                                    "5.8.3" 
                                                 libdeflate 
                                                     "1.25" 
                                                       PCRE 
                                         "10.47 2025-10-21" 
                                                        ICU 
                                                     "78.3" 
                                                        TRE 
                                  "TRE 0.8.0 R_fixes (BSD)" 
                                                      iconv 
              "GNU libiconv 1.18 /usr/lib/libiconv.2.dylib" 
                                                   readline 
                                                      "8.3" 
                                                       BLAS 
"/opt/anaconda3/envs/spatialge-dev/lib/libblas.3.9.0.dylib" 
```

## BLAS / LAPACK

```
/opt/anaconda3/envs/spatialge-dev/lib/liblapack.3.9.0.dylib 
/opt/anaconda3/envs/spatialge-dev/lib/libblas.3.9.0.dylib 
```
