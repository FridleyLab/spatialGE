# Prepare stage: builds the inputs STdiff() hands off to the fit stage.

##
# stdiff_ensure_vst
# Ensure vst.variance.standardized is present in gene_meta for every sample in
# the STlist, computing and joining it if missing.
# @param x an STlist.
# @return the STlist with gene_meta updated.
#
stdiff_ensure_vst = function(x){
  for(i in 1:length(x@tr_counts)){
    # Check that vst is not already calculated
    if(length(grep('vst.variance.standardized', colnames(x@gene_meta[[i]]))) == 0){
      #x@gene_meta[[i]] = Seurat::FindVariableFeatures(x@counts[[i]], verbose=F) %>%
      x@gene_meta[[i]] = Seurat_FindVariableFeatures(sample_counts(x, i)) %>%
        tibble::rownames_to_column(var='gene') %>%
        dplyr::select('gene', 'vst.variance.standardized') %>%
        dplyr::left_join(x@gene_meta[[i]], ., by='gene')
    }
  }
  return(x)
}


##
# stdiff_build_combo_table
# Build the table of unique gene x cluster combinations to test, along with
# the original-to-coded annotation dictionary.
# @param x an STlist.
# @param samples the samples to test.
# @param annot the annotation column in x@spatial_meta being tested.
# @param topgenes the number of top variable genes to keep per sample (or NULL for all).
# @param clusters user-specified clusters to test (passed through to `prepare_stdiff_combo`).
# @param pairwise whether tests are pairwise or one-vs-rest.
# @param verbose either logical or an integer (0, 1, or 2) to increase verbosity.
# @return a list with `combo_df` and `meta_dict`.
#
stdiff_build_combo_table = function(x, samples, annot, topgenes, clusters, pairwise, verbose){
  # Extract spot/cell annotations for each sample
  metas = tibble::tibble()
  for(sample_name in samples){
    meta_tmp = as.character(unique(x@spatial_meta[[sample_name]][[annot]]))
    metas = dplyr::bind_rows(metas, tibble::tibble(samplename=sample_name, orig_annot=meta_tmp))
    rm(meta_tmp) # Clean env
  }

  # Find variable genes for each sample
  genes = tibble::tibble()
  for(sample_name in samples){
    # Get top variable genes
    genes_tmp = x@gene_meta[[sample_name]] %>% dplyr::arrange(dplyr::desc(vst.variance.standardized))
    if(!is.null(topgenes)){
      genes_tmp = genes_tmp %>% dplyr::slice_head(n=topgenes)
    }
    genes_tmp = as.vector(genes_tmp[['gene']])
    genes = dplyr::bind_rows(genes, tibble::tibble(samplename=sample_name, gene=genes_tmp))
    rm(genes_tmp) # Clean env
  }

  # Merge sample names, genes, and annotations in the same data frame.
  # This is deliberately many-to-many: every gene of a sample must combine with
  # every annotation of that sample to form the gene x cluster combinations
  # tested below. Declaring it silences dplyr's relationship warning, which
  # `multiple='all'` does not (multiple and relationship are separate checks).
  gene_and_meta = dplyr::full_join(genes, metas, by='samplename', relationship='many-to-many')
  rm(genes, metas) # Clean env

  # Create "dictionary" with coded annotations (to avoid potentially problematic characters)
  meta_dict = tibble::tibble(orig_annot=unique(gene_and_meta[['orig_annot']]),
                             coded_annot=paste0('c', seq(1, length(unique(gene_and_meta[['orig_annot']])))))

  # Append recoded annotations in gene_and_meta
  gene_and_meta = gene_and_meta %>% dplyr::left_join(., meta_dict, by='orig_annot') %>% dplyr::rename(meta=coded_annot)

  # Create a table with unique combinations of genes and annotations to test using parallelization
  combo_df = prepare_stdiff_combo(to_expand=gene_and_meta, user_clusters=clusters, pairwise=pairwise, verbose=verbose)
  rm(gene_and_meta) # Clean env

  return(list(combo_df=combo_df, meta_dict=meta_dict))
}


##
# prepare_stdiff_combo
# @return data frame with combinations
#
#prepare_stdiff_combo = function(to_expand=NULL, clusters=NULL, annot_dict=NULL, pairwise=NULL){
prepare_stdiff_combo = function(to_expand=NULL, user_clusters=NULL, pairwise=NULL, verbose=NULL){
  combo_df = tibble::tibble()
  for(sample_name in unique(to_expand[['samplename']])){
    # Extract annotations for a given sample
    to_expand_subset = to_expand %>% dplyr::filter(samplename == sample_name) #%>%
    #  dplyr::select('meta') %>% unlist() %>% unique()

    # Define clusters to test if NULL, or subset to those requested by user
    if(is.null(user_clusters)){
      annots_tmp = to_expand_subset %>% dplyr::select('meta') %>% unlist() %>% unique()
    } else{
      annots_tmp = to_expand_subset %>% dplyr::filter(orig_annot %in% user_clusters) %>% dplyr::select('meta') %>% unlist() %>% unique()
    }
    #annots_tmp = annots_tmp[annots_tmp %in% unlist(annot_dict[['coded_annot']][ annot_dict[['orig_annot']] %in% clusters_tmp ]) ]
    if(length(annots_tmp) >= 2){
      if(pairwise){
        # Get unique combinations of each two annotation, without repeating (e.g., c1 vs c2 and c2 vs c1)
        # Get combinations for each cluster
        combo_meta = tibble::tibble()
        annots_tmp2 = annots_tmp
        for(cl_tmp in annots_tmp){
          combo_meta = dplyr::bind_rows(combo_meta,
                                        expand.grid(cl_tmp,
                                                    grep(cl_tmp, annots_tmp2, value=T, invert=T),
                                                    stringsAsFactors=F))
          annots_tmp2 = grep(cl_tmp, annots_tmp2, value=T, invert=T)
        }
        names(combo_meta) = c('meta1', 'meta2')
        rm(cl_tmp, annots_tmp2) #Clean env

        # Add genes to each combination
        combo_meta_gene = tibble::tibble()
        for(comb in 1:nrow(combo_meta)){
          combo_meta_gene = dplyr::bind_rows(combo_meta_gene,
                                             tibble::tibble(combo_meta[comb, 1],
                                                            combo_meta[comb, 2],
                                                            gene=to_expand %>% dplyr::filter(samplename == sample_name) %>%
                                                              dplyr::select(gene) %>% unlist() %>% unique()))
        }
        combo_meta = combo_meta_gene
        rm(combo_meta_gene, comb) # Clean env
      } else{
        combo_meta = expand.grid(meta1=annots_tmp, meta2='other',
                                 gene=to_expand %>% dplyr::filter(samplename == sample_name) %>%
                                   dplyr::select(gene) %>% unlist() %>% unique(),
                                 stringsAsFactors=F)
      }
      rm(to_expand_subset, annots_tmp) # Clean env
      combo_df = dplyr::bind_rows(combo_df,
                                  combo_meta %>%
                                    tibble::add_column(samplename=sample_name, .before=1))
    } else{
      if(verbose >= 1L){
        cat(paste0('\t\tSkipping sample ', sample_name, '. Less than two clusters to compare.\n'))
      }
    }
  }
  combo_df = as.data.frame(combo_df) %>%
    dplyr::arrange(samplename, meta1, meta2, gene)

  return(combo_df)
}
