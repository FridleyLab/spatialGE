##
# This funtion renames duplicate genes in a gene count data frame. The gene
# retaining the original name has the most counts. Others are named by adding a
# suffix (_g1, _g2, _g3...)
#
#
genedups = function(df=NULL){

  # Get column with gene names.
  gene_names = df[[1]]

  # Find duplicate gene names.
  dup_gene_mask <- duplicated(gene_names)
  dup_gene_names <- gene_names[dup_gene_mask]

  df_mod = df
  # Loop through gene names to identify duplicate genes to remove.
  for(j in 1:length(dup_gene_names)){
    # Find row IDs matching each duplicate gene.
    dup_gene_idx <- grep(paste0("^", dup_gene_names[j]), unlist(df[1]))
#print(dup_gene_names[j])
    # Get sums of counts of the duplicate genes and create mask to flag
    # the gene with the highest count.
    dup_gene_countsums <- rowSums(df[dup_gene_idx, ][-1])
#print(dup_gene_countsums)
    dup_gene_countsums_mask <- dup_gene_countsums >= max(dup_gene_countsums)

    # Safety check in case there are two duplicate genes with the same count number.
    if(sum(dup_gene_countsums_mask) > 1){
      dup_gene_countsums_mask[dup_gene_countsums_mask][1] = FALSE
    }
#print(dup_gene_countsums_mask)

    for(gene in 1:length(dup_gene_idx)){
      if(!dup_gene_countsums_mask[gene]){
        df_mod[[1]][dup_gene_idx[gene]] = paste0(df_mod[[1]][dup_gene_idx[gene]], '_g', gene)
#print(df_mod[[1]][dup_gene_idx[gene]])
      }
    }

  }

  return(df_mod)
}

