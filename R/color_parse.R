##
# @title color_parse: Creates a color palette
# @description Uses Khroma or RColorBrewer to return the colors of a palette name.
# @details
# This function takes a character string and uses either khroma or RColoBrewer to
# create a color palette. The function first looks in khroma, then RColoBrewer. In
# other words, khroma colors have priority.
#
# @param color_pal A name of a Khroma or RColorBrewer. Alternatively, a vector with
# colors of lenght equal or larger than number of categories.
# @paran n_cats The number of colors to produce. If NULL, assumes 5 colors.
# @return cat_cols A color palette
#
#
color_parse = function(color_pal=NULL, n_cats=NULL){
 
  # Get color palette and number of colors needed.
  # Get names of Khroma colors.
  khroma_cols = khroma::info()
  khroma_cols = khroma_cols$palette
  
  # Assume 5 categories if n_cats not provided (for kriging/quilts).
  if(is.null(n_cats)){
    n_cats = 5
  }
  
  # Test if input is a Khroma name or RColorBrewer.
  # If so, create palette.
  if(color_pal[1] %in% khroma_cols){
    p_palette = khroma::colour(color_pal[1])
    cat_cols = as.vector(p_palette(n_cats))
  }else if(color_pal[1] %in% rownames(RColorBrewer::brewer.pal.info)){
    cat_cols = RColorBrewer::brewer.pal(n_cats, color_pal[1])
  }else{ # Test if user provided a vector of colors.
    if(length(color_pal) >= n_cats){
      cat_cols = color_pal[1:n_cats]
    } else{
      stop('Provide enough colors to plot or an appropriate Khroma/RColorBrewer palette.')
    }
  }
  
  return(cat_cols)
}
