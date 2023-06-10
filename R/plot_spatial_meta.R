##
# @title plot_spatial_meta: Plot cluster memberships of ST spots
# @description Generates a plot of the location of spots within an spatial array,
# and colors them according to spatially-weighted hierarchical clustering assignments.
# @details
# The function takes an STList with cluster memberships and plots the spots with
# colors indicating the cluster they belong to. Optionally, the user can annotate
# tumor/stroma compartments if ESTIMATE scores are available.
#
# @param x an STList with hierarchical cluster memberships.
# @param samples an integer indicating the spatial array(s) to be plotted. Numbers
# follow the order of `names(x@counts)`. If NULL, will plot all spatial arrays.
# @param ks the k values to plot
# @param ws the spatial weights to plot
# @param plot_meta a column name in `x@spatial_meta` to plot
# @param color_pal a string of a color palette from khroma or RColorBrewer, or a
# vector with colors with enough elements to plot categories.
# @param visium whether or not to reverse axes for Visium slides.
# @param ptsize, a number specifying the size of the points. Passed to `size` aesthetic.
# @return a list with the requested plots.
#
# @examples
# # In this example, melanoma is an STList.
# # cluster_p <- plot_STclusters(melanoma, samples=c(2,3), visium=F)
#
# @export
#
#' @import ggplot2
#' @importFrom magrittr %>%
#
#
plot_spatial_meta = function(x, samples=NULL, ks='dtc', ws=NULL, deepSplit=NULL, plot_meta=NULL, color_pal=NULL, visium=T, ptsize=NULL){

  # Define which samples to plot
  if(is.null(samples)){
    samples = names(x@spatial_meta)
  } else{
    if(is.numeric(samples)){
      samples = names(x@spatial_meta)[samples]
    }
    if(length(grep(paste0(samples, collapse='|'), names(x@spatial_meta))) == 0){
      stop('The requested samples are not present in the STList spatial metadata.')
    }
  }

  # Define columns to plot
  if(is.null(plot_meta)){
    plot_meta = grep('^stclust_spw', colnames(x@spatial_meta[[1]]), value=T)

    if(!is.null(ws)){
      if(any(ws == 0)){ # To avoid zero ('0') matching other weights
        ws_tmp = ws[ws != 0]
        plot_meta_tmp = grep('stclust_spw0_|stclust_spw0$', plot_meta, value=T)
        if(length(ws_tmp) > 0){ # In case there are other ws in addition to '0'
          plot_meta = c(plot_meta_tmp,
                        grep(paste0('stclust_spw', ws_tmp, collapse='|'), plot_meta, value=T))
        } else{
          plot_meta = plot_meta_tmp
        }
        rm(ws_tmp, plot_meta_tmp)
      } else{
        plot_meta = grep(paste0('stclust_spw', ws, collapse='|'), plot_meta, value=T)
      }
    }

    if(ks[1] != 'dtc'){
      plot_meta = grep(paste0('_k', ks,'$', collapse='|'), plot_meta, value=T)
    } else if(ks[1] == 'dtc'){
      if(!is.null(deepSplit)){
        plot_meta = grep(paste0('_dspl', stringr::str_to_title(as.character(deepSplit)), '$', collapse='|'), plot_meta, value=T)
      } else{
        plot_meta = grep('_dspl', plot_meta, value=T)
      }
    } else{
      stop('Specify one or several k values to plot, or use ks=\'dtc\' (default).')
    }
  }

  # Check that the meta data column exists
  if(length(grep(paste0('^', plot_meta, '$', collapse='|'),  colnames(x@spatial_meta[[1]]))) == 0){
    stop('No metadata column or clustering parameters were specified. Or specified parameters do not exist in metadata.')
  }

  plot_list = list()
  for(s in samples){
    # Extract metadata for specific sample
    df_tmp = x@spatial_meta[[s]]

    # Define size of points
    if(is.null(ptsize)){
      ptsize = 0.5
    }

    for(metacol in plot_meta){
      # Set default color if NULL input
      if(is.null(color_pal)){
        color_pal = 'light'
        if(is.numeric(x@spatial_meta[[s]][[metacol]])){
          color_pal = 'sunset'
        }
      }

      df_tmp2 = df_tmp %>%
        dplyr::select(libname, ypos, xpos, meta:=!!metacol)

      if(!is.numeric(df_tmp2[['meta']]) & length(unique(df_tmp2[['meta']])) < 100){
        # Convert meta data to factor
        df_tmp2 = df_tmp2 %>%
          dplyr::mutate(meta=tidyr::replace_na(as.character(meta), 'No_Data')) %>%
          dplyr::mutate(meta=as.factor(meta))

        # Create color palette.
        meta_cols = color_parse(color_pal, n_cats=length(unique(df_tmp2[['meta']])))
        names(meta_cols) = unique(df_tmp2[['meta']])
        if(any(names(meta_cols) == 'No_Data')){
          meta_cols[names(meta_cols) == 'No_Data'] = 'gray50'
        }
      }

      # Prepare titles for plots
      if(grepl('^stclust_', metacol)){
        title_w = as.character(stringr::str_extract(metacol, paste0("spw0\\.?[0-9]*"))) %>% gsub('spw', '', .)
        if(grepl('_k[0-9]+$', metacol)){
          title_k = as.character(stringr::str_extract(metacol, paste0("_k[0-9]+"))) %>% gsub('_k', '', .)
          title_p = paste0("STclust k=", title_k, "\nspatial weight=", title_w, '\nsample: ')
        } else if(grepl('_dspl[\\.0-9TrueFalse]+$', metacol)){
          title_dspl = stringr::str_extract(metacol, '[\\.0-9]+$|True$|False$')
          title_p = paste0("STclust (dtc; deepSplit=", title_dspl, ")\nSpatial weight=", title_w, '\nSample: ')
        }
        title_p = paste0(title_p, s)
        title_leg = 'Clusters'
      } else{
        title_p = paste0('Sample: ', s)
        title_leg = plot_meta
      }

      # Create plot
      p = ggplot2::ggplot(df_tmp2) +
        ggplot2::geom_point(ggplot2::aes(x=xpos, y=ypos, color=meta), size=ptsize)
      # Assign color palette to plot for categorical or numerical variable
      if(is.factor(df_tmp2[['meta']])){
        p = p + ggplot2::scale_color_manual(values=c(meta_cols))
      } else{
        if(!is.null(color_pal)){
          # Get color palette and number of colors needed.
          # Get names of Khroma colors.
          khroma_cols = khroma::info()
          khroma_cols = khroma_cols$palette
          if(color_pal %in% khroma_cols){
            p = p + khroma::scale_color_picker(palette=color_pal)
          } else{
            p = p + ggplot2::scale_color_distiller(palette=color_pal)
          }
        }
      }

      #p = p + ggplot2::ggtitle(title_p) + ggplot2::theme_void() # MAY 09, 2023 PUT META DATA NAME ON LEGEND TITLE, NOT PLOT TITLE

      if(!is.numeric(df_tmp2[['meta']])){ # Test if it's not numeric and make legend spots/dots larger
        p = p +
          ggplot2::guides(color=guide_legend(override.aes=list(size=ptsize+1)))
      }
      p = p +
        labs(color=title_leg, title=title_p) + ggplot2::theme_void()

      if(visium){
        p = p + ggplot2::scale_y_reverse() + ggplot2::coord_fixed(ratio=1)
      } else{
        p = p + ggplot2::coord_fixed(ratio=1)
      }
      #p = p + ggplot2::theme(legend.title=ggplot2::element_blank()) # MAY 09, 2023 PUT META DATA NAME ON LEGEND TITLE, NOT PLOT TITLE

      plot_list[[paste0(s, '_', metacol)]] = p
    }
  }
  return(plot_list)
}

