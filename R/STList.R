##
#' @title STList: Creation of STList objects
#' @description Creates an STList object from one or several spatial transcriptomic experiments.
#' @details
#' Objects of the class STList can be created from two sources:
#' \itemize{
#' \item Raw RNA counts (genes in rows, spots in columns) and spot coordinates (spots
#' in rows, with three columns: spot IDs, y position, and x position).
#' \item Visium output folders from Space Ranger. The folder should have the
#' structure resulting from spaceranger count.
#'}
#' Optionally, the user can input a system path to a table with data associated to
#' each spatial array (e.g., clinical variables). This sample file should contain
#' sample IDs in the first column matching at least partially the file names of the
#' count/coordinate files or Visium directories. The function will read data in
#' parallel if unix system available (Windows users will experience longer times
#' depending on the size of the data set)
#'
#' @param rnacounts The count data which can be provided in one of these formats:
#' \itemize{
#' \item File paths to comma- or tab-delimited files containing raw RNA counts, one
#' for each spatial array. Thye first column contain gene names with no duplicates allowed.
#' Isoforms can be provided as part of the gene names. Subsequent columns contain data for
#' each spot.
#' \item File paths to Visium output folders (one per spatial array) can be provided.
#' \item One named list of data frames with raw RNA counts (one data frame per spatial array).
#' \item One data frame with raw RNA counts. This option does not support sample data.
#' }
#' @param spotcoords The spot coordinates which can be provided in one of these formats:
#' \itemize{
#' \item File paths to comma- or tab-delimited files containing spot coordinates, one
#' for each spatial array. The files must contain three columns: spot IDs, y positions, and
#' x positions. The spot IDs must match the column names for each spot (column) in the RNA count files.
#' \item One named list of data frames with spot coordinates. The list names must match list
#' names of the RNA counts list.
#' \item One data frame with spot coordinates.
#' }
#' @param samples The data associated to each spatial array. This file can also include
#' system paths to RNA count and spot coordinate without the need to specify `rnacounts`
#' and `spotcoords`. One of the following options must be entered to create an STList:
#' \itemize{
#' \item A vector with sample names, which will be used to partial match RNA count and
#' spot coordinates file paths. A sample name must not match two file paths.
#' \item A sample file containing meta data for each spatial array. This file is a comma-
#' or tab-separated file with sample names in the first column. Paths to RNA count and
#' coordinate files can be in the second and third column respectively, and omiting the
#' `rnacounts` and `spotcoords` arguments. If Visium, only the second column with paths to
#' output folders is expected. Subsequent columns can contain variables associated to each
#' spatial array.
#' }
#' @param filterMT Logical. Whether or not to filtet mtDNA genes in Visium outputs.
#' @return x, the STList object containing the counts and coordinates, and optionally
#' the sample metadata.
#'
#' @examples
#' # Using included melanoma example (Thrane et al.)
#' # data_files <- system.file("extdata", package="spatialGEdev")
#' # count_files <- grep("genes", data_files, value=T)
#' # coord_files <- grep("mapping", data_files, value=T)
#' # clin_file <- grep("clinical", data_files, value=T)
#' # melanoma <- STList(rnacounts=count_files, spotcoords=coord_files, samples=clin_file)
#' # melanoma
#'
#' @export
#
#
STList = function(rnacounts=NULL, spotcoords=NULL, samples=NULL, filterMT=T) {

  # Check input type.
  input_check = detect_input(rnacounts=rnacounts, spotcoords=spotcoords, samples=samples)

  # Output error if input_check is an empty list (input format not recognized).
  if(rlang::is_empty(input_check)){
    stop('Input not recognized. Please refer to documentation.')
  }

  # Create lists to store count matrices from each sample.
  counts_df_list = list()
  coords_df_list = list()

  # CASE: SINGLE ARRAY, ONE RNA AND ONE COORDINATES.
  if(!is.null(input_check$rna)){
    if(input_check$rna == 'df' && input_check$coords == 'df'){
      counts_df_list[['sp_array']] = rnacounts
      coords_df_list[['sp_array']] = spotcoords
    }
  }

  # CASE: NAMED LIST OF DATAFRAMES WITH COUNTS AND NAMED LIST OF DATA FRAMES WITH COORDINATES.
  # METADATA INFO OPTIONAL.
  if(!is.null(input_check$rna)){
    if(input_check$rna == 'list_dfs' && input_check$coords == 'list_dfs'){
      # Check that RNA and coord lists have the same number of elements.
      if(length(rnacounts) == length(spotcoords)){
        # Test that names within the lists are the same.
        if(length(setdiff(names(rnacounts), names(spotcoords))) != 0){
          stop('The RNA count and coordinate lists do not have the same names.')
        }
        # Sort lists' names to order both.
        sorted_names = sort(names(rnacounts))
        counts_df_list = rnacounts[sorted_names]
        coords_df_list = spotcoords[sorted_names]
      } else{
        stop('The RNA count and coordinate do not have the same number of elements.')
      }
    }
  }

  # CASE: SAMPLE FIlE ONLY WITH FILE PATH(S) TO COUNT COORDINATE MATRICE WITHIN.
  if(is.null(input_check$rna) && is.null(input_check$coords) && !is.null(input_check$samples) && input_check$sample != 'samplesfile_visium'){
    # Check that all system paths provided exist.
    if(!file.exists(samples)){
      stop('Could not find samples file.')
    }

    # Read sample file and get sample names.
    delsample = input_check$samples[2]
    sample_names = readr::read_delim(samples, delim=delsample, col_types=readr::cols(), progress=F)
    sample_names = sample_names[[1]]

    # Create objects to store sample names found. Will be used to assess if expected samples exist.
    count_found=c()
    coord_found=c()
    for(i in sample_names){
      sample_i = readLines(samples)
      sample_i = grep(i, sample_i, value=T)
      sample_i = strsplit(sample_i, input_check$sample[2])
      #sample_i = sample_i[[1]][2]
      count_found = append(count_found, sample_i[[1]][2])
      coord_found = append(coord_found, sample_i[[1]][3])
    }

    # Use parallelization (if possible) to read data.
    require("parallel")

    # Define number of available cores to use.
    cores = 1
    if(.Platform$OS.type == 'unix'){
      avail_cores = parallel::detectCores()
      if(avail_cores > (length(rnacounts) + 1)){
        cores = (length(rnacounts) + 1)
      } else if( (avail_cores <= (length(rnacounts) + 1)) && avail_cores > 1){
        cores = avail_cores - 1
      }
    }

    # We have reached at least here...
    cat("Creating STList...\n")

    # Get delimiter of RNA counts and coordinates.
    count_testline = readLines(count_found[1], n = 1)
    coord_testline = readLines(coord_found[1], n = 1)
    if(grepl('\t', count_testline)){
      delrna = '\t'
    } else if(grepl(',', count_testline)){
      delrna = ','
    } else{
      stop('RNA counts are not comma or tab- delimited')
    }

    if(grepl('\t', coord_testline)){
      delcoords = '\t'
    } else if(grepl(',', coord_testline)){
      delcoords = ','
    } else{
      stop('Sport coordinates are not comma or tab- delimited')
    }

    # Use parallelization to read count data if possible.
    counts_df_list = mclapply(seq_along(sample_names), function(i){
      # Get file name as given by sample_names.
      sample_i = grep(sample_names[i], count_found, value=T)
      # Read filepaths.
      counts_df = readr::read_delim(sample_i, delim=delrna, col_types=readr::cols(), progress=F)
      return(counts_df)
    }, mc.cores=cores, mc.preschedule=T)
    # Name list elements.
    names(counts_df_list) = sample_names

    # Use parallelization to read coordinate data if possible.
    coords_df_list = mclapply(seq_along(sample_names), function(i){
      # Get file name as given by sample_names.
      sample_i = grep(sample_names[i], coord_found, value=T)
      # Read filepaths.
      coords_df = readr::read_delim(sample_i, delim=delcoords, col_types=readr::cols(), progress=F)
      return(coords_df)
    }, mc.cores=cores, mc.preschedule=T)
    # Name list elements.
    names(coords_df_list) = sample_names
  }

  # CASE: FILE PATH(S) TO COUNT, COORDINATE MATRICES, AND SAMPLE NAMES (FILE OR VECTOR).
  if(!is.null(input_check$rna)){
    if(input_check$rna == 'rnapath' && input_check$coords == 'coordpath'){
      # Check that all system paths provided exist.
      if(any(!file.exists(rnacounts)) || any(!file.exists(spotcoords))){
        stop('Could not find at least one of the count or coordinate files.')
      }

      # Check what type of sample names were provided.
      if(input_check$sample[1] == 'sample_names'){
        sample_names = samples
      } else if(input_check$sample[1] == 'samplesfile' || input_check$sample[1] == 'samplesfile_matrices'){
        delsample = input_check$samples[2]
        sample_names = readr::read_delim(samples, delim=delsample, col_types=readr::cols(), progress=F)
        sample_names = sample_names[[1]]
      } else{
        stop('Could not get sample names.')
      }

      # Check that sample_names are found once, and only once in RNA and coordinate file paths.
      # Create objects to store sample names found. Will be used to assess if expected samples exist.
      count_found=c()
      coord_found=c()
      for(i in sample_names){
        if(input_check$sample[1] == 'sample_names'){
          if(any(grepl(i, rnacounts)) && any(grepl(i, spotcoords))){
            count_found = append(count_found, i)
            coord_found = append(coord_found, i)
          }
        } else if(input_check$sample[1] == 'samplesfile' || input_check$sample[1] == 'samplesfile_matrices'){
          sample_i = readLines(samples)
          sample_i = grep(i, sample_i, value=T)
          sample_i = strsplit(sample_i, input_check$sample[2])
          count_found = append(count_found, sample_i[[1]][2])
          coord_found = append(coord_found, sample_i[[1]][3])
        }
      }
      # Output error if sample names provided do not match to filepaths.
      if(length(count_found) != length(rnacounts) || length(coord_found) != length(spotcoords)){
        stop(paste0('At least one of the provided sample names do not match with any of the RNA or coordinate filepaths.'))
      }

      # Use parallelization (if possible) to read data.
      require("parallel")

      # Define number of available cores to use.
      cores = 1
      if(.Platform$OS.type == 'unix'){
        avail_cores = parallel::detectCores()
        if(avail_cores > (length(rnacounts) + 1)){
          cores = (length(rnacounts) + 1)
        } else if( (avail_cores <= (length(rnacounts) + 1)) && avail_cores > 1){
          cores = avail_cores - 1
        }
      }

      # We have reached at least here...
      cat("Creating STList...\n")

      # Get delimiter of RNA counts.
      delrna = input_check$rna[2]
      # Use parallelization to read count data if possible.
      counts_df_list = mclapply(seq_along(sample_names), function(i){
        # Get file name as given by sample_names.
        sample_i = grep(sample_names[i], rnacounts, value=T)
        # Read filepaths.
        counts_df = readr::read_delim(sample_i, delim=delrna, col_types=readr::cols(), progress=F)
        return(counts_df)
      }, mc.cores=cores, mc.preschedule=T)
      # Name list elements.
      names(counts_df_list) = sample_names

      # Get delimiter of RNA counts.
      delcoords = input_check$coords[2]
      # Use parallelization to read coordinate data if possible.
      coords_df_list = mclapply(seq_along(sample_names), function(i){
        # Get file name as given by sample_names.
        sample_i = grep(sample_names[i], spotcoords, value=T)
        # Read filepaths.
        coords_df = readr::read_delim(sample_i, delim=delcoords, col_types=readr::cols(), progress=F)
        return(coords_df)
      }, mc.cores=cores, mc.preschedule=T)
      # Name list elements.
      names(coords_df_list) = sample_names
    }
  }

  # CASE: Visium output folders, WITH OR WITHOUT SAMPLE FILE.
  if(input_check$samples[1] != 'names_from_list_or_df'){
    if(input_check$samples[1] != 'samplesfile_matrices'){
      if(input_check$rna == 'visium_out' || input_check$samples == 'samplesfile_visium'){

        # Check what type of sample names were provided.
        if(input_check$sample[1] == 'sample_names'){
          sample_names = samples
        } else if(input_check$sample[1] == 'samplesfile' || input_check$sample[1] == 'samplesfile_visium'){
          delsample = input_check$samples[2]
          sample_names = readr::read_delim(samples, delim=delsample, col_types=readr::cols(), progress=F)
          sample_names = sample_names[[1]]
        } else{
          stop('Could not get sample names.')
        }

        fp_list = list()
        for(i in 1:length(sample_names)){
          # Get file name as given by sample_names or samples_file.
          if(input_check$sample[1] == 'samplesfile' || input_check$sample[1] == 'sample_names'){
            sample_i = grep(sample_names[i], rnacounts, value=T)
          }else if(input_check$sample[1] == 'samplesfile_visium'){
            sample_i = readLines(samples)
            sample_i = grep(sample_names[i], sample_i, value=T)
            sample_i = strsplit(sample_i, input_check$sample[2])
            sample_i = sample_i[[1]][2]
          }
          # Create list to store system paths.
          fp_list[[i]] = list()
          # Get all system paths within output folder.
          temp_fps = list.files(sample_i, recursive=T, include.dirs=T)
          # Find matches of relevant files wihtin output folders.
          vfeatures = grep('filtered_feature_bc_matrix/features.tsv.gz',  temp_fps, value=T)
          vbarcodes = grep('filtered_feature_bc_matrix/barcodes.tsv.gz', temp_fps, value=T)
          vcounts = grep('filtered_feature_bc_matrix/matrix.mtx.gz', temp_fps, value=T)
          vcoords  = grep('spatial/tissue_positions_list.csv', temp_fps, value=T)

          # Get complete system paths.
          vfeatures = paste0(sample_i, '/', vfeatures)
          vbarcodes = paste0(sample_i, '/', vbarcodes)
          vcounts = paste0(sample_i, '/', vcounts)
          vcoords  = paste0(sample_i, '/', vcoords)

          # Filter out 'SPATIAL_RNA_COUNTER' folders (intermediate files from Space Ranger?).
          vcoords = vcoords[!grepl('SPATIAL_RNA_COUNTER', vcoords)]

          # Test that all files have been found.
          needed_fps_test = c(!grepl('gz', vfeatures), !grepl('gz', vbarcodes), !grepl('gz', vcounts), !grepl('csv', vcoords))
          if(any(needed_fps_test)){
            stop(paste0('Visium output folder #', i, ' does not have all necessary files OR sample names in samples file and system paths do not match.'))
          }

          fp_list[[i]]$features = vfeatures
          fp_list[[i]]$barcodes = vbarcodes
          fp_list[[i]]$counts = vcounts
          fp_list[[i]]$coords = vcoords
          fp_list[[i]]$runname = sample_names[i]
        }

        # Use parallelization (if possible) to read data.
        require("parallel")

        # Define number of available cores to use.
        cores = 1
        if(.Platform$OS.type == 'unix'){
          avail_cores = parallel::detectCores()
          if(avail_cores > (length(rnacounts) + 1)){
            cores = (length(rnacounts) + 1)
          } else if( (avail_cores <= (length(rnacounts) + 1)) && avail_cores > 1){
            cores = avail_cores - 1
          }
        }

        # We have reached at least here...
        cat("Creating STList...\n")

        # Use parallelization to read count data if possible.
        output_temp = mclapply(seq_along(sample_names), function(i){
          # Get file name as given by sample_names.
          sample_i = grep(sample_names[i], rnacounts, value=T)
          # Process Visium folder.
          visium_processed = import_Visium(features_fp=fp_list[[i]]$features,
                                               barcodes_fp=fp_list[[i]]$barcodes,
                                               counts_fp=fp_list[[i]]$counts,
                                               coords_fp=fp_list[[i]]$coords,
                                               filterMT=filterMT)
          return(visium_processed)
        }, mc.cores=cores, mc.preschedule=T)

        # Organize the paralellized output into corresponding lists.
        for(i in 1:length(output_temp)){
          counts_df_list[[fp_list[[i]]$runname]] = output_temp[[i]]$rawcounts
          coords_df_list[[fp_list[[i]]$runname]] = output_temp[[i]]$coords
        }

      }
    }
  }
  # Process the count and coordinate list.
  for(i in 1:length(names(counts_df_list))){

    # Get loop current name to be processed.
    name_i = names(counts_df_list)[i]

    # Column names of the count data are simplified using the clean_names() function.
    # NOTE: May need to reconsider use later.
    counts_df_list[[name_i]] = janitor::clean_names(counts_df_list[[name_i]])

    # Clean sample names in coordinates data using the clean_names(), so that they
    # mirror the column names in the count data frame.
    # NOTE: May need to reconsider use later.
    coords_df_list[[name_i]][, 1] = janitor::make_clean_names(unlist(coords_df_list[[name_i]][, 1]))

    # Test that spot names are the same in both count and coordinate data frames.
    if(length(setdiff(colnames(counts_df_list[[name_i]])[-1], unlist(coords_df_list[[name_i]][, 1]))) != 0){
      stop(paste0('The spots in the count  data (columns) and coordinate data (rows) do not match in spatial array (\"', name_i, '\").'))
    }

    # Sort coordinate data according to third column in the coordinate data frame.
    coords_df_list[[name_i]] = coords_df_list[[name_i]][order(coords_df_list[[name_i]][, 3]), ]

    # Order column names in count data frame according to sorted coordinate data.
    counts_df_list[[name_i]] = counts_df_list[[name_i]][, c(colnames(counts_df_list[[name_i]][1]), unlist(coords_df_list[[name_i]][, 1]))]

    # Test for duplicated gene names in count data.
    gene_names = counts_df_list[[name_i]][, 1]
    dup_genes_mask = duplicated(gene_names)
    if(sum(dup_genes_mask) > 0){
      stop(paste0('There are duplicated feature/gene names in spatial array (\"', name_i, '\").'))
    }

    # Put column names to coordinate data.
    colnames(coords_df_list[[name_i]] ) = c('libname', 'ypos', 'xpos')

  }

  # Test if clinical data is available.
  if(!is.null(input_check$samples)){
    if(input_check$samples[1] == 'samplesfile' || input_check$samples[1] == 'samplesfile_visium' || input_check$samples[1] == 'samplesfile_matrices'){
      # Get delimiter of samples file and read file.
      del = input_check$samples[2]
      samples_df = readr::read_delim(samples, delim=del, col_types=readr::cols(), progress=F)

      # Get sorted list names and fetch corresponding rows from the samplefile.
      sorted_names = sort(names(counts_df_list))
      samples_df = samples_df[samples_df[[1]] %in% sorted_names, ]

      if(input_check$samples[1] == 'samplesfile_visium'){
        samples_df = samples_df[, -2]
      } else if(input_check$samples[1] == 'samplesfile_matrices'){
        samples_df = samples_df[, -c(2:3)]
      }

      # Test that number of rows in sample file correspond to number of elements in list.
      if(nrow(samples_df) != length(counts_df_list)){
        stop("The number of rows in the clinical data is not the same as the number of spatial arrays.")
      }
    }else{
      samples_df = tibble::tibble()
    }
  }else{
    samples_df = tibble::tibble()
  }

  # Creates STList object from both count and coordinates data.
  STList_obj <- new("STList",
                    counts=counts_df_list,
                    coords=coords_df_list,
                    clinical=samples_df,
                    voom_counts=list(),
                    gene_stdev=list(),
                    gene_het=list(),
                    gene_krige=list(),
                    cell_deconv=list(),
                    cell_het=list(),
                    cell_krige=list(),
                    gene_krige_data=list(),
                    deconv_krige_data=list(),
                    st_clusters=list(),
                    pheno_plots=list()
  )

  return(STList_obj)

}
