##
# Get the number of cores to use in parallel as a function of the number of
# data to process. It will not yield a higher number of cores than half of the
# total available cores. Will default to 1 core if not Unix
# @param n, an integer representing the number of units to process
# @return cores, the number of cores to use
#
count_cores = function(n){
  cores = 1
  if(.Platform$OS.type == 'unix'){
    # Use parallelization (if possible) to read data.
    require("parallel")
    avail_cores = (parallel::detectCores()) / 2
    if(avail_cores <= n){
      cores = avail_cores
    } else{
      cores = n
    }
  }
  return(cores)
}
