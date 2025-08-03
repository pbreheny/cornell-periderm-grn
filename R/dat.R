#' Check for discrepancies between local data and remote data
#' 
#' @param remote   Path to base of remote directory
#' 
#' @examples
#' if (dir.exists('~/argon_home/GRN')) dat('~/argon_home/GRN')

dat <- function(remote) {
  dats <- list.files(recursive=TRUE, all.files=TRUE) %>% extract(basename(.) == '.dat-config')
  for (d in dats) {
    
    # Setup
    loc_dir <- dirname(d)
    rem_dir <- paste0(remote, '/', loc_dir)
    loc_f <- setdiff(list.files(loc_dir), 'old')
    rem_f <- setdiff(list.files(rem_dir), 'old')
    loc_md5 <- vapply(paste0(loc_dir, '/', loc_f), tools::md5sum, '')
    rem_md5 <- vapply(paste0(rem_dir, '/', rem_f), tools::md5sum, '')
    names(loc_md5) <- loc_f
    names(rem_md5) <- rem_f
    
    cat('In directory ', loc_dir, '/ \n', sep='')
    
    # Local not in remote
    prb <- setdiff(loc_f, rem_f)
    if (length(prb)) {
      cat(crayon::red('  The following local files do not appear in remote:\n    '))
      cat(prb, sep='\n    ')
    } else {
      cat(crayon::green('  All local files appear in remote.\n'))
    }
    
    # Remote not in local
    prb <- setdiff(rem_f, loc_f)
    if (length(prb)) {
      cat(crayon::red('  The following remote files do not appear in local:\n    '))
      cat(prb, sep='\n    ')
    } else {
      cat(crayon::green('  All remote files appear in local.\n'))
    }
    
    # Remote and local don't agree
    x <- intersect(loc_f, rem_f)
    prb <- which(loc_md5[x] != rem_md5[x])
    if (length(prb)) {
      cat(crayon::red('  The following files are different in local and remote:\n    '))
      cat(x[prb], sep='\n    ')
    } else {
      cat(crayon::green('  No differences between remote/local files.\n'))
    }
    cat('\n')
  }
}
