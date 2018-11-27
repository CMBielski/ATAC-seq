#!/usr/bin/env Rscript

'%!in%' <- function(x,y)!('%in%'(x,y))

catverbose <- function(...) {
  cat(format(Sys.time(), "%Y%m%d %H:%M:%S |"), ..., "\n")
}

write_batch_script <- function(peaks, bams, max, out, dir, genome, win.size, squish=T) {
  sink(out)
  peaks <- peaks[1:max]
  cat('new')
  cat(paste0('\ngenome ', genome))
  dat <- peaks
  bams_to_load <- paste(bams$path,
                        collapse = ",")
  cat(paste0('\nload ', bams_to_load))
  cat(paste0('\nsnapshotDirectory ', dir))
  cat('\nmaxPanelHeight 200\n')
  for (i in 1:nrow(dat)) {
    cat(paste0('\ngoto ', dat$Chr[i], ':', dat$Start[i]-win.size, '-', dat$End[i]+win.size))
    cat('\nsort base')
    if (squish) {
      cat('\nsquish')
    }
    cat(paste0('\nsnapshot ', paste0("Peak", i), "_", dat$Chr[i], "_", dat$Start[i], '-', dat$End[i], '.png'))
    cat("\n")
  }
  sink()
}

if( ! interactive() ) {
  
  pkgs = c('data.table', 'argparse', 'stringr')
  tmp <- lapply(pkgs, require, character.only = T)
  rm(tmp)
  
  parser=ArgumentParser()
  parser$add_argument('-p', '--peaks', type='character', help='Differential peaks (BED format)')
  parser$add_argument('-b', '--bams', type='character', help='List of bam files')
  parser$add_argument('-m', '--max', type='character', default='50', help='Max number of peaks')
  parser$add_argument('-o', '--out', type='character', help='Output batch script file')
  parser$add_argument('-d', '--dir', type='character', default = NULL, help='Screenshot directory')
  parser$add_argument('-g', '--genome', type='character', default='hg38', help='Genome version')
  parser$add_argument('-w', '--window', type='character', default='50', help='IGV window size')
  parser$add_argument('-sq', '--squish', action="store_true", default=T, help='Squish reads')
  args=parser$parse_args()
  
  peaks <- fread(args$peaks)
  bams <- fread(args$bams)
  max = as.numeric(args$max)
  out <- args$out
  if (is.null(dir)) {
    dir.create(file.path(getwd(), "igv"))
    dir = "igv"
  } else {
    dir <- args$dir
  }
  genome <- args$genome
  win.size <- as.numeric(args$window)
  squish <- args$squish

  write_batch_script(peaks, bams, out, dir, 
                     genome, win.size, squish)
  
}
