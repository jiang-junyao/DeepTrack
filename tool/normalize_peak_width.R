subdivideGRanges <- function (x, subsize = 500,core = NULL) 
{
  subdivide.range <- function (start.pos, end.pos, subsize = 100) 
  {
    width <- end.pos - start.pos + 1
    if (width < 2 * subsize) {
      stop("Width is less than 2 times subsize")
    }
    nchunks <- floor(width/subsize)
    relative.starts <- round(0:(nchunks - 1) * width/nchunks)
    relative.ends <- round(1:nchunks * width/nchunks) - 1
    return(list(start.pos = start.pos + relative.starts, end.pos = start.pos + 
                  relative.ends, length = nchunks))
  }
  
  subdivideIRanges <- function (x, subsize = 100) 
  {
    if (length(x) == 0) {
      return(x)
    }
    start.pos <- start(x)
    end.pos <- end(x)
    widths <- width(x)
    nsubranges <- pmax(1, floor(widths/subsize))
    out.start.pos <- numeric(sum(nsubranges))
    out.end.pos <- numeric(sum(nsubranges))
    out.idx <- 1
    for (i in 1:length(x)) {
      if (widths[i] < 2 * subsize) {
        out.start.pos[out.idx] <- start.pos[i]
        out.end.pos[out.idx] <- end.pos[i]
        out.idx <- out.idx + 1
      }
      else {
        sr <- subdivide.range(start.pos[i], end.pos[i], 
                              subsize)
        out.start.pos[out.idx:(out.idx + sr$length - 1)] <- sr$start.pos
        out.end.pos[out.idx:(out.idx + sr$length - 1)] <- sr$end.pos
        out.idx <- out.idx + sr$length
      }
    }
    IRanges(start = out.start.pos, end = out.end.pos)
  }
  
  
  if (length(x) == 0) {
    return(x)
  }
  if (length(subsize) > 1) {
    stop("The subsize argument should be a single number: the desired width of the subdivided ranges")
  }
  if (is.null(core)) {
    x <- map(seq_along(x),function(i){
      fold = floor(width(x[i])/500) 
      if (fold < 1 ) {
        fold <- 1
      }
      resize(x[i],width = 500 * fold ,fix = 'center')   
    },.progress = T) 
    x <- do.call(c,x)
  }else{
    require(furrr)
    plan(multisession,workers = core)
    options(future.globals.maxSize= 100*1024*1024^2) 
    x <- furrr::future_map(seq_along(x),function(i){
      fold = floor(width(x[i])/500) 
      if (fold < 1 ) {
        fold <- 1
      }
      resize(x[i],width = 500 * fold ,fix = 'center')   
    },.progress = T) 
    x <- do.call(c,x)
  } 
  x <- sort(GenomicRanges::reduce(x))
  gr_list <- lapply(levels(seqnames(x)), function(seqlvl) {
    if (!any(seqnames(x) == seqlvl)) {
      return(GRanges())
    }
    rg <- ranges(x[seqnames(x) == seqlvl])
    GRanges(seqnames = seqlvl, ranges = subdivideIRanges(rg,subsize), seqlengths = seqlengths(x))
  })
  do.call(c, gr_list) 
  
}