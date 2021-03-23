#' @export
binary <-
function(q,n=10){
  out <- NULL
  while (q > 0) {
    out <- c(q%%2, out)
    l   <- length(out)
    q <- q %/% 2
  }
  if(l>n) print("Error: n must be greater than length q in binary representation")
  return(c(rep(0, n-l), out))
}
