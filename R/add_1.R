#' @export
add_1 <-
function(j) {
  m <- length(j)
  if (!j[1]) 	j[1] <- 1
  else {
    j[1] <- 0
    j[2:m] <- add_1(j[2:m])
  }
  return(j)
}
