#' @importFrom Rfast rowVars
#' @export
CELLiD = function(input.data, ref.data) {

  genes = intersect(rownames(input.data), rownames(ref.data))
  if (length(genes) <= 2000) {
    stop("Too few overlapped gene between input data and reference!")
  } else {

    input.data = input.data[genes,,drop = F]
    ref.data = ref.data[genes,,drop = F]
    predicted.cell = pbmclapply(
      1:ncol(input.data),
      function(j) {
        predicted = as.numeric(
          apply(ref.data, 2, function(i) {
            cor(as.numeric(i), as.numeric(input.data[, j]), method = "spearman", use="complete.obs")
          })
        )
        return(predicted)
      }, mc.cores = 10
    )

    predicted.cell = do.call(cbind, predicted.cell)
    rownames(predicted.cell) = colnames(ref.data)

    ct = apply(predicted.cell, 2, function(x) {
      return(as.numeric(which(rank(-x) <=5)))
    }, simplify = F)

    predicted.cell = pbmclapply(
      1:ncol(input.data), function(i) {
        ref = ref.data[,ct[[i]]]
        g = which(rank(-rowVars(as.matrix(ref))) <= 2000)
        ref = ref[g,]
        input = input.data[g,i]
        predict = apply(ref, 2, function(i) {
          cor(as.numeric(i), input, method = "spearman", use="complete.obs")
        })
        return(c(names(sort(predict, decreasing = T)[1:2]), as.numeric(sort(predict, decreasing = T)[1:2])))
      }, mc.cores = 10
    )

    predicted.cell = do.call(rbind, predicted.cell)
    predicted.cell[,3] = round(as.numeric(predicted.cell[,3]), 3)
    predicted.cell[,4] = round(as.numeric(predicted.cell[,4]), 3)
    return(predicted.cell)
  }
}

