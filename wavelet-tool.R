
multi.wt <- function(xx, mat){
  midp <- function(x)  {    (x[-1] + x[-length(x)]) / 2  }
  Morlet <- function(lag)  {    exp(-lag ^ 2 / 2 + 2i * pi * lag) / pi ^ 4  }
  scales = 2 ^ midp(seq(log2(2 * median(diff(xx))), log2(diff(range(xx))/2), length.out = length(xx) + 1))
  lmat = outer(seq(min(xx), max(xx), length = length(xx)), xx, "-")
  w = foreach(s = scales) %dopar%  { Conj(Morlet(lmat / s) / s ^ 0.5) %*% as.matrix(mat) }
  w = array(unlist(w), dim = c(length(xx), ncol(mat), length(scales)))
  w = aperm(w, perm = c(1, 3, 2))
  return(list(x = seq(min(xx), max(xx), length.out = length(xx)), y = scales, z = w))
}

mb.com <- function( wtlist, method = 'wmr', smoothing = 1 ){
  if(method == 'wmr'){
    with(wtlist, {
      mods = Mod(rowSums(z, dims = 2))
      smod = rowSums(Mod(z), dims = 2)
      lmat = outer(seq(min(x), max(x), length.out = length(x)), x, "-")
      Gauss <- function(lag) {  exp(-lag ^ 2 / 2) / sqrt(2 * pi)  }
      modrat = foreach(i = 1:length(y), .combine = cbind, .export = 'smoothing') %dopar% {
        kern = Gauss(lmat / y[i] / smoothing)
        modsv = kern %*% mods[,i]
        smodv = kern %*% smod[,i]
        modsv / smodv
      }
      resmat= modrat ; dim(resmat) = c(length(x), length(y))
      return(resmat)
    })
  } else
    if(method == 'wmf'){
      with(wtlist, {
        index <- c(1:length(y)) ; for(m in 1:length(y)){ index[m] <- sqrt( sum(Mod( z[,m,] )^2) / (dim(z)[3] * length(x)) ) }
        wmfmat <- matrix(NA, length(x), length(y)) ; for(m in 1:length(x)){ wmfmat[m,] <- rowMeans( z[m,,] / index ) }
        # index <- foreach(i = 1:length(y), .combine = c) %dopar% { sqrt( sum(Mod( z[,i,] )^2) / (dim(z)[3] * length(x)) ) }
        # wmfmat <- foreach(i = 1:length(x), .combine = rbind) %dopar% { if(length(y)>1){rowMeans( z[i,,] / index )}else{( z[i,,] / index ) } }
        resmat <- Mod(wmfmat) ; dim(resmat) = c(length(x), length(y))
        return(resmat)
      })
    }
}

mb.boot <- function(wtlist, reps = 20, method = 'wmr', smoothing = 1 ){
  mb.com <- match.fun(mb.com)
  mb.obs = mb.com(wtlist, method = method)
  with(wtlist, {
    nloc = length(x)
    nvars = dim(z)[3]
    nscales = length(y)
    exports = c("reps", "method", "mb.com", "mb.obs", "nscales", "smoothing")
    if(method == 'wmr' ){
      z.boot = foreach(i = 1:nscales, .combine = c, .export = exports, .packages = c('doParallel','foreach')) %dopar% {
        boot = foreach(j = 1:reps, .combine = cbind, .inorder = FALSE) %dopar% {
          rphase = t(array(runif(nvars, -pi, pi), dim = c(nvars, nloc)))
          zp = z[, i,, drop = FALSE] * complex(argument = rphase)
          as.vector(mb.com(list(x = x, y = y[i], z = zp), method = 'wmr', smoothing = smoothing))
        }
        res = foreach(j = 1:nloc, .combine = c) %dopar% {  ecdf(boot[j, ])(mb.obs[j, i])  }
        res
      }
    }else
      if(method == 'wmf'){
        boot = foreach(i = 1:reps, .combine = rbind, .export = exports, .packages = c('doParallel','foreach')) %dopar% {
          rphase = foreach(j = 1:nvars, .inorder = F) %dopar% { t(array(runif(nvars, -pi, pi), dim = c(nscales, nloc))) }
          rphase <- array(unlist(rphase), dim = c(nloc, nscales, nvars))
          zp = z * complex(argument = rphase)
          as.vector(mb.com(list(x = x, y = y, z = zp), method = 'wmf', smoothing = smoothing))
        }
        # z.boot = foreach(j = 1:(nloc*nscales), .combine = c) %dopar% {  ecdf(boot[ , j ])(mb.obs[ j ])  }
        z.boot <- num <- c(1:(nloc*nscales))
        for( j in num[num%%3==0]){z.boot[j] = ecdf(boot[ , j ])(mb.obs[ j ]) }
        for( j in num[num%%3==1]){z.boot[j] = ecdf(boot[ , j ])(mb.obs[ j ]) }
        for( j in num[num%%3==2]){z.boot[j] = ecdf(boot[ , j ])(mb.obs[ j ]) }
      }
    dim(z.boot) = c(length(x), length(y))
    return(z.boot)
  })
}

b.boot <- function(wtlist, reps = 20, method = 'wmr', smoothing = 1 ){
  mb.com <- match.fun(mb.com)
  mb.obs = mb.com(wtlist, method = method)
  with(wtlist, {
    nloc = length(x)
    nvars = dim(z)[3]
    nscales = length(y)
    exports = c("reps", "method", "mb.com", "mb.obs", "nscales", "smoothing")
    if(method == 'wmr' ){
      z.boot = foreach(i = 1:nscales, .combine = rbind, .export = exports, .packages = c('doParallel','foreach')) %dopar% {
        foreach(j = 1:reps, .combine = c, .inorder = FALSE) %dopar% {
          rphase = t(array(runif(nvars, -pi, pi), dim = c(nvars, nloc)))
          zp = z[, i,, drop = FALSE] * complex(argument = rphase)
          colSums(abs(apply(zp, c(1, 2), sum)))/colSums(apply(abs(zp), c(1, 2), sum))
        }
      }
    }else
      if(method == 'wmf'){
        z.boot = foreach(i = 1:reps, .combine = cbind, .export = exports, .packages = c('doParallel','foreach')) %dopar% {
          rphase = foreach(j = 1:nvars, .inorder = F) %dopar% { t(array(runif(nvars, -pi, pi), dim = c(nscales, nloc))) }
          rphase <- array(unlist(rphase), dim = c(nloc, nscales, nvars))
          zp = z * complex(argument = rphase)
          boot <- mb.com(list(x = x, y = y, z = zp), method = 'wmf', smoothing = smoothing)
          colMeans(boot)
        }
      }
    res = matrix(NA, nscales, 2)
    # res[,3] <- colSums(abs(apply(z, c(1, 2), sum)))/colSums(apply(abs(z), c(1, 2), sum))
    for (i in 1:nscales) {
      res[i,1] = mean(z.boot[i,]) + qnorm(0.025,mean=0.5,sd=1) * sd(z.boot[i,])
      res[i,2] = mean(z.boot[i,]) + qnorm(0.975,mean=0.5,sd=1) * sd(z.boot[i,])
    }
    return(res)
  })
}

maxmin <- function(x){(2*x-(max(x)+min(x)))/(max(x)-min(x))}

