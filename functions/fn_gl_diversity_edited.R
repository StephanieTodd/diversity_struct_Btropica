
gl.diversity.edited <- function (gl, spectrumplot = TRUE, confiplot = FALSE, probar = TRUE, 
          table = "DH") 
{
  if (is.null(pop(gl))) 
    pop(gl) <- factor(rep("all", nInd(gl)))
  pops <- seppop(gl)
  if (probar) {
    pb <- txtProgressBar(0, 8, style = 3, width = 20)
    cat(" Counting missing loci...           ")
  }
  nlocpop <- lapply(pops, function(x) sum(!is.na(colMeans(as.matrix(x), 
                                                          na.rm = T))))
  if (probar) {
    setTxtProgressBar(pb, 1)
    cat(" Calculating zero_H/D_alpha ...           ")
  }
  zero_H_alpha_es <- lapply(pops, function(x) {
    dummys <- ((colMeans(as.matrix(x), na.rm = T)%%2) > 0) + 
      1 - 1
    return(list(estH = mean(dummys, na.rm = T), sdH = sd(dummys, 
                                                         na.rm = T), estD = mean(dummys, na.rm = T) + 1, sdD = sd(dummys, 
                                                                                                                  na.rm = T)))
  })
  zero_H_alpha <- unlist(lapply(zero_H_alpha_es, function(x) x[[1]]))
  zero_H_alpha_sd <- unlist(lapply(zero_H_alpha_es, function(x) x[[2]]))
  zero_D_alpha <- unlist(lapply(zero_H_alpha_es, function(x) x[[3]]))
  zero_D_alpha_sd <- unlist(lapply(zero_H_alpha_es, function(x) x[[4]]))
  if (probar) {
    setTxtProgressBar(pb, 2)
    cat(" Calculating one_H/D_alpha ...          ")
  }
  one_H_alpha_es <- lapply(pops, function(x) {
    p <- colMeans(as.matrix(x), na.rm = T)/2
    p <- p[!is.na(p)]
    logp <- ifelse(!is.finite(log(p)), 0, log(p))
    log1_p <- ifelse(!is.finite(log(1 - p)), 0, log(1 - p))
    dummys <- -(p * logp + (1 - p) * log1_p)
    return(list(estH = mean(dummys), sdH = sd(dummys), estD = mean(exp(dummys)), 
                sdD = sd(exp(dummys)), dummys = dummys))
  })
  one_H_alpha <- unlist(lapply(one_H_alpha_es, function(x) x[[1]]))
  one_H_alpha_sd <- unlist(lapply(one_H_alpha_es, function(x) x[[2]]))
  one_D_alpha <- unlist(lapply(one_H_alpha_es, function(x) x[[3]]))
  one_D_alpha_sd <- unlist(lapply(one_H_alpha_es, function(x) x[[4]]))
  if (probar) {
    setTxtProgressBar(pb, 3)
    cat(" Calculating two_H/D_alpha ...           ")
  }
  two_H_alpha_es <- lapply(pops, function(x) {
    p <- colMeans(as.matrix(x)==1, na.rm = T)
    p <- p[!is.na(p)]
    dummys <- p
    return(list(estH = mean(dummys), sdH = sd(dummys), estD = mean(1/(1 - 
                                                                        dummys)), sdD = sd(1/(1 - dummys)), dummys = dummys))
  })
  two_H_alpha <- unlist(lapply(two_H_alpha_es, function(x) x[[1]]))
  two_H_alpha_sd <- unlist(lapply(two_H_alpha_es, function(x) x[[2]]))
  two_D_alpha <- unlist(lapply(two_H_alpha_es, function(x) x[[3]]))
  two_D_alpha_sd <- unlist(lapply(two_H_alpha_es, function(x) x[[4]]))
  mat_zero_H_beta <- NA
  mat_one_H_beta <- NA
  mat_two_H_beta <- NA
  npops <- length(pops)
  if (npops > 1) {
    if (probar) {
      setTxtProgressBar(pb, 4)
      cat(" Counting pairwise missing loci...")
    }
    pairs <- t(combn(npops, 2))
    nlocpairpop <- apply(pairs, 1, function(x) {
      pop1 <- pops[[x[1]]]
      pop2 <- pops[[x[2]]]
      pp1 <- colMeans(as.matrix(pop1), na.rm = T)/2
      pp2 <- colMeans(as.matrix(pop2), na.rm = T)/2
      index <- !is.na(pp1) & !is.na(pp2)
      return(sum(index))
    })
    mat_nloc_pops <- matrix(NA, nrow = npops, ncol = npops)
    mat_nloc_pops[lower.tri(mat_nloc_pops)] <- nlocpairpop
    colnames(mat_nloc_pops) <- rownames(mat_nloc_pops) <- names(pops)
    if (probar) {
      setTxtProgressBar(pb, 5)
      cat(" Calculating zero_H/D_beta ...         ")
    }
    zero_H_beta_es <- apply(pairs, 1, function(x) {
      pop1 <- pops[[x[1]]]
      pop2 <- pops[[x[2]]]
      pp1 <- colMeans(as.matrix(pop1), na.rm = T)/2
      pp1 <- ifelse(pp1 > 0 & pp1 < 1, 0.5, pp1)
      pp2 <- colMeans(as.matrix(pop2), na.rm = T)/2
      pp2 <- ifelse(pp2 > 0 & pp2 < 1, 0.5, pp2)
      index <- !is.na(pp1) & !is.na(pp2)
      pp1 <- pp1[index]
      pp2 <- pp2[index]
      dummys <- abs(pp1 - pp2)
      return(list(estH = mean(dummys), sdH = sd(dummys), 
                  estD = mean(dummys) + 1, sdD = sd(dummys)))
    })
    zero_H_beta <- unlist(lapply(zero_H_beta_es, function(x) x[[1]]))
    zero_H_beta_sd <- unlist(lapply(zero_H_beta_es, function(x) x[[2]]))
    zero_D_beta <- unlist(lapply(zero_H_beta_es, function(x) x[[3]]))
    zero_D_beta_sd <- unlist(lapply(zero_H_beta_es, function(x) x[[4]]))
    mat_zero_H_beta <- matrix(NA, nrow = npops, ncol = npops)
    mat_zero_H_beta[lower.tri(mat_zero_H_beta)] <- zero_H_beta
    mat_zero_H_beta[pairs] <- zero_H_beta_sd
    colnames(mat_zero_H_beta) <- rownames(mat_zero_H_beta) <- names(pops)
    mat_zero_D_beta <- matrix(NA, nrow = npops, ncol = npops)
    mat_zero_D_beta[lower.tri(mat_zero_D_beta)] <- zero_D_beta
    mat_zero_D_beta[pairs] <- zero_D_beta_sd
    colnames(mat_zero_D_beta) <- rownames(mat_zero_D_beta) <- names(pops)
    if (probar) {
      setTxtProgressBar(pb, 6)
      cat(" Calculating one_H/D_beta ...    ")
    }
    p <- colMeans(as.matrix(gl), na.rm = T)/2
    i0 <- which(!is.na(p))
    logp <- ifelse(!is.finite(log(p)), 0, log(p))
    log1_p <- ifelse(!is.finite(log(1 - p)), 0, log(1 - p))
    one_H_alpha_all <- -(p * logp + (1 - p) * log1_p)
    one_H_beta_es <- apply(pairs, 1, function(x) {
      i1 <- which(!is.na(colMeans(as.matrix(pops[[x[1]]]), 
                                  na.rm = T)/2))
      i2 <- which(!is.na(colMeans(as.matrix(pops[[x[2]]]), 
                                  na.rm = T)/2))
      tt <- table(c(i0, i1, i2))
      index <- as.numeric(names(tt)[tt == 3])
      dummys <- one_H_alpha_all[i0 %in% index] - (one_H_alpha_es[[x[1]]]$dummys[i1 %in% 
                                                                                  index] + one_H_alpha_es[[x[2]]]$dummys[i2 %in% 
                                                                                                                           index])/2
      return(list(estH = mean(dummys), sdH = sd(dummys), 
                  estD = mean(exp(dummys)), sdD = sd(exp(dummys))))
    })
    one_H_beta <- unlist(lapply(one_H_beta_es, function(x) x[[1]]))
    one_H_beta_sd <- unlist(lapply(one_H_beta_es, function(x) x[[2]]))
    one_D_beta <- unlist(lapply(one_H_beta_es, function(x) x[[3]]))
    one_D_beta_sd <- unlist(lapply(one_H_beta_es, function(x) x[[4]]))
    mat_one_H_beta <- matrix(NA, nrow = npops, ncol = npops)
    mat_one_H_beta[lower.tri(mat_one_H_beta)] <- one_H_beta
    mat_one_H_beta[pairs] <- one_H_beta_sd
    colnames(mat_one_H_beta) <- rownames(mat_one_H_beta) <- names(pops)
    mat_one_D_beta <- matrix(NA, nrow = npops, ncol = npops)
    mat_one_D_beta[lower.tri(mat_one_D_beta)] <- one_D_beta
    mat_one_D_beta[pairs] <- one_D_beta_sd
    colnames(mat_one_D_beta) <- rownames(mat_one_D_beta) <- names(pops)
    if (probar) {
      setTxtProgressBar(pb, 7)
      cat(" Calculating two_H/D_beta...    ")
    }
    p <- colMeans(as.matrix(gl), na.rm = T)/2
    i0 <- which(!is.na(p))
    two_H_alpha_all <- (1 - (p * p + (1 - p) * (1 - p)))
    two_H_beta_es <- apply(pairs, 1, function(x) {
      i1 <- which(!is.na(colMeans(as.matrix(pops[[x[1]]]), 
                                  na.rm = T)/2))
      i2 <- which(!is.na(colMeans(as.matrix(pops[[x[2]]]), 
                                  na.rm = T)/2))
      tt <- table(c(i0, i1, i2))
      index <- as.numeric(names(tt)[tt == 3])
      m2Ha <- (two_H_alpha_es[[x[1]]]$dummys[i1 %in% index] + 
                 two_H_alpha_es[[x[2]]]$dummys[i2 %in% index])/2
      dummys <- ((two_H_alpha_all[i0 %in% index] - m2Ha)/(1 - 
                                                            m2Ha)) * (npops/(npops - 1))
      return(list(estH = mean(dummys), sdH = sd(dummys), 
                  estD = mean(exp(dummys)), sdD = sd(exp(dummys))))
    })
    two_H_beta <- unlist(lapply(two_H_beta_es, function(x) x[[1]]))
    two_H_beta_sd <- unlist(lapply(two_H_beta_es, function(x) x[[2]]))
    two_D_beta <- unlist(lapply(two_H_beta_es, function(x) x[[3]]))
    two_D_beta_sd <- unlist(lapply(two_H_beta_es, function(x) x[[4]]))
    mat_two_H_beta <- matrix(NA, nrow = npops, ncol = npops)
    mat_two_H_beta[lower.tri(mat_two_H_beta)] <- two_H_beta
    mat_two_H_beta[pairs] <- two_H_beta_sd
    colnames(mat_two_H_beta) <- rownames(mat_two_H_beta) <- names(pops)
    mat_two_D_beta <- matrix(NA, nrow = npops, ncol = npops)
    mat_two_D_beta[lower.tri(mat_two_D_beta)] <- two_D_beta
    mat_two_D_beta[pairs] <- two_D_beta_sd
    colnames(mat_two_D_beta) <- rownames(mat_two_D_beta) <- names(pops)
  }
  if (probar) {
    setTxtProgressBar(pb, 8)
    cat(" Done.                               ")
  }
  if (spectrumplot) {
    fs <- cbind(zero_D_alpha, one_D_alpha, two_D_alpha)
    cx <- max(1 - (max(-12 + nrow(fs), 0) * 0.025), 0.5)
    bb <- barplot(fs, beside = T, names.arg = rep(rownames(fs), 
                                                  3), ylim = c(1, 2.15), main = "q-profile", 
                  col = rainbow(npops), las = 2, xpd = FALSE, cex.names = cx)
    text(colMeans(bb), rep(2.1, 3), labels = c("q=0", 
                                               "q=1", "q=2"))
    sds <- cbind(zero_D_alpha_sd, one_D_alpha_sd, two_D_alpha_sd)
    up <- fs + sds
    low <- fs - sds
    if (confiplot) 
      for (i in 1:ncol(bb)) {
        for (ii in 1:nrow(bb)) {
          arrows(x0 = bb[ii, i], y0 = up[ii, i], x1 = bb[ii, 
                                                         i], y1 = low[ii, i], angle = 90, code = 3, 
                 length = 0)
        }
      }
  }
  if (!is.na(match(table, c("H", "DH", "HD")))) {
    tt <- data.frame(nloci = unlist(nlocpop), m_0Ha = zero_H_alpha, 
                     sd_0Ha = zero_H_alpha_sd, m_1Ha = one_H_alpha, sd_1Ha = one_H_alpha_sd, 
                     m_2Ha = two_H_alpha, sd_2Ha = two_H_alpha_sd)
    print(knitr::kable(tt, digits = 3))
    if (npops > 1) {
      cat("\n\npairwise non-missing loci")
      print(knitr::kable(mat_nloc_pops, digits = 3))
      cat("\n\n0_H_beta")
      print(knitr::kable(mat_zero_H_beta, digits = 3), 
      )
      cat("\n\n1_H_beta")
      print(knitr::kable(mat_one_H_beta, digits = 3))
      cat("\n\n2_H_beta")
      print(knitr::kable(mat_two_H_beta, digits = 3))
    }
  }
  if (!is.na(match(table, c("D", "DH", "HD")))) {
    tt <- data.frame(nloci = unlist(nlocpop), m_0Da = zero_D_alpha, 
                     sd_0Da = zero_D_alpha_sd, m_1Da = one_D_alpha, sd_1Da = one_D_alpha_sd, 
                     m_2Da = two_D_alpha, sd_2Da = two_D_alpha_sd)
    print(knitr::kable(tt, digits = 3))
    if (npops > 1) {
      cat("\n\npairwise non-missing loci")
      print(knitr::kable(mat_nloc_pops, digits = 3))
      cat("\n\n0_D_beta")
      print(knitr::kable(mat_zero_D_beta, digits = 3))
      cat("\n\n1_D_beta")
      print(knitr::kable(mat_one_D_beta, digits = 3))
      cat("\n\n2_D_beta")
      print(knitr::kable(mat_two_D_beta, digits = 3))
    }
  }
  if (npops > 1) 
    out <- list(nlocpop = unlist(nlocpop), nlocpairpop = mat_nloc_pops, 
                zero_H_alpha = zero_H_alpha, zero_H_alpha_sd = zero_H_alpha_sd, 
                one_H_alpha = one_H_alpha, one_H_alpha_sd = one_H_alpha_sd, 
                two_H_alpha = two_H_alpha, two_H_alpha_sd = two_H_alpha_sd, 
                zero_D_alpha = zero_D_alpha, zero_D_alpha_sd = zero_D_alpha_sd, 
                one_D_alpha = one_D_alpha, one_D_alpha_sd = one_D_alpha_sd, 
                two_D_alpha = two_D_alpha, two_D_alpha_sd = two_D_alpha_sd, 
                zero_H_beta = mat_zero_H_beta, one_H_beta = mat_one_H_beta, 
                two_H_beta = mat_two_H_beta, zero_D_beta = mat_zero_D_beta, 
                one_D_beta = mat_one_D_beta, two_D_beta = mat_two_D_beta)
  else out <- list(nlocpop = unlist(nlocpop), zero_H_alpha = zero_H_alpha, 
                   zero_H_alpha_sd = zero_H_alpha_sd, one_H_alpha = one_H_alpha, 
                   one_H_alpha_sd = one_H_alpha_sd, two_H_alpha = two_H_alpha, 
                   two_H_alpha_sd = two_H_alpha_sd, zero_D_alpha = zero_D_alpha, 
                   zero_D_alpha_sd = zero_D_alpha_sd, one_D_alpha = one_D_alpha, 
                   one_D_alpha_sd = one_D_alpha_sd, two_D_alpha = two_D_alpha, 
                   two_D_alpha_sd = two_D_alpha_sd)
  return(out)
}
