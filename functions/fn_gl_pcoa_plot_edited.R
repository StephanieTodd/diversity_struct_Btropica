# glPca = shell.pca
# x = shell.gl
# pop.labels ="ind"
# xaxis=1
# yaxis=2 
# scale = T
# ellipse = F
# save2tmp = F
# plevel = 0.95
# interactive = FALSE
# as.pop = NULL
# hadjust = 1.5
# vadjust = 1
# xaxis = 1
# yaxis = 2
# zaxis = NULL
# pt.size = 2
# pt.colors = NULL
# pt.shapes = NULL
# label.size = 1
# axis.label.size = 1.5
# save2tmp = FALSE
# verbose = NULL

gl.pcoa.plot_edited <- function (glPca, x, scale = FALSE, ellipse = NULL, plevel = 0.95, pop.labels = "pop", 
                                 interactive = FALSE, as.pop = NULL, hadjust = 1.5, vadjust = 1, xaxis = 1, yaxis = 2, zaxis = NULL, pt.size = 2, 
                                 pt.colors = NULL, pt.shapes = NULL, label.size = 1, axis.label.size = 1.5, save2tmp = FALSE, verbose = NULL) {
  hold_x <- x
  hold_glPca <- glPca
  verbose <- 0 # gl.check.verbosity(verbose)
  funname <- match.call()[[1]]
 # utils.flag.start(func = funname, build = "Josh", verbosity = verbose)
  datatype1 <- utils.check.datatype(glPca, accept = c("glPca", 
                                                      "list"), verbose = verbose)
  datatype2 <- utils.check.datatype(x, accept = c("SNP", "SilicoDArT", 
                                                  "fd", "dist", "list"), verbose = verbose)
  if (interactive | !is.null(zaxis)) {
    pkg <- "plotly"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      stop(error("Package ", pkg, " needed for this function to work. Please install it."))
    }
  }
  if (datatype1 == "list") {
    pkg <- "gganimate"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      stop(error("Package ", pkg, " needed for this function to work. Please install it."))
    }
    pkg <- "tibble"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      stop(error("Package ", pkg, " needed for this function to work. Please install it."))
    }
    x <- x[[1]]
    glPca <- glPca[[1]]
  }
  if (pop.labels != "none" && pop.labels != "ind" && pop.labels != 
      "pop" && pop.labels != "legend") {
    cat(warn("  Warning: Parameter 'pop.labels' must be one of none|ind|pop|legend, set to 'pop'\n"))
    pop.labels <- "pop"
  }
  if (plevel < 0 | plevel > 1) {
    cat(warn("  Warning: Parameter 'plevel' must fall between 0 and 1, set to 0.95\n"))
    plevel <- 0.95
  }
  if (hadjust < 0 | hadjust > 3) {
    cat(warn("  Warning: Parameter 'hadjust' must fall between 0 and 3, set to 1.5\n"))
    hadjust <- 1.5
  }
  if (vadjust < 0 | hadjust > 3) {
    cat(warn("  Warning: Parameter 'vadjust' must fall between 0 and 3, set to 1.5\n"))
    vadjust <- 1.5
  }
  if (xaxis < 1 | xaxis > ncol(glPca$scores)) {
    cat(warn("  Warning: X-axis must be specified to lie between 1 and the number of retained dimensions of the ordination", 
             ncol(glPca$scores), "; set to 1\n"))
    xaxis <- 1
  }
  if (yaxis < 1 | yaxis > ncol(glPca$scores)) {
    cat(warn("  Warning: Y-axis must be specified to lie between 1 and the number of retained dimensions of the ordination", 
             ncol(glPca$scores), "; set to 2\n"))
    yaxis <- 2
  }
  if (!is.null(zaxis)) {
    if (zaxis < 1 | zaxis > ncol(glPca$scores)) {
      cat(warn("  Warning: Z-axis must be specified to lie between 1 and the number of retained dimensions of the ordination", 
               ncol(glPca$scores), "; set to 3\n"))
      zaxis <- 3
    }
  }
  if (datatype1 %in% c("SNP", "SilicoDArT")) {
    pop.hold <- pop(x)
    if (!is.null(as.pop)) {
      if (as.pop %in% names(x@other$ind.metrics)) {
        pop(x) <- as.matrix(x@other$ind.metrics[as.pop])
        if (verbose >= 2) {
           cat(report("  Temporarily setting population assignments to", 
                     as.pop, "as specified by the as.pop parameter\n"))
        }
      }
      else {
        stop(error("Fatal Error: individual metric assigned to 'pop' does not exist. Check names(gl@other$loc.metrics) and select again\n"))
      }
    }
  }
  if (datatype2 == "fd") {
    x <- x$fd
    datatype2 <- utils.check.datatype(x, verbose = 0)
  }
  axis.label.size <- axis.label.size * 10
  gen <- NULL
  if (datatype1 == "list") {
    gen_number <- length(hold_x)
    df_sim <- as.data.frame(matrix(ncol = 5))
    colnames(df_sim) <- c("PCoAx", "PCoAy", "ind", "pop", 
                          "gen")
    test_pos_neg <- as.data.frame(matrix(nrow = gen_number, 
                                         ncol = 3))
    colnames(test_pos_neg) <- c("gen", "test_x", "test_y")
    ind_x_axis <- which.max(abs(hold_glPca[[1]]$scores[, 
                                                       xaxis]))
    ind_y_axis <- which.max(abs(hold_glPca[[1]]$scores[, 
                                                       yaxis]))
    test_pos_neg[1, "test_x"] <- if (hold_glPca[[1]]$scores[ind_x_axis, 
                                                            xaxis] >= 0) 
      "positive"
    else "negative"
    test_pos_neg[1, "test_y"] <- if (hold_glPca[[1]]$scores[ind_y_axis, 
                                                            yaxis] >= 0) 
      "positive"
    else "negative"
    for (sim_i in 1:gen_number) {
      glPca <- hold_glPca[[sim_i]]
      x <- hold_x[[sim_i]]
      m <- cbind(glPca$scores[, xaxis], glPca$scores[, 
                                                     yaxis])
      df <- data.frame(m)
      xlab <- paste("PCA Axis", xaxis)
      ylab <- paste("PCA Axis", yaxis)
      ind <- indNames(x)
      pop <- factor(pop(x))
      gen <- unique(x$other$sim.vars$generation)
      df <- cbind(df, ind, pop, unique(x$other$sim.vars$generation))
      colnames(df) <- c("PCoAx", "PCoAy", "ind", "pop", 
                        "gen")
      test_pos_neg[sim_i, "test_x"] <- if (hold_glPca[[sim_i]]$scores[ind_x_axis, 
                                                                      xaxis] >= 0) 
        "positive"
      else "negative"
      test_pos_neg[sim_i, "test_y"] <- if (hold_glPca[[sim_i]]$scores[ind_y_axis, 
                                                                      yaxis] >= 0) 
        "positive"
      else "negative"
      if (test_pos_neg[1, "test_x"] != test_pos_neg[sim_i, 
                                                    "test_x"]) {
        df$PCoAx <- df$PCoAx * -1
      }
      if (test_pos_neg[1, "test_y"] != test_pos_neg[sim_i, 
                                                    "test_y"]) {
        df$PCoAy <- df$PCoAy * -1
      }
      df_sim <- rbind(df_sim, df)
    }
    df_sim <- tibble::as_tibble(df_sim)
    df_sim <- df_sim[-1, ]
    p <- ggplot(df_sim, aes(PCoAx, PCoAy, colour = pop)) + 
      geom_point(size = 3) + labs(title = "Generation: {frame_time}", 
                                  x = xlab, y = ylab) + gganimate::transition_time(gen) + 
      gganimate::ease_aes("linear")
    return(p)
  }
  PCoAx <- PCoAy <- NULL
  if (is.null(zaxis)) {
    m <- cbind(glPca$scores[, xaxis], glPca$scores[, yaxis])
  } else {
    m <- cbind(glPca$scores[, xaxis], glPca$scores[, yaxis], 
               glPca$scores[, zaxis])
  }
  df <- data.frame(m)
  s <- sum(glPca$eig[glPca$eig >= 0])
  e <- round(glPca$eig * 100/s, 1)
  if (datatype2 == "SNP" | datatype2 == "SilicoDArT") {
    xlab <- paste("PCA Axis", xaxis, "(", e[xaxis], "%)")
    ylab <- paste("PCA Axis", yaxis, "(", e[yaxis], "%)")
    if (!is.null(zaxis)) {
      zlab <- paste("PCA Axis", zaxis, "(", e[zaxis], "%)")
    }
    ind <- indNames(x)
    pop <- factor(pop(x))
    df <- cbind(df, ind, pop)
    if (is.null(zaxis)) {
      colnames(df) <- c("PCoAx", "PCoAy", "ind", "pop")
    }
    else {
      colnames(df) <- c("PCoAx", "PCoAy", "PCoAz", "ind", 
                        "pop")
    }
  } else {
    xlab <- paste("PCoA Axis", xaxis, "(", e[xaxis], "%)")
    ylab <- paste("PCoA Axis", yaxis, "(", e[yaxis], "%)")
    if (!is.null(zaxis)) {
      zlab <- paste("PCA Axis", zaxis, "(", e[zaxis], "%)")
    }
    ind <- rownames(as.matrix(x))
    pop <- ind
    df <- cbind(df, ind, pop)
    if (is.null(zaxis)) {
      colnames(df) <- c("PCoAx", "PCoAy", "ind", "pop")
    } else {
      colnames(df) <- c("PCoAx", "PCoAy", "PCoAz", "ind", 
                        "pop")
    }
    if (interactive) {
      cat(warn("  Sorry, interactive labels are not available for an ordination generated from a Distance Matrix\n"))
      cat(warn("  Labelling the plot with names taken from the Distance Matrix\n"))
    }
    pop.labels <- "pop"
  }
  if (is.null(zaxis)) {
    if (pop.labels == "pop") {
      if (datatype2 == "SNP") {
        if (verbose >= 2) 
          cat(report("  Plotting populations in a space defined by the SNPs\n"))
      }
      else if (datatype2 == "SilicoDArT") {
        if (verbose >= 2) 
          cat(report("  Plotting populations in a space defined by the presence/absence data\n"))
      }
      else {
        if (verbose >= 2) 
         cat(report("  Plotting entities from the Distance Matrix\n"))
      }
        plott <- ggplot(df, aes(x = PCoAx, y = PCoAy, 
                                group = pop, color = pop))
      }
      else {
        plott <- ggplot(df, aes(x = PCoAx, y = PCoAy, 
                                group = pop, color = pop, shape = pop))
      }
      plott <- plott + geom_point(size = pt.size, aes(color = pop)) + 
        directlabels::geom_dl(aes(label = pop), method = list("smart.grid", 
                                                              cex = label.size)) + theme(axis.title = element_text(face = "bold.italic", 
                                                                                                                   size = axis.label.size, color = "black"), axis.text.x = element_text(face = "bold", 
                                                                                                                                                                                        angle = 0, vjust = 0.5, size = axis.label.size), 
                                                                                         axis.text.y = element_text(face = "bold", angle = 0, 
                                                                                                                    vjust = 0.5, size = axis.label.size)) + labs(x = xlab, 
                                                                                                                                                                 y = ylab)
      if (!is.null(pt.shapes)) {
        plott <- plott + scale_shape_manual(values = pt.shapes)
      }
      if (!is.null(pt.colors)) {
        plott <- plott + scale_color_manual(values = pt.colors)
      }
      plott <- plott + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
        theme(legend.position = "none")
      if (scale == TRUE) {
        plott <- plott + coord_fixed(ratio = 1)
      }
      if (ellipse == TRUE) {
        plott <- plott + stat_ellipse(type = "norm", 
                                      level = plevel)
      }
    }
    if (interactive) {
     cat(report("  Displaying an interactive plot\n"))
      cat(warn("  NOTE: Returning the ordination scores, not a ggplot2 compatable object\n"))
      plott <- ggplot(df, aes(x = PCoAx, y = PCoAy, label = ind)) + 
        geom_point(size = pt.size, aes(color = pop)) + 
        theme(axis.title = element_text(face = "bold.italic", 
                                        size = axis.label.size, color = "black"), axis.text.x = element_text(face = "bold", 
                                                                                                             angle = 0, vjust = 0.5, size = axis.label.size), 
              axis.text.y = element_text(face = "bold", angle = 0, 
                                         vjust = 0.5, size = axis.label.size), legend.title = element_text(color = "black", 
                                                                                                           size = axis.label.size, face = "bold"), legend.text = element_text(color = "black", 
                                                                                                                                                                              size = axis.label.size, face = "bold")) + 
        labs(x = xlab, y = ylab) + geom_hline(yintercept = 0) + 
        geom_vline(xintercept = 0) + theme(legend.position = "none")
      if (scale == TRUE) {
        plott <- plott + coord_fixed(ratio = 1)
      }
      if (ellipse == TRUE) {
        plott <- plott + stat_ellipse(aes(color = pop), 
                                      type = "norm", level = plevel)
      }
      cat(warn("  Ignore any warning on the number of shape categories\n"))
    }
    if (pop.labels == "legend") {
      if (verbose >= 2) 
       # cat(report("  Plotting populations identified by a legend\n"))
      Population <- pop
      if (is.null(pt.shapes)) {
        plott <- ggplot(df, aes(x = PCoAx, y = PCoAy, 
                                group = Population, color = Population))
      }
      else {
        plott <- ggplot(df, aes(x = PCoAx, y = PCoAy, 
                                group = pop, color = Population, shape = Population))
      }
      plott <- plott + geom_point(size = pt.size, aes(color = pop)) + 
        theme(axis.title = element_text(face = "bold.italic", 
                                        size = axis.label.size, color = "black"), axis.text.x = element_text(face = "bold", 
                                                                                                             angle = 0, vjust = 0.5, size = axis.label.size), 
              axis.text.y = element_text(face = "bold", angle = 0, 
                                         vjust = 0.5, size = axis.label.size), legend.title = element_text(color = "black", 
                                                                                                           size = axis.label.size, face = "bold"), legend.text = element_text(color = "black", 
                                                                                                                                                                              size = axis.label.size, face = "bold")) + 
        labs(x = xlab, y = ylab)
      if (!is.null(pt.shapes)) {
        plott <- plott + scale_shape_manual(values = pt.shapes)
      }
      if (!is.null(pt.colors)) {
        plott <- plott + scale_color_manual(values = pt.colors)
      }
      plott <- plott + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
      if (scale == TRUE) {
        plott <- plott + coord_fixed(ratio = 1)
      }
      if (ellipse == TRUE) {
        plott <- plott + stat_ellipse(type = "norm", 
                                      level = plevel)
      }
    }
    if (pop.labels == "none" | pop.labels == FALSE) {
      if (verbose >= 0) 
       # cat(report("  Plotting points with no labels\n"))
      if (is.null(pt.shapes)) {
        plott <- ggplot(df, aes(x = PCoAx, y = PCoAy, 
                                color = pop))
      }
      else {
        plott <- ggplot(df, aes(x = PCoAx, y = PCoAy, 
                                color = pop, shape = pop))
      }
      plott <- plott + geom_point(size = pt.size, aes(color = pop)) + 
        theme(axis.title = element_text(face = "bold.italic", 
                                        size = axis.label.size, color = "black"), axis.text.x = element_text(face = "bold", 
                                                                                                             angle = 0, vjust = 0.5, size = axis.label.size), 
              axis.text.y = element_text(face = "bold", angle = 0, 
                                         vjust = 0.5, size = axis.label.size)) + labs(x = xlab, 
                                                                                      y = ylab)
      if (!is.null(pt.shapes)) {
        plott <- plott + scale_shape_manual(values = pt.shapes)
      }
      if (!is.null(pt.colors)) {
        plott <- plott + scale_color_manual(values = pt.colors)
      }
      plott <- plott + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
        theme(legend.position = "none")
      if (scale == TRUE) {
        plott <- plott + coord_fixed(ratio = 1)
      }
      if (ellipse == TRUE) {
        plott <- plott + stat_ellipse(type = "norm", 
                                      level = plevel)
      }
    }
    if (verbose >= 2) {
     # cat(report("  Preparing plot .... please wait\n"))
    }
    if (interactive) {
      plott <- plotly::ggplotly(plott)
      show(plott)
    }
    else {
      show(plott)
    }
  if (!is.null(zaxis)) {
    if (verbose >= 2) {
     cat(report("  Displaying a three dimensional plot, mouse over for details for each point\n"))
    }
    plott <- plotly::plot_ly(df, x = ~PCoAx, y = ~PCoAy, 
                             z = ~PCoAz, marker = list(size = pt.size * 2), colors = pt.colors, 
                             text = ind) %>% plotly::add_markers(color = ~pop) %>% 
      plotly::layout(legend = list(title = list(text = "Populations")), 
                     scene = list(xaxis = list(title = xlab, titlefont = list(size = axis.label.size/2)), 
                                  yaxis = list(title = ylab, titlefont = list(size = axis.label.size/2)), 
                                  zaxis = list(title = zlab, titlefont = list(size = axis.label.size/2))))
    show(plott)
    if (verbose >= 2) {
      cat(warn("  May need to zoom out to place 3D plot within bounds\n"))
    }
  }
  if (save2tmp) {
    temp_plot <- tempfile(pattern = "Plot_")
    match_call <- paste0(names(match.call()), "_", as.character(match.call()), 
                         collapse = "_")
    saveRDS(list(match_call, plott), file = temp_plot)
    if (verbose >= 2) {
      cat(report("  Saving the ggplot to the session tempfile\n"))
    }
    temp_table <- tempfile(pattern = "Table_")
    saveRDS(list(match_call, df), file = temp_table)
    if (verbose >= 2) {
      cat(report("  Saving tabulation to the session tempfile\n"))
    }
  }
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  invisible(NULL)
  return(plott)  # added this
}
