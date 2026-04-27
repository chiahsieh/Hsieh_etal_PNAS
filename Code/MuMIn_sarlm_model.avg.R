
#trace(model.avg.defult, edit = T)
`model.avg.defult` <-
  function(object, ..., beta = c("none", "sd", "partial.sd"),
           rank = NULL, rank.args = NULL, revised.var = TRUE,
           dispersion = NULL, ct.args = NULL) {
    
    if (is.object(object)) {
      models <- list(object, ...)
      rank <- .getRank(rank, rank.args = rank.args, object = object) 
    } else {
      if(length(object) == 0L) stop("'object' is an empty list")
      models <- object
      object <- object[[1L]]
      if (!is.null(rank) || is.null(rank <- attr(models, "rank"))) {
        rank <- .getRank(rank, rank.args = rank.args, object = object)
      }
    }
    
    strbeta <- betaMode <- NULL
    eval(.expr_beta_arg)
    
    nModels <- length(models)
    if(nModels == 1L) stop("only one model supplied. Nothing to do")
    checkIsModelDataIdentical(models)
    
    testSmoothKConsistency(models) # for gam, if any
    
    ICname <- asChar(attr(rank, "call")[[1L]])
    
    allterms1 <- lapply(models, getAllTerms)
    all.terms <- unique(unlist(allterms1, use.names = FALSE))
    
    # sort by level (main effects first)
    all.terms <- all.terms[order(vapply(gregexpr(":", all.terms),
                                        function(x) if(x[1L] == -1L) 0L else length(x), 1L), all.terms)]
    
    # allmodelnames <- modelNames(models, asNumeric = FALSE,
    # withRandomTerms = FALSE, withFamily = FALSE)
    allmodelnames <- .modelNames(allTerms = allterms1, uqTerms = all.terms)
    
    #if(is.null(names(models))) names(models) <- allmodelnames
    
    coefTableCall <- as.call(c(alist(coefTable, models[[i]],
                                     dispersion = dispersion[i]), ct.args))
    if(is.null(dispersion)) coefTableCall$dispersion <- NULL
    
    if(betaMode == 2L) {
      coefTableCall[[1L]] <- as.name("std.coef")
      coefTableCall[['partial.sd']] <- TRUE
    }
    
    .DebugPrint(coefTableCall)
    
    # check if models are unique:
    if(!is.null(dispersion)) dispersion <- rep(dispersion, length.out = nModels)
    coefTables <- vector(nModels, mode = "list")
    for(i in seq_len(nModels))
      coefTables[[i]] <-  eval(coefTableCall)
    
    mcoeffs <- lapply(coefTables, "[", , 1L)
    dup <- unique(sapply(mcoeffs, function(i) which(sapply(mcoeffs, identical, i))))
    dup <- dup[sapply(dup, length) > 1L]
    if (length(dup) > 0L) stop("models are not unique. Duplicates: ",
                               prettyEnumStr(sapply(dup, paste0, collapse = " = "),
                                             quote = "'"))
    
    LL <- .getLik(object)
    logLik <- LL$logLik
    lLName <- LL$name
    
    ic <- vapply(models, rank, 0)
    logLiks <- lapply(models, logLik)
    delta <- ic - min(ic)
    weight <- exp(-delta / 2) / sum(exp(-delta / 2))
    model.order <- order(weight, decreasing = TRUE)
    
    # ----!!! From now on, everything MUST BE ORDERED by 'weight' !!!-----------
    mstab <- cbind(df = vapply(logLiks, attr, 0, "df"),
                   logLik = as.numeric(logLiks), IC = ic, delta = delta, weight = weight,
                   deparse.level = 0L)
    if(!is.null(dispersion)) mstab <- cbind(mstab, Dispersion = dispersion)
    rownames(mstab) <- allmodelnames
    mstab <- mstab[model.order, ]
    weight <- mstab[, "weight"] # has been sorted in table
    models <- models[model.order]
    coefTables <- coefTables[model.order]
    
    if (betaMode == 1L) {
      response.sd <- sd(models[[i]]$y) #change for sarlm
      #response.sd <- sd(model.response(model.frame(object))).     #original code 
      for(i in seq_along(coefTables)) {
        X <- models[[i]]$X #change for sarlm
        #X <- model.matrix(models[[i]])   #original code 
        coefTables[[i]][, 1L:2L] <-
          coefTables[[i]][, 1L:2L] *
          #apply(model.matrix(models[[i]]), 2L, sd) / response.sd.   #original unused code 
          apply(X[, match(rownames(coefTables[[i]]), colnames(X)),
                  drop = FALSE], 2L, sd) / response.sd
      }
    }
    
    cfarr <- coefArray(coefTables)
    cfmat <- array(cfarr[, 1L, ], dim = dim(cfarr)[-2L], dimnames = dimnames(cfarr)[-2L])
    cfmat[is.na(cfmat)]<- 0
    coefMat <- array(NA_real_, dim = c(2L, ncol(cfmat)),
                     dimnames = list(c("full", "subset"), colnames(cfmat)))
    coefMat[1L, ] <- drop(weight %*% cfmat)
    coefMat[2L, ] <- coefMat[1L, ] / colSums(array(weight *
                                                     as.numeric(!is.na(cfarr[, 1L, ])), dim = dim(cfmat)))
    coefMat[is.nan(coefMat)] <- NA_real_
    
    names(all.terms) <- seq_along(all.terms)
    colnames(mstab)[3L] <- ICname
    
    # Benchmark: 3.7x faster
    #system.time(for(i in 1:10000) t(array(unlist(p), dim=c(length(all.terms),length(models)))))
    #system.time(for(i in 1:10000) do.call("rbind", p))
    
    vpresent <- do.call("rbind", lapply(models, function(x)
      all.terms %in% getAllTerms(x)))
    
    if(all(dim(vpresent) > 0L)) {
      sw <- apply(weight * vpresent, 2L, sum)
      names(sw) <- all.terms
      o <- order(sw, decreasing = TRUE)
      sw <- sw[o]
      attr(sw, "n.models") <- structure(colSums(vpresent)[o], names = all.terms)
      class(sw) <- c("sw", "numeric")
    } else {
      sw <- structure(integer(0L), n.models = integer(0L), class = c("sw", "numeric"))
    }
    
    mmxs <- tryCatch(cbindDataFrameList(lapply(models, function(x) x$X)), #change for sarlm
                     error = return_null, warning = return_null)
    #mmxs <- tryCatch(cbindDataFrameList(lapply(models, model.matrix)),  #origical code 
    #                 error = return_null, warning = return_null)
    
    # Far less efficient:
    #mmxs <- lapply(models, model.matrix)
    #mx <- mmxs[[1]];
    #for (i in mmxs[-1])
    #	mx <- cbind(mx, i[,!(colnames(i) %in% colnames(mx)), drop=FALSE])
    
    # residuals averaged (with brute force)
    #rsd <- tryCatch(apply(vapply(models, residuals, residuals(object)), 1L,
    #weighted.mean, w = weight), error = return_null)
    #rsd <- NULL
    ## XXX: how to calc residuals ?
    
    modelClasses <- lapply(models, class)
    frm <-
      if(all(vapply(modelClasses[-1L], identical, FALSE, modelClasses[[1L]]))) {
        trm <- tryCatch(terms(models[[1L]]),
                        error = function(e) terms(formula(models[[1L]])))
        response <- attr(trm, "response")
        m1 <- models[[1L]]
        makeArgs(m1, all.terms, opt = list(
          response = if(response > 0L) attr(trm, "variables")[[response + 1L]] else NULL,
          gmFormulaEnv = environment(formula(m1)),
          intercept = ! identical(unique(unlist(lapply(allterms1, attr, "intercept"))), 0),
          interceptLabel = unique(unlist(lapply(allterms1, attr, "interceptLabel"))),
          #	random = attr(allTerms0, "random"),
          gmCall = get_call(m1),
          gmEnv = parent.frame(),
          allTerms = all.terms,
          random = . ~ .
        ))[[1L]]
      } else NA
    
    #.Debug(.Generic <- "model.avg") #original code
    .Generic <- "model.avg"
    
    ret <- list(
      msTable = structure(as.data.frame(mstab, stringsAsFactors = TRUE),
                          term.codes = attr(allmodelnames, "variables")),
      coefficients = coefMat,
      coefArray = cfarr,
      sw = sw,
      x = mmxs,
      residuals = NULL, # no residuals, as they can be calculated in several ways
      formula = frm,
      call = {
        cl <- match.call()
        cl[[1L]] <- as.name(.Generic)
        cl
      }
    )
    
    attr(ret, "rank") <- rank
    attr(ret, "modelList") <- models
    attr(ret, "beta") <- strbeta
    attr(ret, "nobs") <- nobs(object)
    attr(ret, "revised.var") <- revised.var
    class(ret) <- "averaging"
    return(ret)
  }


##============ Not run ========================
## Other functions to help testing
# testSmoothKConsistency ------------------------------------------------------------------------------------------
testSmoothKConsistency <-
  function(models) {
    
    # method 2: guess from coefficient names:
    if(inherits(models, "model.selection")) {
      
      x <- lapply(attr(models, "coefTables"), rownames)
      # XXX: add label 'te(x1,x2)'
      
      res <- unlist(unname(lapply(x, function(x) {
        s <- grep("^(s|t[ei2])\\(.+\\)\\.\\d+$", x, perl = TRUE)
        if(length(s) != 0L) {
          m <- regexpr("^(?:s|t[ei2])\\((.+)\\)\\.\\d+$",  x[s], perl = TRUE)
          cst <- attr(m, "capture.start")[, 1L]
          y <- substring(x[s], cst, cst + attr(m, "capture.length")[, 1L] - 1L)
          tapply(y, y, length)
        } else NULL
      })), recursive = FALSE)
      names(res) <- sapply(lapply(expr.split(names(res), ","), sort),
                           paste0, collapse = ",")
    } else {
      # use information stored in gam objects:
      .getSmoothK <-
        function(x) {
          if(inherits(x, "gamm") || (is.list(x) &&
                                     (length(x) >= 2L) && identical(names(x)[2L], "gam"))) {
            x <- x[[2L]]
          } else if(!inherits(x, "gam")) return(NULL)
          n <- length(x$smooth)
          rval <- vector("list", n)
          nmv <- character(n)
          for(i in seq_len(n)) {
            y <- x$smooth[[i]]
            if(is.null(y$margin)) {
              nmv[i] <- y$term
              rval[[i]] <- y$bs.dim
            } else {
              nm1 <- vapply(y$margin, `[[`, "", "term")
              o <- order(nm1)
              nmv[i] <- paste0(nm1[o], collapse = ",")
              rval[[i]] <- sapply(y$margin, `[[`, "bs.dim")[o]
            }
            nmv[i] <- paste0(sub("\\(.*", "", y$label), "(", nmv[i], ")")
          }
          names(rval) <- nmv
          rval
        }
      res <- unlist(unname(lapply(models, .getSmoothK)), recursive = FALSE)
    }
    
    if(!is.null(res)) {
      res <- vapply(split(res, names(res)), function(x) {
        k1 <- x[[1L]]
        for(i in 1L:length(x)) if(!identical(x[[i]], k1)) return(TRUE)
        return(FALSE)
      }, FALSE)
      
      
      if(any(res))
        warning("smooth term dimensions differ between models for variables ",
                prettyEnumStr(names(res)[res], quote = "'"),
                ". Related coefficients are incomparable."
        )
    }
    invisible()
  }


# .modelNames ---------------------------------------------------------------------------------
.modelNames <-
  function(models = NULL, allTerms, uqTerms, use.letters = FALSE, ...) {
    if(missing(allTerms)) allTerms <- lapply(models, getAllTerms)
    if(missing(uqTerms) || is.null(uqTerms))
      uqTerms <- unique(unlist(allTerms, use.names = FALSE))
    
    n <- length(uqTerms)
    
    if(use.letters && n > length(LETTERS)) 
      stop("more terms than there are letters")
    sep <- if(!use.letters && n > 9L) "+" else ""
    
    labels <- if (use.letters) LETTERS[seq_len(n)] else as.character(seq_len(n))
    ret <- sapply(allTerms, function(x) paste(labels[sort(match(x, uqTerms))],
                                              collapse = sep))
    
    dup <- table(ret)
    dup <- dup[dup > 1L]
    
    if(length(dup) > 0L) {
      idup <- which(ret %in% names(dup))
      ret[idup] <- sapply(idup, function(i) paste0(ret[i],
                                                   letters[sum(ret[seq.int(i)] == ret[i])]))
    }
    ret[ret == ""] <- "(Null)"
    attr(ret, "variables") <- structure(seq_along(uqTerms), names = uqTerms)
    ret
  }


#coefArray --------------------------------------------------------------------------
`coefArray` <- 
  function(object) {
    coefNames <- fixCoefNames(unique(unlist(lapply(object, rownames),
                                            use.names = FALSE)))
    nCoef <- length(coefNames)
    nModels <- length(object)
    rval <- array(NA_real_, dim = c(nModels, 3L, nCoef),
                  dimnames = list(names(object), c("Estimate", "Std. Error", "df"), coefNames))
    for(i in seq_along(object)) {
      z <- object[[i]]
      rval[i, seq_len(ncol(z)), ] <- t(z[match(coefNames, fixCoefNames(rownames(z))), ])
    }
    rval
  }

#makeArgs--------------------------------------------
makeArgs <- 
  function(obj, termNames, opt, ...) {
    reportProblems <- character(0L)
    termNames[termNames %in% opt$interceptLabel] <- "1"
    ## XXX: what if length(opt$intercept) > 1 ???
    f <- reformulate(c(if(!opt$intercept) "0" else if (!length(termNames)) "1", termNames), response = opt$response)
    
    environment(f) <- opt$gmFormulaEnv
    ret <- list(formula = f)
    if(!is.null(opt$gmCall$start)) {
      coefNames <- fixCoefNames(.getCoefNames(f, opt$gmDataHead,
                                              opt$gmCall$contrasts, envir = opt$gmEnv))
      idx <- match(coefNames, opt$gmCoefNames)
      if(anyNA(idx)) reportProblems <-
        append(reportProblems, "cannot subset 'start' argument. Coefficients in the model do not exist in 'global.model'")
      else ret$start <- substitute(start[idx], list(start = opt$gmCall$start,
                                                    idx = idx))
    }
    #attr(ret, "formulaList") <- list(f)
    attr(ret, "problems") <- reportProblems
    ret
  }


# .Debug----------------------------------------
`.Debug` <- 
  function(expr) {
    if(isTRUE(getOption("debug.MuMIn"))) {
      eval.parent(substitute(expr))
    }
  }

