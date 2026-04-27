## To let MuMIn can successfully run dredge and model.avg function with the format of sarlm, 
## the following functions should be override by the costumed codes for current working environment 
## every time you load the package. 
## PS. Becuase the overided functions only work for sarlm format, reloading the pakage/restart 
##     R to reset the default functions in the pakage
## 
## 2023.02.06
##
## I. Overrided functions: Using  trace(function.name, edit = T) to override
##   1. nobs   <- measure how many response variables in the model to generate tables
##   2. coeffs <- remove lambda/tho estimate to generate only predictors' table 
##               like liner regression models (i.e., (intercept), covariates, ....)
##   3. coefTable <- generate coefficient table. Overide function remove lambda/tho estimate
##                  to generate only predictors' table like liner regression models.
##                  It has a nested function *.makeCoefTable*, might need to define
##                  the function before overriding the coefTable
##   4. std.coef <- generate a two-column coefficient table of Estimates and Std.Error.  
##
## 
## II. Replaced codes in dredge() and model.avg()
##   1. model.matrix(model) <- use to extract response variable values in regression model 
##                             this can be simply extract by replacing with *model$X*
##                             Because this function is basic and broadly used by regression functions, don't override
##                             Line 93, 135 in model.avg 
##   2. response.sd <- sd(model.response(model.frame(x)))  
##                         <- model.response() and model.frame() don't work for sarlm, replacing by 
##                           *response.sd <- sd(x$y)*
##                             Line 91 in model.avg
## 
##  

#I. Overrided functions ==============================================================================================================
## 1. nobs/methods-nobs.R ---------------------------------------------------------------------------------------------
#### https://rdrr.io/cran/MuMIn/src/R/methods-nobs.R
#### Used to extract number of response to generate tables 
trace(nobs, edit = T)
`nobs` <-
  function(object, ...) NROW(fitted(object))
#Manual setup code and replace it in dredge/model.avg, doesn't work 
#`nobs.sarlm` <-
#  function(object, ...) NROW(fitted(object))

## 2. coeffs ----------------------------------------------------------------------------------------------------------
trace(coeffs, edit = T)
`coeffs` <-
  function (model) coef(model)[-1]
#Manual setup code and replace it in dredge/model.avg, doesn't work 
#`coeffs.sarlm` <- function (model) coef(model)[-1]

## 3. coefTable ----------------------------------------------------------------------------------------------------------
#### Used to extract model's output as two-column coefficient table: Est. & Std. Error
#### Extension model format here https://rdrr.io/cran/MuMIn/src/R/coefTable.R
trace(coefTable, edit = T)
`coefTable` <- 
  function (model, ...) {
    coef. <- coef(model)[-1] #remove lambda/tho in first row
    .makeCoefTable(coef., sqrt(diag(summary(model, ...)$resvar))[names(coef.)]) #code from MuMIn
  }
#Nested function for coefTable
.makeCoefTable <- 
  #x = model
  function(x, se, df = NA_real_, coefNames = names(x)) {
    if(n <- length(x)) {
      xdefined <- !is.na(x)
      ndef <- sum(xdefined)
      if(ndef < n) {
        if(length(se) == ndef) {
          y <- rep(NA_real_, n); y[xdefined] <- se; se <- y
        }
        if(length(df) == ndef) {
          y <- rep(NA_real_, n); y[xdefined] <- df; df <- y
        }
      }
    }
    if(n && n != length(se)) stop("length(x) is not equal to length(se)")
    ret <- matrix(NA_real_, ncol = 3L, nrow = length(x),
                  dimnames = list(coefNames, c("Estimate", "Std. Error", "df")))
    if(n) ret[, ] <- cbind(x, se, rep(if(is.null(df)) NA_real_ else df,
                                      length.out = n), deparse.level = 0L)
    class(ret) <- c("coefTable", "matrix")
    ret
  }
#Manual setup code and replace it in dredge/model.avg, doesn't work 
#`coefTable.sarlm` <- 
#  function (model, ...) {
#    coef. <- coef(model)[-1]
#  .makeCoefTable(coef., sqrt(diag(summary(model, ...)$resvar))[names(coef.)])
#}

## 4. std.coef -----------------------------------------------------------------------------
##trace.edit override
trace(std.coef, edit = T)
`std.coef` <-
  function(x, partial.sd, ...) {
    #b <- coefTable(x, ...)[, 1L:2L, drop = FALSE]
    b <- coefTable(x, ...)
    #mm <- model.matrix(x) 
    mm <- x$X #update with sarlm format: extracting covariates
    mm <- mm[, match(rownames(b), colnames(mm)), drop = FALSE]
    colnames(mm) <- names(b)
    #b <- b[colnames(mm), ]
    if(partial.sd) {
      bx <- .partialsd(b[, 1L], apply(mm, 2L, sd),
                       .vif(x), nobs(x), sum(attr(mm, "assign") != 0))
    } else {
      #response.sd <- sd(model.response(model.frame(x)))  
      response.sd <- sd(x$y) #update with sarlm format: extracting response
      bx <- apply(mm, 2L, sd) / response.sd
    }
    b[, 1L:2L] <- b[, 1L:2L] * bx
    colnames(b)[1L:2L] <- c("Estimate*", "Std. Error*")
    return(b)
  }

##============ Not run ========================
## Other functions to help testing
# makeArgs ----------
### https://rdrr.io/cran/MuMIn/src/R/makeArgs.R
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


# .expr_beta_arg/dredge.R -----
.expr_beta_arg <- 
  expression({
    if(is.logical(beta) && beta) {
      betaMode <- as.integer(beta)
      strbeta <- if(beta) "sd" else "none"
    } else if(is.character(beta)) {
      strbeta <- match.arg(beta)
      beta <- strbeta != "none"
      betaMode <- (strbeta != "none") + (strbeta == "partial.sd")
    } else {
      cry(, "invalid value for 'beta' : the argument is taken to be \"none\"",
          warn = TRUE)
      betaMode <- 0L
      strbeta <- "none"
    }
  })


# matchCoef/matchCoef.R -------
#trace(matchCoef, edit = T)
`matchCoef` <-
  function(m1, m2,
           all.terms = getAllTerms(m2, intercept = TRUE),
           beta = 0L,
           terms1 = getAllTerms(m1, intercept = TRUE),
           coef1 = NULL,
           allCoef = FALSE,
           ...
  ) {
    
    if(is.null(coef1)) {
      ct <- if (beta != 0L) std.coef(m1, beta == 2L, ...) else coefTable(m1, ...)
      coef1 <- ct[, 1L]
      names(coef1) <- rownames(ct)
    } else if(allCoef) stop("'coef1' is given and 'allCoef' is not FALSE")
    
    if(any((terms1 %in% all.terms) == FALSE)) stop("'m1' is not nested within 'm2'")
    row <- structure(rep(NA_real_, length(all.terms)), names = all.terms)
    
    fxdCoefNames <- fixCoefNames(names(coef1))
    row[terms1] <- NaN
    pos <- match(terms1, fxdCoefNames, nomatch = 0L)
    row[fxdCoefNames[pos]] <- coef1[pos]
    if(allCoef) {
      i <- match(names(coef1), rownames(ct))
      j <- !is.na(i)
      rownames(ct)[i[j]] <- fxdCoefNames[j]
      attr(row, "coefTable") <- ct
    }
    row
  }

# .getLik(), isGEE/REML /utils-models.R -------------
isGEE <- 
  function(object) 
    inherits(object, c("geeglm", "geese", "gee", "geem", "wgee", "yagsResult"))

.isREMLFit <-
  isREML <- 
  function (x) {
    if (inherits(x, "merMod")) 
      return(lme4::isREML(x))
    if (inherits(x, c("lme", "gls", "gam")) && is.character(x$method)) 
      return(x$method[1L] %in% c("lme.REML", "REML"))
    return(FALSE)
  }

`.getLik` <-
  function(x) {
    if(isGEE(x)) {
      list(logLik = quasiLik, name = "qLik")
    } else {
      list(logLik = logLik, name = "logLik")
    }
  }




# .getRank() /utils-models.R ---------------
`.getRank` <-
  function(rank = NULL, rank.args = NULL, object = NULL, ...) {
    rank.args <- c(rank.args, list(...))
    
    if(is.null(rank)) {
      x <- NULL # just not to annoy R check
      IC <- as.function(c(alist(x =, do.call("AICc", list(x)))))
      attr(IC, "call") <- call("AICc", as.name("x"))
      class(IC) <- c("function", "rankFunction")
      return(IC)
    } else if(inherits(rank, "rankFunction") && length(rank.args) == 0L) {
      return(rank)
    } else if(is.list(rank) && length(rank) == 1L && is.function(rank[[1L]])) {
      srank <- names(rank)[1L]
      rank <- rank[[1L]]
    } else {
      srank <- substitute(rank, parent.frame())
      if(srank == "rank") srank <- substitute(rank)
    }
    
    rank <- match.fun(rank)
    ICName <- switch(mode(srank), call = as.name("IC"), character = as.name(srank), name=, srank)
    ICarg <- c(list(as.name("x")), rank.args)
    ICCall <- as.call(c(ICName, ICarg))
    IC <- as.function(c(alist(x =), list(substitute(do.call("rank", ICarg),
                                                    list(ICarg = ICarg)))))
    
    if(!is.null(object)) {
      test <- IC(object)
      if (!is.numeric(test) || length(test) != 1L)
        stop("'rank' should return numeric vector of length 1")
    }
    
    attr(IC, "call") <- ICCall
    class(IC) <- c("function", "rankFunction")
    IC
  }


# checkIsModelDataIdentical/utils-models.R -------
checkIsModelDataIdentical <-
  function(models, error = TRUE) {
    
    cl <- sys.call(sys.parent())
    err <-  if (error) 	function(x) stop(simpleError(x, cl))
    else function(x) warning(simpleWarning(x, cl))
    res <- TRUE
    
    responses <- lapply(models, function(x) getResponseFormula(formula(x)))
    
    if(!all(vapply(responses[-1L], "==", FALSE, responses[[1L]]))) {
      err("response differs between models")
      res <- FALSE
    }
    
    #datas <- lapply(models, function(x) get_call(x)$data)
    # XXX: need to compare deparse'd 'datas' due to ..1 bug(?) in which dotted
    #  arguments (..1 etc) passed by lapply are not "identical"
    datas <- vapply(lapply(models, function(x) get_call(x)$data), asChar, "")
    
    # XXX: when using only 'nobs' - seems to be evaluated first outside of MuMIn
    # namespace which e.g. gives an error in glmmML - the glmmML::nobs method
    # is faulty.
    nresid <- vapply(models, function(x) nobs(x), 1) # , nall=TRUE
    
    if(!all(sapply(datas[-1L], identical, datas[[1L]])) ||
       !all(nresid == nresid[[1L]])) { # better than 'nresid[-1L] == nresid[[1L]]'
      # XXX: na.action checking here
      err("models are not all fitted to the same data")
      res <- FALSE
    }
    invisible(res)
  }


# getResponseFormula/utils-models.R --------
getResponseFormula <-
  function(f) {
    f <- if(!is.null(tf <- attr(f, "terms"))) {
      formula(tf)
    } else formula(f)
    if((length(f) == 2L) || (is.call(f[[2L]]) && f[[2L]][[1L]] == "~"))
      0 else f[[2L]]
  }


# asChar/substitution.R -----
asChar <- function(x, control = NULL, nlines = 1L, ...)
if(is.character(x)) x[1L:nlines] else
  deparse(x, control = control, nlines = nlines, ...)



#.checkNaAction -----------------
.checkNaAction <-
  function(x, cl = get_call(x),
           naomi = c("na.omit", "na.exclude"), what = "model",
           envir = parent.frame()) {
    naact <- NA_character_
    msg <- NA_character_
    
    # handles strings, symbols and calls (let's naively assume no one tries to pass
    # anything else here)
    .getNAActionString <- function(x) {
      if(is.symbol(x)) {
        x <- as.character(x)
      } else if(is.call(x)) {
        x <- eval(x, envir)
        if(is.symbol(x)) x <- as.character(x)
      }
      return(x)
    }
    # TEST:
    #.checkNaAction(list(call = as.call(alist(fun, na.action = getOption("na.action", default = na.fail)))))
    #.checkNaAction(list(call = as.call(alist(fun, na.action = na.fail))))
    #.checkNaAction(list(call = as.call(alist(fun, na.action = na.omit))))
    
    if (!is.null(cl$na.action)) {
      naact <- .getNAActionString(cl$na.action)
      if (naact %in% naomi)
        msg <- sprintf("%s uses 'na.action' = \"%s\"", what, naact)
    } else {
      naact <- formals(eval(cl[[1L]], envir))$na.action
      if (missing(naact)) {
        naact <- getOption("na.action")
        if(is.function(naact)) {
          statsNs <- getNamespace("stats")
          for(i in naomi) if(identical(get(i, envir = statsNs, inherits = FALSE), naact,
                                       ignore.environment = TRUE)) {
            naact <- i
            break
          }
        }
        
        naact <- .getNAActionString(naact)
        if (is.character(naact) && (naact %in% naomi))
          msg <- sprintf("%s's 'na.action' argument is not set and options('na.action') is \"%s\"",
                         what, naact)
      } else if (!is.null(naact)) {
        naact <- .getNAActionString(naact)
        if (naact %in% naomi)
          msg <- sprintf("%s uses the default 'na.action' = \"%s\"", what, naact)
      }
    }
    res <- is.na(msg)
    attr(res, "na.action") <- naact
    attr(res, "message") <- msg
    res
  }


# return_null-----
return_null <-
  function(...) NULL 

# fixCoefNames--------
`fixCoefNames` <-
  function(x, peel = TRUE) {
    if(!length(x)) return(x)
    ia <- grep(":", x, fixed = TRUE)
    if(!length(ia)) return(structure(x, order = rep.int(1L, length(x))))
    
    ixi <- x[ia]
    if(peel) {
      # peel only when ALL items are prefixed/wrapped
      # for pscl::hurdle. Cf are prefixed with count_|zero_
      if(peel <- all(startsWith(ixi, c("count_", "zero_")))) {
        pos <- regexpr("_", ixi, fixed = TRUE)
        peelpfx <- substring(ixi, 1L, pos)
        peelsfx <- ""
        ixi <- substring(ixi, pos + 1L)
      } else {
        # unmarkedFit with its phi(...), lambda(...) etc...
        if(peel <- all(endsWith(ixi, ")"))) {
          # only if 'XXX(...)', i.e. exclude 'XXX():YYY()' or such
          m <- regexpr("^(([^()]*)\\(((?:[^()]*|(?1))*)\\))$", ixi, perl = TRUE, useBytes = TRUE)
          cptgrps <- .matches(ixi, m)
          if(peel <- all(cptgrps[, 2L] != "")) {
            peelpfx <- paste0(cptgrps[, 2L], "(")
            peelsfx <- ")"
            ixi <- cptgrps[, 3L]
          }
        }
      }
    }
    # replace {...}, [...], (...), ::, and ::: with placeholders
    m <- gregexpr("(?:\\{(?:[^\\{\\}]*|(?0))*\\}|\\[(?:[^\\[\\]]*|(?0))*\\]|\\((?:[^()]*|(?0))*\\)|(?>:::?))", ixi, perl = TRUE)
    xtpl <- ixi
    regmatches(xtpl, m) <- lapply(m, function(x) {
      if((ml <- attr(x, "match.length"))[1L] == -1L) return(character(0L))
      sapply(ml, function(n) paste0(rep("_", n), collapse = ""))
    })
    
    # split by ':' and sort
    splits <- gregexpr(":", xtpl, fixed = TRUE)
    ixi <- mapply(function(x, p) {
      if(p[1L] == -1) return(x)
      paste0(base::sort(substring(x, c(1L, p + 1L), c(p - 1L, nchar(x)))), collapse = ":")
    }, ixi, splits, USE.NAMES = FALSE, SIMPLIFY = TRUE)
    
    if(peel) ixi <- paste0(peelpfx, ixi, peelsfx)
    
    x[ia] <- ixi
    ord <- rep.int(1L, length(x))
    ord[ia] <- vapply(splits, length, 0L) + 1L
    attr(x, "order") <- ord
    x
  }



# prettyEnumStr ------------
`prettyEnumStr` <- function(x, sep = ", ", sep.last = gettext(" and "), quote = TRUE) {
  n <- length(x)
  if(is.function(quote))
    x <- quote(x) else {
      if(identical(quote, TRUE)) quote <- '"'
      if(is.character(quote)) x <- paste0(quote, x, quote)
    }
  paste0(x, if(n > 1L) c(rep(sep, n - 2L), sep.last, "") else NULL,
         collapse = "")
}


# formula_margin_check------
formula_margin_check <- function(j, m) {
  stopifnot(is.logical(j))
  !any(m[!j, j], na.rm = TRUE)
}




# `.makeListNames` ------------
`.makeListNames` <- function(x) {
  nm <- names(x)
  lapply(seq_along(x), function(i) {
    if(is.null(nm) || nm[i] == "") {
      switch(mode(x[[i]]),
             call = {
               v <- asChar(x[[i]], width.cutoff = 20L)
               if(length(v) != 1L) v <- sprintf("%s...", v[1L])
               v },
             symbol =, name = as.character(x[[i]]),
             NULL =, logical =, numeric =, complex =, character = x[[i]], i
      )
    } else nm[i]
  })
}


# evalExprInEnv-----------
evalExprInEnv <- function(expr, env, enclos, ...) {
  list2env(list(...), env)
  eval(expr, envir = env, enclos = enclos)
}


# .DebugPrint-----
`.DebugPrint` <-
  function (x) {
    if (isTRUE(getOption("debug.MuMIn"))) {
      fun <- asChar(sys.call(sys.parent())[[1L]])
      name <- substitute(x)
      cat(sprintf("<%s> ~ ", fun))
      if(is.language(name)) cat(asChar(name), "= \n")
      print(x)
    }
  }


# exprapply0------
`exprapply0` <- function(e, name, func, ...)
  exprApply(e, name, func, ..., symbols = FALSE)



#.subst.with----------------
.subst.with <- function (x, fac, allTerms, vName, envir = parent.frame()) {
  if (length(x) > 4L) cry(x, "too many arguments [%d]", length(x) - 1L)
  if (length(x[[2L]]) == 2L && x[[2L]][[1L]] == "+") {
    fun <- "all"
    sx <- asChar(x[[2L]][[2L]], backtick = FALSE)
  } else {
    fun <- "any"
    sx <- asChar(x[[2L]], backtick = FALSE)
  }
  dn <- dimnames(fac)
  if (!(sx %in% dn[[2L]])) cry(x, "unknown variable name '%s'", sx)
  xorder <- if(length(x) >= 3L) as.integer(eval(x[[3L]], envir))
  else unique(rowSums(fac))
  i <- which(fac[, sx])
  j <- which(is.element(rowSums(fac[i, , drop = FALSE]), xorder))
  if(length(j) == 0L) cry(x, "no terms match the criteria")    
  as.call(c(as.name(fun), call("[", vName, as.call(c(as.name("c"), match(dn[[1L]][i[j]], allTerms))))))
}

#.subst.term-----
.subst.term <- function(x) {
  if(length(x) < 2L) cry(x, "'Term' needs one argument")
  as.name(asChar(x[[2L]]))
}


# updateDeps ------------
updateDeps <-
  function(expr, deps) {
    ret <- list()
    env <- sys.frame(sys.nframe())
    expr <- exprapply0(expr, "dc", function(z) {
      v <- vapply(as.list(z[-1L]), asChar, "")
      n <- length(v)
      k <- match(v, colnames(deps))
      for(i in 2L:n) deps[k[1L:(i - 1L)], k[i]] <- TRUE
      assign("deps", deps, envir = env, inherits = FALSE)
      TRUE
    })
    list(deps = deps, expr = expr)
  }


#subst ----------
subst <-
  function(expr, envir = NULL, ...) {
    eval.parent(call("substitute", expr, c(envir, list(...))))
  }



#.subst.vars.for.args--------
.subst.vars.for.args <- function(e) {
  for(i in 2L:length(e))
    if(!is.name(e[[i]]))
      e[[i]] <- as.name(asChar(e[[i]]))
  e
}


#.subst.vars.for.args--------
`.subst.names.for.items` <-
  function(expr, names, varName, n = length(names), fun = "[") {
    exprApply(expr, names, symbols = TRUE,
              function(x, v, fun, varName, parent) {
                if(is.call(parent) && any(parent[[1L]] == c("I", "$", "@")))
                  return(x)
                if(length(x) == 1L)
                  return(call(fun, varName, match(asChar(x), v)))
                x
              }, v = names, fun = fun, varName = as.name(varName))
  }

#.subst.v---------
.subst.v <- function(x, cVar, fn) {
  if(length(x) > 2L) cry(x, "discarding extra arguments", warn = TRUE)
  i <- which(fn == x[[2L]])[1L]
  if(is.na(i)) cry(x, "'%s' is not a valid name of 'varying' element",
                   as.character(x[[2L]]), warn = TRUE)
  call("[[", cVar, i)
}


# `cry`-------
`cry` <-
  function(Call = NA, Message, ..., warn = FALSE, domain = paste0("R-", .packageName)) {
    if (is.character(Call)) {
      Call <- call(Call[1L], sys.call(-1L)[[1L]])
    } else if(is.numeric(Call)) {
      Call <- sys.call(Call - 1L)
    } else if (!is.call(Call) && !is.null(Call))
      Call <- sys.call(-1L)
    if(warn) warning(simpleWarning(gettextf(Message, ..., domain = domain), Call)) else
      stop(simpleError(gettextf(Message, ..., domain = domain), Call))
  }

#manual function ---- not run -----
`matchCoef.sarlm` <-
  function(m1, m2,
           all.terms = getAllTerms(m2, intercept = TRUE),
           beta = 0L,
           terms1 = getAllTerms(m1, intercept = TRUE),
           coef1 = NULL,
           allCoef = FALSE,
           ...
  ) {
    
    if(is.null(coef1)) {
      ct <- if (beta != 0L) std.coef.sarlm(m1, beta == 2L, ...) else coefTable.sarlm(m1, ...)
      coef1 <- ct[, 1L]
      names(coef1) <- rownames(ct)
    } else if(allCoef) stop("'coef1' is given and 'allCoef' is not FALSE")
    
    if(any((terms1 %in% all.terms) == FALSE)) stop("'m1' is not nested within 'm2'")
    row <- structure(rep(NA_real_, length(all.terms)), names = all.terms)
    
    fxdCoefNames <- fixCoefNames(names(coef1))
    row[terms1] <- NaN
    pos <- match(terms1, fxdCoefNames, nomatch = 0L)
    row[fxdCoefNames[pos]] <- coef1[pos]
    if(allCoef) {
      i <- match(names(coef1), rownames(ct))
      j <- !is.na(i)
      rownames(ct)[i[j]] <- fxdCoefNames[j]
      attr(row, "coefTable") <- ct
    }
    row
  }


