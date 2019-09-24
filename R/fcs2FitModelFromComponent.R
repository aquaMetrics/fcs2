.fcs2FitModelFromComponent <-
function(runTotalVars = NULL, allRunsTotalVar = NULL, allRunsRangeVars = NULL, modelMatrix, surveyAreaVar = "SurveyArea", nRunsVar = NULL, 
         muLinearVars = c(), muRW1Vars = c(), muRW2Vars = c(), muSpatialVar = NULL, muAdjacency = NULL,
         rhoLinearVars = c(), rhoRW1Vars = c(), rhoRW2Vars = c(), rhoSpatialVar = NULL, rhoAdjacency = NULL, 
         subset = 1:nrow(modelMatrix), na.action, dataType = factor(rep("run", nrow(modelMatrix)), levels=c("run", "total", "range")),
         rwNoLevels = 10, rwBoundaries = list(), prior.parameters = list(), initial.values = list(), 
         runINLA, runBUGS = TRUE, n.chains = 2, n.iter = 2000, n.burnin = floor(n.iter/2), 
         n.thin = max(1, floor(n.chains * (n.iter - n.burnin)/n.sims)), n.sims = 1000, 
         bugsFilename = "model.bug", bugsProgram = "OpenBUGS", 
         verbose = TRUE, muFormula, rhoFormula, call = match.call(), fit, ...)
{    
    # get default na.action if missing
    if (missing(na.action)) {
        na.action <- getOption("na.action")
        if (is.null(na.action))
            na.action <- na.omit  # use 'na.omit' if option not set
        else
            na.action <- eval(parse(text=na.action))
    }
        
    # check parameter not given twice where not allowed
    if (length(var <- intersect(muLinearVars, union(muRW1Vars, muRW2Vars))) > 0)
        stop(paste("variable '", var[1], "' cannot be used for both linear and random walk terms", sep=''))
    if (length(var <- intersect(muRW1Vars, muRW2Vars)) > 0)
        stop(paste("variable '", var[1], "' cannot be used for both a first and a second order random walk term", sep=''))
    if (length(var <- intersect(rhoLinearVars, union(rhoRW1Vars, rhoRW2Vars))) > 0)
        stop(paste("variable '", var[1], "' cannot be used for both linear and random walk terms", sep=''))
    if (length(var <- intersect(rhoRW1Vars, rhoRW2Vars)) > 0)
        stop(paste("variable '", var[1], "' cannot be used for both a first and a second order random walk term", sep=''))


    # create fit object, unless already provided
    if (missing(fit)) {
        fit <- c(list(call=call), 
                 mget(ls()[na.omit(match(c("runTotalVars", "allRunsTotalVar", "allRunsRangeVars", "modelMatrix", "surveyAreaVar", "nRunsVar", 
                                           "muLinearVars", "muRW1Vars", "muRW2Vars", "muSpatialVar", "muAdjacency",
                                           "rhoLinearVars", "rhoRW1Vars", "rhoRW2Vars", "rhoSpatialVar", "rhoAdjacency",
                                           "dataType", "rwNoLevels", "rwBoundaries", "prior.parameters", "initial.values"), ls()))], environment()))                                                  
        fit$modelMatrix <- fit$modelMatrix[subset, , drop=FALSE]

        # extract modelMatrix containing only required variables
        fit$modelMatrix <- fit$modelMatrix[, unique(c(fit$runTotalVars, fit$allRunsTotalVar, fit$allRunsRangeVars, fit$surveyAreaVar, fit$nRunsVar,
                                                      fit$muLinearVars, fit$muRW1Vars, fit$muRW2Vars, fit$muSpatialVar,
                                                      fit$rhoLinearVars, fit$rhoRW1Vars, fit$rhoRW2Vars, fit$rhoSpatialVar)), drop=FALSE]
        if (class(na.action) == "function") {
            fit$modelMatrix <- na.action(fit$modelMatrix)
            fit$na.action <- attr(fit$modelMatrix, "na.action")
        } else
            fit$na.action <- na.action
          
        # calculate range of covariates
        fit$covariateMin <- apply(fit$modelMatrix[, unique(c(fit$muLinearVars, fit$muRW1Vars, fit$muRW2Vars, 
                                                             fit$rhoLinearVars, fit$rhoRW1Vars, fit$rhoRW2Vars)), drop=FALSE], 2, min, na.rm=TRUE)
        fit$covariateMax <- apply(fit$modelMatrix[, unique(c(fit$muLinearVars, fit$muRW1Vars, fit$muRW2Vars, 
                                                             fit$rhoLinearVars, fit$rhoRW1Vars, fit$rhoRW2Vars)), drop=FALSE], 2, max, na.rm=TRUE)
       
        # check that min/max are not Inf
        ## NOTE: could alternatively remove these with a warning?
        if (length(var <- which(fit$covariateMin == -Inf)))
            stop(paste("variable '", names(var), "' has -Inf values", sep=""))
        if (length(var <- which(fit$covariateMax == Inf)))
            stop(paste("variable '", names(var), "' has Inf values", sep=""))           
        
        # add or create formulae
        if (missing(muFormula))
            muFormula <- formula(paste("~", 
                    paste(collapse=' + ', 
                          c("1",
                            if(length(fit$muLinearVars) > 0)
                                paste(fit$muLinearVars, collapse=" + ", sep=''),
                            if(length(fit$muRW1Vars) > 0)
                                paste("rw1(", fit$muRW1Vars, ")", collapse=" + ", sep=''),
                            if(length(fit$muRW2Vars) > 0)
                                paste("rw2(", fit$muRW2Vars, ")", collapse=" + ", sep=''),
                            if(!is.null(fit$muSpatialVar))
                                paste("spatial(", fit$muSpatialVar, ")", collapse=" + ", sep='')
                           )), sep=''))
        fit$muFormula <- muFormula
        
        if (missing(rhoFormula))
            rhoFormula <- formula(paste("~", 
                    paste(collapse=' + ', 
                          c("1",
                            if(length(fit$rhoLinearVars) > 0)
                                paste(fit$rhoLinearVars, collapse=" + ", sep=''),
                            if(length(fit$rhoRW1Vars) > 0)
                                paste("rw1(", fit$rhoRW1Vars, ")", collapse=" + ", sep=''),
                            if(length(fit$rhoRW2Vars) > 0)
                                paste("rw2(", fit$rhoRW2Vars, ")", collapse=" + ", sep=''),
                            if(!is.null(fit$rhoSpatialVar))
                                paste("spatial(", fit$rhoSpatialVar, ")", collapse=" + ", sep='')
                           )), sep=''))
        fit$rhoFormula <- rhoFormula
    }                                                   
    class(fit) <- "fcs2Fit" 

    
    # Set total number of observations
    ## NOTE: if na.action = na.pass these may not all be useable
    fit$N <- nrow(fit$modelMatrix)
   
    # set whether using multi-run extension
    # (if no runs > 1)
    if (is.null(fit$multiRun))
        fit$multiRun <- !is.null(fit$nRunsVar) && max(fit$modelMatrix[!is.na(fit$dataType), fit$nRunsVar]) > 1

    # create list of all parameters in the model (as appear in BUGS code, ie including rw and spatial)
    parameters <- variable.names(fit, hyperparams="precision")
    
    # check initial.values
    if (length(var <- setdiff(names(fit$initial.values), parameters)) > 0)
        stop(paste("initial value given for unrecognised variable '", var[1], "'", sep=''))
        
    # check rw boundaries
    if (length(var <- setdiff(names(fit$rwBoundaries), variable.names(fit, q=FALSE, r=FALSE, linear=FALSE, hyperparams=FALSE))) > 0)
        stop(paste("random walk boundaries given for unrecognised variable '", var[1], "'", sep=''))
    
    # fill in missing rw boundaries
    if (length(c(fit$muRW1Vars, fit$muRW2Vars, fit$rhoRW1Vars, fit$rhoRW2Vars)) > 0)
        fit$rwBoundaries <- .fcs2SetRWBoundaries(fit)

    # Fill gaps in priors with default values
    fit$prior.parameters <- .fcs2SetDefaultPriors(fit)

    # Calculate prior parameters from mean/sd etc
    fit$prior.parameters <- fcs2Priors(fit$prior.parameters, print=FALSE, plot=FALSE)
    
    # check prior.parameters
    if (length(var <- setdiff(names(fit$prior.parameters), variable.names(fit, hyperparams="prec", rw=FALSE, spatial=FALSE))) > 0)
        stop(paste("prior parameters given for unrecognised variable '", var[1], "'", sep=''))
  
    # Fill in initial values gaps with INLA estimates
    if ((!(missing(runINLA) || is.null(runINLA)) && runINLA != FALSE) || 
        ((min(parameters %in% names(fit$initial.values)) == 0) && (missing(runINLA) || is.null(runINLA) || runINLA != FALSE))) {
        if (missing(runINLA))
            fit$inlaFits <- fcs2RunINLA(fit, TRUE, verbose)
        else
            fit$inlaFits <- fcs2RunINLA(fit, runINLA, verbose)
        fit$initial.values <- .fcs2SetInitialValues(fit)
    }
    
    
    # Runs BUGS
    if (runBUGS) {
        # check that all initial values are given
        if (min(parameters %in% names(fit$initial.values)) == 0)
            warning("initial values are not all set so these will be sampled from priors")

        # Create data for BUGS
        bugsData <- .fcs2CreateBUGSData(fit)
        
        # write BUGS model file
        fcs2WriteModel(bugsFilename, fit$muLinearVars, c(fit$muRW1Vars, fit$muRW2Vars), fit$muSpatialVar,
                       fit$rhoLinearVars, c(fit$rhoRW1Vars, fit$rhoRW2Vars), fit$rhoSpatialVar, fit$multiRun,
                       "catch" %in% names(bugsData), "totalCatchMin" %in% names(bugsData), 
                       qUnif=is.null(bugsData$q.a))  # q uniform if q.a not in bugsData
        
        # Create initial values for each chain (correcting variable names)
        ### NOTE: currently the same for all chains but could consider jittering all but first chain
        inits <- list()
        for (i in 1:n.chains) {
            inits <- c(inits, list(fit$initial.values))
            names(inits[[length(inits)]]) <- make.names(names(inits[[length(inits)]]))
        }
        inits <- inits
        
        # Correct parameter names
        parameters <- make.names(parameters)
    
        # run BUGS
        if (verbose) {
            cat("Fitting model with ", n.chains, " chains of ", n.iter, " iterations using ", bugsProgram, " ... ", sep="")
            if (bugsProgram %in% c("openbugs", "OpenBUGS"))
                cat("\n")
            flush.console()
        }
        # if OpenBUGS, run fixed version of openbugs function
        if (bugsProgram %in% c("openbugs", "OpenBUGS"))
            fit$bugsFit <- .openbugs.fcs2(bugsData, inits, parameters, bugsFilename, n.chains, n.iter, n.burnin, n.thin, n.sims, fit=fit, ...)
        else
            fit$bugsFit <- bugs(bugsData, inits, parameters, bugsFilename, n.chains, n.iter, n.burnin, n.thin, n.sims, program=bugsProgram, ...)
        if (verbose)
            cat("done\n")
        
        # correct BUGS output
        fit$bugsFit <- .fcs2CorrectBUGSOutput(fit)
        
        # issue warning if MCMC not converged
        if (n.chains > 1 && max(fit$bugsFit$summary[, "Rhat"]) > 1.1)
            warning("MCMC chains may not have converged.  Use summary(fit=...) to view 'Rhat' which should be close to 1 for every parameter.")
    }
    
    # return fit object
    fit
}

