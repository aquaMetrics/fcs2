## openbugs_fixed.R : FIXED openbugs function from R2WinBUGS package :
##    Fixes bug that crashes R when BRugs::samplesMonitors fails
##

.openbugs.fcs2 <- function (data, inits, parameters.to.save, model.file = "model.txt", 
    n.chains = 3, n.iter = 2000, n.burnin = floor(n.iter/2), 
    n.thin = max(1, floor(n.chains * (n.iter - n.burnin)/n.sims)), 
    n.sims = 1000, DIC = TRUE, bugs.directory = "c:/Program Files/OpenBUGS/", 
    working.directory = NULL, digits = 5, over.relax = FALSE, 
    bugs.seed = NULL, fit) 
{
    if (!is.R()) 
        stop("OpenBUGS is not yet available in S-PLUS")
    if (!require("BRugs")) 
        stop("BRugs is required")
    modelFile <- model.file
    numChains <- n.chains
    nBurnin <- n.burnin
    nIter <- n.iter - n.burnin
    nThin <- n.thin
    if (DIC) 
        parameters.to.save <- c(parameters.to.save, "deviance")
    parametersToSave <- parameters.to.save
    if (is.null(working.directory)) {
        working.directory <- tempdir()
        inTempDir <- TRUE
    }
    savedWD <- getwd()
    setwd(working.directory)
    on.exit(setwd(savedWD))
    if (inTempDir && basename(model.file) == model.file) 
        try(file.copy(file.path(savedWD, model.file), model.file, 
            overwrite = TRUE))
    if (!file.exists(modelFile)) {
        stop(modelFile, " does not exist")
    }
    if (file.info(modelFile)$isdir) {
        stop(modelFile, " is a directory, but a file is required")
    }
    if (!length(grep("\r\n", readChar(modelFile, 10^3)))) {
        message("Carriage returns added to model file ", modelFile)
        model <- readLines(modelFile)
        try(writeLines(model, modelFile))
    }
    BRugs::modelCheck(modelFile)
    if (!(is.vector(data) && is.character(data) && all(file.exists(data)))) {
        data <- BRugs::bugsData(data, digits = digits)
    }
    if (inTempDir && all(basename(data) == data)) 
        try(file.copy(file.path(savedWD, data), data, overwrite = TRUE))
    BRugs::modelData(data)
    BRugs::modelCompile(numChains)
    if (!is.null(bugs.seed)) 
        BRugs::modelSetSeed(newSeed = bugs.seed)
    if (missing(inits) || is.null(inits)) {
        BRugs::modelGenInits()
    }
    else {
        if (is.list(inits) || is.function(inits) || (is.character(inits) && 
            !any(file.exists(inits)))) {
            inits <- BRugs::bugsInits(inits = inits, numChains = numChains, 
                digits = digits)
        }
        if (inTempDir && all(basename(inits) == inits)) 
            try(file.copy(file.path(savedWD, inits), inits, overwrite = TRUE))
        BRugs::modelInits(inits)
        BRugs::modelGenInits()
    }
    BRugs::samplesSetThin(nThin)
    if (getOption("BRugsVerbose")) {
        cat("Sampling has been started ...\n")
        flush.console()
    }
    BRugs::modelUpdate(nBurnin, overRelax = over.relax)
    if (DIC) {
        BRugs::dicSet()
        on.exit(BRugs::dicClear(), add = TRUE)
    }
    BRugs::samplesSet(parametersToSave)
    BRugs::modelUpdate(nIter, overRelax = over.relax)
    
    
    ### BUG FIX: avoid call to samplesMonitors as can fail ###
    
    ## OLD: params <- sort.name(BRugs::samplesMonitors("*"), parametersToSave)
    
    # extract rw and spatial vars
    groupPars <- variable.names.fcs2Fit(fit, r=FALSE, q=FALSE, linear=FALSE, hyper=FALSE)
    groupPars <- make.names(groupPars)  # convert to names
    
    # create new list of pars without these
    pars.long <- setdiff(parametersToSave, groupPars)  # removed these names
    
    # Now add long versions of rw and spatial terms:    
    
    # Random walk terms:
    for (var in c(fit$muRW1Vars, fit$muRW2Vars))
        pars.long <- c(pars.long, paste(make.names(paste("beta.", var, sep="")),
                                        "[", 1:length(fit$rwBoundaries[[paste("beta", var, sep=".")]]), "]", sep=""))
    for (var in c(fit$rhoRW1Vars, fit$rhoRW2Vars))
        pars.long <- c(pars.long, paste(make.names(paste("gamma.", var, sep="")),
                                        "[", 1:length(fit$rwBoundaries[[paste("gamma", var, sep=".")]]), "]", sep=""))
                                        
    # Spatial terms:
    if (!is.null(fit$muSpatialVar))
        pars.long <- c(pars.long, paste(make.names(paste("beta.", fit$muSpatialVar, sep="")),
                                        "[", 1:length(fit$muAdjacency$num), "]", sep=""))
    if (!is.null(fit$rhoSpatialVar))
        pars.long <- c(pars.long, paste(make.names(paste("gamma.", fit$rhoSpatialVar, sep="")),
                                        "[", 1:length(fit$rhoAdjacency$num), "]", sep=""))
    
    # create full params list with order or parametersToSave
    params <- .sort.name(pars.long, parametersToSave)

    ### END BUG FIX ###
    
    
    samples <- sapply(params, BRugs::samplesSample)
    n.saved.per.chain <- nrow(samples)/numChains
    samples.array <- array(samples, c(n.saved.per.chain, numChains, 
        ncol(samples)))
    dimnames(samples.array)[[3]] <- dimnames(samples)[[2]]
    if (DIC) {
        DICOutput <- BRugs::dicStats()
    }
    else {
        DICOutput <- NULL
    }
    bugs.output <- as.bugs.array(sims.array = samples.array, 
        model.file = modelFile, program = "OpenBUGS", DIC = DIC, 
        DICOutput = DICOutput, n.iter = n.iter, n.burnin = n.burnin, 
        n.thin = n.thin)
    bugs.output
}

.sort.name <- function(a, b)
{
    bracket.pos <- regexpr("\\[", a)
    a.stem <- substr(a, 1, ifelse(bracket.pos > 0, bracket.pos - 1, nchar(a)))
    a[order(match(a.stem, b))]
}