variable.names.fcs2Fit <-
function(object, q = object$multiRun, r = TRUE, abundance = TRUE, prevalence = TRUE, linear = TRUE, 
                                   rw = "group", spatial = "group", hyperparams = "scale")
{
    attach(object, warn.conflicts=FALSE)
    on.exit(detach(object))
    params <- c()
    
    # Catch probability q
    if (q)
        params <- c(params, "q")
    
    # Shape parameter r
    if (r)
        params <- c(params, "r")
      
    # Abundance terms:
    if (abundance) {
        # Constant and linear terms
        if (linear)
            params <- c(params, paste("beta", c("const", muLinearVars), sep="."))
        
        # Random walk terms
        if (length(c(muRW1Vars, muRW2Vars)) > 0) {
            if (!is.na(pmatch(rw, c("group", TRUE))))
                params <- c(params, paste("beta", c(muRW1Vars, muRW2Vars), sep="."))
            else if (!is.na(pmatch(rw, "singular"))) {
                for (var in c(muRW1Vars, muRW2Vars))
                    params <- c(params, paste("beta.", var, "[", 1:length(object$rwBoundaries[[paste("beta", var, sep=".")]]), "]", sep=""))
            }
        }    
    
        # Spatial term
        if (!is.null(muSpatialVar)) {
            if (!is.na(pmatch(spatial, c("group", TRUE))))
                params <- c(params, paste("beta", muSpatialVar, sep="."))
            else if (!is.na(pmatch(spatial, "singular")))
                params <- c(params, paste("beta.", muSpatialVar, "[", 1:length(object$muAdjacency$num), "]", sep=""))
        }
        
        # Random walk hyperparameters
        if (length(c(muRW1Vars, muRW2Vars)) > 0) {
            if (!is.na(pmatch(hyperparams, "scale")))
                params <- c(params, paste("sigma", c(muRW1Vars, muRW2Vars), sep="."))
            else if (!is.na(pmatch(hyperparams, "precision")))
                params <- c(params, paste("tau", c(muRW1Vars, muRW2Vars), sep="."))
        }    

        # Spatial hyperparameter
        if (!is.null(muSpatialVar)) {
            if (!is.na(pmatch(hyperparams, "scale")))
                params <- c(params, paste("sigma", muSpatialVar, sep="."))
            else if (!is.na(pmatch(hyperparams, "precision")))
                params <- c(params, paste("tau", muSpatialVar, sep="."))
        }
    }
    
    # Prevalence terms:
    if (prevalence) {
        # Constant and linear terms
        if (linear)
            params <- c(params, paste("gamma", c("const", rhoLinearVars), sep="."))
        
        # Random walk terms
        if (length(c(rhoRW1Vars, rhoRW2Vars)) > 0) {
            if (!is.na(pmatch(rw, c("group", TRUE))))
                params <- c(params, paste("gamma", c(rhoRW1Vars, rhoRW2Vars), sep="."))
            else if (!is.na(pmatch(rw, "singular"))) {
                for (var in c(rhoRW1Vars, rhoRW2Vars))
                    params <- c(params, paste("gamma.", var, "[", 1:length(object$rwBoundaries[[paste("gamma", var, sep=".")]]), "]", sep=""))
            }
        }    
    
        # Spatial terms
        if (!is.null(rhoSpatialVar)) {
            if (!is.na(pmatch(spatial, c("group", TRUE))))
                params <- c(params, paste("gamma", rhoSpatialVar, sep="."))
            else if (!is.na(pmatch(spatial, "singular")))
                params <- c(params, paste("gamma.", rhoSpatialVar, "[", 1:length(object$rhoAdjacency$num), "]", sep=""))
        }
        
        # Random walk hyperparameterw
        if (length(c(rhoRW1Vars, rhoRW2Vars)) > 0) {
            if (!is.na(pmatch(hyperparams, "scale")))
                params <- c(params, paste("nu", c(rhoRW1Vars, rhoRW2Vars), sep="."))
            else if (!is.na(pmatch(hyperparams, "precision")))
                params <- c(params, paste("phi", c(rhoRW1Vars, rhoRW2Vars), sep="."))
        }    

        # Spatial hyperparameter
        if (!is.null(rhoSpatialVar)) {
            if (!is.na(pmatch(hyperparams, "scale")))
                params <- c(params, paste("nu", rhoSpatialVar, sep="."))
            else if (!is.na(pmatch(hyperparams, "precision")))
                params <- c(params, paste("phi", rhoSpatialVar, sep="."))
        }
    }
    
    params
}

