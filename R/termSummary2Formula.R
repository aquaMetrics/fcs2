termSummary2Formula <-
function(termSum)
{
    terms <- character(0)
    termName <- names(termSum)
    
    for (i in 1:length(termSum)) {
        if (termSum[i] == "rw2")
            terms <- c(terms, paste('rw2(', termName[i], ')', sep=''))
        else if (termSum[i] == "spatial")
            terms <- c(terms, paste('spatial(', termName[i], ', adjacency)', sep=''))
        else if (termSum[i] != "") { # must be a number
            terms <- c(terms, termName[i])
            
            # add higher order linear terms
            if (termSum[i] != "1") {
                termOrder <- as.numeric(as.matrix(termSum[i]))
                for (j in 2:termOrder)
                    terms <- c(terms, paste('I(', termName[i], '^', j, ')', sep=''))
            }
        }
    }
    if (length(terms) == 0)
        terms <- "1"
    
    # return as formula
    formula(paste("~", paste(terms, collapse=" + ")))
}

