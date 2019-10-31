##>>......................................................................<<##
##>>.... Functions to process SNOMED descriptions into readable text    ..<<## 
##>>......................................................................<<##

##>>.. Makes a name in the form of Name1 (a.k.a. Name2)
makeDoubleName <- function(name1, name2) {
    name_dist_cutoff <- 0.15  
    if(length(name1)==0) {
        return(name2)
    }
    else {
        name_dist <- stringdist(name1, name2, method='jw')
        if(name_dist > name_dist_cutoff) {
            return(paste0(name2, " (a.k.a. ", name1, ")"))
        } else {
            return(name1)
        }
    }
}

##>>.. Makes a name in the form of Name1 (a.k.a. Name2)
makeDoubleNameFromMany <- function(name1, names) {
    name_dist_cutoff <- 0.15  
    if(length(name1)==0 | name1=="") {
        return(pickShortest(names))
    }
    else {
        if(length(names)>0) {
            for(i in 1:length(names)) {  #try to find the best second name
                # 1) find one without parenthesis
                # 2) make sure it's not the same as the first name
                # 3) if it is - go to the next one
                if( length( grep('\\(|\\)', names[i])) == 0) {
                    name_dist <- stringdist(name1, names[i], method='jw')
                    if(name_dist > name_dist_cutoff) {
                        return(paste0(names[i], " (a.k.a. ", name1, ")")) 
                    }
                }
            }
        }
    }
    return(capitalize(name1))  #only need to capitalize name1 - because names come from SNOMED already capitalized
}

##>>.. returns a value in names most different from name1 - or "" if the two strings are exactly the same
returnMostDiff<- function(name1, names) {
    name_dist_cutoff <- 0.15  
    name_dist <- 0
    prev_dist <- 0
    index <- 0
    if(length(name1)==0 | name1=="") {
        return(pickShortest(names))
    }
    else {
        if(length(names)>0) {
            for(i in 1:length(names)) {  #try to find the best second name
                name_dist <- stringdist(tolower(name1), tolower(names[i]), method='jw')
                if(name_dist > prev_dist) {
                    prev_dist <- name_dist
                    index <- i
                }
            }
        }
    }
    if(prev_dist == 0) { return ("") }
    return(names[index])  #only need to capitalize name1 - because names come from SNOMED already capitalized
}

##>>.. Trying to remove spurious redundant verbiage
##>>.. Ideally want to remove "structure" in "Brain structure" etc, but leave in Blood vessel structure 
cleanSnoName <- function(snoname) {
    stopwords <- c("\\(disorder\\)|\\(body structure\\)")
    t <- unlist(strsplit(snoname, stopwords))
    t <- unlist(lapply(t, trimws)) 
    return(paste(t[!(t %in% stopwords)], collapse=' '))
}

##>>.. Removes redundant verbiage for a whole vector of strings
cleanSnoNames <- function(snonames) {
    return(unlist(lapply(snonames, cleanSnoName)))
}

##>>.. Concatenates similar terms in the form of A, B, and C
concatTerms <- function (terms) {
    len <- length(terms)
    if(len<2) { return (terms) }
    if(len == 2) { return (paste0(terms[1], " and ", terms[2])) }
    str <- ""
    str <- paste0(str, terms[1: (len-1)], collapse=", ")
    str <- paste0(str, ", and ", terms[len])
    return(paste0(str,collapse=""))
}

##>>.. Replaces what is known to be bad redundant structure and other concept names with better ones
##>>.. This is manual labor. Automating it would require injecting domain knowledge...
##>>.. Automating was not attempted by OntoVerbal either. 
englishize1 <- function(term) {
    mappings <- c("Cerebrovascular system structure"="the cerebrovascular system",
                  "Blood vessel structure" = "the blood vessels", 
                  "Brain structure" = "the brain")
    return(mappings[term])
}

##>>.. same as above, for a vector of terms
englishize <- function(terms) {
    return(unlist(lapply(terms, englishize1)))
}

##>>.. returns the shortest name from a list of names
pickShortest <- function(terms) {
    if(length(terms)<2) {
        return(terms)
    }
    minlength <- 300
    shortestI <- 0
    for(i in 1:length(terms)) {
        if(length(terms[i]) < minlength) {
            minlength <- length(terms[i])
            shortestI = i
        }
    }
    return(terms[i])
}

##>>.. Picks shorted unique names for all terms given - a term comes with lists $nodeId and $nodeName
pickShortestUnique <- function(terms) {
    #terms is a dataframe with an id
    if(nrow(terms) <2) {
        return(unique(terms$nodeName))
    }
    
    uniqueNodes <- unique(terms$nodeId)
    names <- c()
    for(i in 1:length(uniqueNodes)) {
        names <- append(names, pickShortest(terms$nodeName[terms$nodeId==uniqueNodes[i]]))
    }
    return(names)
}

pickTopN <- function(terms, n) {
    t <- c()
    if(nrow(terms)>0) {
        # terms <- snomed$subjRels 
        # n=3
        t <- unique(terms$nodeName)
        if(length(t)>n) {
            t <- t[1:n]
        }
    }
    cleanSnoNames(t)
    return(t)
}

pickTopNVec <- function(vterms, n) {
    t <- c()
    if(length(vterms)>0) {
        t <- unique(vterms)
        if(length(t)>n) {
            t <- t[1:n]
        }
    }
    cleanSnoNames(t)
    return(t)
}

##>>.. Takes in a paragraph, strips out stop words, returns a vector of words
para2words <- function (paragraph) {
    paragraph <- gsub("[[:punct:]]", "", paragraph) 
    paragraph <- removeWords(paragraph, stopwords("english"))
    words <- unlist(strsplit(paragraph, " "))
    words <- words[words != ""]  
    return(words)
}

##>>......................................................................<<##
##>>....       Helper functions to process XML or HTML strings          ..<<## 
##>>......................................................................<<##

##>>.. Function to clean an HTML string - remove all HTML tags
cleanFun <- function(htmlString) {
    return(gsub("<.*?>", "", htmlString))
}

