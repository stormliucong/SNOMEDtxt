
### ABOUT ##########################################################
#>>.. Program: Sno2Eng                                               
#>>.. Author: Olga Lyudovyk                                          
#>>.. Symbolic Methods, Spring 2018, DBMI, Columbia University       
#>>.. Purpose:                                                       
#>>.. Convert SNOMED disease concept definitions into readable        
#>>.. fluid English text, and evaluate vs. Google descriptions.      
#>>.. This is the main code to look up and "verablise" SNOMED        
#>>.. concepts. The UI is defined in a separate file.                


#>>.. This is a trimmed down version for the rShiny app

source("Sno2Eng_string_process.R")

#setwd("~/desktop/columbia/symb_methods/project/Sno2Eng/rshiny")
### 1. DATA: read (pre-prepared) SNOMED datasets ##############################################

diseases <- fread("dz_desc.csv")   #file made below - only disease concepts + their descriptions
#>>.. tried looking at only diseases and their relationships, but was missing viruses (HPV, Herpes...)
#>>.. so reverting to the full dataset
relps <- fread("relps.csv")     #only 1'st degree relationships of diseases (disease = source or destination)
concepts <- fread("concepts.csv")  #all concepts, no descriptions - needed for non-disease attributes of diseases
conceptterms <- fread("terms.csv") #concept term descriptions
relp_types <- fread("relp_types.csv")
relp_types$relTypeId <- as.integer(relp_types$relTypeId)


### 2. The Main function ##########################################################

getSno2Eng <- function (searchterm)  {
    snomed <- NULL
    #test# searchterm <- "HPV"
    snoID <- findSNOMED(searchterm)    #>> this can actually return multiple matches
    if (length(snoID) >0 ) {
        snomed <- makeSnomed (snoID[1], searchterm)       #>> current version: only deal with the first match
    }
    return(snomed)
}

#test# 
# res <- getSno2Eng("flu")
# res$story
# 
# getSnoOriginal(res$id)$text

getSno2EngByID <- function (snoID)  {
    snomed <- NULL
    if (length(snoID) >0 ) {
        snomed <- makeSnomed (snoID[1])       #>> current version: only deal with the first match
    }
    return(snomed)
}


### 3. SNOMED concept search   ###########################################################

findSNOMED <- function(searchterm) {
    #>>.. first, in the database
    #>>.. TODO LATER: do a better job with database look up:
    #>>..     1) try looking in UMLS concept names or synonyms
    #>>..     2) do a partial string lookup in SNOMED with TRIE lookup (maybe package ??)
    #>>..     3) worst case: try partial word search and return the best partial match
    searchterm <- capitalize(searchterm)
    snos <- diseases$conceptId[diseases$concept_term==searchterm]
    if(is.null(snos)  | length(snos)<1 ) {
        snos <- searchSnomedCT(searchterm)    #>>.. if not found, via API
    }
    return(snos[1])
}


searchSnomedCT <- function (searchterm) {
    disease_types <- c("mobd", "dsyn", "neop", "virs", "cgab", "inpo", "bact", "inpo", "cgab")
    dat <- getSnomedsViaAPI(searchterm)

    # snomedid <- c()   # if we want multiple IDs, uncomment this and the line starting with #snomedid <- ....
    if(length(dat)>0 & length(dat$Phrases)>0) { # if we got any results
        for(i in 1: length(dat$Phrases$Phrase)) {
            if(length(dat$Phrases$Phrase[[i]]$Mappings$Mapping)>0) {
                for (j in 1:length(dat$Phrases$Phrase[[i]]$Mappings$Mapping)) {
                    #searchterm <- gsub(",", "" , dat$Phrases$Phrase[[i]]$PhraseText[j])  # forgot why I did this
                    mappingDF <- dat$Phrases$Phrase[[i]]$Mappings$Mapping[[j]]
                    if(length(mappingDF)>0) { # found mappings
                        for(k in 1:length(mappingDF)) {  # iter through all mappings - can 1 mapping have >1 candidate?
                            candidate <- mappingDF$MappingCandidates$Candidate[[k]]
                            semtypes <- unlist(candidate$SemTypes$SemType)
                            for(l in 1:length(semtypes)) {
                                if(semtypes[l] %in% disease_types) {
                                    return(candidate$snomedId[l])  # return the first match
                                    # snomedid <- append(snomedid, candidate$snomedId[l]) #this was for returning multiple ids
                                    break
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return (NULL)
}


#' getSnomedsViaAPI: Helper function to call SnomedCT API
#' @param searchterm :  a string
#' @return text from JSON object returned by the API (no processing done here)
#' @example dat <- getSnomedsViaAPI(searchterm)
getSnomedsViaAPI <- function(searchterm) {
    result <- httr::POST("http://snomedct.t3as.org/snomed-coder-web/rest/v1.0/snomedctCodes",
                         body=searchterm, encode="form")
    text <- httr::content(result, "text", encoding="UTF-8")
    data <- jsonlite::fromJSON(text)

    return(data)
}


### 4. SNOMED content organization + verbalization functions  ######################################

#>>.. to do: make 2 snomeds: verbose and normal
#>>.. to do: handle multiple groups
#>>.. to do: make do more with relationships where our concept is the destination, not the source?

#' makeSnomed: The main SNOMED object construction function
#' Gets SNOMED relationships via getSNOMED, organizes into a snomed object, calls the verbalization function (snoStory)
#' uses: getSNOMED, snoStory, <<from Sno2Eng_string_process.R:>> makeDoubleName, cleanSnoNames, pickShortestUnique
#' @param snomedid : SNOMED concept ID (a string of integers or an integer)
#' @param searchterm : string: initial search term
#' @return snomed object (snomed$story is the paragraph we want as the end goal)
#' @example snomed <- makeSnomed (53084003)
makeSnomed <- function(snomedid, searchterm="") {
    snomed<-getSNOMED(snomedid)
    #test#
    #snomed<-getSNOMED(396331005 )  #test  53084003; 733417008: has multiple groups; 9482002: HPV (not dz); 6142004
    #searchterm <-""

    if(nrow(snomed$names) == 0) {  # it's the case when SNOMEDCT API returned a "dead" node - so we didn't find a disease
        snomed <- NULL
        return(snomed)
    }
    snomed$searchterm<-searchterm
    if(length(snomed$searchterm)>0 & snomed$searchterm!="") {
        snomed$longName <- makeDoubleNameFromMany(snomed$searchterm, snomed$names$concept_term)
    } else {
        snomed$longName <- snomed$names$concept_term[1]
    }
    snomed$parents <- cleanSnoNames(unique(snomed$parents$parentName))
    snomed$children <- cleanSnoNames(unique(snomed$children$childName))
    snomed$name <- ifelse(length(snomed$searchterm)>1, snomed$searchterm, snomed$names$concept_term[1])
    rels <- snomed$rels

    #>>.. to do: handle everything with a groupID in a special way - these qualifiers are specific to the group
    # if(length(rels$relGroup[rels$relGroup>1])) { ......  }

    snomed$affectedSites <-cleanSnoNames(pickShortestUnique(rels[rels$relpTypeId == 363698007,]))
    snomed$causedby <-cleanSnoNames(pickShortestUnique(rels[rels$relpTypeId == 246075003,]))
    snomed$assocWith <-cleanSnoNames(pickShortestUnique(rels[rels$relpTypeId == 47429007,]))
    snomed$occur <-cleanSnoNames(pickShortestUnique(rels[rels$relpTypeId == 246454002,]))
    snomed$manifests <-cleanSnoNames(pickShortestUnique(rels[rels$relpTypeId == 363705008,]))
    snomed$pathProcess <-cleanSnoNames(pickShortestUnique(rels[rels$relpTypeId == 370135005,]))
    snomed$tempRelTo <-cleanSnoNames(pickShortestUnique(rels[rels$relpTypeId == 726633004,]))
    snomed$after <-cleanSnoNames(pickShortestUnique(rels[rels$relpTypeId == 255234002,]))
    snomed$morph <-cleanSnoNames(pickShortestUnique(rels[rels$relpTypeId == 116676008,]))
    snomed$clinCourse <-cleanSnoNames(pickShortestUnique(rels[rels$relpTypeId == 263502005,]))
    snomed$discoveredBy <-cleanSnoNames(pickShortestUnique(rels[rels$relpTypeId == 419066007,]))
    snomed$discoverMethod <-cleanSnoNames(pickShortestUnique(rels[rels$relpTypeId == 418775008,]))
    snomed$interprets <-cleanSnoNames(pickShortestUnique(rels[rels$relpTypeId == 363714003,]))
    snomed$hasInterpret <-cleanSnoNames(pickShortestUnique(rels[rels$relpTypeId == 363713009,]))
    snomed$dueTo <-cleanSnoNames(pickShortestUnique(rels[rels$relpTypeId == 42752001,]))
    snomed$severity <-cleanSnoNames(pickShortestUnique(rels[rels$relpTypeId == 246112005,]))
    snomed$discoveredBy <-cleanSnoNames(pickShortestUnique(rels[rels$relpTypeId == 418775008,]))
    snomed$following <-cleanSnoNames(pickShortestUnique(rels[rels$relpTypeId == 255234002,]))
    snomed$during <-cleanSnoNames(pickShortestUnique(rels[rels$relpTypeId == 371881003,]))
    snomed$episodicity <-cleanSnoNames(pickShortestUnique(rels[rels$relpTypeId == 246456000,]))   #deprecated?

    snomed$otherRelated <- pickTopN(snomed$subjRels, 3)
    #>>.. And when we are done - build the actual SNOMED story (i.e. translate into English :)
    snomed$story <- snoStory(snomed)
    return(snomed)
}

#.... todo: ensure that all relationship/concept names returned are of type 900000000000013009, and not 900000000000003001

#' getSNOMED: The database operations function: input is a snomed concept ID (a string of digits)
#' Looks up all relevant info for a snomedID in the DB, and builds an object to contain it
#' Expects global datasets available - see part 1.
#' @param snomedid : SNOMED concept ID (a string of integers or an integer)
#' @return snomed object with id, children, parents, rels, and subjRels
#' @example snomed<-getSNOMED(snomedid)
getSNOMED <- function(snomedid) {
    r_parent <- 116680003  #>>..is-a relationship type id
    syn_typeid <- c("900000000000013009")  #>>..the other type 900000000000003001 is the less readable "fully specified"
    sno <- list()
    sno$id <- snomedid
    #test# sno$id <- 9482002  ##>>.. HPV

    ### There are some strictly not diseases which would be found - like HPV
    ### I will unrestrict from looking only at diseases - for name, parents, and children

    #sno$names <- diseases %>% filter(conceptId == sno$id)
    sno$names <- concepts %>% filter(id==sno$id) %>% inner_join(conceptterms, by=c("id"="conceptId")) %>%
        filter(typeId==syn_typeid) %>% dplyr::select(id, term, typeId) %>% distinct() %>%
        rename(conceptId=id, concept_term=term, concept_typeId=typeId)

    sno$parents <- relps %>% filter(sourceId==sno$id & typeId==r_parent) %>% dplyr::select(destinationId, typeId) %>%
        distinct() %>% inner_join(conceptterms, by=c("destinationId"="conceptId")) %>%
        filter(typeId.y==syn_typeid) %>% dplyr::select(V1, destinationId, term) %>% #arrange(V1)
        group_by(destinationId) %>% slice(1) %>% dplyr::select(destinationId, term)  %>% distinct()
    colnames(sno$parents) <- c("parentId", "parentName")

    sno$children <- relps %>% filter(destinationId==sno$id & typeId == r_parent) %>%
        dplyr::select(sourceId, typeId) %>% distinct() %>%
        inner_join(conceptterms, by=c("sourceId"="conceptId")) %>%
        filter(typeId.y==syn_typeid) %>% dplyr::select(V1, sourceId, term) %>% #arrange(V1)
        group_by(sourceId) %>% slice(1) %>% dplyr::select(sourceId, term)  %>% distinct()
    colnames(sno$children) <- c("childId", "childName")

    #>>.. relationships where the concept at hand is the source, so <our concept> is associated with <something else>
    sno$rels <- relps %>% filter(sourceId==sno$id & typeId != r_parent) %>%
        dplyr::select(destinationId, typeId, relationshipGroup) %>%
        distinct() %>% inner_join(conceptterms, by=c("destinationId"="conceptId")) %>%
        filter(typeId.y==syn_typeid) %>%
        dplyr::select(V1, destinationId, term, typeId.x, relationshipGroup) %>% #arrange(V1)
        group_by(destinationId) %>% slice(1) %>% dplyr::select(destinationId, term, typeId.x, relationshipGroup) %>%
        distinct() %>% rename(nodeId=destinationId, nodeName=term, relpTypeId=typeId.x, relGroup=relationshipGroup)


    #>>.. relationships where the concept at hand is the destination, so <something else> is associated with <our concept>
    #>>.. this query could probably be combined with the one above, but I didn't figure out how to do an "or" in the join
    sno$subjRels <- relps %>% filter(destinationId==sno$id & typeId != r_parent) %>%
        dplyr::select(sourceId, typeId) %>%
        distinct() %>% inner_join(conceptterms, by=c("sourceId"="conceptId")) %>%
        filter(typeId.y==syn_typeid) %>%
        dplyr::select(V1, sourceId, term, typeId.x) %>% #arrange(V1)
        group_by(sourceId) %>% slice(1) %>% dplyr::select(sourceId, term, typeId.x) %>%
        distinct() %>% rename(nodeId=sourceId, nodeName=term, relpTypeId=typeId.x)
    return(sno)
}


#>>.. to do: try different "translation" options for relationships: initially randomize, later select based on user eval
#>>.. to do: order of relationships is also hard to decide, can also randomize and pick the best based on user eval
#>>.. to do: concepts associated with - and other relationships where this concept is a destination, not a source

#' snoStory: The main Verbalization function: Builds the "story":
#' processes all info about the concept, groups similar concepts, applies functions to improve readability
#' input is a snomed object; is called by makeSnomed, which passes this object to it
#' @param sno : snomed object
#' @return snomed object with the "story"
#' @example snomed$story <- snoStory(snomed)
snoStory <- function(sno) {
    #>>.. relationship definitions taken directly from SNOMED documentation:
    #>>.. https://confluence.ihtsdotools.org/display/DOCSTART/6.+SNOMED+CT+Concept+Model

    #>>.. For now phrases are hardcoded. Note different grammatial structure - to reduce the feeling of redundancy.
    #>>.. The order of relationships is hard to decide... Can also be randomized and best options later selected based on
    #>>.. user evaluations.
    snoOut <- paste0(capitalize(sno$longName), " is a kind of ", concatTerms(sno$parents))

    #>>.. |Finding site| specifies the body site affected by a condition.
    if (length(sno$affectedSites)>0 ) {
        ### options for finding site: "is the disease of the "... instead of "that affects"
        snoOut <- paste0(snoOut, " that affects the ", concatTerms(sno$affectedSites), ". ")
    } else {
        snoOut <- paste0(snoOut, ". ")
    }
    
    #>>.. |Has definitional manifestation| links disorders to the manifestations (observations) that define them.
    if(length(sno$manifests)>0 ) {
        snoOut <- paste0(snoOut, "It manifests itself in ", concatTerms(sno$manifests), ". ")
    }

    #>>.. |Associated morphology| specifies the morphologic changes seen at the tissue or cellular level that are characteristic features of a disease.
    if(length(sno$morph)>0 ) {
        snoOut <- paste0(snoOut, "The associated morphology is ", concatTerms(sno$morph), ". ")
    }
    
    #>>.. |Pathological process| provides information about the underlying pathological process for a disorder, but only when the results of that process are not structural and cannot be represented by the |associated morphology| attribute.
    #>>.. TODO the name of the pathological process often already includes "pathological process" - this needs to be rephrased
    if(length(sno$pathProcess)==1 ) {
        snoOut <- paste0(snoOut, "Pathological processes associated with ", sno$name, " is ", concatTerms(sno$pathProcess), ". ")
    } else if(length(sno$pathProcess)>1 ) {
        snoOut <- paste0(snoOut, "Pathological processeses associated with ", sno$name, " are ", concatTerms(sno$pathProcess), ". ")
    }
    
    #>>.. Based on descriptions in online references, subtypes (i.e. examples) are useful in definitions.
    if(length(sno$children)>0 ) {
        #>>.. check if any of the terms are significantly different from the search term - use string distance
        # candidate <- returnMostDiff (sno$name, sno$children)
        # if (candidate != "") {
        #     snoOut <- paste0(snoOut, "An example of ", sno$name, " is ", candidate, ". ")
        # }
        if (length(sno$children)==1) {
            snoOut <- paste0(snoOut, "An example of ", sno$name, " is ", sno$children, ". ")
        } else if (length(sno$children) < 4) {
            snoOut <- paste0(snoOut, "Examples of ", sno$name, " are ", concatTerms(sno$children), ". ")
        } else {
            snoOut <- paste0(snoOut, "Some examples of ", sno$name, " are ", concatTerms(pickTopNVec(sno$children, 3)), ". ")
        }
    }

    #>>.. |Causative agent| identifies the direct causative agent of a disease such as an |organism|, |substance| or |physical force|. (Note: This attribute is not used for vectors, such as mosquitos transmitting malaria).
    if(length(sno$causedby)>0 ) {
        snoOut <- paste0(snoOut, capitalize(sno$name), " is caused by ", concatTerms(sno$causedby), ". ")
    }

    #>>.. |Due to| relates a |clinical finding| directly to a cause such as another |clinical finding| or a |procedure|.
    if(length(sno$dueTo)>0 ) {
        snoOut <- paste0(snoOut, capitalize(sno$name), " occurs due to ", concatTerms(sno$dueTo), ". ")
    }

    #>>.. |Associated with| represents a clinically relevant association between concepts without either asserting or excluding a causal or sequential relationship between the two.
    if(length(sno$assocWith)>0 ) {
        snoOut <- paste0(snoOut, capitalize(sno$name), " is also associated with ", concatTerms(sno$assocWith), ". ")
    }

    #>>.. |Occurrence| refers to a specific period of life during which a condition first presents.
    if(length(sno$occur)>0 ) {
        snoOut <- paste0(snoOut, capitalize(sno$name), " presents in ", concatTerms(sno$occur), ". ")
    }
    
    #>>.. need to combine temporal concepts....
    #>>.. |After| represents a sequence of events where a clinical finding occurs after another |clinical finding| or a |procedure|.
    if(length(sno$during)>0 | length(sno$following)>0  | length(sno$after)>0 ) {
        snoOut <- paste0(snoOut, capitalize(sno$name), " can occur ")
        conj = ""
        if(length(sno$during)>0) {
            snoOut <- paste0(snoOut, " during ", concatTerms(sno$during))
            conj=", "
        }
        if(length(sno$following)>0) {
            snoOut <- paste0(snoOut, conj, " following ", concatTerms(sno$following))
            if(conj!="") {
                conj=", or "
            }
        }
        if(length(sno$after)>0) {
            snoOut <- paste0(snoOut, conj, " after ", concatTerms(sno$after))
        }
        snoOut <- paste0(snoOut, ".")
    }

    if(length(sno$tempRelTo)>0 ) {
        snoOut <- paste0(snoOut, capitalize(sno$name), " can be temporarily related to ", concatTerms(sno$tempRelTo), ". ")
    }
    
    #>>.. |Finding method| specifies the means by which a clinical finding was determined. This attribute is frequently used in conjunction with |finding informer|.
    #>>.. |Finding informer| specifies the person (by role) or other entity (e.g. a monitoring device) from which the clinical finding information was obtained. This attribute is frequently used in conjunction with |finding method|.
    if(length(sno$discoverMethod)>0 | length(sno$discoveredBy)>0) {
        snoOut <- paste0(snoOut, capitalize(sno$name), " is discovered ")
        if(length(sno$discoveredBy)>0) {
            snoOut <- paste0(snoOut, " by ", concatTerms(sno$discoveredBy))
        }
        if(length(sno$discoverMethod)>0) {
            snoOut <- paste0(snoOut, " through ", concatTerms(sno$discoverMethod))
        }
        snoOut <- paste0(snoOut, ".")
    }
    
    #>>.. |Clinical course| represents both the onset and course of a disease.
    if(length(sno$clinCourse)>0 ) {
        snoOut <- paste0(snoOut, "Clinical course is ", concatTerms(sno$clinCourse), ". ")
    }

    #>>.. |Severity| used to sub-class a |clinical finding| concept according to its relative severity.
    if(length(sno$severity)>0 ) {
        snoOut <- paste0(snoOut, "The severity of ", sno$name, " is ", concatTerms(sno$severity), ". ")
    }

    #>>.. |Episodicity| represents episodes of care provided by a physician or other care provider, such as a general practitioner. This attribute is not used to represent episodes of disease experienced by the patient.
    if(length(sno$episod)>0 ) {
        snoOut <- paste0(snoOut, "The episodicity of ", sno$searchterm, " is ", concatTerms(sno$episod), ". ")
    }

    #>>.. |Interprets| refers to the entity being evaluated or interpreted, when an evaluation, interpretation or judgment is intrinsic to the meaning of a concept.
    #>>.. |Has interpretation|, when grouped with the attribute |interprets|, designates the judgment aspect being evaluated or interpreted for a concept (e.g. presence, absence etc.)

    if(length(sno$interprets)>0 | length(sno$hasInterpret)>0) {
        snoOut <- paste0(snoOut, capitalize(sno$name))
        if(length(sno$interprets)>0) {
            snoOut <- paste0(snoOut, " interprets or evaluates ", concatTerms(sno$interprets))
            if(length(sno$hasInterpret)>0) {
                snoOut <- paste0(snoOut, " as ", concatTerms(sno$hasInterpret))
            }
        }
        snoOut <- paste0(snoOut, ".")
    }

    if(nrow(sno$subjRels)>0) {
        snoOut <- paste0(snoOut, " Other related concepts include ", concatTerms(sno$otherRelated), ".")
    }

    #>> this may be convenient for debugging: strings get too long for R display in concole
    # write.csv(snoOut, "output")
    return (snoOut)
}


#' getSnoStory: Get the SNOMED concept description from the snomed object
#' @param snomed : snomed object
#' @return sno$story (a string)
#' @example getSnoStory(eng)   or   print(getSnoStory(sno))
getSnoStory <- function(sno) {
    return(sno$story)
}



### 5. Functions for evaluation #############################################

#>>.. Getting a reference text for the searchterm

#' getMedline: retrieves a disease description from National Library of Medicine XML service
#' @param searchterm : a string
#' @return a string (medline search results, somewhat cleaned of XML tags)
#' @example getMedline("flu")
getMedline <- function (searchterm) {
    #test# searchterm <- "asthma"
    medline <- NULL
    searchterm2 <- gsub(" ", "_", searchterm, fixed=TRUE)  #>> to handle spaces...
    url_lookup <- paste0("https://wsearch.nlm.nih.gov/ws/query?db=healthTopics&term=", searchterm2, "&rettype=brief")
    data =  getURL(url = url_lookup, verbose=F)    #postfields=xml.request, httpheader=myheader ??
    xmltext  <- xmlTreeParse(data, useInternalNodes=T)
    #str <- unlist(xpathApply(xmltext,'//nlmSearchResult/list[1]/document[1]/content[4]',xmlValue))
    top = xmlRoot(xmltext)
    nodes = getNodeSet(top,"//nlmSearchResult/list[1]/document[1]/content[@name='FullSummary']")[[1]]
    if(length(nodes)>0) {
        str <- xmlSApply(nodes, xmlValue)
        medline <- unlist(cleanFun(str))
    }
    return(medline[[1]])  # not sure why unlist doesn't take care of this...
}


#>>.. Making a reference (un-processed) SNOMED dump


#' getSnoOriginal: retrieves all SNOMED relationships (preferred term) and all SNOMED names for a concept as they are
#' @param snoID : a string of digits or an integer that is a SNOMED concept ID
#' @return snoOr object: a list
#' @example getSnoOriginal(408856003)
getSnoOriginal <- function (snoID) {
    if(is.null(snoID) | snoID < 1) {
        return (NULL)
    }
    snoOr <- list()
    #test#  snoID <- 408856003

    prefType <- c("900000000000003001")  ##<< "preferred term" type only -  more verbose but yet more correct
    snoOr$id <- snoID
    snoOr$names <- conceptterms %>% filter(conceptId==snoID) %>% dplyr::select(typeId, term)

    snoOr$rels <- relps %>% filter(sourceId==snoID) %>%
        inner_join(conceptterms, by=c("destinationId"="conceptId"))  %>%
        filter(typeId.y==prefType) %>%
        dplyr::select(relationshipGroup, destinationId, typeId.x, typeId.y, term) %>% distinct()  %>%
        inner_join(relp_types, by=c("typeId.x"="relTypeId")) %>% filter(nameType==prefType) %>%
        select(relationshipGroup,destinationId,typeId.x, relType,term) %>% distinct %>%
        rename(relGroup=relationshipGroup, destinationId=destinationId, relTypeId=typeId.x, relName=relType, term=term)

    snoOr$subjRels <- relps %>% filter(destinationId==snoID) %>%
        inner_join(conceptterms, by=c("sourceId"="conceptId"))  %>%
        filter(typeId.y==prefType) %>%
        dplyr::select(relationshipGroup, sourceId, typeId.x, typeId.y, term) %>% distinct() %>%
        inner_join(relp_types, by=c("typeId.x"="relTypeId")) %>% filter(nameType==prefType) %>%
        select(relationshipGroup,sourceId,typeId.x, relType,term) %>% distinct %>%
        rename(relGroup=relationshipGroup, sourceId=sourceId, relTypeId=typeId.x, relName=relType, term=term)

    snoOr$text <- makeSnoOriginalText(snoOr)
    return(snoOr)
}

#' makeSnoOriginalText: concatenates all "preferred" original SNOMED terms/relationships concatenated for a concept
#' used by getSnoOriginal
#' @param snoOr object: a list
#' @return snoOut : string
#' @example makeSnoOriginalText(snoOr)
makeSnoOriginalText <- function (snoOr) {
    snoOut <- ""
    #test# snoID <- 408856003    snoOr <- getSnoOriginal(snoID)

    prefType <- c("900000000000003001")
    snoOut <- paste0("ConceptID: ", snoOr$id, ". Terms: ", paste(snoOr$names$term, collapse=", "), ". ")

    if(nrow(snoOr$rels)>0) {
        snoOut <- paste0(snoOut, "Relationships: ")
        for(i in 1:length(snoOr$rels$relName)) {
            snoOut <- paste0(snoOut, snoOr$rels$term[i], " = ", snoOr$rels$relName[i], ". ")
        }
        #snoOut <- paste0(snoOut, ".")
    }
    if(nrow(snoOr$subjRels)>0) {
        snoOut <- paste0(snoOut, " Related concepts: ")
        for(i in 1:length(snoOr$rels$relName)) {
            snoOut <- paste0(snoOut, snoOr$subjRels$term[i], " - ", snoOr$subjRels$relName[i], ". ")
        }
        #snoOut <- paste0(snoOut, ".")
    }

    return(snoOut)
}


