

#' Setup prediction scenarios
#'
#' @param myPheno data frame containing at least the following columns: gid (genotype identifier), env, value (trait phenotype)
#' @param scenario one of knEnv, knLoc.knYr, knLoc.nYr, nLoc.knYr, nLoc.nYr.
#' knEnv is known environment, including a sample of observations from the testing environment, as well as all other environments, knLoc.knYr is known location, known year, knLoc.nYr is known location / unknown year, nLoc.knYr is unknown location, known year, nLoc.nYr is unknown location / unknown year.
#' @param envs.train (optional) restriction of environments to be in training set, if not provided all environments are considered
#' @param envs.pred (optional) restriction of environments to be in testing set, if not provided all environments are considered
#' @param ignore.genos (optional) character vector of genotypes to exclude from testing set, while including them in training set, typically check genotypes
#' @param traits (optional) character vector for trait column(s)
#' @param genos (optional) character vector of restricted genotype list, typically genotyped lines.
#' @param prop.CVS1 numeric, proportion of genotypes in testing set for scenario knEnv, default is 0.8.
#'
#' @return list with one lement per testing (environment) set.
#' Within each testing set, list of training environment, testing environment, genotypes in training and testing sets, and phenotypic data in training and testing sets.
#' @author Charlotte Brault
#' @export

setup_scenarios <- function(myPheno, scenario, envs.train=NULL, envs.pred=NULL,
                            ignore.genos=NULL,traits=NULL, genos=NULL,prop.CVS1=0.8){
  # envs.train=NULL
  # envs.pred <- envs2pred
  # myPheno = pheno.env
  # scenario ="knLoc.knYr"
  # ignore.genos=c("Wheaton","Oslo","Bacup")
  # genos=rownames(GRM)

  out <- list()

  ## if list of genotypes is provided, restrict the training and testing sets to that list
  if(!is.null(genos)){
    cmon.geno <- intersect(genos, unique(as.character(myPheno$gid)))
    myPheno <- dplyr::filter(myPheno, gid %in% cmon.geno) %>% droplevels()
  }

  ## 1. determine environments in training and testing sets for all scenarios
  all_env <- unique(myPheno$env)
  if(is.null(envs.train)) envs.train <- all_env else envs.train <- intersect(envs.train, all_env)
  if(is.null(envs.pred)) envs.pred <- all_env else envs.pred <- intersect(envs.pred, all_env)
  stopifnot(scenario %in% c("knLoc.knYr","knLoc.nYr","nLoc.knYr","nLoc.nYr","knEnv"))
  all_env <- unique(c(envs.train, envs.pred))
  ## create a character vector with all environments in testing set
  tsn <- sapply(envs.pred, function(x) stringr::str_subset(all_env, x))

  if(scenario %in% c("knLoc.knYr","knEnv")){
    ## remove environment to be predicted
    trn <- lapply(tsn, function(x) intersect(setdiff(all_env, x), envs.train))
  } else if(scenario %in% "knLoc.nYr"){
    ## remove environments from same year
    trn <- lapply(tsn, function(x){
      envs <- stringr::str_subset(all_env, unlist(lapply(strsplit(x,"_",fixed=T),`[`,2)), negate=TRUE)
      return(intersect(envs, envs.train))
    })

  } else if(scenario %in% "nLoc.knYr"){
    ## remove environments from same location
    trn <- lapply(tsn, function(x) {
      envs <- stringr::str_subset(all_env, unlist(lapply(strsplit(x,"_",fixed=T),`[`,1)), negate=TRUE)
      return(intersect(envs, envs.train))
    })

  }else if(scenario %in% "nLoc.nYr"){
    trn <- lapply(tsn, function(x) {
      ## remove environments from same location and year
      sub1 <- stringr::str_subset(all_env, unlist(lapply(strsplit(x,"_",fixed=T),`[`,1)), negate=TRUE)
      sub1 <- stringr::str_subset(sub1, unlist(lapply(strsplit(x,"_",fixed=T),`[`,2)), negate=TRUE)
      return(intersect(sub1, envs.train))
    })
  }

  ## determine gid in training and testing sets for all scenarios
  ### genotypes in testing set = genotypes in environment to be predicted
  genos.tsn <- purrr::map(names(tsn), function(x)
    unique(as.character(myPheno$gid[myPheno$env %in% x])))
  ## for check genotypes, don't predict them but keep them in training set
  if(!is.null(ignore.genos)){
    genos.tsn <- purrr::map(genos.tsn, function(x) setdiff(x, ignore.genos))
  }
  ## for scenario knEnv, take some lines in testing environment to be in training set
  if(scenario == "knEnv"){
    set.seed(54446)
    genos.tsn <- purrr::map(genos.tsn, function(x)
      sample(x, size=round(prop.CVS1*length(x))))
  }

  ### genotypes in training set=genotypes in training set environments
  genos.trn <- purrr::map2(.x=trn, .y=genos.tsn, \(x, y){
    all.genos <- unique(as.character(myPheno$gid[myPheno$env %in% x]))
    setdiff(all.genos, y)
  })

  ## output: subset phenotypic data for training and testing set
  if(is.null(traits)) {traits <- colnames(myPheno)[sapply(myPheno, is.numeric)]}
  #myPheno <- myPheno %>% dplyr::select(tidyselect::all_of(c("env", "gid", traits)))
  ### for training set
  pheno.train <- purrr::map2(.x=trn, .y=genos.tsn, \(x, y){
    ## select gid for training environments
    sel <- filter(myPheno, env %in% x)
    ## set as NA phenotypic values for genotypes in testing set
    sel[sel$gid %in% y, traits] <- NA
    return(sel)
  })

  ### for testing set
  pheno.test <- purrr::map2(.x=tsn, .y=genos.tsn, \(x, y){
    ## phenotypic data from testing environments and gids
    filter(myPheno, env %in% x & gid %in% y) %>% droplevels()
  })

  ## return all outputs in one list
  out.list=list(envs.trn=trn, envs.tsn=tsn,
                genos.trn=genos.trn, genos.tsn=genos.tsn,
                pheno.trn=pheno.train, pheno.tsn=pheno.test)
  out.list <- purrr::map(out.list, setNames, unlist(tsn))
  return(out.list)

}







#' Get cross-validation folds
#'
#' @param index length of total set to select from, numeric vector of length 1
#' @param nb.folds number of fold to divide, default is 5
#' @param seed (optional) for reproducibility, set a seed.
#'
#' @return list of testing set indexes for each fold.
#' @author Charlotte Brault
#' @export

getFolds <- function(index, nb.folds=5, seed=NULL){
  stopifnot(nb.folds <= index,
            ifelse(! is.null(seed), is.numeric(seed), TRUE))
  if(!is.null(seed)){set.seed(seed, kind="Mersenne-Twister")}
  split <- split(sample(1:index), rep(1:nb.folds, length=index))
  names(split) <- paste0("Fold", seq(1:nb.folds))
  return(split)
}


#' Format accession names
#'
#' @param names character vector of accession names
#' @param toupper logical, set to upper case, default is TRUE
#' @param dashToUnderscore logical, transform dash to underscore, default is FALSE
#' @param pedigree is character vector a pedigree? if TRUE, don't remove slash strings
#' @param verbose logical, default to FALSE
#'
#' @importFrom stringr str_remove_all str_replace_all
#' @return list of initial and corrected names
#' @export

format_names <- function(names, toupper = TRUE, dashToUnderscore = FALSE,
                         pedigree = FALSE, verbose = FALSE) {
  names_init <- names

  # Define special Unicode characters explicitly
  u2019 <- intToUtf8(0x2019) # ’
  u00A7 <- intToUtf8(0x00A7) # §
  u03C6 <- intToUtf8(0x03C6) # φ
  u00B6 <- intToUtf8(0x00B6) # ¶
  u2014 <- intToUtf8(0x2014) # —

  ## patterns to remove
  pat <- paste0(" \\(S check\\)| \\(MR check\\)| = ALVORADA|'S\"|",
                u2019, "S", u2019, "|'S'|", u2019, "s", u2019,
                "|\\*|", u00A7)
  names1 <- stringr::str_remove_all(names_init, pattern = pat)

  # Remove trailing special symbols
  pat_end <- paste0("[_", u00A7, u03C6, "\\*", u00B6, "]+$")
  names2 <- stringr::str_remove_all(names1, pattern = pat_end)

  # Remove single quotes
  names3 <- stringr::str_remove_all(names2, pattern = "'")

  # Remove all spaces
  names4 <- stringr::str_remove_all(names3, pattern = " ")

  # Trim whitespace just in case
  names5 <- stringr::str_trim(names4)

  # Remove leading underscore
  names6 <- stringr::str_remove(names5, "^_")

  # Remove trailing underscores
  names7 <- stringr::str_remove(names6, "_$")
  names8 <- stringr::str_remove(names7, "_$")

  ## Replace dash / long dash
  if (dashToUnderscore) {
    names9 <- stringr::str_replace_all(names8, u2014, "_")
    names10 <- stringr::str_replace_all(names9, "-", "_")
  } else {
    names10 <- stringr::str_replace_all(names8, u2014, "-")
  }

  ## Replace slash by underscore unless pedigree
  if (!pedigree) {
    names11 <- stringr::str_replace_all(names10, "/", "_")
  } else {
    names11 <- names10
  }

  ## Set to upper case
  if (toupper) {
    names11 <- toupper(names11)
  }

  ## Output initial and corrected names
  out <- data.frame(initial.name = names_init, corrected.names = names11)
  if (verbose) {
    cat(length(unique(names_init)), "initial unique names and",
        length(unique(names11)), "final unique names\n")
  }

  return(list(init = names_init, corrected = names11, df = out))
}


#' Concordance of names between reference and match data
#' Matching of names will be done based on the spelling of the names, without punctuation.
#'
#' @param match.name character vector of names to match (format)
#' @param ref.name character vector of reference names. Default is NULL but it needs to be provided if corresp.list is NULL.
#' @param allowNoMatch logical, keep all values from match.name even when there is no match, default is TRUE
#' @param corresp.list named list of correspondence between names, with correct spelling as name and possible synonym as character vector, default is NULL.
#' If provided, add new correspondence to the list. If there is a conflict between ref.names and corresp.list, the name in corresp.list will be used.
#' @param verbose integer, level of verbosity, default is 1
#'
#' @return a list of 3 elements: the correspondence data frame, the updated correspondence list, and a vector of new names to be applied to match name.
#' @author Charlotte Brault
#' @export

concordance_match_name <- function(match.name,ref.name=NULL, allowNoMatch=TRUE,
                                   corresp.list=NULL, verbose=1){
  stopifnot(!all(is.null(ref.name) & is.null(corresp.list)))
  ref.name <- unique(ref.name)
  input.name <- match.name ## for output of the same length
  match.name <- unique(match.name)
  if(!is.null(ref.name)){
    df.ref <- data.frame(ref.name=ref.name,
                         noPunct=gsub("[[:punct:]]","",ref.name),
                         stringsAsFactors = F)
  } else {
    ## create df.ref based on corresp.list
    df.ref <- data.frame(ref.name=unique(names(corresp.list)),
                         noPunct=gsub("[[:punct:]]","",unique(names(corresp.list))))
  }
  df.ref <- dplyr::distinct(df.ref)

  df.match <- data.frame(match.name=match.name,
                         format.name=format_names(match.name)$corrected,
                         stringsAsFactors = F)
  df.match$noPunct <- gsub("[[:punct:]]","",df.match$format.name)
  df.match <- dplyr::distinct(df.match)

  ## merge match and ref based on name without punctuation
  df.merged <- merge(df.ref, df.match, by.x="noPunct", by.y="noPunct",
                     all.x=F, all.y=allowNoMatch)
  if(verbose > 0){
    print(paste0("Number of unique genotype names in reference data: ", nrow(df.ref)))
    print(paste0("Number of unique genotype names in match data: ", nrow(df.match)))
    #print(paste0("Number of genotype names in reference data and match data: ", nrow(df.merged)))
    print(paste0("Number of matches: ",sum(!is.na(df.merged$ref.name))))
  }

  ## add a column with new name to be applied to match
  df.merged$new.name <- dplyr::coalesce(df.merged$ref.name, df.merged$format.name)

  ## verify that there is one name for a given spelling

  # tmp <- distinct(df.merged[,c("new.name","noPunct")])
  # if(length(unique(tmp$new.name)) != length(unique(tmp$noPunct))){
  #   dup <- tmp$noPunct[duplicated(tmp$noPunct)]
  #   warning(print(paste0("Error several names for the same spelling: ", dup,": ",
  #                        paste(tmp$new.name[tmp$noPunct %in% dup], collapse=", "))))
  #
  # }

  ## create a list of to store genotype names correspondence
  if(is.null(corresp.list)){
    corresp.list <- list(df.merged$ref.names)
  } else {
    ## for existing correspondence, add to new name the corresp name
    for(nam in df.merged$new.name){
      ## look in corresp.list if this genotype exists
      match <- names(which(purrr::map_lgl(corresp.list, \(x)
                                          length(intersect(x,gsub("[[:punct:]]","",nam))) !=0)))
      if(length(match) ==1){
        df.merged$new.name[df.merged$new.name == nam] <- match
      } else if(length(match) > 1){
        warning(paste0("Several names in corresp.list for ",nam,": ",match))
      }
    } ## end for nam
  } ## end if/else


  # ## remove from corresp.list elements that are similar to the name
  # corresp.list <- Filter(function(a) any(!is.na(a)),
  #                        purrr::map2(.x=corresp.list, .y=names(corresp.list), \(x,y) x=x[x!=y]))


  for(nam in df.merged$new.name){
    corresp.list[[nam]] <- unique(c(corresp.list[[nam]],
                                    unlist(df.merged[df.merged$new.name == nam,
                                                     c("noPunct","format.name","match.name")])))

  }

  ## remove from corresp.list elements that are similar to the name
  # corresp.list <- Filter(function(a) any(!is.na(a)),
  #                        purrr::map2(.x=corresp.list, .y=names(corresp.list), \(x,y) x=x[x!=y]))

  # if(any(duplicated(unlist(corresp.list)))){
  #   warning(paste0("Duplicated names in corresp.list, please verify: ",
  #                        unlist(corresp.list)[duplicated(unlist(corresp.list))]))
  # }

  corresp.list <- corresp.list[order(names(corresp.list))]
  new.names <- plyr::mapvalues(df.merged$match.name[match(input.name,df.merged$match.name)],
                               from=df.merged$match.name,
                               to=df.merged$new.name)
  df.changes <- df.merged %>% dplyr::filter(new.name != match.name) %>%
    dplyr::select(tidyselect::all_of(c("match.name","new.name")))

  return(list(df.corresp=df.merged,
              df.changes=df.changes,
              corresp.list=corresp.list,
              new.names=new.names))
}



#' Format phenotypic data from GrainGenes (Excel tables)
#'
#' Load Excel files containing phenotypic data from GrainGenes, from multiple locations and years.
#' Combine them into one data frame and separate genotype information with phenotypic data.
#'
#' @param p2d path to directory where tables are saved
#' @param years numeric vector of years to look for
#' @param locs character vector of tab names to look for (including Entry) or to location names to identify the trial.
#' If several names are corresponding to one trial, repeat the different versions in the vector and add the final name as vector name.
#' @param traits character vector of trait names to look for.
#' If several names are corresponding to one trait, name the different versions and use as vector name the sought version.
#' for example: traits=c("VSK","Heading","FDK"); names(traits)=c("VSK","HD","VSK")
#' @param cols2rem character vector of column names to remove, to avoid bad matching
#' @param distMatchTrait numeric value, distance for string matching. Default is 8.
#' Increased distance would lead to more matching and is more prone to errors.
#'
#' @seealso [format_names()]
#'
#' @return list of 4 components:
#' - var.match.info: data frame of variable matching
#' - sheet.match.info: data frame of sheet matching (finding the relevant tabs)
#' - phenot: data frame of combined phenotypic data for all years and locations
#' - entry.info: data frame of combined genotype information
#'
#' @author Charlotte Brault
#' @export

format_phenot <- function(p2d, years, locs, traits, cols2rem=NULL,distMatchTrait=8){

  if(is.null(names(traits))) names(traits) = traits
  if(is.null(names(locs))) names(locs) = locs

  df.match.all <- list()

  list.files <- list.files(p2d)
  ## restrict to tables of phenotypic data by searching "table" word in file name
  list.files <- list.files[grep("table",x=list.files,ignore.case = TRUE)]
  if(length(list.files) == 0){
    print("Not files found")
  }

  ## restrict to the years under study
  list.files2 <- c()

  for(yr in year_ranges){
    idx <- grep(as.character(yr),list.files)
    if(length(idx) < 1){
      print(paste0("No file found for year ",yr))
    } else if(length(idx) > 1){
      print(paste0("Several files found for year ",yr))
    } else {
      list.files2 <- c(list.files2,list.files[idx])
      names(list.files2)[length(list.files2)] <- yr
    }
  }

  yrs <- names(list.files2)
  ## retrieve the whole file paths
  list.files2 <- paste0(p2d,"/",list.files2)
  stopifnot(all(file.exists(list.files2)))
  names(list.files2) <- yrs

  load.yrs <- list() ;entry.list <- list()
  clean.phenot <- list() ; df.match.sheet <- list()

  ## -- For loop over years -- ##
  for(yr in names(list.files2)){
    print(paste0("Working on year ",yr))
    ## extract sheet names
    sheets <- readxl::excel_sheets(path=list.files2[yr])

    ## handle case where sheet name doesn't contain the location name
    ### find tables with phenotypic data or entry info
    if(all(grepl("Table",sheets))){
      print("Need to associate each table with location")
      sheets2locs <- c() ; sheets2keep <- c()
      for(sh in sheets){
        suppressMessages(tmp <- readxl::read_excel(path=list.files2[yr], sheet=sh))
        ## check if a location name is present in the first row
        detect <- stringr::str_detect(string=colnames(tmp)[1],pattern=locs)
        ## second test to detect tables for MAS, summaries (to exclude)
        detect2 <- stringr::str_detect(
          string=colnames(tmp)[1],
          pattern=c("USDA-ARS","rust","Molecular markers",
                    "correlation","summary","mean"))
        if(any(detect) & !any(detect2)) {
          sheets2keep <- c(sheets2keep,sh)
          sheets2locs <- c(sheets2locs,names(locs)[which(detect)])
        }
      }
      ### case where sheet have location names
    } else {
      print("Take sheet name for location information")
      ## Select sheets based on name
      sheets2keep <- c() ; sheets2locs <- c()
      for(sh in sheets){
        idx <- which(!is.na(charmatch(c(locs), sh)))
        if(length(idx)> 0){
          sheets2keep <- c(sheets2keep,sh)
          sheets2locs <- c(sheets2locs, c(names(locs)[idx]))
        }
      }
    }

    stopifnot(length(sheets2keep) == length(sheets2locs))
    names(sheets2keep) <- sheets2locs
    df.match.sheet[[yr]] <- data.frame(year=yr,sheet=sheets2keep,loc=sheets2locs)
    # print(paste0("Establish correspondance between",
    #              " Excel sheets and location:\n"))
    # print(data.frame(sheet=sheets2keep,
    #                  loc=sheets2locs))

    ## Load all selected sheets in one element per year and location
    suppressMessages(load.all <- lapply(sheets2keep, function(x)
      readxl::read_excel(path=list.files2[yr], sheet=x, skip=2, na="-")))
    load.yrs[[yr]] <- load.all


    ## --- Now load and curate data within tabs --- ##

    ## clean each year

    if(!("Entry" %in% names(load.yrs[[yr]]))){
      print(paste0("Entry tab not found for ",yr))
    }
    ## for each location
    for(tab in names(load.yrs[[yr]])){

      dat <- as.data.frame(load.yrs[[yr]][[tab]])

      ## handle genotype info
      if(tab == "Entry"){
        ## detect first row
        first_row <- which(dat[,1] == 1)
        dat <- dat[first_row:nrow(dat),]


        ### try to guess order of columns first year, market class, submitter and organization

        #### organization
        # define pattern to look for
        pat.orga <- "\\bMN\\b|\\bND\\b|\\bSD\\b|\\bUMN\\b|\\bSDSU\\b|\\bNDSU\\b|\\bAAFC\\b"
        # find the column number corresponding to the most matched
        idx.orga <- which.max(apply(dat, 2, function(x)
          sum(stringr::str_detect(stats::na.omit(x), pattern=pat.orga))))

        #### first year
        # define pattern to look for
        pat.yr <- "\\b[0-9]{4}\\b"
        # find the column number corresponding to the most matched
        app <- apply(dat, 2, function(x) sum(grepl(pattern="^19|20",x)))
        idx.yr <- ifelse(max(app) < 10, NA, which.max(app))

        #### submitter
        # define pattern to look for
        pat.sub <- "Anderson|Glover|Mergoum|Smith|Cooper"
        # find the column number corresponding to the most matched
        app <- apply(dat, 2, function(x) sum(stringr::str_detect(stats::na.omit(x), pattern=pat.sub)))
        idx.sub <- ifelse(max(app) < 8, NA, which.max(app))

        #### market class
        # define pattern to look for
        pat.class <- "HRS|HRSW|Durum|CWRS"
        # find the column number corresponding to the most matched
        app <- apply(dat, 2, function(x) sum(stringr::str_detect(stats::na.omit(x), pattern=pat.class)))
        idx.class <- ifelse(max(app) < 8, NA, which.max(app))

        colOrder <- c(1:3,idx.orga)
        if(!is.na(idx.yr)) colOrder <- c(colOrder, idx.yr)
        if(!is.na(idx.sub)) colOrder <- c(colOrder, idx.sub)
        if(!is.na(idx.class)) colOrder <- c(colOrder, idx.class)
        dat <- dat[,colOrder]
        colnames(dat)[1:4] <- c("Entry","Line","Pedigree","Organization")
        if(ncol(dat) > 4) colnames(dat)[5] = "FirstYear"
        if(ncol(dat) > 5) colnames(dat)[6] = "Submitter"
        if(ncol(dat) == 7) colnames(dat)[7] = "MarketClass"
        ## remove lines with missing line name
        dat <- dat[!is.na(dat$Line),]

        # ## correct column names
        # if(ncol(dat) == 7){
        #   colnames(dat) <- c("Entry","Line","Pedigree","MarketClass",
        #                      "FirstYear","Submitter","Organization")
        # }else if(ncol(dat) == 6){
        #   colnames(dat) <- c("Entry","Line","Pedigree","FirstYear",
        #                      "Submitter","Organization")
        # } else if(ncol(dat) ==5){
        #   colnames(dat) <- c("Entry","Line","Pedigree","FirstYear","Organization")
        # }


        entry.list[[yr]] <- dat
        entry.list[[yr]]$Line <- format_names(entry.list[[yr]]$Line)$corrected

        ## curate phenotypic data
      } else{
        ## select first column
        tmp <- suppressWarnings(stats::na.omit(as.numeric(dat[,1])))

        ## detect first row with genotype names
        if(any(grepl("Line",dat[,1]))){
          ### clean column names
          colN <- which(dat[,1] == "Line")
          colnames(dat) <- stringr::str_remove(paste0(colnames(dat),dat[colN,]),"NA")
          colnames(dat) <- stringr::str_remove(colnames(dat),"ppm|g")
          colnames(dat) <- janitor::make_clean_names(colnames(dat),case="none",
                                                     replace=c("%"=""))
          dat <- dat[(colN+1):nrow(dat),]

          ## case when first column is composed of entry number

        } else if(length(tmp) > 10){

          #print(paste0(yr,", ",tab))
          dat[,1] <- NULL
          ### clean column names
          if(!"Line" %in% colnames(dat)){
            #print(paste0("No column names with Line for tab ",tab," for year ",yr))
            colN <- which(dat[,1] == "Line")
            colnames(dat) <- stringr::str_remove(paste0(colnames(dat),dat[colN,]),"NA")
            colnames(dat) <- stringr::str_remove(colnames(dat),"ppm|g")
            colnames(dat) <- janitor::make_clean_names(colnames(dat),case="none",
                                                       replace=c("%"=""))
            dat <- dat[(colN+1):nrow(dat),]
          }

        }

        colnames(dat)[c(1)] <- c("Line")
        ## handle the case when X + number is prefixed to column name
        colnames(dat)[2:ncol(dat)] <- gsub(pattern="X[1-9]",
                                           x=colnames(dat)[2:ncol(dat)],
                                           replacement="")
        ## remove summary lines at the bottom
        dat <- dat[!is.na(dat$Line),]
        # remove first and summary rows
        pat.sum <- "CV|C.V|LSD|LINE|MEAN|Mean|MAXIMUM|Maximum|MINIMUM|Minimum|AVERAGE|Average"
        idx <- grep(pattern = pat.sum, dat$Line)
        if(length(idx) > 0) {dat <- dat[-idx,]}

        ## add year and location information
        dat$Year <- as.character(yr)
        dat$Location <- as.character(tab)
        dat <- dat %>% dplyr::relocate(c("Year","Location"))

        ## set as numeric trait columns
        dat = suppressWarnings(dat %>%
                                 ## first 3 columns should be year,loc, entry
                                 dplyr::mutate_at(c(1:3), as.factor) %>%
                                 dplyr::mutate_at(c(4:ncol(dat)), as.numeric))
        dat <- janitor::remove_empty(dat, which="cols")
        ### print(head(dat))
        ### print(dim(dat))
        ## remove traits with bad matching
        if(!is.null(cols2rem)){
          dat[,cols2rem] <- NULL
        }

        ## match trait names with dat table
        traits_tmp = colnames(dat)[4:ncol(dat)]

        ## amatch is a function from stringdist package to find the closest match
        idx = stringdist::amatch(x=janitor::make_clean_names(traits_tmp,use_make_names = FALSE),
                                 table=janitor::make_clean_names(traits),
                                 method="lcs", maxDist=distMatchTrait)

        ## verify no duplicated traits
        if(length(unique(stats::na.omit(idx))) != length(stats::na.omit(idx))){
          print(paste0("Problem with trait matching"))
        }
        #stopifnot(length(unique(idx)) == length(idx))
        # print(paste0(length(idx[!is.na(idx)])," traits found in ",tab," ",yr,
        #              "\nCheck matching:"))
        ## verify the matching
        df.match <- data.frame(year=yr, location=tab,
                               initial_name=traits_tmp[which(!is.na(idx))],
                               match=traits[stats::na.omit(idx)],
                               final_name=names(traits)[stats::na.omit(idx)],
                               match_problem=ifelse(length(unique(stats::na.omit(idx))) !=
                                                      length(stats::na.omit(idx)),TRUE,FALSE))
        #print(df.match)
        df.match.all[[yr]][[tab]] <- as.data.frame(df.match)
        ## missing matches
        if(any(is.na(idx))){
          print(paste("The following traits won't be matched and be discarded:",
                      paste(traits_tmp[which(is.na(idx))],collapse=", "),
                      " for ",tab," and year ",yr))
          ## remove traits not matched
          dat[,traits_tmp[which(is.na(idx))]] <- NULL
        }

        # replace new trait names
        colnames(dat) <- plyr::mapvalues(colnames(dat),
                                         from=df.match$initial_name,
                                         to=df.match$final_name,warn_missing = FALSE)

        ## detect rows filled with NA values and remove them
        rowsNA <- which(apply(dat[,4:ncol(dat), drop=FALSE],1,function(x) all(is.na(x))))
        if(length(rowsNA) >0){
          dat <- dat[-rowsNA,]
        }

        ## check genotype name with custom R function
        dat$Line <- format_names(dat$Line)$corrected
        ## if lines have been phenotyped but are not in the entry table, add them
        sd <- setdiff(as.character(dat$line),  entry.list[[yr]]$Line)
        if(length(sd) > 0){
          print("New genotypes added to the entry list")
          entry.list[[yr]] <- dplyr::bind_rows(entry.list[[yr]],Line=sd)
        }

        ## check trait boundaries?
        ## save cleaned table in a list for final combining
        clean.phenot[[yr]][[tab]] <- as.data.frame(dat)

      } # end if/else tab = Entry

    } # end for tab (loc)
    df.match.all[[yr]] <- do.call(plyr::rbind.fill,df.match.all[[yr]])
    clean.phenot[[yr]] <- do.call(plyr::rbind.fill,clean.phenot[[yr]])
  } ## end for year

  ## format outputs
  match.info <- do.call(plyr::rbind.fill, df.match.all)
  match.info.sheets <- do.call(plyr::rbind.fill, df.match.sheet)
  clean.phenot <- do.call(plyr::rbind.fill, clean.phenot)
  entry.info <- do.call(plyr::rbind.fill, entry.list)
  return(list(var.match.info=match.info,
              sheet.match.info=match.info.sheets,
              phenot=clean.phenot,
              entry.info=entry.info))
}


#' Format and curate vcf file
#'
#' @param vcf.p2f path to the vcf file
#' @param matrix.gt matrix of genotypes, must contains CHROM and POS columns, default is NULL
#' @param p2f.export.vcf path to export curated vcf file, default is NULL (no export), useful for further Beagle imputation.
#' @param IDnum logical, remove the prefix idxx. from the genotype name, default is FALSE
#' @param mrk.info data frame with marker information pulled from the vcf, with columns CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, default is NULL. Needed if `vcf.p2f` is not provided and `p2f.export.vcf` is not null.
#' @param corresp.geno.name table of correspondence between genotype names, with the first column being the name to be matched and the second one being the updated name.
#' @param remove.chrUkn logical, remove markers from unknown chromosome, default is TRUE
#' @param check.mrk.dups logical, check for duplicated markers, default is TRUE
#' @param check.geno.dups logical, check for duplicated genotypes, default is TRUE
#' @param remove.nonPolyMrk logical, remove non-segregating markers, default is TRUE
#' @param tresh.heterozygous numeric, threshold of the maximum proportion of heterozygous genotypes, default is NULL
#' @param thresh.NA.ind numeric, threshold of the maximum proportion missing values for individuals, default is 0.5
#' @param thresh.NA.mrk numeric, threshold of the maximum proportion of missing values for markers, default is 0.5
#' @param format.names logical, use a custom function to transform genotype names, default is FALSE
#' @param concordance.function logical, use a custom function to match genotype names, default is FALSE
#' @param corresp.geno.list named list of correspondence between genotype names, with correct spelling as name and possible synonym as character vector, default is NULL.
#' @param imputation character, impute missing values with kNNI or Beagle or don't impute, default is NULL (no imputation)
#' @param p2f.beagle path to export file for Beagle imputation, default is NULL
#' @param thresh.MAF numeric, threshold of minor allele frequency, default is NULL
#' @param verbose numeric, level of verbosity, default is 1
#'
#' @return data frame with 0/1/2 values from the curated vcf file, with genotypes in row and markers in columns.
#' @author Charlotte Brault
#' @seealso [format_names()], [getGenoTas_to_DF()], [concordance_match_name()]
#' @importFrom vcfR read.vcfR extract.gt getFIX
#' @importFrom methods new
#' @importFrom dplyr arrange filter mutate select distinct
#' @importFrom tibble as_tibble

#' @export

format_curate_vcf <- function(vcf.p2f=NULL,
                              matrix.gt=NULL,
                              mrk.info=NULL,
                              corresp.geno.name=NULL,
                              p2f.export.vcf=NULL,
                              IDnum=FALSE,
                              remove.chrUkn=TRUE,
                              check.mrk.dups=TRUE,
                              remove.nonPolyMrk=TRUE,
                              tresh.heterozygous=NULL,
                              check.geno.dups=TRUE,
                              thresh.NA.ind=0.5,
                              thresh.NA.mrk=0.5,
                              format.names=FALSE,
                              concordance.function=FALSE,
                              corresp.geno.list=NULL,
                              imputation=NULL,
                              p2f.beagle=NULL,
                              thresh.MAF=NULL,
                              verbose=1){

  ### Verifications
  stopifnot(thresh.NA.mrk >= 0, thresh.NA.mrk <= 1,
            thresh.NA.ind >=0, thresh.NA.ind <= 1,
            is.logical(check.geno.dups), is.logical(check.mrk.dups),
            is.logical(remove.chrUkn), is.logical(format.names),
            is.logical(IDnum),is.logical(remove.nonPolyMrk),
            ## check if vcf.p2f or matrix.gt is provided
            xor(is.null(vcf.p2f), is.null(matrix.gt)),
            imputation %in% c(NULL,"kNNI","Beagle"),
            length(imputation) <=1)
  #vcf.file <- as.data.frame(vcf.file)


  ### 1. Load VCF file
  if(!is.null(vcf.p2f)){
    stopifnot(file.exists(vcf.p2f))

    if(verbose > 0) print("Load the vcf file...")
    ## load vcf file and convert to gt
    vcf <- vcfR::read.vcfR(vcf.p2f, verbose = FALSE)
    gt <- vcfR::extract.gt(vcf, element = "GT", as.numeric=F, IDtoRowNames = TRUE)
    meta <- vcf@meta
    mrk.info <- tibble::as_tibble(vcfR::getFIX(vcf))
    mrk.info$POS <- as.numeric(mrk.info$POS)
    mrk.info <- dplyr::arrange(mrk.info, CHROM,POS)
    ## combine marker information with genotype
    vcf.file <- cbind(mrk.info[,c("CHROM","POS")],gt[mrk.info$ID,])
    vcf.file$POS <- as.numeric(vcf.file$POS)

    rm(vcf, gt)

  } else {
    stopifnot(!is.null(matrix.gt), all(c("CHROM","POS") %in% colnames(matrix.gt)[c(1,2)]))
    vcf.file <- matrix.gt
    meta <- NULL
  }

  if(verbose >0){
    print(paste0("Initial dimensions: ", nrow(vcf.file)," markers and ",
                 (ncol(vcf.file) - 2), " genotypes"))
  }

  ## remove the prefix idxx. from the genotype name
  if(IDnum){
    colnames(vcf.file) <- gsub(x=colnames(vcf.file),pattern="^id[0-9]+\\.",replacement="")
  }
  ## remove duplicated rows
  ##vcf.file <- distinct(vcf.file) #

  vcf.file[,3:ncol(vcf.file)] <- as.data.frame(apply(vcf.file[,3:ncol(vcf.file)],2, function(x)
    gsub(pattern="\\.",x=x,replacement=NA)))
  ## replace patterns "./." by NA and '2/2" by "1/1"
  vcf.file[,3:ncol(vcf.file)] <- as.data.frame(apply(vcf.file[,3:ncol(vcf.file)],2, function(x)
    gsub(pattern="2/2",x=x,replacement="1/1")))

  ## set position as numeric
  vcf.file$POS <- as.numeric(stringr::str_trim(vcf.file$POS))

  ## update genotype names
  if(!is.null(corresp.geno.name)){
    stopifnot(ncol(corresp.geno.name) >=2)
    if(verbose >0){
      print("Update genotype names..")
    }

    colnames(vcf.file) <- plyr::mapvalues(colnames(vcf.file),
                                          from=corresp.geno.name[,1],
                                          to=corresp.geno.name[,2], warn_missing=F)
  }

  ## handle unknown chromosome, remove markers if few
  if(remove.chrUkn){
    chr.unk <- c("CHRUnknown","U","ChrUnknown","Unknown","Un","ChrUn")
    if(verbose >0){
      print(paste0("Removed ",nrow(vcf.file[vcf.file$CHROM %in% chr.unk,]),
                   " markers from unknown chromosome"))
    }
    vcf.file <- vcf.file[!vcf.file$CHROM %in% chr.unk,]

  }

  ## check for duplicated markers
  # extract map information
  vcf.file$chr_pos <- paste0(vcf.file$CHROM, "_",vcf.file$POS)
  if(check.mrk.dups){
    dup_snp <- unique(vcf.file$chr_pos[duplicated(vcf.file$chr_pos)])
    if(verbose >0){
      print(paste0("Found ",length(dup_snp)," duplicated markers"))
    }
    ## for duplicated loci, keep the one with the lowest missing value
    idx2rem <- c()
    for(mrk in dup_snp){
      ## get row id for duplicated markers
      idx.mrkdup <- grep(mrk,vcf.file$chr_pos)
      ## detect rows with less missing values
      idx2rem <- c(idx2rem,idx.mrkdup[-which.min(rowSums(is.na(vcf.file[vcf.file$chr_pos == mrk,])))])
    }
    ## remove duplicated columns with more missing values
    if(length(idx2rem)>0){
      vcf.file <- vcf.file[-idx2rem,]
    }
  }
  rownames(vcf.file) <- vcf.file$chr_pos
  ## if there is no chr before chr_pos, add it
  if(all(grepl("[1-7]|U",substr(rownames(vcf.file),1,1)))){
    rownames(vcf.file) <- paste0("chr",rownames(vcf.file))
  }
  cols2rem <- c("CHROM","POS", "chr_pos")
  vcf.file[,cols2rem]<- NULL

  ## 2. remove duplicated genotypes

  ## format genotype name with custom function format_names
  if(format.names){
    colnames(vcf.file) <- format_names(colnames(vcf.file))$corrected
  } else if(any(grep(" ", colnames(vcf.file))) & !is.null(p2f.export.vcf)){
    ## need to remove white spaces for synbreed export, replace it by underscore
    colnames(vcf.file) <- gsub(x=colnames(vcf.file),pattern=" ",replacement = "_")
  }
  if(concordance.function){
    stopifnot(!is.null(corresp.geno.list))
    conc <- concordance_match_name(match.name=colnames(vcf.file), corresp.list=corresp.geno.list)
    colnames(vcf.file) <- conc$new.names
  }

  if(check.geno.dups){
    ## classic detection of duplicated genotypes
    dup.genos1 <- unique(colnames(vcf.file)[duplicated(colnames(vcf.file))])

    ### detect duplicated genotypes with .1/.2/... appended at the end of the name
    dup.genos2 <- grep(pattern="\\.[1-9]{1,2}$", colnames(vcf.file), value=T)
    dup.genos2 <- unique(stringr::str_remove(dup.genos2, pattern="\\.[1-9]{1,2}$"))
    dup.genos <- unique(c(dup.genos1,dup.genos2))
    if(verbose >0){
      print(paste0("Found ",length(dup.genos)," duplicated genotypes"))
    }
    ## for duplicated genotypes, fill if missing data and keep the one with less missing data
    for(gen in dup.genos){
      idx.genodup <- grep(paste0("^",gen,"[.1-9]{0,2}"),colnames(vcf.file))
      #colnames(vcf.file.sel)[idx.genodup]
      if(length(idx.genodup) > 1){
        tmp <- vcf.file[,idx.genodup]

        ## find the column with the least NA values
        idxLessNa <- which.min(colSums(is.na(tmp)))
        ## fill NA values with data from other columns if available
        idxna <- which(is.na(tmp[,idxLessNa]))
        for(i in 1:ncol(tmp)){
          tmp[idxna,idxLessNa] <- dplyr::coalesce(tmp[idxna,idxLessNa],tmp[idxna,i])
        }
        ## suppress duplicated columns, keep the one with less NA values and filled
        vcf.file[,idx.genodup[idxLessNa]] <- tmp[,idxLessNa]
        colnames(vcf.file)[idx.genodup[idxLessNa]] <- gen
        idx2rem <- setdiff(idx.genodup,idx.genodup[idxLessNa])
        vcf.file = vcf.file[,-idx2rem]
      }
    }
    stopifnot(length(unique(colnames(vcf.file))) == ncol(vcf.file))
  }


  ## 3. Curate: remove SNPs/ genos with too many missing data

  # remove individuals with too many NA
  gen.miss <- apply(vcf.file, 2, function(x) sum(is.na(x)))
  geno2rem <- names(gen.miss)[gen.miss > thresh.NA.ind*nrow(vcf.file)]
  if(verbose > 0){
    print(paste0("Removed ",length(geno2rem)," individuals with too many missing data"))
  }
  vcf.file <- vcf.file[,!colnames(vcf.file) %in% geno2rem]

  ### remove markers with too many NA
  snp.miss <- apply(vcf.file, 1, function(x) sum(is.na(x)))
  snp2rem <- names(snp.miss)[snp.miss > thresh.NA.mrk*ncol(vcf.file)]
  if(verbose > 0){
    print(paste0("Removed ",length(snp2rem)," markers with too many missing data"))
  }
  vcf.file <- vcf.file[!rownames(vcf.file) %in% snp2rem,]

  ## Remove non-polymorphic markers
  if(remove.nonPolyMrk){
    nb.allele <- apply(vcf.file, 1, function(x) length(table(t(x))))
    mrk2rem <- names(nb.allele)[nb.allele <= 1]
    if(verbose > 0){
      print(paste0("Removed ",length(mrk2rem)," non-polymorphic markers"))
    }
    if(length(mrk2rem) >0) vcf.file <- vcf.file[!rownames(vcf.file) %in% mrk2rem,]
    rm(mrk2rem)
  }

  ## Remove markers with too many heterozygous genotypes
  if(!is.null(tresh.heterozygous)){
    hetero <- apply(vcf.file, 2, function(x) {sum(x %in% c("0/1","1/0","0|1","1|0"), na.rm = T)})
    geno2rem <- names(hetero)[hetero > tresh.heterozygous * nrow(vcf.file)]
    if(verbose > 0){
      print(paste0("Removed ",length(geno2rem)," genotypes with too many heterozygous markers"))
    }
    vcf.file <- vcf.file[,!colnames(vcf.file) %in% geno2rem,]
  }

  rownames(vcf.file) <- gsub(x = rownames(vcf.file),pattern="Chr",replacement = "chr")
  ## recode numerically: replace 0/0=0, 0/1 or 1/0=1 and 1/1=2
  rowN <- rownames(vcf.file)

  vcf.num <- apply(vcf.file, 2, function(x) {
    as.numeric(stringr::str_replace_all(gsub("[[:punct:]]", "", x),
                                        c("00"="0","01"="1","10"="1","11"="2")))})
  vcf.num <- as.data.frame(vcf.num)
  rownames(vcf.num) <- rowN

  ### 4. (optional) export to a formatted vcf file for Beagle imputation

  if(!is.null(p2f.export.vcf)){
    ## if no extension or vcf extension provided, extend to vcf.gz extension
    if(!grepl("vcf.gz$",p2f.export.vcf)){
      p2f.export.vcf <- paste0(gsub(".vcf","",p2f.export.vcf),".vcf.gz")
    }

    ## recreate vcfR object
    vcf2exp <- methods::new(Class = "vcfR")
    ### meta element
    if(!is.null(meta)){
      vcf2exp@meta <- meta
    } else {
      vcf2exp@meta <- c("##fileformat=VCFv4.2",
                        "##reference=xxx",
                        "##comment=Generated from format_curate_vcf function",
                        #"##filedate= 20065",
                        "##source='write.vcf of vcfR'",
                        "##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>")
    }

    ### fix element: add marker information
    #### define a default mrk.info data frame if not provided
    if(is.null(vcf.p2f) & is.null(mrk.info)){
      chr_pos <- rownames(vcf.num)
      mrk.info <- data.frame(CHROM=unlist(lapply(strsplit(chr_pos,"_"),function(x) x[1])),
                             POS=unlist(lapply(strsplit(chr_pos,"_"),function(x) x[2])),
                             ID=chr_pos,
                             REF="A",ALT="T",
                             QUAL=".",FILTER="PASS",INFO=".",
                             FORMAT="GT")
      mrk.info$CHROM <- gsub("Chr","chr",mrk.info$CHROM)

    } else {
      ## use the mrk.info from the vcf file, subset selected markers
      mrk.info$ID <- paste0(gsub("Chr","chr",mrk.info$CHROM),"_",mrk.info$POS)
      mrk.info <- mrk.info[match(rownames(vcf.num), mrk.info$ID),]
    }
    if(!all(c("INFO","FORMAT") %in% colnames(mrk.info))){
      mrk.info <- cbind(mrk.info,
                        #QUAL=".",FILTER="PASS",
                        INFO=".",FORMAT="GT")
      mrk.info <- mrk.info[,c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")]
    }
    mrk.info <- mrk.info[match(rowN, mrk.info$ID),]
    stopifnot(nrow(mrk.info) == nrow(vcf.num))
    mrk.info <- as.matrix(mrk.info)
    mrk.info[,"POS"] <- gsub(" ","",mrk.info[,"POS"])
    vcf2exp@fix <- mrk.info

    ### gt element, coded in 0/0, 0/1, 1/1
    vcf2exp@gt <- as.matrix(vcf.file)

    ### export
    vcfR::write.vcf(x=vcf2exp, file=p2f.export.vcf, mask=FALSE, APPEND=FALSE)

  } ## end if !is.null(p2f.export.vcf)


  ### Imputation of missing markers with kNNI
  if(length(imputation) >0 && imputation == "kNNI"){
    if(is.null(p2f.export.vcf)){
      stop("Need to export vcf file to perform imputation")
    }
    ## load rTASSEL and rJava packages
    if(!requireNamespace("rTASSEL", quietly = TRUE)){
      stop("Please install rTASSEL package to perform kNNI imputation")
    }
    if(!requireNamespace("rJava", quietly = TRUE)){
      stop("Please install rJava package to perform kNNI imputation")
    }

    tasGenoVCF <- rTASSEL::readGenotypeTableFromPath(path = p2f.export.vcf)
    ## use KNNI as imputation method, using default parameters
    tasGenoImp <- rTASSEL::imputeLDKNNi(tasObj=tasGenoVCF, highLDSSites = 30,
                                        knnTaxa = 10, maxDistance = 1e+07)
    ## use an internal function to convert this object to a data frame
    vcf.file <- getGenoTas_to_DF(tasGenoImp)

    ## format it to export it as a vcf file
    vcf2exp@gt <- as.matrix(vcf.file)
    # gp$geno <- t(vcf.file[rownames(gp$map),rownames(gp$geno)])
    # stopifnot(identical(colnames(gp$geno),rownames(gp$map)))

    ## detect and warn if missing values left

    if(verbose > 0) {
      ### proportion of missing value per column
      print("Summary of missing data percentage per marker: ")
      print(summary(colMeans(is.na(gp$geno)) *100))

      ### proportion of missing value per line
      print("Summary of missing data percentage per genotype: ")
      print(summary(rowMeans(is.na(gp$geno)) *100))

    }
    # # export VCF
    # synbreed::write.vcf(gp,file=p2f.export.vcf)
    ### export
    vcfR::write.vcf(x=vcf2exp, file=p2f.export.vcf, mask=FALSE, APPEND=FALSE)

  } else if(length(imputation) >0 && imputation == "Beagle"){
    if(is.null(p2f.export.vcf)){
      stop("Need to export vcf file to perform imputation")
    }
    if(!file.exists(p2f.beagle)){
      stop("Beagle jar file not found")
    }
    ## run Beagle imputation
    p2f.imp <- gsub(".vcf.gz","_imputed",p2f.export.vcf)
    system(paste0("java -Xmx2g -jar ",p2f.beagle," gt=",p2f.export.vcf,
                  " out=",p2f.imp,
                  " nthreads=4 ne=100000"))

    ## load imputed vcf
    vcf.imp <- vcfR::read.vcfR(paste0(p2f.imp,".vcf.gz"), verbose = FALSE)
    gt.imp <- vcfR::extract.gt(vcf.imp, element = "GT", as.numeric=F, IDtoRowNames = FALSE)
    mrk.info.imp <- tibble::as_tibble(vcfR::getFIX(vcf.imp))
    vcf.file <- cbind(mrk.info.imp[,c("CHROM","POS")],gt.imp)
    vcf.file$POS <- as.numeric(vcf.file$POS)
    vcf.file <- dplyr::arrange(vcf.file, CHROM,POS)
    rownames(vcf.file) <- paste0(vcf.file$CHROM, "_",vcf.file$POS)
    cols2rem <- c("CHROM","POS")
    vcf.file[,cols2rem]<- NULL
    rowN <- rownames(vcf.file)
    ## recode numerically: replace 0/0=0, 0/1 or 1/0=1 and 1/1=2
    vcf.num <- apply(vcf.file, 2, function(x) {
      as.numeric(stringr::str_replace_all(gsub("[[:punct:]]", "", x),
                                          c("00"="0","01"="1","10"="1","11"="2")))})
    vcf.num <- as.data.frame(vcf.num)
    rownames(vcf.num) <- rowN
  } ## end if Beagle imput


  ## Remove markers with MAF below threshold after imputation
  if(!is.null(thresh.MAF)){
    thresh.MAF <- as.numeric(thresh.MAF)
    if(thresh.MAF < 0 | thresh.MAF > 0.5){
      stop("MAF should be between 0 and 0.5")
    }

    maf <- apply(vcf.num, 2, function(x) {
      t = table(x) ; t = t[names(t) %in% c(0,2)]
      if(length(t) < 2) maf = 0
      else maf = sort(t)[1]/length(x)
      return(maf)
    })

    mrk2rem <- names(maf)[maf < thresh.MAF]

    if(verbose > 0){
      print(paste0("Removed ",length(mrk2rem)," markers with MAF below ",thresh.MAF))
    }
    vcf.num <- vcf.num[!rownames(vcf.num) %in% mrk2rem,]
  }


  if(verbose >0){
    print(paste0("Final dimensions: ",ncol(vcf.num), " genotypes and ", nrow(vcf.num)," markers"))
  }
  return(t(vcf.num)) ## output in dimension geno x marker
}


## author: Vishnu Ramasubramanian
#' Convert TASSEL object to genotypic data frame
#'
#' @param tasGeno Genotypic object from TASSEL
#'
#' @importFrom SummarizedExperiment assays colData
#' @returns matrix of genotypic data in 0/1/2 data frame format
#' @keywords internal
#' @author Vishnu Ramasubramanian
getGenoTas_to_DF <- function(tasGeno){
  ## load rTASSEL and rJava packages
  if(!requireNamespace("rTASSEL", quietly = TRUE)){
    stop("Please install rTASSEL package to perform kNNI imputation")
  }
  if(!requireNamespace("rJava", quietly = TRUE)){
    stop("Please install rJava package to perform kNNI imputation")
  }


  tasSumExp <- rTASSEL::getSumExpFromGenotypeTable(tasObj=tasGeno)
  tasGenoDF <- (SummarizedExperiment::assays(tasSumExp)[[1]])
  colnames(tasGenoDF) <- SummarizedExperiment::colData(tasSumExp)[,"Sample"]

  ### Extract Table Report DF
  tableReport <- rJava::new(
    rJava::J("net.maizegenetics.dna.map.PositionListTableReport"),
    tasGeno %>% rTASSEL::positionList()) %>%
    rTASSEL::tableReport() %>% as.data.frame()

  varSplit <- strsplit(tableReport[,"VARIANT"],"/")
  varSplitTab <- cbind.data.frame(unlist(lapply(varSplit,function(x) x[1])),
                                  unlist(lapply(varSplit,function(x) x[2])))

  gt2d_tasGeno <- tasGenoDF
  rownames(gt2d_tasGeno) <- tableReport$Name
  return(gt2d_tasGeno)
}



#' Compute genomic prediction different methods
#'
#' @param geno genomic data with genotypes in row (GID in rownames) and marker in columns. Values should be column centered and scaled.
#' @param pheno phenotypic data with genotypes in row (in GID column) and traits in columns. Phenotypic value should be corrected for year and location effects beforehand.
#' @param traits character vector of trait names
#' @param GP.method character vector of length one of genomic prediction methods to use. Must be one of "rrBLUP", "GBLUP", "RKHS", "RKHS-KA", "RandomForest", "BayesA", "BayesB" or "LASSO".
#' @param nreps number of repetitions for cross-validation, default is 10
#' @param nfolds number of folds for cross-validation, default is 10
#' @param h bandwith parameter for RKHS, default is 1. If multiple values are provided, method will be RKHS Kernel Averaging.
#' @param nIter number of iterations for RKHS, default is 6000
#' @param burnIn number of burn-in iterations for RKHS, default is 1000
#' @param ntree number of trees for RandomForest, default is 100
#' @param nb.mtry number of mtry for RandomForest, default is 10
#' @param p2d.temp path to directory to export temporary genomic prediction results, default is NULL (could cause error in parallelization if NULL).
#' @param nb.cores number of cores to parallelize the computation, default is 1 (no parallelization)
#' @param p2f.stats path to file to export genomic prediction results, default is NULL
#'
#' @return a list with the following elements: `obspred` with observed vs. predicted genotypic values, and `gp.stats` with genomic prediction statistics
#' @author Charlotte Brault
#' @seealso [getFolds()]
#' @importFrom BGLR BGLR
#' @importFrom parallel makeCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach %dopar%
#' @importFrom rrBLUP kinship.BLUP
#' @importFrom caret train predict.train
#' @importFrom randomForest randomForest
#' @importFrom dplyr bind_rows
#' @importFrom stats na.omit cor
#' @importFrom utils write.csv write.table
#' @importFrom glmnet cv.glmnet
#' @importFrom purrr map_dfr
#' @export

compute_GP_methods <- function(geno, pheno, traits, GP.method, nreps=10,
                               nfolds=10, h=1, nb.mtry=10, nIter=6000,burnIn=1000,
                               ntree=100,p2d.temp=NULL,
                               nb.cores=1, p2f.stats=NULL){


  stopifnot(GP.method %in% c("rrBLUP","GBLUP","RKHS","RKHS-KA","RandomForest",
                             "BayesA","BayesB","LASSO"),
            all(is.numeric(h)), length(GP.method) == 1,
            "GID" %in% colnames(pheno))
  geno <- as.matrix(geno)
  # prepare outputs
  ## create a data frame to store accuracy
  gp.stats <- expand.grid(trait=traits, rep=seq(nreps), corP=NA,
                          GP.method=GP.method, nb.genos.TS=NA,
                          nb.genos.VS=NA, nb.snps=ncol(geno))
  ## list to store observed / predicted breeding values
  pred.list <- list()

  ## for loop across traits
  for(tr in traits){
    write(paste0("Working on trait ",tr), stderr())
    ## format geno and pheno for the trait, after removing missing values
    pheno_tr <- stats::na.omit(pheno[,c("GID",tr)])
    inds <- intersect(pheno_tr$GID, rownames(geno))
    pheno_tr <- pheno_tr[match(inds, pheno_tr$GID),]
    geno_tr <- geno[inds,]
    ## vector of phenotypic values
    y <- pheno_tr[,tr] ; names(y) <- pheno_tr$GID

    ## genomic relationship matrix
    G <- as.matrix(crossprod(geno_tr)/ncol(geno_tr))


    ## create cluster
    cl <- parallel::makeCluster(type="SOCK",spec=nb.cores)
    doSNOW::registerDoSNOW(cl)

    for(r in 1:nreps){
      write(r, stderr())
      ## get cross-validation folds
      folds <- getFolds(nrow(pheno_tr), nb.folds=nfolds, seed=r+57425)


      if(GP.method %in% "rrBLUP"){
        ## parallel computation of GP for rrBLUP
        out <- foreach::foreach(f=1:nfolds) %dopar%{
          ## estimate marker effects on training set
          res <- rrBLUP::kinship.BLUP(y=y[-folds[[f]]],
                                      G.train=geno_tr[-folds[[f]],],
                                      G.pred=geno_tr[folds[[f]],], K.method="RR")$g.pred
          return(res)
        }

      } else if(GP.method %in% "GBLUP"){
        ## parallel computation of GP for GBLUP
        out <- foreach::foreach(f=1:nfolds) %dopar%{
          ## compute the genomic relationship matrix
          K <- tcrossprod(geno_tr)/ncol(geno_tr)

          ## function inputs
          data <- pheno_tr
          ## set phenotypes to NA for the training set
          data[[tr]][folds[[f]]] <- NA
          ## estimate marker effects on training set
          res <- rrBLUP::kin.blup(data=data,geno="GID",pheno=tr,K=K)
          return(res$pred[folds[[f]]])
        }


      } else if(GP.method %in% "RKHS" & length(h) ==1){
        ## define distance matrix for RKHS
        D <- as.matrix(dist(geno_tr,method="euclidean"))^2
        D <- D/mean(D)
        ### compute kernel
        K <- exp(-h*D)

        ## parallel processing
        out <- foreach::foreach(f=1:nfolds,.errorhandling = "pass") %dopar% {
          ## estimate marker effects on training set
          yNA <- y
          yNA[folds[[f]]] <- NA
          if(!is.null(p2d.temp)){
            p2f.temp <- paste0(p2d.temp,"/",GP.method,"_",tr,
                               "_",r,"_",f)
          }
          fit <- BGLR::BGLR(y=yNA,ETA=list(list(K=K,model='RKHS')),
                            nIter=nIter,burnIn=burnIn,
                            saveAt=p2f.temp,verbose=FALSE)
          return(as.data.frame(fit$yHat[folds[[f]]]))
        }

      } else if(GP.method %in% "RKHS-KA" | length(h) >1){
        ## define distance matrix for RKHS
        D <- as.matrix(dist(geno_tr,method="euclidean"))^2
        D <- D/mean(D)
        KList <- list()
        for(i in 1:length(h)){
          KList[[i]]<-list(K=exp(-h[i]*D),model='RKHS')
        }

        ## parallel processing
        out <- foreach::foreach(f=1:nfolds,.errorhandling = "pass") %dopar% {
          ## estimate marker effects on training set
          yNA <- y
          yNA[folds[[f]]] <- NA
          if(!is.null(p2d.temp)){
            p2f.temp <- paste0(p2d.temp,"/",GP.method,"_",
                               tr,"_",r,"_",f)
          }
          fit <- BGLR::BGLR(y=yNA,ETA=KList,
                            nIter=nIter,burnIn=burnIn,
                            saveAt=p2f.temp,verbose=FALSE)
          return(as.data.frame(fit$yHat[folds[[f]]]))
        }


      } else if(GP.method %in% "RandomForest"){

        ## optimize mtry = number of randomly selected variables at each split
        #tunegridrf <- expand.grid(.mtry=seq(1,ncol(geno_tr)/3,length.out=nb.mtry))
        tunegridrf <- expand.grid(.mtry=c(100,500,1000,2000,5000))
        ## parallel processing
        out <- foreach::foreach(f=1:nfolds,.errorhandling="pass") %dopar% {
          ## estimate marker effects on training set

          fit <- caret::train(y=y[-folds[[f]]],
                              x=geno_tr[-folds[[f]],],
                              method = "rf",tuneGrid = tunegridrf, ntree=ntree)
          return(as.data.frame(caret::predict.train(fit,geno_tr[folds[[f]],])))
        }
      }else if(GP.method %in% "LASSO"){
        out <- foreach::foreach(f=1:nfolds) %dopar% {
          ### Fit the model with glmnet package
          cv <- glmnet::cv.glmnet(y=y[-folds[[f]]],x=geno_tr[-folds[[f]],], alpha=1)
          fit <- glmnet::glmnet(y=y[-folds[[f]]],x=geno_tr[-folds[[f]],],
                                alpha=1,lambda=cv$lambda.min)
          return(geno_tr[folds[[f]],] %*% as.matrix(fit$beta))
        }

      } else if(GP.method %in% "BayesA"){
        out <- foreach::foreach(f=1:nfolds) %dopar% {

          if(!is.null(p2d.temp)){
            p2f.temp <- paste0(p2d.temp,"/",GP.method,"_",tr,"_",r,"_",f)
          }
          ## estimate marker effects on training set
          fit <- BGLR::BGLR(y=y[-folds[[f]]],
                            ETA=list(list(X=geno_tr[-folds[[f]],],model='BayesA')),
                            nIter=nIter,burnIn=burnIn,
                            saveAt=p2f.temp,verbose=FALSE)
          return(geno_tr[folds[[f]],] %*% as.matrix(fit$ETA[[1]]$b))
        }

      } else if(GP.method %in% "BayesB"){
        out <- foreach::foreach(f=1:nfolds) %dopar% {
          if(!is.null(p2d.temp)){
            p2f.temp <- paste0(p2d.temp,"/",GP.method,"_",tr,"_",r,"_",f)
          }
          ## estimate marker effects on training set
          fit <- BGLR::BGLR(y=y[-folds[[f]]],
                            ETA=list(list(X=geno_tr[-folds[[f]],],model='BayesB')),
                            nIter=nIter,burnIn=burnIn,
                            saveAt=p2f.temp,verbose=FALSE)
          ## return predicted value on validation set for this fold
          return(geno_tr[folds[[f]],] %*% as.matrix(fit$ETA[[1]]$b))
        }

      } ## end if method


      ## combine all predicted values across folds
      ### (too few individuals to calculate predictive ability for each fold)
      pred.obs.all <- purrr::map_dfr(out,as.data.frame)
      pred.obs.all <- data.frame(id_geno=rownames(pred.obs.all),
                                 obs=pheno[match(rownames(pred.obs.all),pheno$GID),tr],
                                 pred=pred.obs.all[,1],
                                 trait=tr,rep=r, GP.method=GP.method)
      pred.list[[tr]][[r]] <- pred.obs.all


      ## estimnate predictive ability as Pearson's correlation between observed and predicted
      idx <- which(gp.stats$trait %in% tr & gp.stats$rep %in% r)
      gp.stats$corP[idx] <- round(cor(pred.obs.all$obs,pred.obs.all$pred),3)
      gp.stats$nb.genos.TS[idx] <- length(y)-mean(sapply(folds, length))
      gp.stats$nb.genos.VS[idx] <- mean(sapply(folds, length))

    } ## end for rep

    parallel::stopCluster(cl)

  } ## end for trait

  ## if file path provided, export results.
  if(!is.null(p2f.stats)){
    ext <- tools::file_ext(p2f.stats)
    if(ext %in% c("csv")){
      utils::write.csv(gp.stats, file=p2f.stats, row.names=FALSE)
    } else if(ext %in% c("txt","tsv")){
      utils::write.table(file=p2f.stats, gp.stats, sep="\t", row.names=FALSE)
    } else if(ext %in% "rds"){
      saveRDS(gp.stats, file=p2f.stats)
    }

  }

  return(list(obspred=pred.list, gp.stats=gp.stats))

}


#' Run genomic prediction on all genotypes, output predicted values
#'
#' @param geno genomic data with genotypes in row (GID in rownames) and marker in columns. Values should be column centered and scaled.
#' @param pheno phenotypic data with genotypes in row (in GID column) and traits in columns
#' @param traits character vector of trait names, must correspond to `pheno` column names
#' @param GP.method character vector of genomic prediction methods to use, must be one of "rrBLUP", "RKHS", "BayesA", "BayesB"
#' @param runCV logical, if TRUE, run cross-validation on the common genotypes in `pheno` and `geno`, default is FALSE
#' @param testSetGID GID of the test set genotypes, if NULL, will provide predicted genotypic values for all genotypes in `geno`
#' @param nreps number of repetitions for cross-validation, default is 10
#' @param nfolds number of folds for cross-validation, default is 10
#' @param h bandwith parameter for RKHS, default is 1.
#' @param nb.mtry number of randomly selected variables at each split for RandomForest, default is 10
#' @param nIter number of iterations for RKHS, default is 6000
#' @param burnIn number of burn-in iterations for RKHS, default is 1000
#' @param ntree number of trees for RandomForest, default is 100
#' @param p2d.temp path to directory to export temporary genomic prediction results, default is NULL (could cause error in parallelization if NULL).
#' @param nb.cores number of cores to parallelize the computation, default is 1 (no parallelization)
#' @param p2f.stats.cv path to file to export cross-validation genomic prediction results, default is NULL
#' @param p2f.pred path to file to export genotypic values from genomic prediction, default is NULL
#' @param verbose integer, level of verbosity, default is 1
#'
#' @return a data frame with the predicted values for all genotypes and traits
#' @seealso [compute_GP_methods()]
#' @author Charlotte Brault
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach %dopar% foreach
#' @importFrom rrBLUP kinship.BLUP
#' @importFrom tidyr pivot_wider
#' @importFrom tidyselect all_of
#' @export

compute_GP_allGeno <- function(geno, pheno, traits, GP.method,
                               runCV=FALSE, testSetGID=NULL,
                               nreps=10,nfolds=10, h=1, nb.mtry=10,
                               nIter=6000,burnIn=1000,
                               ntree=100,p2d.temp=NULL,
                               nb.cores=1, p2f.stats.cv=NULL,
                               p2f.pred=NULL, verbose=1){

  ## initial verifications
  stopifnot(GP.method %in% c("rrBLUP","RKHS","BayesA","BayesB","BayesC"),
            #"LASSO","RKHS-KA","RandomForest"),
            all(is.numeric(h)), length(GP.method) == 1,
            "GID" %in% colnames(pheno),
            all(traits %in% colnames(pheno)))

  ## formatting
  geno <- as.matrix(geno)
  ntraits <- length(traits)
  pheno.traits <- pheno[,c("GID",traits)]

  ## cross-validation
  if(runCV){
    cmonGID <- intersect(rownames(geno), pheno$GID)
    if(verbose > 0){
      print(paste0("Compute genomic prediction in cross-validation for ",
                   length(cmonGID)," common genotypes"))
    }
    ## use custom function
    gp.cv <- compute_GP_methods(
      geno = geno[cmonGID,],
      pheno = pheno[match(cmonGID, pheno$GID),],
      traits = traits,
      nfolds = nfolds,
      nreps = nreps,
      GP.method = GP.method,
      h = h,
      nIter = nIter,
      burnIn = burnIn,
      p2d.temp=p2d.temp,
      nb.cores=nb.cores,
      p2f.stats=p2f.stats.cv)
  }


  ## create cluster
  cl <- parallel::makeCluster(type="SOCK",spec=nb.cores)
  doSNOW::registerDoSNOW(cl)
  ## run genomic prediction on all genotypes

  if(GP.method %in% "rrBLUP"){
    ## parallel computation of GP for rrBLUP
    out <- foreach::foreach(f=1:ntraits,.errorhandling = "pass",
                            .combine = "rbind",.packages = "rrBLUP") %dopar% {
                              ## estimate marker effects on all available data
                              yNA <- merge(data.frame(GID=c(rownames(geno))),
                                           pheno.traits[,c("GID",traits[f])], by="GID", all.x=TRUE)
                              yGID <- yNA$GID
                              ## for test set geno if provided, set to NA (remove from training set)
                              if(!is.null(testSetGID)){
                                yNA[yNA$GID %in% testSetGID,tr] <- NA
                              }
                              fit <- rrBLUP::kinship.BLUP(y=yNA[[traits[f]]],
                                                          G.train=geno,
                                                          G.pred=geno, K.method="RR")$g.pred
                              ## output predicted values
                              tmp <- data.frame(GID=yGID, GP.method=GP.method,trait=traits[f],yHat=fit)
                              return(tmp)
                            }
  } else if(GP.method == "RKHS"){
    D <- as.matrix(dist(geno,method="euclidean"))^2
    D <- D/mean(D, na.rm=T)
    ### compute kernel
    K <- exp(-h*D)
    rm(D)

    if(is.null(p2d.temp)){
      print(warning("No temporary directory provided, this may cause issues."))
    }

    ## parallel processing over traits
    out <- foreach::foreach(
      f=1:ntraits,.errorhandling = "pass",
      .combine = "rbind",.packages = "BGLR") %dopar% {
        yNA <- merge(data.frame(GID=c(rownames(geno))),
                     pheno.traits[,c("GID",traits[f])],
                     by="GID", all.x=TRUE)
        yNA <- yNA[match(rownames(K),yNA$GID),]
        ## if testSetGID provided, set to NA (remove from training set)
        if(!is.null(testSetGID)){
          yNA[yNA$GID %in% testSetGID,traits[f]] <- NA
        }
        stopifnot(identical(yNA$GID, rownames(geno)))

        p2f.temp <- paste0(p2d.temp,"/RKHS_PredAll_testSet_",traits[f])
        ## fit model
        fit <- BGLR::BGLR(y=yNA[[traits[f]]],
                          ETA=list(list(K=K,model='RKHS')),
                          nIter=nIter,burnIn=burnIn,
                          saveAt=p2f.temp,verbose=FALSE)
        ## output predicted values
        tmp <- data.frame(GID=yNA$GID, GP.method=GP.method,
                          trait=traits[f],yHat=fit$yHat)
        return(tmp)
      }

  } else if(GP.method %in% c("BayesA","BayesB","BayesC")){
    if(is.null(p2d.temp)){
      print(warning("No temporary directory provided, this may cause issues."))
    }
    ## parallel processing over traits
    out <- foreach::foreach(
      f=1:ntraits,.errorhandling = "pass",
      .combine = "rbind",.packages = "BGLR") %dopar% {
        yNA <- merge(data.frame(GID=c(rownames(geno))),
                     pheno.traits[,c("GID",traits[f])], by="GID", all.x=TRUE)

        yNA <- yNA[match(rownames(K),yNA$GID),]
        ## if testSetGID provided, set to NA (remove from training set)
        if(!is.null(testSetGID)){
          yNA[yNA$GID %in% testSetGID,traits[f]] <- NA
        }
        stopifnot(identical(yNA$GID, rownames(geno)))

        p2f.temp <- paste0(p2d.temp,"/",GP.method,"_PredAll_testSet_",traits[f])
        ## fit model
        fit <- BGLR::BGLR(y=yNA[[traits[f]]],
                          ETA=list(list(X=geno,model=GP.method)),
                          nIter=nIter,burnIn=burnIn,
                          saveAt=p2f.temp,verbose=FALSE)
        ## output predicted values
        tmp <- data.frame(GID=yNA$GID, GP.method=GP.method,
                          trait=traits[f],yHat=fit$yHat)
        return(tmp)
      }
  }

  parallel::stopCluster(cl)
  ## set the predicted values in a wide format, with one column per trait
  all.pred <- out %>%
    dplyr::select(tidyselect::all_of(c("GID","GP.method","trait","yHat"))) %>%
    tidyr::pivot_wider(names_from="trait", values_from="yHat")

  if(!is.null(testSetGID)){
    ## if testSetID provided, keep only the genotypes in the test set
    all.pred <- all.pred[all.pred$GID %in% testSetGID,]
  }


  ## export predicted values
  if(!is.null(p2f.pred)){
    ## if file path provided, export results.
    if(!is.null(p2f.stats)){
      ext <- tools::file_ext(p2f.pred)
      if(ext %in% c("csv")){
        utils::write.csv(all.pred, file=p2f.pred, row.names=FALSE)
      } else if(ext %in% c("txt","tsv")){
        utils::write.table(file=p2f.pred, all.pred, sep="\t",
                           row.names=FALSE)
      } else if(ext %in% "rds"){
        saveRDS(all.pred, file=p2f.pred)
      }
    }

  } ## end if p2f.pred

  if(runCV){
    return(list(obspred=gp.cv$obspred, gp.stats=gp.cv$gp.stats,
                all.pred=all.pred))
  } else {
    return(all.pred)
  }
}


#' Estimate genetic gain for selected traits
#'
#' @param data data frame with columns for traits and first year of evaluation
#' @param traits character vector of trait names
#' @param first_year character vector of column names for first year of evaluation.
#' The first year will be subtracted by the minimal first year for the estimation of the intercept and percentage of change
#'
#' @returns data frame with columns for trait, intercept, slope, R2.adj, pvalue and percentage of change
#' @author Charlotte Brault
#' @export

GGAIN <- function(data, traits, first_year){
  stopifnot(all(traits %in% colnames(data)),
            all(first_year %in% colnames(data)))
  ## subtract the minimum first year to the calculation of the percentage of change
  data$FIRST_YEAR0 <- data[[first_year]] - min(data[[first_year]], na.rm=T) +1

  df.fit <- data.frame(trait=traits, intercept=NA, slope=NA,
                       slope.min=NA, slope.max=NA, R2.adj=NA, pvalue=NA)
  for(tr in traits){
    fit <- lm(paste0(tr, "~ 1 + FIRST_YEAR0"), data = data)
    idx <- which(df.fit$trait %in% tr)
    ## add intercept and slope
    df.fit[idx,c("intercept","slope")] <- round(c(coef(fit)[c(1,2)]),3)
    ## get confidence interval
    df.fit[idx,c("slope.min","slope.max")] <- round(confint(fit, "FIRST_YEAR0", level = 0.95)[,1:2],3)
    ## get R2
    df.fit[idx,"R2.adj"] <- round(summary(fit)$adj.r.squared,3)
    ## get p-value
    df.fit[idx,"pvalue"] <-  formatC(summary(fit)$coefficients[2,4], format="e", digits=2)
  }
  ## estimate the percentage of change
  df.fit$percChange <- (df.fit$slope/df.fit$intercept)*100


  return(df.fit)
}


## ----- Print / plot functions -----






## from R metan package
# Function to make HTML tables
#' Print nice table using DT
#'
#' @param table data frame to print
#' @param rownames logical, if TRUE, print row names, default is FALSE
#' @param digits integer, number of digits to print for numeric columns, default is 3
#' @param ... additional arguments to pass to DT::datatable
#'
#' @returns a DT datatable object
#' @export
#' @examples
#' print_table(mtcars)

print_table <- function(table, rownames = FALSE, digits = 3, ...){
  df <- DT::datatable(table, rownames = rownames, extensions = 'Buttons',
                      options = list(scrollX = TRUE,
                                     dom = '<<t>Bp>',
                                     buttons = c('copy','csv', 'excel', 'pdf', 'print')), ...)
  num_cols <- c(as.numeric(which(sapply(table, class) == "numeric")))
  if(length(num_cols) > 0){
    DT::formatSignif(df, columns = num_cols, digits = digits)
  } else{
    df
  }
}




#' Print data summary for numeric variables using gtsummary
#'
#' @param data data frame with numeric variables
#' @param variables character vector of numeric variable names to include in the summary
#' @param by character vector of variable names to group by
#' @param add.pval logical, if TRUE, add p-values to the summary table
#'
#' @returns a gtsummary object
#' @author Charlotte Brault
#' @importFrom gtsummary tbl_summary add_p
#' @export

data_summary <- function(data, variables=NULL, by=NULL, add.pval=T){

  ## if not provided, select numeric variables
  if(is.null(variables)){
    variables <- names(data)[sapply(data, is.numeric)]
  }
  data[,variables] <- apply(data[,variables], 2, as.numeric)


  gt <- gtsummary::tbl_summary(data,
                               include=tidyselect::all_of(variables),by = by,
                               type = all_continuous() ~ "continuous2",
                               statistic = all_continuous() ~ c(
                                 "{mean} ({sd})",
                                 "{min}, {max}",
                                 "{N_nonmiss}"))
  ## add p-value of comparison between groups if groups provided
  if(add.pval & !is.null(by)) gt <- gt %>% gtsummary::add_p()

  return(gt)

}



#' Plot distribution of breeding values
#' Plot distribution of breeding values for each trait and add position of selected genotypes
#'
#' @param BV data frame with breeding values, with GID in rows and traits in columns
#' @param trait character vector of trait name (length 1)
#' @param genoCol character vector of column name for GID (default is "GID")
#' @param selGen character vector of selected genotypes to highlight in the plot
#' @param colorCol character vector of column name for color (default is NULL)
#'
#' @importFrom ggplot2 aes geom_density geom_area geom_vline labs ggplot geom_rug scale_fill_viridis_d scale_color_viridis_d
#' @importFrom ggrepel geom_label_repel
#' @importFrom bigstatsr theme_bigstatsr
#' @returns a ggplot object
#' @export
plotDistrib_selGen <- function(BV, trait, genoCol="GID", selGen=NULL, colorCol=NULL){

  ## add verifications
  if(!trait %in% colnames(BV)){
    stop(paste0("Trait ", trait, " not found in the dataset"))
    stopifnot(length(trait) == 1)
  }
  if(!genoCol %in% colnames(BV)){
    stop(paste0("Genotype column ", genoCol, " not found in the dataset"))
  }
  if(!is.null(selGen) && !all(selGen %in% BV[[genoCol]])){
    stop(paste0("Selected genotypes not found in the dataset"))
  }
  if(!is.null(colorCol) && !colorCol %in% colnames(BV)){
    stop(paste0("Color column ", colorCol, " not found in the dataset"))
  }
  ## split the distribution into densities
  x <- BV[[trait]]
  q15.9 <- quantile(x, .159,na.rm=TRUE) # 1 Std 68.2%
  q84.1 <- quantile(x, .841,na.rm=TRUE)
  q2.3  <- quantile(x, .023,na.rm=TRUE) # 2 Std 95.4%
  q97.7 <- quantile(x, .977,na.rm=TRUE)
  q0.01 <- quantile(x, .001,na.rm=TRUE) # 3 Std 99.8%
  q99.9 <- quantile(x, .999,na.rm=TRUE)
  meanx <- mean(x,na.rm=TRUE)
  medx  <- median(x,na.rm=TRUE)
  x.dens  <- density(x,na.rm=TRUE)
  df.dens <- data.frame(x=x.dens$x, y=x.dens$y)


  ## compose the plot
  p <- ggplot2::ggplot(BV, aes(x=.data[[trait]]))+
    ggplot2::geom_density(color = 'skyblue') +
    ggplot2::geom_area(data = subset(df.dens, x >= q15.9 & x <= q84.1), # 1 Std 68.2%
                       aes(x=x,y=y), fill='skyblue', alpha=0.8) +
    ggplot2::geom_area(data = subset(df.dens, x >= q2.3 & x <= q97.7), # 2 Std 95.4%
                       aes(x=x,y=y), fill='skyblue', alpha=0.6) +
    ggplot2::geom_area(data = subset(df.dens, x >= q0.01 & x <= q99.9), # 3 Std 99.8%
                       aes(x=x,y=y), fill='skyblue', alpha=0.3) +
    ggplot2::geom_vline(xintercept=meanx, color="grey60", linewidth=1.5, linetype="dashed") +
    ggplot2::geom_vline(xintercept=medx, color='#FFFFFF',linewidth=1.5, linetype="dashed") +
    ggplot2::labs(title=trait,x=trait, y="Density") +
    bigstatsr::theme_bigstatsr(size.rel=0.8)+
    ggplot2::geom_rug(alpha=0.8, color="skyblue")

  if(!is.null(selGen)){
    ## subset the dataset for the selected genotypes
    dat.sel <- BV[BV[[genoCol]] %in% selGen,]
    p <- p+ ggplot2::geom_rug(data=dat.sel, mapping= aes(x=.data[[trait]], color=.data[[colorCol]]),
                              linewidth=1, sides="b", inherit.aes = FALSE)

    if(is.numeric(dat.sel[[colorCol]])){
      p <- p +
        ggrepel::geom_label_repel(data=dat.sel,
                                  mapping = aes(x=.data[[trait]],y=0,
                                                label=.data[[genoCol]],fill=.data[[colorCol]]),
                                  color="black",direction="y",point.padding = 0.01,
                                  force_pull = 0.1,size=3,
                                  vjust=0.6, max.overlaps = Inf,
                                  fontface="bold",alpha=0.7) +
        ggplot2::scale_fill_viridis_d(begin=0.35, direction=1,option="H")+
        ggplot2::scale_color_viridis_d(begin=0.35, direction=1,option="H")
    } else {
      p <- p +
        ggrepel::geom_label_repel(data=dat.sel,
                                  mapping = aes(x=.data[[trait]],y=0,
                                                label=.data[[genoCol]],fill=.data[[colorCol]]),
                                  color="black", direction="y",point.padding = 0.01,
                                  force_pull = 0.1,size=3,
                                  vjust=0.6, max.overlaps = Inf,
                                  fontface="bold",alpha=0.7)+
        ggplot2::scale_fill_viridis_d(begin=0.35, direction=1,option="H")+
        ggplot2::scale_color_viridis_d(begin=0.35, direction=1,option="H")
    } ## end if colorCol
  } ## end if selGen
  plot(p)
  return(p)
}


#' Interactively select rows from a table with filters and download
#'
#' @param data data frame with numeric variables
#' @param vars.inc character vector of variable names sought to be increased to include in the color gradient
#' @param vars.dec character vector of variable names sought to be decreased to include in the color gradient
#' @param output.cols character vector of variable names to include in the output file, if not provided, all columns will be included
#'
#' @returns a shiny app with a table and a download button
#' @author Charlotte Brault
#' @importFrom shiny shinyApp fluidPage titlePanel downloadButton reactive
#' @importFrom DT datatable formatSignif formatStyle styleInterval DTOutput renderDT
#' @importFrom htmlwidgets JS
#' @export

SelectFromTable <- function(data, vars.inc=NULL, vars.dec=NULL, output.cols=NULL){
  num_cols <- c(as.numeric(which(sapply(data, class) == "numeric")))

  if(is.null(output.cols)) output.cols <- colnames(data)
  vars <- c(vars.inc, vars.dec)
  stopifnot(all(vars %in% colnames(data)))
  data[,vars] <- apply(data[,vars],2, as.numeric)

  ## set character columns to factor
  data <- data %>% dplyr::mutate(dplyr::across(dplyr::where(is.character), as.factor))

  ## if not defined, select all numeric columns
  if(is.null(vars.inc) & is.null(vars.dec)){
    vars.inc <- sapply(data, function(x) is.numeric(x) & !is.integer(x))
    vars.dec <- NULL
  }

  # Create selectable table
  dt <- data %>%
    as.data.frame() %>%
    DT::datatable(data, rownames = FALSE, extensions = 'Buttons',
                  filter=list(position="top",clear=F,selection = "multiple"),
                  options = list(scrollX = TRUE,
                                 autoWidth = TRUE,
                                 pageLength=7,
                                 dom = '<<t>Bp>',
                                 buttons = c('copy', 'excel', 'pdf', 'print'),

                                 columnDefs = list(list(className = 'dt-center', targets = "_all")),
                                 render=htmlwidgets::JS(
                                   "function(data, type, row, meta) {",
                                   "return type === 'display' && data.length > 10 ?",
                                   "'<span title=\"' + data + '\">' + data.substr(0, 10) + '...</span>' : data;",
                                   "}")
                  )
    )

  if(length(num_cols) > 0){
    dt <- DT::formatSignif(dt, columns = num_cols, digits = 3)
  }

  #print_table(digits = 3,pageLength=10, filter=list(position="top",clear=F,selection = "multiple"))

  pal <- "ggthemes::Temperature Diverging" ## "grDevices::Temps"
  #"ggthemes::Red-Green-Gold Diverging"

  # Loop through the columns formatting according to that column's distribution
  for (tr in vars) {
    # Create breaks for shading column values high to low
    brks <- stats::quantile(x <- data[[tr]], probs = seq(.05, .95, .01), na.rm = TRUE)
    # Create shades of green for backgrounds
    y <- round(seq(255, 40, length.out = length(brks) + 1), 0)
    #clrs <- paste0("rgb(", y, ", 255,", y, ")")
    cols <- paletteer::paletteer_c(pal,n=length(brks)+1, direction=ifelse(tr %in% vars.inc,-1,1))
    # Format cells in j-th column
    dt <- DT::formatStyle(dt, tr, backgroundColor = DT::styleInterval(brks, cols))
    rm(brks,cols)
  }

  ## create the shinyApp part
  shiny::shinyApp(
    ui=shiny::fluidPage(
      shiny::titlePanel("Filter and Select Rows in Table"),
      DT::DTOutput("table"),
      shiny::downloadButton("download", "Download Selected genotypes")
    ),
    server=function(input, output, session) {
      output$table <- DT::renderDT({dt} ,server = FALSE)

      selected_ids <- shiny::reactive({
        input$table_rows_selected
      })

      output$download <- downloadHandler(
        filename = function() { "selected_ids.csv" },
        content = function(file) {
          utils::write.csv(data[selected_ids(), output.cols, drop = FALSE], file, row.names = FALSE)
        }
      )
    }
  )

}


## ----- Format T3 -----

#' Create table and file to add accession to T3
#'
#' @param dat data frame with genotype information
#' @param checkDB logical, whether to check the T3 database with `BrAPI` package to exclude accessions already in the database, default is TRUE.
#' @param p2f path to file to write the output of the function, need to be with an xlsx or csv extension
#' @param return_table logical, whether to return the table or not, default is TRUE.
#' @param accession_name character string, name of the column in `dat` that contains the accession name
#' @param population_name character string, name of the column in `dat` that contains the population name
#' @param organization_name charac6er string, name of the column in `dat` that contains the organization name
#' @param synonym character string, name of the column in `dat` that contains the synonym name
#' @param variety character string, name of the column in `dat` that contains the variety name
#' @param country_of_origin character string, name of the column in `dat` that contains the country of origin
#' @param notes character string, name of the column in `dat` that contains the notes
#' @param accession_number character string, name of the column in `dat` that contains the accession number
#' @param purdy_pedigree character string, name of the column in `dat` that contains the purdy pedigree
#' @param filial_generation character string, name of the column in `dat` that contains the filial generation
#' @param species_name character string, name of the species, default is "Triticum aestivum" (fixed).
#'
#' @seealso \code{\link{create_trials_T3}}, \code{\link{create_phenot_T3}}
#'
#' @return a data frame with the accessions to add to T3 if return_table is TRUE, otherwise a file is written to the path `p2f`
#' @importFrom dplyr distinct relocate
#' @importFrom utils write.csv
#' @importFrom tools file_ext
#' @author Charlotte Brault
#' @export
create_accessions_T3 <- function(dat=NULL, checkDB=TRUE, p2f=NULL,return_table=TRUE,
                                 accession_name=NULL, population_name=NULL,
                                 organization_name=NULL, synonym=NULL,
                                 variety= NULL,country_of_origin=NULL,
                                 notes= NULL, accession_number=NULL,
                                 purdy_pedigree=NULL,filial_generation=NULL,
                                 species_name="Triticum aestivum"){

  ## Verifications
  stopifnot(!is.null(dat), !is.null(p2f) |!return_table,
            !is.null(accession_name), !is.null(population_name),
            length(accession_name) == 1, length(population_name) == 1,
            is.data.frame(dat))


  colsFormat <- c(accession_name,population_name ,
                  organization_name,synonym,
                  variety,country_of_origin,
                  notes,accession_number ,purdy_pedigree ,
                  filial_generation)
  names(colsFormat) <- c("accession_name","population_name" ,
                         "organization_name(s)","synonym(s)",
                         "variety(s)","country_of_origin(s)","notes(s)",
                         "accession_number(s)","purdy_pedigree","filial_generation")

  colsGeno <- stats::na.omit(colsFormat)
  ## verify that all column written exists in the table
  stopifnot(all(colsGeno %in% colnames(dat)))
  dat.geno <- dat[,colsGeno] %>% dplyr::distinct()
  ## change column names to what is expected in T3
  colnames(dat.geno) <- plyr::mapvalues(colnames(dat.geno), from=stats::na.omit(colsFormat),
                                        to=names(stats::na.omit(colsFormat)))

  ## add other empty columns
  for(col in names(colsFormat)){
    if(is.na(colsFormat[col])){
      dat.geno[[col]] <- ""
    }
  }
  #head(dat.geno)

  T3.accession <- cbind(dat.geno,species_name=species_name)
  T3.accession <- dplyr::relocate(T3.accession,"species_name", .after="accession_name")

  ## check and subset to accessions that are not in the database
  if(checkDB){
    print("Check if accessions are already in the database")
    conn <- BrAPI::getBrAPIConnection("T3/Wheat")

    searchGeno = conn$post("search/germplasm",
                           body=list(germplasmNames=as.list(T3.accession$accession_name)))

    # extract the search id from the first response
    id = searchGeno$content$result$searchResultsDbId
    # get the search results
    searchGeno = conn$get(paste0("search/germplasm/", id), pageSize=10000)$data
    ## check if the accessions are already in the database
    identical(purrr::map_chr(searchGeno, "germplasmName"),
              purrr::map_chr(searchGeno, "defaultDisplayName"))
    ## find the accessions that are in accession_list and not in T3
    toAdd <- setdiff(T3.accession$accession_name,purrr::map_chr(searchGeno, "germplasmName"))
    length(toAdd) ## number of accessions to add
    print(paste0("Number of accessions to add: ",length(toAdd)))
    ## keep only those accessions to the final table
    T3.accession <- T3.accession[T3.accession$accession_name %in% toAdd,]

  }
  ## write to file
  if(!is.null(p2f)){
    if(tools::file_ext(p2f) == "xlsx"){
      openxlsx::write.xlsx(x=T3.accession, file=p2f,colNames = TRUE)
    } else if(tools::file_ext(p2f) %in% c("csv")){
      utils::write.csv(file=p2f, T3.accession,
                       fileEncoding = "UTF-8", row.names=FALSE)

    } else {
      stop("File extension not recognized. Use xlsx or csv.")
    }
    print(paste0("Accessions successfully written to ",p2f))
  }

  if(return_table){
    return(T3.accession)
  }
}




#' Create table and file to add trials to T3
#'
#' @param dat data frame with trial information
#' @param checkDB logical, whether to check the T3 database with `BrAPI` package to exclude trials already in the database, default is TRUE.
#' @param p2f path to file (with xlsx extension) to write the output of the function
#' @param return_table logical, whether to return the table or not, default is TRUE.
#' @param location_state_corresp list, a named list with the correspondence between state (or Canadian province) and location, default is NULL. If not provided, a default list is used.
#' @param keepID character string, ID column to keep in the table (returned not exported)
#' @param year character string, name of the column in `dat` that contains the year, mandatory.
#' @param location character string, name of the column in `dat` that contains the location, mandatory.
#' @param accession_name character string, name of the column in `dat` that contains the accession name, mandatory.
#' @param trial_name character string, name of the column in `dat` that contains the trial name, if not present, it is computed as `paste0(prefix_BP,"_",year,"_",location)`.
#' @param transplanting_date character string, name of the column in `dat` that contains the transplanting date.
#' @param plot_number character string, name of the column in `dat` that contains the plot number. If not present, increasing number from 1000.
#' @param planting_date character string, name of the column in `dat` that contains the planting date.
#' @param harvest_date character string, name of the column in `dat` that contains the harvest date.
#' @param block_number character string, name of the column in `dat` that contains the block number, if not present, filled with 1.
#' @param is_a_control character string, name of the column in `dat` that contains the control status, filled with 0/1.
#' @param rep_number character string, name of the column in `dat` that contains the rep number, if not present, filled with 1.
#' @param range_number character string, name of the column in `dat` that contains the range number.
#' @param row_number character string, name of the column in `dat` that contains the row number.
#' @param col_number character string, name of the column in `dat` that contains the column number.
#' @param plot_width character string, name of the column in `dat` that contains the plot width in meters.
#' @param plot_length character string, name of the column in `dat` that contains the plot length in meters.
#' @param field_size character string, name of the column in `dat` that contains the field size.
#' @param seedlot_name character string, name of the column in `dat` that contains the seedlot name.
#' @param num_seed_per_plot character string, name of the column in `dat` that contains the number of seed per plot.
#' @param prefix_BP character string, prefix for the breeding program, default is "URSN".
#' @param weight_gram_seed_per_plot character string, name of the column in `dat` that contains the weight of seed per plot in grams.
#' @param is_private character string, name of the column in `dat` that contains the private status, set as 1 if private.
#' @param design_type character string, name of the column in `dat` that contains the design type, default is "RCBD".
#' @param trial_type character string, name of the column in `dat` that contains the trial type, default is "phenotyping_trial". Should be one of: Seedling Nursery,
#'  phenotyping_trial, Advanced Yield Trial, Preliminary Yield Trial, Uniform Yield Trial,
#'  Variety Release Trial, Clonal Evaluation, genetic_gain_trial, storage_trial, heterosis_trial,
#' health_status_trial, grafting_trial, Screen House, Seed Multiplication, crossing_block_trial, Specialty Trial
#' @param description character string, name of the column in `dat` that contains the description of the trial, default is "Cooperative nursery of scab trials".
#' @param breeding_program character string, name of the column in `dat` that contains the breeding program, default is "Regional Scab Nursery Cooperative".
#'
#' @return a data frame with the trials to add to T3 if return_table is TRUE, otherwise a file is written to the path `p2f`
#' @seealso \code{\link{create_accessions_T3}}, \code{\link{create_phenot_T3}}
#' @importFrom dplyr distinct relocate select
#' @importFrom utils write.csv
#' @importFrom tools file_ext
#' @author Charlotte Brault
#' @export

create_trials_T3 <- function(dat=NULL,
                             checkDB=TRUE,
                             p2f=NULL,
                             return_table=TRUE,
                             location_state_corresp=NULL,
                             keepID=NULL,
                             year=NULL,
                             location=NULL,
                             accession_name=NULL,
                             trial_name=NULL,
                             transplanting_date=NA,
                             plot_number=NA,
                             planting_date=NA,
                             harvest_date=NA,
                             block_number=NA,
                             is_a_control=NA,
                             rep_number=NA,
                             range_number=NA,
                             row_number=NA,
                             col_number=NA,
                             plot_width=NA,
                             plot_length=NA,
                             field_size=NA,
                             seedlot_name=NA,
                             num_seed_per_plot=NA,
                             weight_gram_seed_per_plot=NA,
                             prefix_BP="URSN",
                             is_private=NA,
                             design_type="RCBD",
                             trial_type="phenotyping_trial",
                             description="Cooperative nursery of scab trials",
                             breeding_program="Regional Scab Nursery Cooperative"){


  ## Verifications
  stopifnot(!is.null(dat), !is.null(p2f) |!return_table,
            !is.null(year), !is.null(location), !is.null(accession_name),
            is.data.frame(dat))

  if(is.null(location_state_corresp)){
    location_state_corresp <-
      list("MN"=c("St. Paul","Crookston", "Morris","Fergus Falls","Sabin","Waseca"),
           "ND"=c("Fargo","Prosper","Langdon","Carrington", "Minot","Casselton",
                  "Williston","Hettinger","Thompson","Berthold","Dickinson","Forman"),
           "SD"=c("Brookings", "Selby","Groton","Aberdeen","Madison","Highmore",
                  "Letcher","Faulkton",
                  "Watertown","Redfield","Miller","Agar","Claire City","Selby","Dakota Lakes"),
           "MB, Canada"=c("Morden","Glenlea", "Winnipeg"),
           "NE"=c("Sidney_NE","Mead")
      )
  }

  ## format the columns that correspond to columns in dat.trial
  colsFormat <- c(year,location,accession_name,trial_name,plot_number,
                  planting_date, harvest_date,transplanting_date,
                  block_number, is_a_control, rep_number, range_number,
                  row_number, col_number, plot_width, plot_length, field_size,
                  seedlot_name,num_seed_per_plot,
                  weight_gram_seed_per_plot,is_private,description,trial_type)
  names(colsFormat) <- c("year","location","accession_name","trial_name","plot_number",
                         "planting_date", "harvest_date","transplanting_date",
                         "block_number", "is_a_control", "rep_number",
                         "range_number","row_number", "col_number",
                         "plot_width", "plot_length", "field_size",
                         "seedlot_name","num_seed_per_plot",
                         "weight_gram_seed_per_plot","is_private","description","trial_type")



  colsTrial <- stats::na.omit(colsFormat)

  ## verify that all columns provided exist in the table
  if(!all(colsTrial %in% colnames(dat))){
    sdf <- setdiff(colsTrial, colnames(dat))
    print(paste0("Columns not found in the data: ", paste(sdf, collapse=", ")))
  }

  dat.trial <- dat %>%
    dplyr::select(tidyselect::all_of(c(keepID,colsTrial))) %>%
    dplyr::distinct()


  colnames(dat.trial) <- plyr::mapvalues(colnames(dat.trial), from=stats::na.omit(colsFormat),
                                         to=names(stats::na.omit(colsFormat)), warn_missing = FALSE)


  ## define a unique id for trial name: paste location and year
  if(is.null(trial_name)){
    dat.trial$trial_name <- paste0(prefix_BP,"_",dat.trial$year,"_",dat.trial$location)
    ### Replace string St. Paul by StPaul
    dat.trial$trial_name <- gsub(pattern="St. ", replacement="St", x=dat.trial$trial_name)
  }


  ## to execute only if plot number is not available: compute it as sequential number
  if(is.na(plot_number)){
    dat.trial <-  purrr::map_dfr(split(dat.trial, dat.trial$trial_name),
                                 function(x) cbind(x, plot_number=seq(1000,length.out=nrow(x))))
  }

  ## add plot name by pasting trial name and plot number
  dat.trial$plot_name <- paste0(dat.trial$trial_name,"_PLOT",
                                dat.trial$plot_number)
  stopifnot(!any(duplicated(dat.trial$plot_name)))

  ## associate locations with states
  for(st in names(location_state_corresp)){
    locs <- location_state_corresp[[st]]
    if(length(locs) > 0){
      idx <- which(dat.trial$location %in% locs)
      if(length(idx) > 0){
        dat.trial$location[idx] <- paste0(dat.trial$location[idx],", ",st)
      }
    }
  }

  idxMiss <- grep(pattern=", ", dat.trial$location, invert = TRUE)
  if(length(idxMiss) > 0){
    print(paste("Missing location states: ",unique(dat.trial$location[idxMiss])))
  }


  if(checkDB){
    selected_breeding_program <- breeding_program
    conn <- BrAPI::getBrAPIConnection("T3/Wheat")
    resp <- conn$get("/studies", query=list(programName=selected_breeding_program), page="all")
    ## extract relevant information, using purrr package
    df.T3.trial <- data.frame(trialName=purrr::map_chr(list_flatten(resp$data),"studyName"),
                              trialId=purrr::map_chr(list_flatten(resp$data),"studyDbId"),
                              locName=purrr::map_chr(list_flatten(resp$data),"locationName"),
                              year=unlist(purrr::list_flatten(
                                purrr::map(purrr::list_flatten(resp$data),"seasons"))))
    # plantDate=purrr::map_chr(list_flatten(resp$data),"startDate"))

    dat.trial <- dat.trial[!dat.trial$trial_name %in% df.T3.trial$trialName,]

    ## remove the name of state or country
    df.T3.trial$locName2 <- stringr::str_remove(string = df.T3.trial$locName, pattern="(, [a-zA-Z]+)+$")
    ## check if there is already trials for the specific year and location that has not been removed
    idx <- which(dat.trial$location %in% df.T3.trial$locName2 &
                   dat.trial$year %in% df.T3.trial$year)
    if(length(idx) > 0){
      print(paste("Trials already in the database: ",
                  paste(unique(dat.trial$trial_name[idx]),", ")))
      dat.trial <- dat.trial[-idx,]
    }


    ## test if trial is already present in T3
    stopifnot(!any(dat.trial$trial_name %in% df.T3.trial$trialName))
  } # end if checkDB



  ## format planting and harvest date columns in date format
  ### date should be in format yyyy-mm-dd or empty
  if(!is.na(planting_date)){
    dat.trial$planting_date <- as.character(as.Date(dat.trial$planting_date,format="%Y-%m-%d"))
    dat.trial$planting_date[is.na(dat.trial$planting_date)] <- ""
  } else{
    dat.trial$planting_date <- ""
  }
  if(!is.na(harvest_date)){
    dat.trial$harvest_date <- as.character(as.Date(dat.trial$harvest_date,format="%Y-%m-%d"))
    dat.trial$harvest_date[is.na(dat.trial$harvest_date)] <- ""
  } else{
    dat.trial$harvest_date <- ""
  }

  ## set as blank if column is not present, set as 1 for block_number and rep_number
  for(col in names(colsFormat)){
    if(col %in% c("plot_number","plot_name")){
      next
    } else{
      if(is.na(colsFormat[col])){
        dat.trial[[col]] <- ""
        if(col %in% "block_number"){
          dat.trial[[col]] <- 1
        }
        if(col %in% "rep_number"){
          dat.trial[[col]] <- 1
        }
      }
    }
  }

  ## for control column, replace TRUE/FALSE with 1/0
  if(!is.na(colsFormat["is_a_control"])){
    if(is.logical(dat.trial$is_a_control)){
      dat.trial$is_a_control <- dplyr::case_match(
        dat.trial$is_a_control, TRUE ~ 1,FALSE ~ 0)
      # dat.trial$is_a_control <- plyr::mapvalues(dat.trial$is_a_control,
      #                                           from=c(TRUE,FALSE),to=c(1,0))
    }
  }



  colN <- c("trial_name","breeding_program","location","year",
            "transplanting_date","design_type",
            "description","trial_type","plot_width","plot_length","field_size",
            "planting_date","harvest_date","plot_name","accession_name",
            "plot_number","block_number","is_a_control","rep_number",
            "range_number","row_number","col_number","seedlot_name",
            "num_seed_per_plot","weight_gram_seed_per_plot","is_private")
  if(!is.null(keepID)){colN <- c(keepID,colN)}

  T3.trial <- cbind(dat.trial,#[,!colnames(dat.trial) %in% "ID"],
                    breeding_program=breeding_program,
                    design_type=design_type#,#description=description,
                    # trial_type=trial_type
  )

  ## reorder columns
  stopifnot(all(colnames(T3.trial) %in% colN),
            all(colN %in% colnames(T3.trial)))
  T3.trial <- T3.trial[,colN]
  head(T3.trial)

  stopifnot(T3.trial$design_type %in% c("","CRD","RCBD","RRC","DRRC","ARC","Alpha",
                                        "Lattice","Augmented","MAD","greenhouse",
                                        "splitplot","p-rep","Westcott"),
            is.numeric(T3.trial$plot_number),
            T3.trial$block_number != "")

  if(!is.null(p2f)){
    if(tools::file_ext(p2f) == "xlsx"){
      openxlsx::write.xlsx(x=T3.trial[,!colnames(T3.trial) %in% keepID],
                           file=p2f,colNames = TRUE)
    } else if(tools::file_ext(p2f) == "csv"){
      utils::write.csv(x=T3.trial[,!colnames(T3.trial) %in% keepID],
                       file=p2f, fileEncoding="UTF-8", row.names = FALSE)
    } else {
      stop("File extension must be either xlsx or csv")
    }
    print(paste0("Trials successfully written to ",p2f))
  }
  if(return_table){
    return(T3.trial)
  }
}


#' Create table and file to add phenotype data in T3
#'
#' @param dat data frame with phenotypic data, containing at least the plot name and the traits
#' @param df_corresp_trait data frame with the correspondence between the column names in `dat` and the traits in T3, must countain the columns `T3_trait` and `corresp_col`.
#' @param p2f path to file to write the output of the function, need to be with an xlsx or csv extension.
#' @param return_table logical, whether to return the table or not, default is TRUE.
#' @param plot_name_col character string, name of the column in `dat` that contains the plot name, default is "plot_name"
#' @seealso \code{\link{create_accessions_T3}}, \code{\link{create_trials_T3}}
#' @return a data frame with the phenotypes to add to T3 if return_table is TRUE, otherwise a file is written to the path `p2f`.
#' @importFrom stats na.omit
#' @importFrom openxlsx write.xlsx
#' @importFrom dplyr select distinct
#' @importFrom utils write.csv
#' @importFrom tools file_ext
#' @author Charlotte Brault
#' @export

create_phenot_T3 <- function(dat=NULL, df_corresp_trait=NULL, p2f=NULL,
                             return_table=TRUE, plot_name_col=NULL){

  df_corresp_trait <- stats::na.omit(df_corresp_trait)

  ## Verifications
  stopifnot(!is.null(df_corresp_trait), !is.null(p2f) |!return_table,
            is.data.frame(df_corresp_trait), is.data.frame(dat),
            all(df_corresp_trait$corresp.col %in% colnames(dat)),
            all(c("T3_trait","corresp_col") %in% colnames(df_corresp_trait)))

  ## select ID column and traits that are in the correspondence table
  T3.phenot <- dat %>%
    dplyr::select(tidyselect::all_of(c(plot_name_col,
                                       intersect(df_corresp_trait$corresp_col, colnames(dat))))) %>%
    dplyr::distinct()
  ## set column names to what is expected in T3
  colnames(T3.phenot) <- plyr::mapvalues(colnames(T3.phenot),
                                         from=df_corresp_trait$corresp_col,
                                         to=df_corresp_trait$T3_trait)
  traits <- intersect(colnames(T3.phenot),df_corresp_trait$T3_trait)

  ## set traits as numeric
  T3.phenot[,traits] <- apply(T3.phenot[,traits],2, as.numeric)

  ## verify some units
  idxPerc <- grep("%",colnames(T3.phenot))

  if(length(idxPerc) > 0){
    stopifnot(all(apply(T3.phenot[,idxPerc,drop=FALSE], 2, min, na.rm=T) >= 0),
              all(apply(T3.phenot[,idxPerc,drop=FALSE], 2, max, na.rm=T) <= 100),
              all(apply(T3.phenot[,traits,drop=FALSE],2, is.numeric)))
  }

  ## export to Excel file
  if(!is.null(p2f)){
    if(tools::file_ext(p2f) == "xlsx"){
      openxlsx::write.xlsx(x=T3.phenot, file=p2f,colNames = TRUE)
    } else if(tools::file_ext(p2f) == "csv"){
      utils::write.csv(x=T3.phenot, file=p2f, fileEncoding="UTF-8", row.names = FALSE)
    } else {
      stop("File extension must be either xlsx or csv")
    }
    print(paste0("Phenotypes successfully written to ",p2f," for ",
                 length(traits)," traits"))
  }
  if(return_table){
    return(T3.phenot)
  }

}







## ----- Plot palette ----


#' Color palette for locations
#'
#' @returns a vector of colors for the locations
#' @export
#' @examples
#' palette=pal_loc()[c("GLL","MRD","LGD","CRK","CRT","PRP","FRG","BRK","STP")]
pal_loc <- function(){
  return(c("GLL"="#AA3939",
           "MRD"="#FFA91B",
           "LGD"="#FFE400",
           "CRK"="#2FD21E",
           "CRT"="#0C2576",
           "PRP"="#D4F74B",
           "FRG"="#FD242A",
           "BRK"="#79DCDC",
           "STP"="#7829E0",
           ## locations specific to URSN
           "SWC"="#C2AD00",
           "BRD"="#9C0058",
           "WNP"="#77B500",
           "WLT"="#3610E1",
           "SDN"="#130AAF",
           "PLM"="#FB7200",
           "BZM"="#4977E2",
           "PWL"="#02BADC",
           "HTG"="#017E95",
           "SLB"="#ED6E94",
           "MNT"="#723900",
           "TPS"="#00941A",
           "FRM"="#FFF709",
           "GRT"="#BA0038",
           "MRS"="khaki4"
  ))
}


#' Color palette for random colors
#'
#' @param n number of colors to return, default is 25
#' @export
#' @returns a vector of colors
#'
#' @examples
#' palette=pal_rand(n=5)
pal_rand <- function(n=25){
  pal <-  c("#36BA55","#FAB9AC","#57D2D9","#F5586A","#F0BF41",
            "#BFC240","#3A2DE0","#6F8CA8","#8B1F1D","#DF9F58",
            "#80944E","#DCD6B2","#1962A0","#1A1C23","#735bcc",
            "#826C37","#565654","#6A79B0","#465946","#d7f96e",
            "dodgerblue","ivory3","magenta","#61393C","cyan")
  return(pal[1:n])
}

#' Color palette for breeding programs
#'
#' @returns a vector of colors for the breeding programs
#' @export
#' @examples
#' palette=pal_orga()[c("UMN","NDSU","SDSU","MSU","Check","Unknown")]
pal_orga <- function(){
  return(c("UMN"="#ffcc33",#"#7A0018CE",#"#7a0019",
           "UMN_BP"="#ffcc33",
           "NDSU"="#2E8B57",#127A5DF8",
           "NDSU_BP"="#36BA55",
           "SDSU"="#3A5FCD",#"#0642CC",
           "MSU"="#bd9b04",#"#BAA15BE2",#"#B9975B",
           "AgriPro"="#8A42C9",
           "Syngenta"="#F5586A",
           "Trigen"="#F54954",
           "AAFC_CRC"="#46688F",
           "Limagrain"="#61393C",
           "AAFC_AAC"="#43AABA",
           "AAFC_SPARC"="#79ADF2",
           "AAFC"="dodgerblue",
           "WestBred"="lightpink2",
           "UIDO"="#d7f96e",#"#f1b300",
           "Ukn_URSN"="ivory2",#"#EEAEEE",
           "all"="grey67","Ukn"="grey67","Unknown"="grey67",
           "Check"="#1A1C23")
  )
}

#' Color palette for the 21 wheat chromosomes
#'
#' @returns a vector of colors for the 21 wheat chromosomes
#' @export
#' @examples
#' palette=pal_21wheatchr()[c(1:21)]
pal_21wheatchr <- function(){
  return(c("#8B0000", "#FF4500", "#FA8072", "#27408B", "#4169E1", "#1C86EE",
           "#008B45", "#00CD66", "#00FF7F", "goldenrod4", "goldenrod2", "gold",
           "mediumpurple3", "mediumpurple1", "mediumorchid2", "#008B8B", "#00CDCD", "#98F5FF",
           "#4D4D4D", "#8C8C8C", "#B5B5B5"))
}
