########################################
#File: makeBaldGelData.R
#Desc: as in Baldassarri and Gelman AJS paper using ANES, 
#constructs data set of year-specific correlations
#between item pairs nested within issue
########################################

##########
##########      LIBRARIES
##########

library(DBI)
library(RSQLite)
library(data.table)

##########
##########      CONSTANTS
##########

#names of variables in original data and their new names by category
vec_ovars_core <- c("vcf0004", "vcf0803", "vcf0301")
vec_newnames_core <- c("year","ideoid","partyid")

vec_ovars_design <- c("vcf0017","vcf0009x", "vcf0009y", "vcf0009z")
vec_newnames_design <- c("mode","faceWT","webWT","pooledWT")

vec_ovars_econ <- c("vcf0839", "vcf0886", "vcf0809", "vcf0806", "vcf0887",
                  "vcf9050",  "vcf0893", "vcf0894", "vcf9046", "vcf0890",
                  "vcf9049", "vcf0889",  "vcf0891", "vcf9047")

vec_newnames_econ <- c("fedspend_govserv","fedspend_poor","gov_jobguar",
                     "healthinsur", "fedspend_childcare","fedspend_asstblacks",
                     "fedspend_homeless", "fedspend_welfare","fedspend_foodstamp",
                     "fedspend_schools", "fedspend_socsec","fedspend_aids",
                     "fedspend_college","fedspend_envi")

vec_ovars_civrights <- c("vcf9016", "vcf0867a", "vcf9018", "vcf9042", "vcf9037",
                       "vcf0830",  "vcf9014","vcf9040", "vcf9013","vcf9017",
                       "vcf0817", "vcf0816",  "vcf9039", "vcf9041", "vcf0814",
                       "vcf9015","vcf0813")


vec_newnames_civrights <- c("more_equal", "affirmact", "equal_treat",
                            "blacks_deservemore", "blacks_fair", "aid2blacks",
                            "rights2much", "blacks_specfav", "equal_opp",
                            "worryLessEqual", "school_bus", "school_integ",
                            "hard4blacks", "blacks_shldTry", "civrights2fast",
                            "chancesNEQ", "blacks_poschange")

vec_ovars_moral <- c("vcf0877", "vcf0876a", "vcf0852","vcf0854", "vcf0878",
                   "vcf0853", "vcf0851", "vcf0834", "vcf0837","vcf0838", "vcf9043")

vec_newnames_moral <- c("gayer_military", "protect_homo", "adjust_morals",
                      "tolerate_difvals", "gays_adoptkids", "tradvals",
                      "newLstylesBad", "women_equals", "abortion1", "abortion2",
                      "school_prayer")

vec_ovars_frgnp <- c("vcf0843", "vcf0811", "vcf0841", "vcf0892", "vcf0888","vcf9048")
vec_newnames_frgnp <- c("defense_spend","urban_unrest","ruskies","fedspend_frgnaid",
                      "fedspend_crime","fedspend_space")

vec_all_ovars <- c(vec_ovars_core, vec_ovars_design, vec_ovars_econ,
                   vec_ovars_civrights, vec_ovars_moral, vec_ovars_frgnp)

vec_all_newnames <- c(vec_newnames_core, vec_newnames_design, vec_newnames_econ,
                      vec_newnames_civrights, vec_newnames_moral, vec_newnames_frgnp)

#final version of issue variables that will be used for correlations
#code dummies for issue areas after recodes
vec_econissues  <- c("fedspend_govservRC", "fedspend_poor", "gov_jobguar",
                   "healthinsur", "fedspend_childcare", "fedspend_asstblacks",
                   "fedspend_homeless", "fedspend_welfare", "fedspend_foodstamp",
                   "fedspend_schools", "fedspend_socsec", "fedspend_aids",
                   "fedspend_college", "fedspend_envi")

vec_civrissues <- c("more_equalRC", "affirmact", "equal_treat", "blacks_deservemore",
                    "blacks_fair", "aid2blacks", "rights2muchRC", "blacks_specfavRC",
                    "equal_opp", "worryLessEqualRC", "school_bus", "school_integ",
                    "hard4blacks", "blacks_shldTryRC", "civrights2fast", "chancesNEQ",
                    "blacks_poschange")

vec_moralissues <- c("gayer_military", "protect_homo", "adjust_morals",
                     "tolerate_difvals", "gays_adoptkids", "tradvalsRC",
                     "newLstylesBadRC", "women_equals", "abortion4RC",
                     "school_prayer")

vec_frgnpissues <- c("defense_spend","urban_unrest","ruskies","fedspend_frgnaid",
                   "fedspend_crimeRC","fedspend_space")

vec_issueVars <- c(vec_econissues, vec_civrissues, vec_moralissues, vec_frgnpissues)

##########
##########      FUNCTIONS
##########

#functions that will take a correlation matrix
#and turn it into a data frame for a single year

#=======================================
#f_extract_cell
#extracts cell value if it's on upper triangle of correlation matrix
#Params: p_irow   =   row number
#        p_icol   =   column number
#        p_cormat =   correlation matrix
#        p_boomat =   matrix of logical values indicating whether
#                     the row-column entry should be extracted
#=======================================
f_extract_cell <- function(p_irow, p_icol, p_cormat, p_boomat) {

  #conditional to determine if cell matches 
  if(p_boomat[p_irow, p_icol]==TRUE & is.na(p_cormat[p_irow, p_icol])==FALSE) {
  
    val <- as.numeric(p_cormat[p_irow, p_icol])
    rname <- rownames(p_cormat)[p_irow]
    cname <- colnames(p_cormat)[p_icol]
    
    vec <- c(rname, cname, as.numeric(val))
    invisible(vec)
  }
}

#=======================================
#f_loop_over_row
#Params:  p_irow    =   row of p_cormat to take values from
#         p_cormat  =   correlation matrix that contains values to extract
#         p_boomat  =   matrix of logical values saying whether to extract cell
#=======================================
#extractsall desired cell values from 1 row of correlation matrix
f_loop_over_row <- function(p_irow, p_cormat, p_boomat) {

  l_vals1row <- invisible(lapply(1:ncol(p_cormat),
                          FUN=f_extract_cell,
                          p_irow=p_irow, p_cormat=p_cormat, p_boomat=p_boomat))
  
  #combine those cells together
  mat_fromRow <- do.call(rbind, l_vals1row)
  
}

#=======================================
#f_df_from_cormat
#Params:  p_DT  =   data table containing individual-level item responses
#         p_year  =   year of p_DT to use
#         p_noInternet  =   logical; TRUE excludes Internet respondents
#         p_weighted    =   logical; TRUE applies weights appropriate to the sample
#=======================================
#embeds last two functions to create data frame for correlation matrix
f_df_from_cormat <- function(p_year, p_DT, p_noInternet, p_weighted) {
  
  vec_its2cor <- c("ideoid", "partyid", vec_issueVars)
  
  #weight and case selection based on inclusion/exclusion of Internet sample
  if(p_noInternet==FALSE){
  
    #if keeping Internet respondents, choose pooled weight for potential weighting
    dt_1year <- p_DT[year==p_year, c("pooledWT",vec_its2cor)]
    
    vec_weights <- dt_1year$pooledWT
    
    dt_1year[ , pooledWT := NULL]

  } else{

      #excclude Internet mode
      vec_netkeep <- c("faceWT", vec_its2cor)
      dt_1year <- p_DT[year==p_year & mode != 4, ..vec_netkeep]
      
      #keep face-to-face weight for potential weighting
      vec_weights <- dt_1year$faceWT
      dt_1year[ , faceWT := NULL]
  }

  #weighted or unweighted
  if(p_weighted==TRUE) {
    
    mat_cor <- WGCNA::cor(x=dt_1year[, ..vec_its2cor], y= dt_1year[, ..vec_its2cor],
                            weights.x=vec_weights, weights.y=vec_weights, 
                            use="pairwise.complete.obs")
  } else{
    
    mat_cor <- cor(dt_1year[, ..vec_its2cor], use="pairwise.complete.obs")
  }

  #create a matrix of logical values indicating if entry is on upper triangle
  #done so not to replicate pairs
  mat_isupdi <- upper.tri(mat_cor, diag=FALSE)
  
  l_fromAllRows <- invisible(lapply(1:nrow(mat_cor),
                             FUN= f_loop_over_row,
                             p_cormat=mat_cor,  p_boomat=mat_isupdi))
  
  #stack together; add year and group label before converting to data frame
  mat_fullyear <- do.call(rbind, l_fromAllRows)
  mat_fullyear <- cbind(mat_fullyear, p_year)
  df_fullyear <- as.data.frame(mat_fullyear, stringsAsFactors=FALSE)
  
  invisible(df_fullyear)
}

##########
##########			MAIN CODE
##########

dt_orig <- fread("../data/anes_timeseries_cdf_rawdata.txt", header=TRUE, sep="|")

setnames(dt_orig, old = names(dt_orig), new = tolower(names(dt_orig)))

#check if duplicates in keeper vars
if(length(unique(vec_all_ovars)) != length(vec_all_ovars)) {
  stop("ERROR: DUPLICATE VARIABLES ENTERED")
}

#the second condition excludes web respondents in 2016
dt_sub72 <- dt_orig[vcf0004 >= 1972, ..vec_all_ovars]

rm(dt_orig)

setnames(dt_sub72, vec_all_ovars, vec_all_newnames)

#convert numeric missing codes to NAs
vec_logrec <- !(vec_all_newnames %in% c("year", vec_newnames_design))
vec_rec1vars <- vec_all_newnames[vec_logrec]

#convert NA
f_conv2na <- function(p_var) {
	dt_sub72[[p_var]] <<- ifelse(dt_sub72[[p_var]] %in% 1:7,
                               yes=dt_sub72[[p_var]], no=NA)
}

invisible(lapply(vec_rec1vars, f_conv2na))

#additional recoding of variables where 7 is not a valid (or is questionable)
dt_sub72[ , `:=` (affirmact = ifelse(affirmact %in% 1:5, affirmact, NA) )]

f_conv2na_gt3 <- function(p_x){
  dt_sub72[[p_x]] <<- ifelse(dt_sub72[[p_x]] %in% 1:3, yes=dt_sub72[[p_x]], no=NA)
}

vec_naGT3 <- c("fedspend_asstblacks", "fedspend_foodstamp", "fedspend_socsec",
               "fedspend_space")
invisible(lapply(vec_naGT3, f_conv2na_gt3))


#reverse code that what needs reversing so that higer values 
#reflect more conservative attitudes
dt_sub72[ , `:=` (fedspend_govservRC = -1*fedspend_govserv + 8,
                  more_equalRC = -1*more_equal + 6,
                  rights2muchRC = -1*rights2much + 6,
                  blacks_specfavRC = -1*blacks_specfav + 6,
                  worryLessEqualRC = -1*worryLessEqual + 6,
                  blacks_shldTryRC = -1*blacks_shldTry + 6,
                  tradvalsRC = -1*tradvals + 6,
                  newLstylesBadRC = -1*newLstylesBad + 6,
                  abortion1RC = -1*abortion1 + 5,
                  abortion2RC = -1*abortion2 + 5,
                  fedspend_crimeRC = -1*fedspend_crime+ 4) ]

#new abortion variable combining original two
#abortion1 is from before 1980, abortion2 is after 1980
#but they overlapped in 1980 as two highly similar wordings
#two different combined versions: 
#FIRST (3RC) uses post-1980 version unless missing since it is the wording
#for most of the time series. If that's missing use the prior version
#SECOND (4RC) uses the pre-1980 version unless it's missing since Austin
#thinks it gets results closer to B&G's original
dt_sub72[ , `:=` (abortion3RC = ifelse(is.na(abortion2RC)==TRUE,
                                       yes=abortion1RC,  no=abortion2RC),
                  abortion4RC = ifelse(is.na(abortion1RC)==TRUE,
                                       yes=abortion2RC,  no=abortion1RC) ) ]

#remove old variables once recode done to avoid confusion
vec_vars2remove <- c("abortion2", "abortion1", "abortion2RC", "abortion1RC",
                     "newLstylesBad", "tradvals", "blacks_shldTry",
                     "worryLessEqual", "blacks_specfav", "more_equal",
                     "rights2much", "fedspend_govserv")
dt_sub72[ , (vec_vars2remove) := NULL ]

#print out descriptives of the individual-level data set so that it
#can be checked for any weirdness
sink("output/makeBaldGelData.out")
print(Hmisc::describe(dt_sub72))
sink()

#run function over all years to create data frame
#face-to-face only
dt_face72 <- dt_sub72[ mode == 0, ]
vec_years <- unique(dt_face72$year)

l_allyears <- lapply(vec_years, FUN= f_df_from_cormat, p_DT = dt_face72,
                     p_noInternet=TRUE, p_weighted=TRUE)

#stack together
dt_full <- rbindlist(l_allyears)
setnames(dt_full, c("var1","var2","pair_cor","year"))

dt_full[ , c("pair_cor","year") := lapply(.SD, as.numeric),
        .SDcols = c("pair_cor","year")]

dt_full[ , jointName := paste0(var1, "BY", var2)]

#=====    subset into different types of correlations

#issue partisanship via ideology
dt_issAlignIdeo <- dt_full[(var1=="ideoid" | var2=="ideoid")
                          & var1 != "partyid" & var2 != "partyid" , ]
                          
#issue partisanship via party identification
dt_issAlignParty <- dt_full[(var1=="partyid" | var2=="partyid")
                            & var1 != "ideoid" & var2 != "ideoid", ]
#issue alignments
dt_issconst <- dt_full[!(var1 %in% c("partyid","ideoid")
                         | var2 %in% c("partyid","ideoid") ) , ]

#correlation between ideology and partisanship
dt_partyAlignIdeo <- dt_full[ var1 %in% c("partyid","ideoid") 
                              & var2 %in% c("partyid","ideoid"), ]

#=====    DUMMIES FOR ISSUE DOMAINS 

dt_issconst[ , `:=` (dBothEcon = ifelse(var1 %in% vec_econissues 
                                        & var2 %in% vec_econissues, yes=1, no=0),
                     dBothFrgn = ifelse(var1 %in% vec_frgnpissues
                                        & var2 %in% vec_frgnpissues, yes=1, no=0),
                     dBothMoral = ifelse(var1 %in% vec_moralissues
                                          & var2 %in% vec_moralissues, yes=1, no=0),
                     dBothCivr = ifelse(var1 %in% vec_civrissues
                                        & var2 %in% vec_civrissues, yes=1, no=0)), ]
dt_issconst[ , dSameDomain := dBothEcon + dBothFrgn + dBothMoral + dBothCivr ]
dt_issconst[ , dMixedDomain := -1*dSameDomain + 1 ]

#=====    FOR ISSUE ALIGNMENT WITH PARTY ID AND IDEOLOGY

#align with ideology id
dt_issAlignIdeo[ , `:=` ( dFrgnpDomain = ifelse(var1 %in% vec_frgnpissues 
                                                  | var2 %in% vec_frgnpissues,
                                                yes=1, no=0),
                        
                          dCivRightsDomain = ifelse(var1 %in% vec_civrissues 
                                                      | var2 %in% vec_civrissues,
                                                    yes=1, no=0),

                          dEconDomain = ifelse(var1 %in% vec_econissues 
                                                | var2 %in% vec_econissues,
                                              yes=1, no=0),
                        
                          dMoralDomain = ifelse(var1 %in% vec_moralissues 
                                                  | var2 %in% vec_moralissues,
                                               yes=1, no=0) ) ]
#with party id
dt_issAlignParty[ , `:=` ( dFrgnpDomain = ifelse(var1 %in% vec_frgnpissues 
                                                | var2 %in% vec_frgnpissues,
                                                yes=1, no=0),
                          
                          dCivRightsDomain = ifelse(var1 %in% vec_civrissues 
                                                    | var2 %in% vec_civrissues,
                                                    yes=1, no=0),
                          
                          dEconDomain = ifelse(var1 %in% vec_econissues 
                                               | var2 %in% vec_econissues,
                                               yes=1, no=0),
                          
                          dMoralDomain = ifelse(var1 %in% vec_moralissues 
                                                | var2 %in% vec_moralissues,
                                                yes=1, no=0) ) ]

#====     WRITE TO DB

db_project <- dbConnect(SQLite(), "../data/baldiGel_project.db")

dbWriteTable(db_project, value=dt_issAlignIdeo, name="issue_ideo_align",
             row.names=FALSE, overwrite=TRUE, append=FALSE)

dbWriteTable(db_project, value=dt_issAlignParty, name="issue_partyID_align",
             row.names=FALSE, overwrite=TRUE, append=FALSE)

dbWriteTable(db_project, value=dt_issconst, name="issue_align", row.names=FALSE,
             overwrite=TRUE, append=FALSE)

dbWriteTable(db_project, value=dt_partyAlignIdeo, name="ideo_partyid_align",
             overwrite=TRUE, append=FALSE)

dbDisconnect(db_project)
