########################################
#File: subgroupData.R
#Desc: makes the correlation data by demographic
#subgroup in Baldassarri and Gelman
########################################
##########
##########      LIBRARIES
##########

library(DBI)
library(RSQLite)
library(data.table)

##########
##########      NOTES
##########

#LINE OF BALDASSARRI's ORIGINAL CODE
#sub.cut <- 4 #interest =3; edu=5; FamInc 1-2=low(1-33 percentile),
#4-5=high (68-100 percentile);
#church.attendance(1,2=high (almost every week or more), 
#3-5=low from once or twice a month to never; 
#pol.activism (camp.part=1 low, else high); PartyID3 (Dem=1 (low); Rep=3(high)) 

#Spending: note an ambiguity with the federal spending on crime variable
#while in original form it is negatively related to conservatism, it is also
#positively correlated with being a republican...
#While it seems like one might want to reverse code it, we do not, as attempts
#to replicate B&G's pre-2004 results were closer to their original using the
#crime spending variable's original ordering.

##########
##########      CONSTANTS
##########

#names of variables in original data by category and new name
vec_ovars_core <- c("vcf0004", "vcf0803", "vcf0301")
vec_newnames_core <- c("year","ideoid","partyid")

vec_ovars_design <- c("vcf0006", "vcf0017","vcf0009x", "vcf0009y", "vcf0009z")
vec_newnames_design <- c("caseid", "mode","faceWT","webWT","pooledWT")

vec_ovars_subgrps <- c("vcf0101", "vcf0104","vcf0106","vcf0112", "vcf0310",
                       "vcf0723", "vcf0140", "vcf0114", "vcf0130","vcf0303",
                       "vcf0050a")
vec_newnames_subgrps <- c("age", "gender", "race", "region", "interestElections",
                          "campact6", "edu7", "faminc5", "religAttend",
                          "partyid3", "polknow_interview")

vec_ovars_econ <- c("vcf0839", "vcf0886", "vcf0809", "vcf0806", "vcf0887",
                    "vcf9050", "vcf0893", "vcf0894", "vcf9046", "vcf0890",
                    "vcf9049", "vcf0889", "vcf0891", "vcf9047")
vec_newnames_econ <- c("fedspend_govserv", "fedspend_poor", "gov_jobguar",
                       "healthinsur", "fedspend_childcare", "fedspend_asstblacks",
                       "fedspend_homeless", "fedspend_welfare", "fedspend_foodstamp",
                       "fedspend_schools", "fedspend_socsec", "fedspend_aids",
                       "fedspend_college", "fedspend_envi")

vec_ovars_civrights <- c("vcf9016", "vcf0867a", "vcf9018", "vcf9042", "vcf9037", 
                       "vcf0830", "vcf9014","vcf9040", "vcf9013","vcf9017",
                       "vcf0817", "vcf0816", "vcf9039", "vcf9041", "vcf0814",
                       "vcf9015","vcf0813")
vec_newnames_civrights <- c("more_equal", "affirmact", "equal_treat",
                          "blacks_deservemore", "blacks_fair", "aid2blacks",
                          "rights2much", "blacks_specfav", "equal_opp",
                          "worryLessEqual", "school_bus", "school_integ",
                          "hard4blacks", "blacks_shldTry", "civrights2fast",
                          "chancesNEQ", "blacks_poschange")

vec_ovars_moral <- c("vcf0877", "vcf0876a", "vcf0852","vcf0854", "vcf0878", "vcf0853", 
                   "vcf0851", "vcf0834", "vcf0837","vcf0838", "vcf9043")
vec_newnames_moral <- c("gayer_military","protect_homo","adjust_morals",
                      "tolerate_difvals","gays_adoptkids","tradvals","newLstylesBad",
                      "women_equals","abortion1","abortion2","school_prayer")

vec_ovars_frgnp <- c("vcf0843", "vcf0811", "vcf0841", "vcf0892", "vcf0888","vcf9048")
vec_newnames_frgnp <- c("defense_spend","urban_unrest","ruskies","fedspend_frgnaid",
                      "fedspend_crime","fedspend_space")

vec_all_ovars <- c(vec_ovars_core, vec_ovars_design, vec_ovars_subgrps, vec_ovars_econ,
                 vec_ovars_civrights, vec_ovars_moral, vec_ovars_frgnp)
vec_all_newnames <- c(vec_newnames_core, vec_newnames_design, vec_newnames_subgrps,
                    vec_newnames_econ,vec_newnames_civrights, vec_newnames_moral,
                    vec_newnames_frgnp)

#final version of issue variables that will be used for correlations
#code dummies for issue areas
vec_econissues  <- c("fedspend_govservRC","fedspend_poor", "gov_jobguar",
                     "healthinsur", "fedspend_childcare","fedspend_asstblacks",
                     "fedspend_homeless", "fedspend_welfare","fedspend_foodstamp",
                     "fedspend_schools", "fedspend_socsec","fedspend_aids",
                     "fedspend_college","fedspend_envi")
vec_civrissues <- c("more_equalRC","affirmact","equal_treat","blacks_deservemore",
                    "blacks_fair","aid2blacks","rights2muchRC","blacks_specfavRC",
                    "equal_opp","worryLessEqualRC","school_bus","school_integ",
                    "hard4blacks","blacks_shldTryRC","civrights2fast","chancesNEQ",
                    "blacks_poschange")
vec_moralissues <- c("gayer_military", "protect_homo", "adjust_morals",
                    "tolerate_difvals", "gays_adoptkids", "tradvalsRC",
                    "newLstylesBadRC", "women_equals", "abortion4RC","school_prayer")
vec_frgnpissues <- c("defense_spend","urban_unrest","ruskies","fedspend_frgnaid",
                   "fedspend_crimeRC","fedspend_space")

vec_issueVars <- c(vec_econissues, vec_civrissues, vec_moralissues, vec_frgnpissues)



##########
##########      FUNCTIONS
##########

#functions that will take a correlation matrix and turn it into a data frame
#for a single year

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
  if(p_boomat[p_irow, p_icol]==TRUE & is.na(p_cormat[p_irow, p_icol])==FALSE) {
    val <- as.numeric(p_cormat[p_irow, p_icol])
    rname <- rownames(p_cormat)[p_irow]
    cname <- colnames(p_cormat)[p_icol]
    
    vec <- c(rname, cname, as.numeric(val))
    return(vec)
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
#embeds last two functions to create data frame for correlation matrix
#Params:  p_DT  =   data table containing individual-level item responses
#         p_year  =   year of p_DT to use
#         p_noInternet  =   logical; TRUE excludes Internet respondents
#         p_weighted    =   logical; TRUE applies weights appropriate to the sample
#=======================================

f_df_from_cormat <- function(p_DT, p_cat, p_year, p_noInternet, p_weighted=TRUE) {
  
  vec_its2cor <- c("ideoid","partyid",vec_issueVars)

  #weight and case selection based on inclusion/exclusion of Internet sample
  if(p_noInternet==FALSE){
  
    #if keeping Internet respondents, choose pooled weight for potential weighting
    dt_1year <- p_DT[year==p_year, c("pooledWT",vec_its2cor)]
    
    vec_weights <- dt_1year$pooledWT
    
    dt_1year[ , pooledWT := NULL]

  } else{
      
      #excclude Internet mode
      vec_netkeep <- c("faceWT", vec_its2cor)
      dt_1year <- p_DT[year == p_year & mode != 4, ..vec_netkeep]
      
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
  
  #stacks results, one row per pairwise correlation for this year
  mat_fullyear <- do.call(rbind, l_fromAllRows)

  #add year and group label before converting to data frame
  mat_fullyear <- cbind(mat_fullyear, p_year, p_cat)
  df_fullyear <- as.data.frame(mat_fullyear, stringsAsFactors=FALSE)
  
  invisible(df_fullyear)
}

#=======================================
#f_corsBySubgroup
#makes and creates database table classifying thingies by subgroup
#p_DT       =   data table containing respondent-level ANES data
#p_groupvar   =   grouping variable  
#p_tablepref  =   prefix to attach to table that will be written to database
#p_noInternet =   logical; if TRUE exclude Internet respondents
#p_weighted   =   logical; if TRUE apply sample-appropriate survey weights 
#=======================================
f_corsBySubgroup <- function(p_DT, p_groupvar, p_tablepref, p_noInternet,
                             p_weighted) {
  #only people who have valid values on the grouping variable
  dt_nonNAs <- p_DT[is.na(get(p_groupvar)) == FALSE, ]

  #vector with all possible categories
  vec_vals <- unique(dt_nonNAs[ , get(p_groupvar)])
  
  #vector of years; passed to f_cors1cat below; derived outside the function to avoid
  #possible errors stemming from years without any respondents with given value
  vec_years <- unique(dt_nonNAs$year)
  
  #get correlations for one category of the grouping variable
  #define function locally here instead of globally, outside of f_corsBySubgroup
  #so that one can plug in other variables defined within f_corsBySubgroup
  #and avoid a parameterization of f_cors1cat that would lead to lots
  #of recursion-induced errors
  f_cors1cat <- function(p_cat) {
    
    dt_1subgroup <- dt_nonNAs[ get(p_groupvar) == p_cat, ]
  
    #get rid of the subgroup variables b/c I don't want to accidentally
    #calculate their correlations with stuff
    dt_1subgroup[, (vec_subgroupVars) := NULL] 
    
    #calculate subgroup
    l_allyears1cat <- lapply(vec_years, FUN= f_df_from_cormat,
                             p_DT = dt_1subgroup, p_noInternet = p_noInternet,
                             p_weighted = p_weighted, p_cat = p_cat)
    
    #stack together
    dt_allyears1cat <- rbindlist(l_allyears1cat)
    
    invisible(dt_allyears1cat)
  }
  
  #will return a list of matrices 
  l_allcats <- lapply(vec_vals, FUN= f_cors1cat)
  
  #put together into a data table
  dt_full <- rbindlist(l_allcats)
  setnames(dt_full, c("var1","var2","pair_cor","year","group"))

  #functions will have turned year and correlation into string
  #convert back to numeric
  dt_full[ , c("pair_cor","year") := lapply(.SD, as.numeric),
        .SDcols = c("pair_cor","year")]

  dt_full[ , jointName := paste0(var1, "BY", var2)]

  #subset into different types of correlations
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
  
  #ideoological ID
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

  #party ID
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
  
  #write to database
  db_project <- dbConnect(SQLite(), "../data/baldiGel_project.db")
  
  dbWriteTable(db_project, value=dt_issAlignIdeo,
               name=paste0(p_tablepref,"_issue_ideo_align"), 
               row.names=FALSE, overwrite=TRUE, append=FALSE)

  dbWriteTable(db_project, value=dt_issAlignParty,
              name=paste0(p_tablepref,"_issue_partyID_align"), 
              row.names=FALSE, overwrite=TRUE, append=FALSE)
  
  dbWriteTable(db_project, value=dt_issconst,
               name=paste0(p_tablepref,"_issue_align"),
               row.names=FALSE, overwrite=TRUE, append=FALSE)

  dbDisconnect(db_project)
}

##########
##########      MAIN CODE
##########


#=====      DATA SET-UP

dt_orig <- fread("../data/anes_timeseries_cdf_rawdata.txt", header=TRUE, sep="|")

setnames(dt_orig, old = names(dt_orig), new = tolower(names(dt_orig)))

#check if duplicates in keeper vars
if(length(unique(vec_all_ovars)) != length(vec_all_ovars)) {
  stop("ERROR: DUPLICATE VARIABLES ENTERED")
}

dt_sub72 <- dt_orig[vcf0004 >= 1972, ..vec_all_ovars]

setnames(dt_sub72, vec_all_newnames)

#everything that's not year or a design variable needs an NA recoding
vec_needs_rec <- !(vec_all_newnames %in% c("year",vec_newnames_design, "age"))
vec_rec1vars <- vec_all_newnames[vec_needs_rec]

f_conv2na <- function(p.var) {
  dt_sub72[[p.var]] <<- ifelse(dt_sub72[[p.var]] %in% 1:7,
                              yes=dt_sub72[[p.var]], no=NA)
}

invisible(lapply(vec_rec1vars, f_conv2na))

#recode of age
dt_sub72[ age < 17, age := NA]

#additional recoding of variables where 7 is not a valid (or is questionable)
dt_sub72[ , `:=` (affirmact = ifelse(affirmact %in% 1:5, affirmact, NA) )]

f_conv2na_gt3 <- function(p_x){
  dt_sub72[[p_x]] <<- ifelse(dt_sub72[[p_x]] %in% 1:3, yes=dt_sub72[[p_x]], no=NA)
}

vec_naGT3 <- c("fedspend_asstblacks", "fedspend_foodstamp", "fedspend_socsec",
               "fedspend_space")
invisible(lapply(vec_naGT3, f_conv2na_gt3))

#reverse code that which needs reversing so that higer values 
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
                  fedspend_crimeRC = -1*fedspend_crime + 4) ]

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

#=====      SETUP SUPPLEMENTARY 2016 MEDIA CONSUMPTION DATA

#import 2016 data
vec_orig16 <- c("V160001", "V161495", "V161008", "V161009")
vec_new16 <- c("caseid", "social_media", "days_news", "interest_news")
dt_16 <- fread("../data/anes2016/anes_timeseries_2016_rawdata.txt")

dt_16 <- dt_16[ , ..vec_orig16]
setnames(dt_16, old = vec_orig16, new = vec_new16)

#define new variables
#days in typical week on social media
#the following categories break roughly into thirds
dt_16[social_media == 0, new_social_media := "None"]
dt_16[social_media %in% c(1:6), new_social_media := "1 to 6"]
dt_16[social_media == 7, new_social_media := "Everyday"]

#days in typical week watch/read/listen news
#60 percent claim to consume news every day. Other splits comprise
#about 20 percent each
dt_16[days_news %in% c(1:3), new_days_news := "0 to 3"]
dt_16[days_news %in% c(4:6), new_days_news := "4 to 6"]
dt_16[days_news == 7, new_days_news := "Everyday"]

#interest in news
dt_16[interest_news %in% c(4:5), new_news_interest := "Little or None"]
dt_16[interest_news == 3, new_news_interest := "Moderate"]
dt_16[interest_news %in% c(1,2), new_news_interest := "A Lot or Great Deal"]

dt_16[ , `:=` (social_media = NULL, days_news = NULL,
               interest_news = NULL, year = 2016)]

#=====      DEFINING SUBGROUPS ON DEMOGRAPHICS

#new categorical variables for different sorts of subgroups
dt_sub72[gender==1, newGender := "Male"]
dt_sub72[gender==2, newGender := "Female"]

dt_sub72[race==1, newRace := "White"]
dt_sub72[race==2, newRace := "Black"]

dt_sub72[, newRegion := ifelse(region == 3, yes = "South", no= "Non-south")]

#two options for interest variable
#one that distinguishes into three groups (B&G originals were High and Low)
#because of small number of low (especially in 2012 and 2016), the second
#option combines mid and low
dt_sub72[interestElections==3, newInterest1 := "Very"]
dt_sub72[interestElections==2, newInterest1 := "Somewhat"]
dt_sub72[interestElections==1, newInterest1 := "Not Much"]

dt_sub72[ , newInterest2 := ifelse(newInterest1=="Very", yes="High", no="Low")]

#division of campaign activities
dt_sub72[campact6==1, newActivist := "None"]
dt_sub72[campact6==2, newActivist := "1"]
dt_sub72[campact6 >= 3, newActivist := "2+"]

dt_sub72[edu7 %in% c(1,2,3,4), newCollege := "No College"]
dt_sub72[edu7 %in% c(5,6,7), newCollege := "College"]

dt_sub72[faminc5 %in% c(4,5), newIncome := "68-100"]
dt_sub72[faminc5==3, newIncome := "34-67"]
dt_sub72[faminc5 %in% c(1,2), newIncome := "0-33"]

dt_sub72[polknow_interview %in% c(1,2), newKnow := "High"]
dt_sub72[polknow_interview == 3, newKnow := "Average"]
dt_sub72[polknow_interview %in% c(4,5), newKnow := "Low"]

dt_sub72[age %in% c(17:39), age_group := "17-39"]
dt_sub72[age %in% c(40:59), age_group := "40-59"]
dt_sub72[age >= 60, age_group := "60+"]

### religion
#high group is "almost every week" or "every week"
#Note that there is no 6
dt_sub72[religAttend %in% c(1,2), newReligAttend := "Almost Weekly or More"]
dt_sub72[religAttend %in% c(3:4), newReligAttend := "Less than Weekly"]
dt_sub72[religAttend %in% c(5, 7), newReligAttend := "Never"]                  

dt_sub72[partyid3==1, newParty := "Dem"]
dt_sub72[partyid3==3, newParty := "GOP"]
dt_sub72[partyid3==2, newParty := "Ind"]

#=====      FINALIZE DATA

#add in supplementary 2016 variables
dt_complete <- merge(dt_sub72, dt_16, by = c("caseid", "year"), all.x = TRUE)

#remove old variables once recode done to avoid confusion
vec_vars2remove <- c("caseid", "abortion2", "abortion1", "abortion2RC",
                     "abortion1RC", "newLstylesBad", "tradvals",
                     "blacks_shldTry", "worryLessEqual", "blacks_specfav",
                     "more_equal", "rights2much",  "fedspend_govserv", "race",
                     "gender", "region", "interestElections", "edu7", "faminc5",
                     "religAttend", "partyid3", "polknow_interview")

dt_complete[ , (vec_vars2remove) := NULL]

#vectors of subgroup variables to use and corresponding tables prefices 
vec_subgroupVars <- c("newGender", "newRace", "newRegion", "newCollege",
                      "newActivist", "newInterest1", "newInterest2", "newIncome",
                      "newReligAttend", "newParty", "newKnow", "age_group",
                      "new_news_interest", "new_days_news", "new_social_media")
vec_tablePrefs <- c("bygender", "byrace", "byregion", "bycollege", "byactivist",
                    "byinterest1", "byinterest2", "byincome", "byrelig", "byparty",
                    "bypolknow", "byage", "bynewsint", "bydaysnews", "bysocial")

#print out descriptives before moving to correlations
sink("output/subgroupData.out")
print(Hmisc::describe(dt_complete))
sink()

#NOTE:  WEIGHTING WILL INCREASE THE NUMBER OF GROUPS WITH VARIANCE (VIRTUALLY) ZERO
#limit to face to face sample only
dt_face2face <- dt_complete[ mode == 0, ]

f_doThings <- function(p_i) {
  f_corsBySubgroup(p_DT = dt_face2face, p_groupvar=vec_subgroupVars[p_i],
                   p_tablepref=vec_tablePrefs[p_i], p_weighted=TRUE,
                   p_noInternet=TRUE )
}

invisible(lapply(1:length(vec_subgroupVars), FUN=f_doThings))
