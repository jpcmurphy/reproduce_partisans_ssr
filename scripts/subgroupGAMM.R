##########################################
#subgroupGAMM.R
#File: GAMMS of for inter-issue alignment, partisan sorting,
#and ideological ID sorting by demographichic subgroups
##########################################

##########
##########		LIBRARIES
##########

library(gamm4)
library(voxel)    #ggplot2 graphing functions for mgcv/gamm4 objects
library(ggplot2)
library(DBI)
library(RSQLite)
library(data.table)


##########
##########      CONSTANTS
##########

#proper orderings for graphing of the different model types and issue domains
vec_modlabs <- c("Inter-issue", "Party ID", "Ideo ID")
vec_domlabs <- c("Civil Rights", "Economic", "Moral", "Cross-Domain")


##########
##########    	FUNCTIONS
##########

#===========================================
#f_fitmod_1group
#function that is embedded in f_gammByGroup
#it actually runs the GAMM models for one subgroup
#of demographic
#p_val:       one value of the grouping variable indicating the subgroup
#             that the model is for (e.g. black, college-educated, etc.)
#p_prtyDF:    data frame containing the information on alignment of issue
#             attitudes with partisan self-identification
#p_ideoDF:    data frame containing information on alignment of issue
#             attitudes with ideological self-identification
#p_issIssDF:  data frame containing inter-issue alignment
#===========================================
f_fitmod_1group <- function(p_val, p_prtyDF, p_ideoDF, p_issIssDF) {
  
  #NOTE: k=10 is what s() would choose automatically by default; set here
  #explicitly for reproducibility in case this changes with later versions
  #of mgcv/gamm4 (actual k=9 because one degree of freedom is lost for
  #identifiability constraint)
  fitgamm_ideosort <- gamm4(pair_cor ~ s(year, by=domain, bs = "tp", k = 10)
                                      + domain,
                            random = ~(1|jointName),
                            data=subset(p_ideoDF, group==p_val))
  fitgamm_issconst <- gamm4(pair_cor ~ s(year, by=domain, bs = "tp", k = 10)
                                      + domain,
                            random=~(1|jointName),
                            data=subset(p_issIssDF, group==p_val))
  #this conditional exists for when I run models by party ID
  #Modeling party ID sorting model separately by party ID would crash
  #since it's meaningless. The conditional allows one to run models
  #for issue constraint and ideological ID sorting fg.inor Dems, Reps, and indies.
  if(is.null(p_prtyDF) == FALSE){
    fitgamm_prtysort <- gamm4(pair_cor ~ s(year, by=domain, bs = "tp", k = 10)
                                        + domain, 
                              random=~(1|jointName),
                              data=subset(p_prtyDF, group==p_val))
    l_mods1value <- list(fitgamm_prtysort, fitgamm_ideosort, fitgamm_issconst)
    names(l_mods1value) <- c("prty","ideo","iss")
  }else{
    l_mods1value <- list(fitgamm_ideosort, fitgamm_issconst)
    names(l_mods1value) <- c("ideo","iss")
  }
  invisible(l_mods1value)
}

#===========================================
#f_gammByGroup
#function that runs GAMM models (via f_fitmod_1group) for all groups
#for a given demographic variable
#p_tablePref:   prefix for database table containing data of interest
#===========================================
f_gammByGroup <- function(p_tablePref, p.selectYears=NULL){

  dt_issconst <- setDT(dbReadTable(conn=db_project,
                                   name=paste0(p_tablePref,"_issue_align")))

	dt_ideoIss <- setDT(dbReadTable(conn=db_project,
	                               name=paste0(p_tablePref,"_issue_ideo_align")))
	
	#creating factor variables for domains
	#ideological self-ID
	dt_ideoIss[ dEconDomain==1, domain := "Economic"]
	dt_ideoIss[ dFrgnpDomain==1, domain := "Security"]
	dt_ideoIss[ dMoralDomain==1, domain := "Moral"]
	dt_ideoIss[ dCivRightsDomain==1, domain := "Civil Rights"]
	dt_ideoIss[ , domain := as.factor(dt_ideoIss$domain)]

	#inter-issue alignment
	dt_issconst[dBothEcon==1, domain := "Economic"]
	dt_issconst[dBothFrgn==1, domain := "Security"]
	dt_issconst[dBothMoral==1, domain := "Moral"]
	dt_issconst[dBothCivr==1, domain := "Civil Rights"]
	dt_issconst[dMixedDomain==1, domain := "Cross-Domain"]
	dt_issconst[ , domain := as.factor(domain)]
	
	#alignment of issues with partyID is meaningless if the subgroupings are
	#partisan self-identifications so only do that if party ID is NOT the grouping factor
	if(p_tablePref != "byparty"){
	  dt_prtyIss <- setDT(dbReadTable(conn=db_project,
	                            name=paste0(p_tablePref,"_issue_partyID_align")))

	  dt_prtyIss[dEconDomain==1 , domain := "Economic"]
	  dt_prtyIss[dFrgnpDomain==1, domain := "Security"]
	  dt_prtyIss[dMoralDomain==1, domain := "Moral"]
	  dt_prtyIss[dCivRightsDomain==1, domain := "Civil Rights"]
	  dt_prtyIss[ , domain := as.factor(domain) ]
	}
	
	#=====    model fitting
	
	vec_group_vals <- unique(dt_issconst$group) #values of demographic variable
	
	if(p_tablePref != "byparty"){
	  l_all_values <- lapply(vec_group_vals, FUN = f_fitmod_1group,
	                         p_prtyDF=dt_prtyIss, p_ideoDF=dt_ideoIss,
	                         p_issIssDF=dt_issconst)
	  names(l_all_values) <- c(vec_group_vals)
	}else{
	  l_all_values <- lapply(vec_group_vals, FUN = f_fitmod_1group,
	                         p_prtyDF=NULL, p_ideoDF=dt_ideoIss,
	                         p_issIssDF=dt_issconst)
	  names(l_all_values) <- c(vec_group_vals)
	}
	
	invisible(l_all_values)
}

#=======================================
#f_fitvals_1group
#creates big ol' data frame of fitted values
#from ONE demographic group
#p_model_list: list containing gamm4 objects
#p_demo_label: label for the demographic subgrouping variable
#=======================================
f_fitvals_1group <- function(p_model_list, p_demo_label){
  
  #use plotGAMM function to get a data frame of fitted values
  #for each outcome
  dt_fit_ideo <- setDT(plotGAMM(p_model_list[["ideo"]], smooth.cov = "year",
                                groupCovs = "domain")$data)
  dt_fit_ideo[ , model_outcome := "Ideo ID"]
  
  dt_fit_iss <- setDT(plotGAMM(p_model_list[["iss"]], smooth.cov = "year",
                               groupCovs = "domain")$data)
  dt_fit_iss[ , model_outcome := "Inter-issue"]
  
  #party one only applies if party contained in the list
  if("prty" %in% names(p_model_list) ){
    
    dt_fit_prty <- setDT(plotGAMM(p_model_list[["prty"]], smooth.cov = "year",
                                  groupCovs = "domain")$data)
    dt_fit_prty[ , model_outcome := "Party ID"]

    #put into one data frame (will subset when actually graphing)
    dt_fit_allmodels <- rbindlist(list(dt_fit_ideo, dt_fit_prty, dt_fit_iss))

  }else{
    dt_fit_allmodels <- rbindlist(list(dt_fit_ideo, dt_fit_iss))
  }
  
  #get rid of securoty since it won't be included in the graphs;
  #small number of items makes it noisy
  dt_fit_allmodels <- dt_fit_allmodels[ domain != "Security", ]
  
  #variables indicating the demo and upper and lower bounds 
  dt_fit_allmodels[ , `:=` (demo = p_demo_label, upper_ci = fit + 1.96*se.fit,
                            lower_ci = fit - 1.96*se.fit)]

  invisible(dt_fit_allmodels)
}

#=======================================
#f_fitvals_1demovar
#uses f_fitvals_1group to create big data frame
#of all models' fitted values for all groups of
#a given demographic variable
#p.list_of_list: list of lists that contain models
#=======================================
f_fitvals_1demovar <- function(p_list_of_lists) {
  
  #vector of the demographic groups contained in the list of lists
  vec_demo_labs <- names(p_list_of_lists)
  
  #loop over all demographic groups of specified grouping variables
  f_makeAllFitted <- function(p_i, p_list_of_lists, p_all_labels) {
    
    l_all_groupdt <- f_fitvals_1group(p_model_list = p_list_of_lists[[p_i]],
                                      p_demo_label = p_all_labels[p_i])
    
    invisible(l_all_groupdt)
  }
  
  l_all_dts <- lapply(c(1:length(p_list_of_lists)), FUN = f_makeAllFitted,
                      p_list_of_lists=p_list_of_lists, p_all_labels=vec_demo_labs)
  
  dt_all <- rbindlist(l_all_dts)
  
  #creates a factor version of the demographic variable that retains original
  #order of demographic categories for labeling purposes instead of defaulting
  #to alphabetical order
  dt_all[ , demo_label := ordered(demo, levels=vec_demo_labs)]
  
  #similarly for domain
  dt_all[ , domain_label := ordered(domain, levels=vec_domlabs)]
  invisible(dt_all)
}

#=======================================
#f_graph_fitted
#creates a graph of item-pair alignment by demographic
#group for a given outcome
#p_DT   =   data frame containing fitted values
#p_legendTitle  =   legend title
#=======================================
f_graph_fitted <- function(p_DT, p_legendTitle){
  
  dt_copy <- copy(p_DT)
  
  #for proper ordering
  #note I'm choosing the model outcome order this way so that
  #the first row of the grid will be inter-issue alignment and
  #domain order so that cross-domain will be last
  #makes things prettier
  vec_modlabs <- c("Inter-issue","Party ID","Ideo ID")

  dt_copy[ , model_outcome := ordered(model_outcome, levels = vec_modlabs)]
  dt_copy[ , domain_label := ordered(domain_label, levels = vec_domlabs)]
  
  #generate colors
  vec_colors <- RColorBrewer::brewer.pal(length(vec_domlabs), "Set2")
  
  ggplot(data= dt_copy,
         aes(year, fit, color=demo_label, linetype=demo_label, fill=demo_label)) +
    geom_ribbon(aes(x=year, ymax= upper_ci, ymin= lower_ci, fill = demo_label),
                alpha=0.3, colour = NA) +
    geom_line() +
    scale_color_manual(name = p_legendTitle, values = vec_colors,
                       aesthetics = c("color", "fill")) +
    scale_x_continuous("Year", breaks = c(1972,1992,2012),
                       labels = c("'72","'92","'12"), minor_breaks = NULL ) +
    facet_grid(rows = vars(model_outcome),  cols = vars(domain_label), switch = "y") +
    labs(linetype=p_legendTitle, color=p_legendTitle, y = "Correlation Coefficient") + 
    annotate("rect",xmin=2004, xmax=Inf, ymin=-Inf, ymax=Inf,
             alpha=0.2, fill="burlywood1") +
    theme_bw() +
    theme(legend.position = "bottom", panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), plot.title = element_text(hjust=0.5),
          strip.background = element_blank(),
          strip.placement = "outside")
}

##########
##########    	MAIN CODE
##########

db_project <- dbConnect(SQLite(), "../data/baldiGel_project.db")

#run models
l_mods4age <- f_gammByGroup("byage")
l_mods4gender <- f_gammByGroup("bygender")
l_mods4race <- f_gammByGroup("byrace")
l_mods4region <- f_gammByGroup("byregion")
l_mods4income <- f_gammByGroup("byincome")
l_mods4college <- f_gammByGroup("bycollege")
l_mods4religAttend <- f_gammByGroup("byrelig")
l_mods4activist <- f_gammByGroup("byactivist")
l_mods4interest1 <- f_gammByGroup("byinterest1")
l_mods4interest2 <- f_gammByGroup("byinterest2")
l_mods4polknow <- f_gammByGroup("bypolknow")
l_mods4party <- f_gammByGroup("byparty")

dbDisconnect(db_project)

#creates data frames of fitted values
dt_age_fitvals <- f_fitvals_1demovar(l_mods4age)
dt_gender_fitvals <- f_fitvals_1demovar(l_mods4gender)
dt_activ_fitvals <- f_fitvals_1demovar(l_mods4activist)
dt_coll_fitvals <- f_fitvals_1demovar(l_mods4college)
dt_inc_fitvals <- f_fitvals_1demovar(l_mods4income)
dt_race_fitvals <- f_fitvals_1demovar(l_mods4race)
dt_region_fitvals <- f_fitvals_1demovar(l_mods4region)
dt_relig_fitvals <- f_fitvals_1demovar(l_mods4religAttend)
dt_interest_fitvals <- f_fitvals_1demovar(l_mods4interest1)
dt_polknow_fitvals <- f_fitvals_1demovar(l_mods4polknow)
dt_byprty_fitvals <- f_fitvals_1demovar(l_mods4party)

#re-order as necessary so higher values are on the left
dt_polknow_fitvals[ , demo_label := factor(demo_label,
                                           levels = c("High", "Average", "Low"))]

dt_interest_fitvals[ , demo_label := factor(demo_label,
                                            levels = c("Very", "Somewhat",
                                                       "Not Much"))]
dt_activ_fitvals[ , demo_label := factor(demo_label,
                                         levels = c("2+", "1", "None") )]
dt_coll_fitvals[ , demo_label := factor(demo_label,
                                        levels = c("College", "No College") )]
#want age progressing lowest to highest
dt_age_fitvals[ , demo_label := factor(demo_label,
                                       levels = c("17-39", "40-59", "60+"))]

#recode and reorder the religion variables to have line breaks
dt_relig_fitvals[ demo_label == "Almost Weekly or More",
                  demo_label := "Almost Weekly \n or More"]
dt_relig_fitvals[ demo_label == "Less than Weekly",
                  demo_label := "Less than \n Weekly"]
dt_relig_fitvals[ demo_label == "Never", demo_label := "Never"]

vec_relig_order <- c("Almost Weekly \n or More", "Less than \n Weekly", "Never")
dt_relig_fitvals[ , demo_label := factor(demo_label, levels = vec_relig_order)]

#graphing
g_age <- f_graph_fitted(dt_age_fitvals, "")
g_gender <- f_graph_fitted(dt_gender_fitvals, "")
g_activ <- f_graph_fitted(dt_activ_fitvals, "Campaign Participation Activities")
g_college <- f_graph_fitted(dt_coll_fitvals, "")
g_income <- f_graph_fitted(dt_inc_fitvals, "Income (percentile)")
g_race <- f_graph_fitted(dt_race_fitvals, "")
g_region <- f_graph_fitted(dt_region_fitvals, "")
g_relig <- f_graph_fitted(dt_relig_fitvals, "Religious Attendance")
g_interest <- f_graph_fitted(dt_interest_fitvals, "Interested in Elections")
g_polknow <- f_graph_fitted(dt_polknow_fitvals,
                            "Political Knowledge (Interviewer-rated)")
  
#for party identification subgrouping, I'm just going to use inter-issue
#because it makes for a cleaner story
#I'm doing the party ID one outside of the function because (1) I want to specify
#specific colors for each party ID and (2) there's just one set of graphs
#so I don't need to panel

dt_byprty_fitvals[ , model_outcome := ordered(model_outcome, levels = vec_modlabs)]

dt_byprty_fitvals[ , domain_label := ordered(domain_label, levels = vec_domlabs) ]

g_prtyID <- ggplot(data= dt_byprty_fitvals[ model_outcome=="Inter-issue", ],
                   aes(x=year, y=fit, color=demo_label, linetype=demo_label) )  +
            geom_ribbon(aes(x=year, ymax= upper_ci, ymin= lower_ci, fill = demo_label),
                        alpha=0.3, colour = NA) +
            geom_line() +
            facet_grid(rows = vars(model_outcome), cols = vars(domain_label),
                       switch = "y") +
            scale_x_continuous("Year", breaks = c(1972,1992,2012),
                              labels = c("'72","'92","'12"), minor_breaks = NULL ) +
            labs(linetype="Party ID", color="Party ID", y="Correlation Coefficient") +
            scale_color_manual(name = "Party ID", aesthetics = c("color", "fill"),
                               values=c("dodgerblue3", "darkorange1", "firebrick")) +
            annotate("rect",xmin=2004, xmax=Inf, ymin=-Inf, ymax=Inf,
                     alpha=0.2, fill="burlywood1") +
            theme_bw() +
            theme(legend.position = "bottom",
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(hjust=0.5),
                  strip.background = element_blank(),
                  strip.placement = "outside")

#writing to files        
l_graphs <- list(g_age, g_relig, g_race, g_region, g_interest, g_activ,
                 g_gender, g_college, g_income, g_prtyID, g_polknow)
vec_suffix <- c("age", "relig","race", "region", "polint", "polact",
                "gender", "college", "income", "partyid", "polknow")

for(i in c(1:length(l_graphs))){
  path_png <- paste0("../graphs/png/subgamm_", vec_suffix[i], ".png")
  ggsave(plot=l_graphs[[i]], path_png, width = 5.5, height = 6.3, units = "in")

  path_svg <- paste0("../graphs/svg/subgamm_", vec_suffix[i], ".svg")
  ggsave(plot=l_graphs[[i]], path_svg, width = 5.5, height = 6.3, units = "in")
}
