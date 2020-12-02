########################################
#File: media2016.R
#Author: Jim Murphy
#Desc: some basic descriptive analyses by media
#consumption for 2016 wave (only wave that the
#media items are available)
########################################

##########
##########      LIBRARIES
##########

library(DBI)
library(RSQLite)
library(data.table)
library(ggplot2)

##########
##########      FUNCTIONS
##########

#=======================================
#f_recode_idsorting
#function for doing the recodes of correlations between
#attitudes and partisan or ideological IDs. Intended
#for embedding in f_prep_data
#=======================================

f_recode_idsorting <- function(p_dt){
  
  temp <- p_dt
  temp[ get("dEconDomain")==1, eval("domain") := "Economic"]
  temp[ get("dFrgnpDomain")==1, eval("domain") := "Security"]
  temp[ get("dMoralDomain")==1, eval("domain") := "Moral"]
  temp[ get("dCivRightsDomain")==1, eval("domain") := "Civil Rights"]
  
  invisible(temp)
}

#=======================================
#f_prep_data
#function for preparing data for graphing for a given demographic subgrouping
#p_table_pref   =   prefix of table containing data in question
#p_group_order  =   vector with levels of subgroups in the order
#                   desired for graphing
#=======================================
f_prep_data <- function(p_table_pref, p_group_order = NULL){
  #read in the data sets
  dt_issalign <- setDT(dbReadTable(conn=db_project,
                                   name=paste0(p_table_pref,"_issue_align")))
  
  dt_ideoIss <- setDT(dbReadTable(conn=db_project,
                                  name=paste0(p_table_pref,"_issue_ideo_align")))
  
  dt_prtyIss <- setDT(dbReadTable(conn=db_project,
                                  name=paste0(p_table_pref,"_issue_partyID_align")))
  
  #recodes of the partisan and ideological ID variables
  dt_ideoIss <- f_recode_idsorting(dt_ideoIss)
  dt_prtyIss <- f_recode_idsorting(dt_prtyIss)
  
  #recodes for inter-issue correlations
  dt_issalign[dBothEcon==1, domain := "Economic"]
  dt_issalign[dBothFrgn==1, domain := "Security"]
  dt_issalign[dBothMoral==1, domain := "Moral"]
  dt_issalign[dBothCivr==1, domain := "Civil Rights"]
  dt_issalign[dMixedDomain==1, domain := "Cross-Domain"]
  
  #subsetting to a consistent set of variables
  vec_subvars <- c("var1", "var2", "pair_cor", "group", "jointName", "domain")
  
  dt_issalign <- dt_issalign[ , ..vec_subvars]
  dt_prtyIss <- dt_prtyIss[ , ..vec_subvars]
  dt_ideoIss <- dt_ideoIss[ , ..vec_subvars]
  
  #variable saying what the outcome variable is
  dt_issalign[ , outcome := "Inter-issue"]
  dt_prtyIss[ , outcome := "Party ID"]
  dt_ideoIss[ , outcome := "Ideo ID"]
  
  dt_combo <- rbindlist(list(dt_issalign, dt_prtyIss, dt_ideoIss), use.names = TRUE)
  
  #ensure order of outcomes
  vec_outcome_order <- c("Inter-issue", "Party ID", "Ideo ID")
  dt_combo[ , eval("outcome") := factor(get("outcome"),
                                        levels = vec_outcome_order)]
  
  #get rid of security issues and re-order the domain variable for graphing
  vec_dom_order <- c("Civil Rights", "Economic", "Moral", "Cross-Domain")
  dt_combo <- dt_combo[get("domain") != "Security", ]
  dt_combo <- dt_combo[ , eval("domain") := factor(get("domain"),
                                                   levels = vec_dom_order)]
  
  #put groups in desired order
  if(is.null(p_group_order) == FALSE){
    
    dt_combo <- dt_combo[ , eval("group") := factor(get("group"),
                                                    levels = p_group_order)]
  }
  
  invisible(dt_combo)
}

#=======================================
#f_boxing
#makes paneled boxplots
#p_dt   =   data table
#p_legend_title   =   legend title
#=======================================
f_boxing <- function(p_dt, p_legend_title){
  g_boxes <- ggplot(p_dt, aes_string(x="group", y="pair_cor", fill= "group")) +
    geom_boxplot(alpha = 0.4) + coord_cartesian(ylim = c(-1, 1)) + 
    labs(x = " ", y= "Correlation Coefficient", fill = p_legend_title) +
    facet_grid(rows = vars(outcome),  cols = vars(domain), switch = "y") +
    theme_bw() +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(hjust=0.5),
          legend.position = "bottom",
          strip.background = element_blank(),
          strip.placement = "outside")
}
##########
##########      MAIN CODE
##########

db_project <- dbConnect(SQLite(), "../data/baldiGel_project.db")


#import and prep the data
#interest in news
vec_newsint_order <- c("Little or None", "Moderate", "A Lot or Great Deal")
dt_news_interest <- f_prep_data("bynewsint", vec_newsint_order)

#social media activity
vec_social_order <- c("None", "1 to 6", "Everyday")
dt_social_media <- f_prep_data("bysocial", vec_social_order)

#daily news consumption
dt_days_news <- f_prep_data("bydaysnews")

dbDisconnect(db_project)

#graphing
#two versions of each, one with ideological ID and one without
#this is because I'm a bit concerned that having all three in the graph
#will look too cramped
g_news_interest3 <- f_boxing(dt_news_interest, "Interest in News")
g_news_days3 <- f_boxing(dt_days_news, "Days News per Week")
g_social_media3 <- f_boxing(dt_social_media, "Days Social Media per Week")

g_news_interest2 <- f_boxing(dt_news_interest[outcome != "Ideo ID", ],
                             "Interest in News")
g_news_days2 <- f_boxing(dt_days_news[outcome != "Ideo ID", ],
                         "Days News per Week")
g_social_media2 <- f_boxing(dt_social_media[outcome != "Ideo ID", ],
                            "Days Social Media per Week")

l_graphs <- list(g_news_interest3, g_news_days3, g_social_media3,
                 g_news_interest2, g_news_days2, g_social_media2)
vec_suffix <- c("news_interest_3row", "news_days_3row", "social_media_3row",
                "news_interest_2row", "news_days_2row", "social_media_2row")

for(i in c(1:length(l_graphs))){
  path_png <- paste0("../graphs/png/boxes_", vec_suffix[i], ".png")
  ggsave(plot=l_graphs[[i]], path_png, width = 5.5, height = 6.3, units = "in")
  
  path_svg <- paste0("../graphs/svg/boxes_", vec_suffix[i], ".svg")
  ggsave(plot=l_graphs[[i]], path_svg, width = 5.5, height = 6.3, units = "in")
}