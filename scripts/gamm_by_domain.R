########################################
#File: gamm_tests.R
#Author: Jim Murphy
#Desc: generalized additive mixed models
#to capture non-linearities in trends over
#time in issue alignment and constraint
########################################

##########
##########    LIBRARIES
##########

library(data.table)
library(gamm4)
library(DBI)
library(RSQLite)
library(voxel)    #ggplot2 graphing functions for mgcv/gamm4 objects
library(ggplot2)

##########
##########    CONSTANTS
##########

#colors to be used
vec_colors <- c("dodgerblue3", "red4", "forestgreen", "darkorange1", "darkorchid4")

##########
##########    FUNCTIONS
##########

#=======================================
#f_recode_domain
#Desc: recodes categorical variable for issue domain
#      Applies only to partisan id and ideological id data
#=======================================
f_recode_domain <- function(p_input){
  dt_output <- p_input
  dt_output[dEconDomain == 1, domain := "Economic"]
  dt_output[dFrgnpDomain == 1, domain := "Security"]
  dt_output[dMoralDomain == 1, domain := "Moral"]
  dt_output[dCivRightsDomain == 1, domain := "Civil Rights"]
  
  dt_output[ , `:=` (jointName = as.factor(jointName),
                     domain = as.factor(domain))]
  
  invisible(dt_output)
}

#=======================================
#f_plot_fit
#Desc: fits GAMM for given data sets and plots fitted values
#Params: p_model  = data table with data for model
#        p_title  = title for graph
#=======================================
f_plot_fit <- function(p_model, p_title){
  
  #get fitted values
  dt_fitvals <- setDT(plotGAMM(p_model, smooth.cov = "year",
                               groupCovs = "domain")$data)
  dt_fitvals[ , `:=` (upper_fit = fit + 1.96*se.fit,
                      lower_fit = fit - 1.96*se.fit)]
  
  n_cols <- length(unique(dt_fitvals$domain))
  vec_col <- vec_colors[1:n_cols]
  
  g_fit <- ggplot(dt_fitvals, aes(x=year, y=fit, color=domain, linetype=domain)) +
    geom_ribbon(aes(x=year, ymax= upper_fit, ymin= lower_fit, fill = domain),
                alpha=0.2, colour = NA) +
    geom_line() +
    scale_x_continuous("Year", breaks = c(1972,1992,2012),
                       labels = c("'72","'92","'12"),
                       minor_breaks = NULL ) +
    coord_cartesian(ylim=c(-0.1,0.5)) +   #prevents cutoff of CI below ylim
    labs(title = p_title, x="Year", y="Correlation Coefficient",
         linetype="", color="", fill = "") +
    scale_color_manual(name = "", aesthetics = c("color", "fill"),
                       values= vec_colors) +
    annotate("rect",xmin=2004, xmax=Inf, ymin=-Inf, ymax=Inf,
             alpha=0.2, fill="burlywood1") +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5),
          legend.position = "right",
          legend.title=element_blank())
  
  invisible(g_fit)            
}

#=======================================
#f_scat_rfx_time
#diagnostic function to check whether for whether random effects
#are correlated with time, which would violate model assumptions
#Params: p_model_obj  =   model object
#        p_origData   =   original data frame on which model was fit
#        p_title      =   title for plot
#=======================================
f_scat_rfx_time <- function(p_model_obj, p_origData, p_title){
  dt_rfx <- random.effects(p_model_obj$mer)$jointName  
  names(dt_rfx) <- "u_intercept"
  dt_rfx$jointName <- rownames(dt_rfx)
  #merging in the old data to get the years the item was asked
  dt_rfx <- merge(dt_rfx, p_origData[ , c("jointName","year")], by="jointName")
  
  g_scatter_rfx <- ggplot(data=dt_rfx, aes(x=year, y=u_intercept)) +
    geom_point() + ggtitle(p_title) + ylab("Random Effect (Intercept)") +
    theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
    geom_smooth()
  invisible(g_scatter_rfx)                    
}

##########
##########    MAIN CODE
##########

db_project <- dbConnect(SQLite(), "../data/baldiGel_project.db")

dt_issconst <- setDT(dbReadTable(db_project, name = "issue_align"))
dt_ideoalign <- setDT(dbReadTable(db_project, name = "issue_ideo_align"))
dt_prtyalign <- setDT(dbReadTable(db_project, name = "issue_partyID_align"))

dbDisconnect(db_project)

#inter-issue constraint
#turn dummy variables into a categorical
dt_issconst[ dBothEcon == 1, domain := "Economic"]
dt_issconst[ dBothFrgn == 1, domain := "Security"]
dt_issconst[ dBothMoral == 1, domain := "Moral"]
dt_issconst[ dBothCivr == 1, domain := "Civil Rights"]
dt_issconst[ dMixedDomain == 1, domain := "Cross-Domain"]

dt_issconst[ , `:=` (jointName = as.factor(jointName),
                     domain = as.factor(domain) )]

#fit the model and plot fit
#NOTE: k=10 is what s() would choose automatically by default; set here
#explicitly for reproducibility in case this changes with later versions
#of mgcv/gamm4 (actual k=9 because one degree of freedom is lost for
#identifiability constraint)
fit_issconst <- gamm4(pair_cor ~ s(year, by=domain, bs = "tp", k = 10) + domain,
                      random=~(1|jointName), data=dt_issconst)

g_issconst <- f_plot_fit(fit_issconst, "Inter-issue Alignment by Issue Domain")

#partisan id sorting
dt_prtyalign <- f_recode_domain(dt_prtyalign)
fit_prtyalign <- gamm4(pair_cor ~ s(year, by=domain, bs = "tp", k = 10) + domain,
                      random=~(1|jointName), data=dt_prtyalign)
g_prtyalign <- f_plot_fit(fit_prtyalign,
                          "Alignment of Issues and Party ID by Issue Domain")

#ideological id sorting
dt_ideoalign <- f_recode_domain(dt_ideoalign)
fit_ideoalign <- gamm4(pair_cor ~ s(year, by=domain, bs = "tp", k = 10) + domain,
                       random=~(1|jointName), data=dt_ideoalign)
g_ideoalign <- f_plot_fit(fit_ideoalign,
                          "Alignment of Issues and Ideological ID by Issue Domain")

ggsave(plot=g_ideoalign, filename="../graphs/png/gamm_ideosort_wghtF2F.png",
       width = 6.3, height = 6.3, units = "in")
ggsave(plot=g_ideoalign, filename="../graphs/svg/gamm_ideosort_wghtF2F.svg",
       width = 6.3, height = 6.3, units = "in")

ggsave(plot=g_prtyalign, filename="../graphs/png/gamm_prtyid_wghtF2F.png",
       width = 6.3, height = 6.3, units = "in")
ggsave(plot=g_prtyalign, filename="../graphs/svg/gamm_prtyid_wghtF2F.svg",
       width = 6.3, height = 6.3, units = "in")

ggsave(plot=g_issconst, filename="../graphs/png/gamm_issalign_wghtF2F.png",
       width = 6.3, height = 6.3, units = "in")
ggsave(plot=g_issconst, filename="../graphs/svg/gamm_issalign_wghtF2F.svg",
       width = 6.3, height = 6.3, units = "in")

#=====    CHECKING CORRELATION OF RANDOM EFFECTS WITH TIME

g_rfx_ideoalign <- f_scat_rfx_time(p_model_obj= fit_ideoalign,
                                  p_origData=dt_ideoalign,
                                  p_title="Ideological Self ID Model")
g_rfx_interiss <- f_scat_rfx_time(fit_issconst, dt_issconst, "Inter-issue Model")
g_rfx_partyid <- f_scat_rfx_time(fit_prtyalign, dt_prtyalign, "Party ID Model")

pdf("../graphs/rfxTime_gamm_bydomain.pdf")
print(g_rfx_interiss)
print(g_rfx_partyid)
print(g_rfx_ideoalign)
dev.off()
