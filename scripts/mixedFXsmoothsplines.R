########################################
#File: spline_tests.R
#Author: Jim Murphy
#Desc: program to run some mixed effects 
#semi-parameteric spline models in search 
#of non-linearities over the time series, 
#which can then inform respecifications of
#the parametric model
########################################

##########
##########    LIBRARIES
##########

library(sme)
library(DBI)
library(RSQLite)
library(ggplot2)

##########
##########    FUNCTIONS
##########

#=======================================
#f_plotSplineFit
#uses elements of mixed effects spline model object
#to make a fitted graph with confidence interval
#Params:  p_fit   =   model object
#         p_graph =   whether to make the graph or just
#                     produce fitted values
#=======================================

f_plotSplineFit <- function(p_fit, p_graph = TRUE) {
  
  vec_Xinputs <- as.numeric(colnames(p_fit$coefficients))
  vec_Yinputs <- p_fit$coefficients[1,] 
  
  #variance
  mu_variance <- diag(vcov(p_fit))
  
  #get values for fitted mean over time
  mu <- spline(x=vec_Xinputs, y=vec_Yinputs, n=500,method="natural")
  
  #get confidence bands
  vec_upperInp <- vec_Yinputs + 1.96 * sqrt(mu_variance)
  vec_lowerInp <- vec_Yinputs - 1.96 * sqrt(mu_variance)
  
  upper_band <- spline(x=vec_Xinputs, y=vec_upperInp, method="natural",n=500)
  lower_band <- spline(x=vec_Xinputs, y=vec_lowerInp, method="natural",n=500)
  
  if(p_graph == TRUE){
    #plot this thing
    g_splineFit <- ggplot() + geom_line(aes(x=mu$x, y=mu$y), color="blue", lwd=1) +
      geom_ribbon(aes(x=mu$x, ymax=upper_band$y, ymin=lower_band$y),
                  color="lightblue", fill="lightblue", alpha=0.4) +
      xlim(1972, 2016) + ylim(0, .45) +
      scale_x_continuous(breaks = c(seq(1976, 2016, by=8), 2016) ) +
      xlab("Year") +ylab("Correlation Coefficient") + theme_bw() +
      annotate("rect",xmin=2004, xmax=Inf, ymin=-Inf, ymax=Inf,
               alpha=0.2, fill="burlywood1") +
      theme(plot.title = element_text(hjust=0.5))
    invisible(g_splineFit)
  }else{
    df_fitted <- as.data.frame(cbind(mu$x, mu$y, upper_band$y, lower_band$y))
    names(df_fitted) <- c("year", "fitted_y", "up95", "low95")
    invisible(df_fitted)
  }    
}

##########
##########    MAIN CODE
##########

#read in that there data from the data base
db_project <- dbConnect(SQLite(), "../data/baldiGel_project.db")
df_issconst <- dbReadTable(db_project, name = "issue_align")  #inter-issue alignment
df_prtyalign <- dbReadTable(db_project, name = "issue_partyID_align") #partisan alignment
df_ideoalign <- dbReadTable(db_project, name = "issue_ideo_align") #ideological alignment
dbDisconnect(db_project)

#since there are very few issue-pairs per year prior to 1984
#I'm going to limit the possible knots to '84 and later
#note that sme automatically includes the start (1974) and
#end (2016) as knots, so this vector should only go up to 2012
vec_potent_knots <- unique(df_issconst$year)
vec_potent_knots <- vec_potent_knots[which(vec_potent_knots %in% c(1984:2012))]

#when running the sme function, the columns of the data frame
#need to be ordered as outcome, year, cluster
fit_issconst <- sme(df_issconst$pair_cor, tme= df_issconst$year,
                    ind = df_issconst$jointName, knots = vec_potent_knots)
g_iss_const <- f_plotSplineFit(fit_issconst)  +
                ggtitle("Inter-issue Alignment, 1972-2016")

###
#repeat for partisan alignment
###
fit_prtyalign <- sme(df_prtyalign$pair_cor, tme= df_prtyalign$year,
                    ind = df_prtyalign$jointName, knots = vec_potent_knots)

g_prty_align <- f_plotSplineFit(fit_prtyalign) +
                ggtitle("Alignment of Issues and Party ID, 1972-2016")

###
#repeat for ideological alignment
###
fit_ideoalign <- sme(df_ideoalign$pair_cor, tme= df_ideoalign$year,
                     ind = df_ideoalign$jointName, knots = vec_potent_knots)

g_ideo_align <- f_plotSplineFit(fit_ideoalign) +
                  ggtitle("Alignment of Issues and Ideological ID, 1972-2016")

#one copy to png, another to svg
ggsave(plot=g_iss_const, filename="../graphs/png/splines_issalign_wghtF2F.png")
ggsave(plot=g_iss_const, filename="../graphs/svg/splines_issalign_wghtF2F.svg")

ggsave(plot=g_prty_align, filename="../graphs/png/splines_prtysort_wghtF2F.png")
ggsave(plot=g_prty_align, filename="../graphs/svg/splines_prtysort_wghtF2F.svg")

ggsave(plot=g_ideo_align, filename="../graphs/png/splines_ideosort_wghtF2F.png")
ggsave(plot=g_ideo_align, filename="../graphs/svg/splines_ideosort_wghtF2F.svg")


#===    version that puts all three on one graph

df_fitted_ideo <- f_plotSplineFit(fit_ideoalign, p_graph = FALSE)
df_fitted_party <- f_plotSplineFit(fit_prtyalign, p_graph = FALSE)
df_fitted_issue <- f_plotSplineFit(fit_issconst, p_graph = FALSE)

df_fitted_ideo$label <- "Ideological Alignment"
df_fitted_party$label <- "Party Alignment"
df_fitted_issue$label <- "Issue Alignment"

df_combo <- rbind(df_fitted_ideo, df_fitted_party, df_fitted_issue)
df_combo$label <- factor(df_combo$label,
                         levels = c("Ideological Alignment", "Party Alignment",
                                    "Issue Alignment") )

#
vec_colors <- c("darkblue", "deepskyblue1", "darkgoldenrod3")
g_combo_spline <- ggplot(data = df_combo,
                         aes(year, fitted_y, color = label,
                             fill = label, linetype = label)) +
                    geom_line(lwd=.7) +
                    geom_ribbon(data = df_combo,
                                aes(x=year, ymax=up95, ymin=low95, fill = label),
                                alpha=0.2, color = NA) +
                    xlim(1972, 2016) + ylim(0, .45) +
                    scale_x_continuous(breaks = c(seq(1976, 2016, by=8), 2016) ) +
                    scale_linetype_manual(values = c("dashed", "solid", "dotted")) +
                    scale_color_manual(values = vec_colors) +
                    scale_fill_manual(values = vec_colors) +
                    xlab("Year") +ylab("Correlation Coefficient") + theme_bw() +
                    annotate("rect",xmin=2004, xmax=Inf, ymin=-Inf,
                             ymax=Inf, alpha=0.2, fill="burlywood1") +
                    theme(plot.title = element_text(hjust=0.5),
                          legend.title = element_blank())

ggsave(plot=g_combo_spline, filename="../graphs/svg/splines_all_wghtF2F.svg")
ggsave(plot=g_combo_spline, filename="../graphs/png/splines_all_wghtF2F.png")