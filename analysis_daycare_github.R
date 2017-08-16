#######################
# Analysis day care density
# Citation: Goldstein ND, Newbern EC, Tabb LP, Welles SL. Density of Day Cares in Relation to Reported Pertussis Incidence in Philadelphia. Public Health. 2017 May;146:126-133.
# 8/20/2014 -- Neal Goldstein
#######################


### FUNCTIONS ###

library(psych) #describe, describeBy
library(gmodels) #CrossTable
library(lme4) #nonlinear mixed effects (glmer)
library(boot) #bootstrapping
library(maptools) #mapping
library(gpclib) #aggregate census tracts
library(spdep) #morans i

#returns the MOR from a glmer model, see: http://www.ncbi.nlm.nih.gov/pubmed/16537344
bootMOR = function(model)
{
  return(exp(0.95*getME(model,"theta")))
}


### READ DATA ###

load("study_population_analysis.RData")
neighborhood_utd = read.csv(file="Neighborhood UTD.csv", header=T, as.is=T, na.strings="")


### SUBSET DATA ###

#drop labels from neighborhood utd and intialize var
neighborhood_utd = neighborhood_utd[-1,]
studypop_final$neighborhood_utd = NA

analysis = NA
matchid = studypop_final$matchid[1]
case_asis = studypop_final$case_asis[1]
case_age = studypop_final$age_weeks[1]

#only include cases that are classified as confirmed or probable per the AS IS definition
for (i in 1:nrow(studypop_final))
{
  cat("\n\n**************** Observation: ",i," ****************\n",sep="")
  
  #neighborhood utd
  studypop_final$neighborhood_utd = ifelse(length(as.numeric(neighborhood_utd$X.3[which(studypop_final$neighborhood[i]==neighborhood_utd$X)]))>0, as.numeric(neighborhood_utd$X.3[which(studypop_final$neighborhood[i]==neighborhood_utd$X)]), NA)
  
  if (case_asis != 0)
    analysis = rbind(analysis, studypop_final[i, ])
  
  if ((matchid != studypop_final$matchid[i+1]) && (!is.na(studypop_final$matchid[i+1])))
  {
    matchid = studypop_final$matchid[i+1]
    case_asis = studypop_final$case_asis[i+1]
    case_age = studypop_final$age_weeks[i+1]
  }
}
rm(i,case_asis,matchid,case_age)
analysis = analysis[-1,]
rownames(analysis) = NULL


### RECODE ###

#dichotomize parity
analysis$mother_parous = ifelse(analysis$mother_parity>=1, 1, 0)

#dichotomize age
analysis$age_1yr = ifelse(analysis$age==0,0,1)

#dichotomize day care density by median
analysis$neighborhood_day_cares_density_high = ifelse(analysis$neighborhood_day_cares_density<median(analysis$neighborhood_day_cares_density,na.rm=T), 0, 1)

#dichotomize neighborhood income by median
analysis$neighborhood_median_income_high = ifelse(analysis$neighborhood_median_income<median(analysis$neighborhood_median_income,na.rm=T), 0, 1)

#alternate specification of day care density
analysis$neighborhood_day_cares_percapita = analysis$neighborhood_day_cares / ((analysis$neighborhood_under5/100)*(analysis$neighborhood_pop_density*analysis$neighborhood_area))
analysis$neighborhood_day_cares_percapita_density = analysis$neighborhood_day_cares_percapita / analysis$neighborhood_area


### IMPUTE mother insurance ###

for (i in 1:nrow(analysis))
{
  cat("\n\n**************** Observation: ",i," ****************\n",sep="")

  #check for missing insurance
  if (is.na(analysis$mother_insurance[i]))
  {
    #try to impute insurance based on neighborhood, where >=50% is cut point
    analysis$mother_insurance[i] = ifelse(mean(analysis$mother_insurance[analysis$neighborhood==analysis$neighborhood[i]],na.rm=T)>=0.5, 1, 0)
  }
}
rm(i)


### GRAND MEAN CENTER ###

#see: http://www.ncbi.nlm.nih.gov/pubmed/17563168

analysis$mother_age_centered = scale(analysis$mother_age, center=T, scale=F)
analysis$neighborhood_median_age_centered = scale(analysis$neighborhood_median_age, center=T, scale=F)
analysis$neighborhood_race_white_centered = scale(analysis$neighborhood_race_white, center=T, scale=F)
analysis$neighborhood_race_black_centered = scale(analysis$neighborhood_race_black, center=T, scale=F)
analysis$neighborhood_race_other_centered = scale(analysis$neighborhood_race_other, center=T, scale=F)
analysis$neighborhood_hispanic_centered = scale(analysis$neighborhood_hispanic, center=T, scale=F)
analysis$neighborhood_foreign_born_centered = scale(analysis$neighborhood_foreign_born, center=T, scale=F)
analysis$neighborhood_single_parent_centered = scale(analysis$neighborhood_single_parent, center=T, scale=F)
analysis$neighborhood_under5_centered = scale(analysis$neighborhood_under5, center=T, scale=F)
analysis$neighborhood_highschool_centered = scale(analysis$neighborhood_highschool, center=T, scale=F)
analysis$neighborhood_bachelors_centered = scale(analysis$neighborhood_bachelors, center=T, scale=F)
analysis$neighborhood_unemployed_centered = scale(analysis$neighborhood_unemployed, center=T, scale=F)
analysis$neighborhood_median_income_log_centered = scale(log(analysis$neighborhood_median_income), center=T, scale=F)
analysis$neighborhood_pop_density_centered = scale(log(analysis$neighborhood_pop_density), center=T, scale=F)
analysis$neighborhood_utd_centered = scale(analysis$neighborhood_utd, center=T, scale=F)
analysis$neighborhood_day_cares_density_centered = scale(analysis$neighborhood_day_cares_density, center=T, scale=F)
analysis$neighborhood_day_cares_percapita_centered = scale(analysis$neighborhood_day_cares_percapita, center=T, scale=F)
analysis$neighborhood_day_cares_percapita_density_centered = scale(analysis$neighborhood_day_cares_percapita_density, center=T, scale=F)
analysis$neighborhood_day_cares_centered = scale(analysis$neighborhood_day_cares, center=T, scale=F)


### ANALYSIS: Overall study population ###

#individual
CrossTable(analysis$age_1yr)
CrossTable(analysis$gender)
CrossTable(analysis$race)
CrossTable(analysis$ethnicity)
CrossTable(analysis$pertussis_vax_utd)
describe(analysis$mother_age)
CrossTable(analysis$mother_married)
CrossTable(analysis$mother_foreign_born)
CrossTable(analysis$mother_education)
CrossTable(analysis$mother_insurance)
CrossTable(analysis$mother_parous)

#neighborhood
describe(analysis$neighborhood_median_age)
describe(analysis$neighborhood_race_white)
describe(analysis$neighborhood_race_black)
describe(analysis$neighborhood_race_other)
describe(analysis$neighborhood_hispanic)
describe(analysis$neighborhood_foreign_born)
describe(analysis$neighborhood_single_parent)
describe(analysis$neighborhood_under5)
describe(analysis$neighborhood_highschool)
describe(analysis$neighborhood_bachelors)
describe(analysis$neighborhood_unemployed)
describe(analysis$neighborhood_median_income)
describe(analysis$neighborhood_pop_density)
describe(analysis$neighborhood_utd)
describe(analysis$neighborhood_day_cares_density)


### ANALYSIS: Comparison of cases and controls, individual variables ###

CrossTable(analysis$age_1yr, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$gender, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$race, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$ethnicity, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$pertussis_vax_utd, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(analysis$mother_age, analysis$case); t.test(analysis$mother_age ~ analysis$case)
CrossTable(analysis$mother_married, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$mother_foreign_born, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$mother_education, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$mother_insurance, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$mother_parous, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(analysis$pertussis_vax_utd, analysis$case, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)


### ANALYSIS: Comparison of cases and controls, neighborhood variables ###

describeBy(analysis$neighborhood_median_age, analysis$case); summary(glmer(case ~ (1 | neighborhood) + neighborhood_median_age_centered, family=binomial(), data=analysis))
describeBy(analysis$neighborhood_race_white, analysis$case); summary(glmer(case ~ (1 | neighborhood) + neighborhood_race_white_centered, family=binomial(), data=analysis))
describeBy(analysis$neighborhood_race_black, analysis$case); summary(glmer(case ~ (1 | neighborhood) + neighborhood_race_black_centered, family=binomial(), data=analysis))
describeBy(analysis$neighborhood_race_other, analysis$case); summary(glmer(case ~ (1 | neighborhood) + neighborhood_race_other_centered, family=binomial(), data=analysis))
describeBy(analysis$neighborhood_hispanic, analysis$case); summary(glmer(case ~ (1 | neighborhood) + neighborhood_hispanic_centered, family=binomial(), data=analysis))
describeBy(analysis$neighborhood_foreign_born, analysis$case); summary(glmer(case ~ (1 | neighborhood) + neighborhood_foreign_born_centered, family=binomial(), data=analysis))
describeBy(analysis$neighborhood_single_parent, analysis$case); summary(glmer(case ~ (1 | neighborhood) + neighborhood_single_parent_centered, family=binomial(), data=analysis))
describeBy(analysis$neighborhood_under5, analysis$case); summary(glmer(case ~ (1 | neighborhood) + neighborhood_under5_centered, family=binomial(), data=analysis))
describeBy(analysis$neighborhood_highschool, analysis$case); summary(glmer(case ~ (1 | neighborhood) + neighborhood_highschool_centered, family=binomial(), data=analysis))
describeBy(analysis$neighborhood_bachelors, analysis$case); summary(glmer(case ~ (1 | neighborhood) + neighborhood_bachelors_centered, family=binomial(), data=analysis))
describeBy(analysis$neighborhood_unemployed, analysis$case); summary(glmer(case ~ (1 | neighborhood) + neighborhood_unemployed_centered, family=binomial(), data=analysis))
describeBy(analysis$neighborhood_median_income, analysis$case); summary(glmer(case ~ (1 | neighborhood) + neighborhood_median_income_log_centered, family=binomial(), data=analysis))
describeBy(analysis$neighborhood_pop_density, analysis$case); summary(glmer(case ~ (1 | neighborhood) + neighborhood_pop_density_centered, family=binomial(), data=analysis))
describeBy(analysis$neighborhood_utd, analysis$case); summary(glmer(case ~ (1 | neighborhood) + neighborhood_utd_centered, family=binomial(), data=analysis))
describeBy(analysis$neighborhood_day_cares_density, analysis$case); summary(glmer(case ~ (1 | neighborhood) + neighborhood_day_cares_density_centered, family=binomial(), data=analysis))


### ANALYSIS: DISTRIBUTION of CASES by NEIGHBORHOOD ###

table(analysis$neighborhood[analysis$case==1])

#read census tracts, 2009 shapefiles (2000 boundaries)
census_tracts = readShapePoly("2009 Census Tract Shapefiles/tl_2009_42101_tract00")

#get ordered list of census tracts
ct_list = round(as.numeric(as.character(slot(census_tracts, "data")$NAME00)),0)

#aggregate into neighborhoods
gpclibPermit() 
neighborhoods = rep(NA,length(ct_list))
neighborhoods[which(is.element(ct_list, c(344,345,355,356)))] = "1"
neighborhoods[which(is.element(ct_list, c(1,2,3,4,5,6,7,8,9,10,11,12,366)))] = "2"
neighborhoods[which(is.element(ct_list, c(224,225,226,227,228,229,230,231,232,233,234,235,236,237,239,257)))] = "3"
neighborhoods[which(is.element(ct_list, c(80,81,82,83,84,85)))] = "4"
neighborhoods[which(is.element(ct_list, c(206,207,208,240,243,244)))] = "5"
neighborhoods[which(is.element(ct_list, c(250,251,252,253,254,255,256)))] = "6"
neighborhoods[which(is.element(ct_list, c(54,55,56,57,58,59,60,61,62,67,68)))] = "7"
neighborhoods[which(is.element(ct_list, c(125,133,134,135,136)))] = "8"
neighborhoods[which(is.element(ct_list, c(293,294,295,296,299,300,301,302)))] = "9"
neighborhoods[which(is.element(ct_list, c(238,241,242,245,246,247,248,249)))] = "10"
neighborhoods[which(is.element(ct_list, c(33,34,35,36,46,47,48,50,51)))] = "11"
neighborhoods[which(is.element(ct_list, c(93,94,95,96,100,101,102,112,113,114,115)))] = "12"
neighborhoods[which(is.element(ct_list, c(174,175,176,194,195,196,197,198,199)))] = "13"
neighborhoods[which(is.element(ct_list, c(187,188,189,190,191,192,193)))] = "14"
neighborhoods[which(is.element(ct_list, c(291,292,303,304,305,306,307,367)))] = "15"
neighborhoods[which(is.element(ct_list, c(280,281,282,283,284)))] = "16"
neighborhoods[which(is.element(ct_list, c(143,158,159,160,161)))] = "17"
neighborhoods[which(is.element(ct_list, c(315,316,317,318,329,330,331,332,333)))] = "18"
neighborhoods[which(is.element(ct_list, c(92,103,104,105,106,107,108,109,110,111,124)))] = "19"
neighborhoods[which(is.element(ct_list, c(170,171,172,173,200,201,202,203,204,205)))] = "20"
neighborhoods[which(is.element(ct_list, c(126,127,128,129,130,142,144,156,157,162,163,164)))] = "21"
neighborhoods[which(is.element(ct_list, c(268,269,270,271,272,276)))] = "22"
neighborhoods[which(is.element(ct_list, c(265,267,277,278,279)))] = "23"
neighborhoods[which(is.element(ct_list, c(273,274,275,285,286,287,288,289,290)))] = "24"
neighborhoods[which(is.element(ct_list, c(97,98,99,116,117,118,119,120,121,122,123)))] = "25"
neighborhoods[which(is.element(ct_list, c(308,309,310,311,312,313,314)))] = "26"
neighborhoods[which(is.element(ct_list, c(63,64,65,66,69,70,71,72,73,74)))] = "27"
neighborhoods[which(is.element(ct_list, c(15,16,17,18,25,26,27)))] = "28"
neighborhoods[which(is.element(ct_list, c(131,132,141,145,146,154,155,165,166)))] = "29"
neighborhoods[which(is.element(ct_list, c(334,335,336,337,338,339,340,341,342,343)))] = "30"
neighborhoods[which(is.element(ct_list, c(180,181,182,183,184,185,186)))] = "31"
neighborhoods[which(is.element(ct_list, c(209,210,211,212,213,214,215,216,217,218,219,220,221,222)))] = "32"
neighborhoods[which(is.element(ct_list, c(13,14,19,20,21,22,31,32)))] = "33"
neighborhoods[which(is.element(ct_list, c(138,139,140,147,148,153,167)))] = "34"
neighborhoods[which(is.element(ct_list, c(41,42,44)))] = "35"
neighborhoods[which(is.element(ct_list, c(357,358,359,360,365)))] = "36"
neighborhoods[which(is.element(ct_list, c(23,24,28,29,30)))] = "37"
neighborhoods[which(is.element(ct_list, c(37,38,39,40,45)))] = "38"
neighborhoods[which(is.element(ct_list, c(137,149,151,152,168,169)))] = "39"
neighborhoods[which(is.element(ct_list, c(353,354,361,362,363,364)))] = "40"
neighborhoods[which(is.element(ct_list, c(328,346,347,348,349,351,352)))] = "41"
neighborhoods[which(is.element(ct_list, c(76,77,78,79,86,87,88,89,90,91)))] = "42"
neighborhoods[which(is.element(ct_list, c(177,178,179)))] = "43"
neighborhoods[which(is.element(ct_list, c(258,259,260,261,262,263,264,266)))] = "44"
neighborhoods[which(is.element(ct_list, c(297,298,319,320,321,322,323,324,325,326,327)))] = "45"
neighborhoods_sp = unionSpatialPolygons(census_tracts, neighborhoods)

#create color palette from neighborhoods, order is from t(sapply(slot(neighborhoods_sp, "polygons"), function(i) slot(i, "ID")))
col_density = data.frame(neighborhood=c(1,10,11,12,13,14,15,16,17,18,19,2,20,21,22,23,24,25,26,27,28,29,3,30,31,32,33,34,35,36,37,38,39,4,40,41,42,43,44,45,5,6,7,8,9), name=NA, density=0, cases=0, daycares=0, utd=0, color="#FFFFFF", stringsAsFactors=F) 

#tally cases, day cares, and vax utd in each neighborhood
for (i in 1:nrow(col_density))
{
  #cases of pertussis by neighborhood
  col_density$cases[i] = sum(analysis$case[analysis$neighborhood==col_density$neighborhood[i]],na.rm=T)

  #day care density by neighborhood
  col_density$daycares[i] = mean(analysis$neighborhood_day_cares_density[analysis$neighborhood==col_density$neighborhood[i]],na.rm=T)
  
  #utd % by neighborhood
  col_density$utd[i] = neighborhood_utd$X.3[neighborhood_utd$X==col_density$neighborhood[i]]
}
rm(i)

#set mapping variable of interest
col_density$density = as.numeric(col_density$cases)

#color palette
plot_palette = c("#DDDDDD","#BBBBBB","#999999")
#plot_palette = heat.colors(3)[3:1]

#grayscale by tertiles of density
col_density$color = ifelse(col_density$density<quantile(col_density$density, probs=c(0,.33,.67,1))[2], plot_palette[1], col_density$color)
col_density$color = ifelse(col_density$density>quantile(col_density$density, probs=c(0,.33,.67,1))[3], plot_palette[3], col_density$color)
col_density$color = ifelse(col_density$density>=quantile(col_density$density, probs=c(0,.33,.67,1))[2] & col_density$density<=quantile(col_density$density, probs=c(0,.33,.67,1))[3], plot_palette[2], col_density$color)

#output to a file (then open with gimp, scale to 1200, export as pdf)
tiff("Philly.tif",height=4,width=5,units='in',res=1200)

#plot cases on map
par(oma=c(0.5, 0.5, 0.5, 0.5))
par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
par(plt=c(0, 1, 0, 1))
plot(neighborhoods_sp, border="#333333", col=col_density$color)
legend("bottomright", fill=plot_palette, legend=c("<6 cases", "6-10 cases",">10 cases"), xpd=T, cex=0.75)
#legend("bottomright", fill=plot_palette, legend=c("Low", "Medium","High"), xpd=T)
#title("\nDistribution of pertussis cases in Philadelphia")
#title("\nDistribution of day care density in Philadelphia")
#title("UTD % by neightborhood in Philadelphia")

#add neighborhood labels
#text(t(sapply(slot(neighborhoods_sp, "polygons"), function(i) slot(i, "labpt"))), cex=0.35, labels=c("Bustleton","Germantown","Grays Ferry-Passyunk","Haddington-Overbrook","Hunting Park-Fairhill","Juniata Park-Harrowgate","Lawndale-Crescentville","Logan","Lower Kensington","Mayfair-Holmesburg","Mill Creek-Parkside","Center City","Nicetown-Tioga","Northern Liberties-West Kensington","Oak Lane-Fernrock","Ogontz","Olney-Feltonville","Overbrook Park-Wynnefield Heights","Paschall-Kingsessing","Pennsport-Queen Village","Poplar-Temple","Chestnut Hill-W. Mt. Airy","Rhawnhurst-Fox Chase","Richmond-Bridesburg","Roxborough-Manayunk","Schuylkill-Point Breeze","Sharswood-Stanton","Oxford Circle","Snyder-Whitman","Somerton","Southwark-Bella Vista","South Broad-Girard Estates","Strawberry Mansion","Cobbs Creek","Torresdale-North","Torresdale-S.-Pennypack Park","University City","Upper Kensington","West Oak Lane-Cedarbrook","Wissinoming-Tacony","East Falls-Westside","East Mt. Airy","Eastwick-Elmwood","Fairmount-Spring Garden","Frankford"))
text(t(sapply(slot(neighborhoods_sp, "polygons"), function(i) slot(i, "labpt"))), cex=0.5, labels=1:45)

#close file
dev.off()
 

# ### ANALYSIS: Comparison of high and low density day care neighborhoods, individual variables ###
# 
# CrossTable(analysis$age_1yr, analysis$neighborhood_day_cares_density_high, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(analysis$gender, analysis$neighborhood_day_cares_density_high, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(analysis$race, analysis$neighborhood_day_cares_density_high, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(analysis$ethnicity, analysis$neighborhood_day_cares_density_high, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(analysis$pertussis_vax_utd, analysis$neighborhood_day_cares_density_high, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# describeBy(analysis$mother_age, analysis$neighborhood_day_cares_density_high); t.test(analysis$mother_age ~ analysis$neighborhood_day_cares_density_high)
# CrossTable(analysis$mother_married, analysis$neighborhood_day_cares_density_high, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(analysis$mother_foreign_born, analysis$neighborhood_day_cares_density_high, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(analysis$mother_education, analysis$neighborhood_day_cares_density_high, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(analysis$mother_insurance, analysis$neighborhood_day_cares_density_high, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(analysis$mother_parous, analysis$neighborhood_day_cares_density_high, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(analysis$pertussis_vax_utd, analysis$neighborhood_day_cares_density_high, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)


# ### ANALYSIS: Comparison high and low density day care neighborhoods, neighborhood variables ###
# 
# describeBy(analysis$neighborhood_median_age, analysis$neighborhood_day_cares_density_high); t.test(analysis$neighborhood_median_age ~ analysis$neighborhood_day_cares_density_high)
# describeBy(analysis$neighborhood_race_white, analysis$neighborhood_day_cares_density_high); wilcox.test(analysis$neighborhood_race_white ~ analysis$neighborhood_day_cares_density_high)
# describeBy(analysis$neighborhood_race_black, analysis$neighborhood_day_cares_density_high); wilcox.test(analysis$neighborhood_race_black ~ analysis$neighborhood_day_cares_density_high)
# describeBy(analysis$neighborhood_race_other, analysis$neighborhood_day_cares_density_high); wilcox.test(analysis$neighborhood_race_other ~ analysis$neighborhood_day_cares_density_high)
# describeBy(analysis$neighborhood_hispanic, analysis$neighborhood_day_cares_density_high); wilcox.test(analysis$neighborhood_hispanic ~ analysis$neighborhood_day_cares_density_high)
# describeBy(analysis$neighborhood_foreign_born, analysis$neighborhood_day_cares_density_high); wilcox.test(analysis$neighborhood_foreign_born ~ analysis$neighborhood_day_cares_density_high)
# describeBy(analysis$neighborhood_single_parent, analysis$neighborhood_day_cares_density_high); wilcox.test(analysis$neighborhood_single_parent ~ analysis$neighborhood_day_cares_density_high)
# describeBy(analysis$neighborhood_under5, analysis$neighborhood_day_cares_density_high); t.test(analysis$neighborhood_under5 ~ analysis$neighborhood_day_cares_density_high)
# describeBy(analysis$neighborhood_highschool, analysis$neighborhood_day_cares_density_high); wilcox.test(analysis$neighborhood_highschool ~ analysis$neighborhood_day_cares_density_high)
# describeBy(analysis$neighborhood_bachelors, analysis$neighborhood_day_cares_density_high); wilcox.test(analysis$neighborhood_bachelors ~ analysis$neighborhood_day_cares_density_high)
# describeBy(analysis$neighborhood_unemployed, analysis$neighborhood_day_cares_density_high); t.test(analysis$neighborhood_unemployed ~ analysis$neighborhood_day_cares_density_high)
# describeBy(analysis$neighborhood_median_income, analysis$neighborhood_day_cares_density_high); t.test(analysis$neighborhood_median_income ~ analysis$neighborhood_day_cares_density_high)
# describeBy(analysis$neighborhood_pop_density, analysis$neighborhood_day_cares_density_high); t.test(analysis$neighborhood_pop_density ~ analysis$neighborhood_day_cares_density_high)
# describeBy(analysis$neighborhood_utd, analysis$neighborhood_day_cares_density_high); wilcox.test(analysis$neighborhood_utd ~ analysis$neighborhood_day_cares_density_high)


### ANALYSIS: Potential confounder selection ###

#age is always included as it was a frequency matched variable

#starting variable list: theory driven and association p<0.15 with case status
#individual: age, race, ethnicity, UTD, maternal age, nativity, education, insurance, parity
#neighborhood: children <5, high school graduate, UTD

#full model with potential confounders
summary(glmer(case ~ (1 + neighborhood_day_cares_density_centered | neighborhood) + neighborhood_day_cares_density_centered+as.factor(age_1yr)+as.factor(race)+as.factor(ethnicity)+mother_age_centered+as.factor(mother_foreign_born)+as.factor(mother_education)+as.factor(mother_insurance)+as.factor(mother_parous)+as.factor(pertussis_vax_utd)+neighborhood_under5_centered+neighborhood_highschool_centered+neighborhood_utd_centered, family=binomial(), data=analysis))

#select individual: backward remove nonsignificant vars and check for change in estimate >10%
summary(glmer(case ~ (1 + neighborhood_day_cares_density_centered | neighborhood) + neighborhood_day_cares_density_centered+as.factor(age_1yr)+as.factor(race)+as.factor(ethnicity)+mother_age_centered+as.factor(mother_foreign_born)+as.factor(mother_education)+as.factor(mother_insurance)+as.factor(mother_parous)+as.factor(pertussis_vax_utd), family=binomial(), data=analysis))
summary(glmer(case ~ (1 + neighborhood_day_cares_density_centered | neighborhood) + neighborhood_day_cares_density_centered+as.factor(age_1yr)+as.factor(race)+as.factor(ethnicity)+as.factor(mother_foreign_born)+as.factor(mother_education)+as.factor(mother_insurance)+as.factor(mother_parous)+as.factor(pertussis_vax_utd), family=binomial(), data=analysis))
#summary(glmer(case ~ (1 + neighborhood_day_cares_density_centered | neighborhood) + neighborhood_day_cares_density_centered+as.factor(age_1yr)+as.factor(race)+as.factor(ethnicity)+as.factor(mother_foreign_born)+as.factor(mother_insurance)+as.factor(mother_parous)+as.factor(pertussis_vax_utd), family=binomial(), data=analysis))
summary(glmer(case ~ (1 + neighborhood_day_cares_density_centered | neighborhood) + neighborhood_day_cares_density_centered+as.factor(age_1yr)+as.factor(race)+as.factor(ethnicity)+as.factor(mother_foreign_born)+as.factor(mother_education)+as.factor(mother_parous)+as.factor(pertussis_vax_utd), family=binomial(), data=analysis))

#select neighborhood: backward remove nonsignificant vars and check for change in estimate >10%
summary(glmer(case ~ (1 + neighborhood_day_cares_density_centered | neighborhood) + neighborhood_day_cares_density_centered+neighborhood_under5_centered+neighborhood_highschool_centered+neighborhood_utd_centered, family=binomial(), data=analysis))
summary(glmer(case ~ (1 + neighborhood_day_cares_density_centered | neighborhood) + neighborhood_day_cares_density_centered+neighborhood_highschool_centered+neighborhood_utd_centered, family=binomial(), data=analysis))

#full model
summary(glmer(case ~ (1 + neighborhood_day_cares_density_centered | neighborhood) + as.factor(age_1yr)+as.factor(race)+as.factor(ethnicity)+as.factor(mother_foreign_born)+as.factor(mother_education)+as.factor(mother_parous)+as.factor(pertussis_vax_utd)+neighborhood_day_cares_density_centered+neighborhood_highschool_centered+neighborhood_utd_centered, family=binomial(), data=analysis))

#final variable list
#individual: age, race, ethnicity, nativity, education, parity, UTD
#neighborhood: high school graduate, UTD


### ANALYSIS: Empty model ###

model = glmer(case ~ (1 | neighborhood), family=binomial(), data=analysis)
summary(model)

#area level variance
getME(model,"theta")

#compute median OR, see: http://www.ncbi.nlm.nih.gov/pubmed/16537344
exp(0.95*getME(model,"theta"))

#bootstrap precision around MOR
boot_model = bootMer(model, bootMOR, nsim=1000, parallel="multicore", ncpus=4)
boot.ci(boot_model, type="norm", index=1)

#check autocorrelation via morans i; see https://cran.r-project.org/web/packages/CARBayes/vignettes/CARBayesvignette.pdf
neighbors_list = nb2listw(poly2nb(neighborhoods_sp), style="B", zero.policy=T)
moran.mc(x=ranef(model)[[1]][[1]], listw=nb2listw(poly2nb(neighborhoods_sp), style="B", zero.policy=T), nsim=1000, zero.policy=T)


### ANALYSIS: Day care model, random slope ###

model = glmer(case ~ (1 + neighborhood_day_cares_density_centered | neighborhood) + neighborhood_day_cares_density_centered, family=binomial(), data=analysis)
summary(model)

#OR and CI estimates for fixed effects
exp(fixef(model))
exp(confint.merMod(model, method="Wald"))

#area level variance
getME(model,"theta")

#compute median OR, see: http://www.ncbi.nlm.nih.gov/pubmed/16537344
exp(0.95*getME(model,"theta"))

#bootstrap precision around MOR
boot_model = bootMer(model, bootMOR, nsim=1000, parallel="multicore", ncpus=4)
boot.ci(boot_model, type="norm", index=1)


### ANALYSIS: Model with individual level ###

model = glmer(case ~ (1 | neighborhood) + as.factor(age_1yr)+as.factor(race)+as.factor(ethnicity)+as.factor(mother_foreign_born)+as.factor(mother_education)+as.factor(mother_parous)+as.factor(pertussis_vax_utd), family=binomial(), data=analysis)
summary(model)

#OR and CI estimates for fixed effects
exp(fixef(model))
exp(confint.merMod(model, method="Wald"))

#area level variance
getME(model,"theta")

#compute median OR, see: http://www.ncbi.nlm.nih.gov/pubmed/16537344
exp(0.95*getME(model,"theta"))

#bootstrap precision around MOR
boot_model = bootMer(model, bootMOR, nsim=1000, parallel="multicore", ncpus=4)
boot.ci(boot_model, type="norm", index=1)


### ANALYSIS: Model with individual+neighborhood level ###

model = glmer(case ~ (1 + neighborhood_day_cares_density_centered | neighborhood) + as.factor(age_1yr)+as.factor(race)+as.factor(ethnicity)+as.factor(mother_foreign_born)+as.factor(mother_education)+as.factor(mother_parous)+as.factor(pertussis_vax_utd)+neighborhood_highschool_centered+neighborhood_utd_centered+neighborhood_day_cares_density_centered, family=binomial(), data=analysis)
summary(model)

#OR and CI estimates for fixed effects
exp(fixef(model))
exp(confint.merMod(model, method="Wald"))

#area level variance
getME(model,"theta")

#compute median OR, see: http://www.ncbi.nlm.nih.gov/pubmed/16537344
exp(0.95*getME(model,"theta"))

#bootstrap precision around MOR
boot_model = bootMer(model, bootMOR, nsim=1000, parallel="multicore", ncpus=4)
boot.ci(boot_model, type="norm", index=1)
