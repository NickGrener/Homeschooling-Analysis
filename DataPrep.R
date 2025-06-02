library(maps)
library(tidycensus)
library(units)
library(tidyr)
library(sf)
library(stringr)
library(dplyr)
library(tigris)
library(mapview)
library(tmap)
library(leafsync)
library(ggspatial)
library(remotes)
library(readxl)
library(ggplot2)
library(readr)
library(car)

### Start with ACS sf file of all 613 school districts
OH.sf.base <- get_acs(
    geography = "school district (unified)",
    state = "OH",
    variables = c("S0101_C01_003", "S0101_C01_004", "S0101_C01_005"),
    survey = "acs5",
    year = 2020,
    geometry = TRUE
)


### Change to wider format then add values from the numeric columns (after
# dropping geometry attributes temporarily) to get total school-age population
OH.acs.df <- OH.sf.base %>% 
    pivot_wider(names_from = variable, values_from = c(estimate, moe))
numeric_cols <- `st_geometry<-`(OH.acs.df, NULL)
OH.acs.df$School_age_pop <- apply(numeric_cols[, 3:5], 1, sum)
OH.acs.df <- dplyr::select(OH.acs.df, -c("estimate_S0101_C01_003","estimate_S0101_C01_004","estimate_S0101_C01_005", "moe_S0101_C01_003","moe_S0101_C01_004","moe_S0101_C01_005"))
# Remove two districts-(North Bass and undefined) with a school age pop of zero
summary(OH.acs.df$School_age_pop)
OH.acs.df <- OH.acs.df[-c(9,176), ]


### Get raw versions of WaPost data for 2017-18 and 2022-23 from Ohio, then full join on IRN
# Source: https://reports.education.ohio.gov/report/nonpublic-data-homeschool-student-counts
homeschool.df <- read.csv("./Data-sources/home_school_counts_22.csv")
homeschool.df <- homeschool.df %>% dplyr::select(-1) %>% rename(Homeschoolers22 = "Homeschool.Students")
homeschool.old.df <- read.csv("./Data-sources/home_school_counts_17.csv")
homeschool.old.df <- homeschool.old.df %>% dplyr::select(-c(1, 2, 4)) %>% rename(Homeschoolers17 = "Homeschool.Students")
response.df <- full_join(homeschool.df, homeschool.old.df, by = "IRN")
response.df <- dplyr::select(response.df, -County)
sum(!complete.cases(response.df))
# 66 of the 636 rows have an NA in one or more column; 
# 32 of these have NA's in every column except for "<10" in'17 count, so delete these rows:
response.df <- response.df %>% filter(!(is.na(District) & Homeschoolers17 == "<10"))
sum(!complete.cases(response.df))


### Most recent District Profile Report (Cupp Report) found at
# https://data.ohio.gov/wps/portal/gov/data/view/ode-cupp-report
# Select just variables of interest and rename to be more workable
DPR.sheet.names <- excel_sheets("./Data-sources/DPR_FY22.xlsx")
predictors.df <- read_excel("./Data-sources/DPR_FY22.xlsx", sheet = "District Data")
predictors.df <- predictors.df %>% dplyr::select(1, 2, 4:13, 16, 19, 33, 34, 42)
predictors.df <- predictors.df %>% rename(Cupp_District = "District", Pupil_Density = "District Pupil Density FY22", Enrollment22 = "Enrolled ADM FY22", Asian = "Asian Students as % of Total FY22", 
                            Pacific_Islander = "Pacific Islanders Students as % of Total FY22", Black = "Black Students as % of Total FY22", American_Indian = "American Indian/Alaskan Native Students as % of Total FY22",
                            Hispanic = "Hispanic Students as % of Total FY22", White = "White Students as % of Total FY22", Multiracial = "Multiracial Students as % of Total FY22",
                            Economically_Disadvantaged = "% of Economically Disadvantaged Students FY22", Avg_Teacher_Salary = "Classroom Teachers' Average Salary FY22", Experienced_Teachers = "% Teachers with 10+ Years Experience FY22",
                            Median_Income = "Ohio Median Income TY20", Mean_Income = "Federal Average Income TY20", Tax_Effort_Index = "Local Tax Effort Index FY22")
str(predictors.df)


### Counts at nonpublic (private) schools; location is given at county level embedded in larger chr string
private.schools.df <- read_excel("./Data-sources/private_school_counts.xlsx")
private.schools.df <- private.schools.df %>% dplyr::select(-1) %>% 
    rename(School_plus_code_plus_county = "Nonpublic School")
# Extract county name (two rows are NA for county so drop these before proceeding)
count_parentheses <- str_count(private.schools.df$School_plus_code_plus_county, "\\)")
summary(count_parentheses)
code_present <- grepl("\\)", private.schools.df$School_plus_code_plus_county)
print(private.schools.df[!code_present,])
private.schools.df <- private.schools.df[code_present, ]
private.schools.df <- separate(private.schools.df, School_plus_code_plus_county, 
                               into = c("School_plus_code", "County"), sep = "\\)")

sum(private.schools.df$Enrollment == "<10")
# Zero out all censored data (negligible number of schools) and get county totals
private.schools.df$Students <- ifelse(private.schools.df$Enrollment == "<10", 0, as.integer(private.schools.df$Enrollment))
private.county.totals.df <- private.schools.df %>% 
    group_by(County) %>%
    summarize(Pri_county_sum = sum(Students)) %>%
    ungroup()
str(private.county.totals.df)
# Fix county names to fit format of receiving df; manually adjust "Van Wert", which is only county with a space in its name
private.county.totals.df$County <- sub("^\\s*-\\s*(\\w+)\\s+\\w+$", "\\1", private.county.totals.df$County)
private.county.totals.df$County[71] <- "Van Wert"


### Republican vote share in 2020 presidential election predictor (at county level)
politics.df <- read_csv("./Data-sources/countypres_2000-2020.csv")
politics.df <- politics.df %>% filter(state_po == "OH", year == 2020, party == "REPUBLICAN")
politics.df <- dplyr::select(politics.df, c(county_name, candidatevotes, totalvotes))
politics.df$County <- paste(substring(politics.df$county_name, 1, 1), tolower(substring(politics.df$county_name, 2)), sep = "")
politics.df$Rep_Vote_2020 <- politics.df$candidatevotes / politics.df$totalvotes
politics.df <- politics.df %>% dplyr::select(County, Rep_Vote_2020)
politics.df$County <- ifelse(politics.df$County == "Van wert", "Van Wert", politics.df$County)
str(politics.df)


### Merge all county-level data: Rep_Vote_2020 & Pri_county_sum; for those 11 counties
# with no private school students, replace NA with zero
county.data.df <- full_join(politics.df, private.county.totals.df, by = "County")
county.data.df$Pri_county_sum[is.na(county.data.df$Pri_county_sum)] <- 0


### 4-year graduation rate predictor 
# Note ID code that will be used to join is embedded in larger chr string
grad.rate.df <- read_csv("./Data-sources/grad_rate.csv")
grad.rate.df <- grad.rate.df %>% dplyr::select(1, 5)
grad.rate.df <- grad.rate.df %>% mutate(Name = sapply(strsplit(grad.rate.df$District, " - "), `[`, 1)) %>% 
    mutate(ID = substring(grad.rate.df$District, regexpr(" - \\d{6} ", grad.rate.df$District) + 3, 
                          regexpr(" - \\d{6} ", grad.rate.df$District) + 7)) %>%
    mutate(county = substring(grad.rate.df$District, regexpr("\\(.*\\)", grad.rate.df$District) + 1, regexpr("\\)", grad.rate.df$District) - 1)) %>% 
    dplyr::select(-1) %>% 
    rename(Grad_Rate_4Yr = "4-Year Graduation Rate")
#Convert character string percentages to numeric proportions
grad.rate.df$Grad_Rate_4Yr <- as.numeric(gsub("%", "", grad.rate.df$Grad_Rate_4Yr)) / 100
str(grad.rate.df)


### Chronic absenteeism rate and avg. daily attendance rate predictors
# Clean using same procedure as graduation rate
absentee.rate.df <- read_csv("./Data-sources/absenteeism_rate.csv")
absentee.rate.df <- absentee.rate.df %>% dplyr::select(1, 4, 5) %>% 
    mutate(a.Name = sapply(strsplit(absentee.rate.df$Distict, " - "), `[`, 1)) %>% 
    mutate(ID = substring(absentee.rate.df$Distict, regexpr(" - \\d{6} ", absentee.rate.df$Distict) + 3, 
                          regexpr(" - \\d{6} ", absentee.rate.df$Distict) + 7)) %>% 
    mutate(a.county = substring(absentee.rate.df$Distict, regexpr("\\(.*\\)", absentee.rate.df$Distict) + 1, regexpr("\\)", absentee.rate.df$Distict) - 1)) %>% 
    dplyr::select(-1) %>% 
    rename(Chronic_Absentee = "Chronic Rate", Attendance_Rate = "Attendance Rate")
# Convert character string percentages to numeric proportions
absentee.rate.df$Chronic_Absentee <- as.numeric(gsub("%", "", absentee.rate.df$Chronic_Absentee)) / 100
absentee.rate.df$Attendance_Rate <- as.numeric(gsub("%", "", absentee.rate.df$Attendance_Rate)) / 100
str(absentee.rate.df)

### Merge these last two df's and identify rows with multiple matches;
# A visual scan shows that they are all individual schools, not districts, so they will get filtered out eventually
merged.df <- merge(grad.rate.df, absentee.rate.df, by = "ID", all.x = TRUE)
multiple_matches <- merged.df[duplicated(merged.df$ID), ]
print(multiple_matches)


### Left join to the predictors data frame after making ID's match
sort(unique(predictors.df$IRN))
sort(unique(merged.df$ID))
predictors.df <- predictors.df %>% 
    mutate(ID = substr(IRN, start = 1, stop = 5)) %>%
    dplyr::select(-2)
sort(unique(predictors.df$ID))
predictors.df <- left_join(predictors.df, merged.df, by = "ID")

### Drop redundant columns (these checked out as all equal once NA's were removed)
# Merge with county level data after fixing column title
# Take care of county name that got entered strangely manually
predictors.df <-predictors.df %>% 
    dplyr::select(-c("a.Name", "a.county")) %>% 
    rename(County = "county")
predictors.df$County <- ifelse(predictors.df$County == "Formerly Berlin-Milan", "Erie", predictors.df$County)
predictors.df <- left_join(predictors.df, county.data.df, by = "County")
str(predictors.df)


### Estimate the private enrollment in each district by obtaining the public county sum, attaching to
# predictors data frame, then multiplying district enrollment by ratio of private to public enrollment in that county.
public.county.totals.df <- predictors.df %>% 
    group_by(County) %>%
    summarize(Pub_county_sum = sum(Enrollment22)) %>%
    ungroup()
predictors.df <- left_join(predictors.df, public.county.totals.df, by = "County")
predictors.df <- predictors.df %>%
    mutate(Est_private = round(Enrollment22 * Pri_county_sum / Pub_county_sum ))
str(predictors.df)


### Get pubic school enrollment counts from 2017 and attach to the response df by IRN 
DPR.sheet.names.17 <- excel_sheets("./Data-sources/DPR_FY17.xlsx")
DPR.df.17 <- read_excel("./Data-sources/DPR_FY17.xlsx", sheet = "District Data")
DPR.df.17 <- DPR.df.17 %>% rename(Enrollment17 = "District Total Average Daily Membership FY17")
DPR.df.17 <- DPR.df.17[-1, c("IRN", "Enrollment17")]
str(DPR.df.17)
response.df <- left_join(response.df, DPR.df.17, by = "IRN")


### Now need to change the ID column in the response df into a GEOID column so that this df can be
# attached to both the response df and the sf df:
response.df$ID <- as.character(response.df$IRN)
sort(unique(response.df$ID))


### Take care of the two cases whose ID strings are wrong length manually
for (i in 1:length(response.df$ID)) {
    if (response.df$ID[i] == "139303") {
        response.df$ID[i] <- "13930"
    } else if (response.df$ID[i] == "442") {
        response.df$ID[i] <- "00044"
    } else {
        response.df$ID[i] <- paste("0", substr(response.df$ID[i], 1, 4), sep = "")
    }
}


### Join all predictors to response variables
master.df <- left_join(predictors.df, response.df, by = "ID")
str(master.df)


### Reorganize columns to get all identifying columns at head for comparison
master.df <- master.df %>% 
    dplyr::select(-Name, -District) %>%
    relocate(ID, IRN, .before = Cupp_District) %>% 
    relocate(Homeschoolers22, Homeschoolers17, .after = Cupp_District) %>% 
    relocate(Enrollment22, Enrollment17, .after = County)
sum(complete.cases(master.df))


### Print the names of columns with NA values
problem.cases <- master.df %>% filter(!complete.cases(master.df))
columns_with_na <- colSums(is.na(problem.cases)) > 0
print(names(problem.cases)[columns_with_na])
# Examine rows with NA in the IRN slot
for (j in 1:nrow(master.df)) {
    if (is.na(master.df[j, "IRN"])) {
        print (master.df[j, "Cupp_District"])
        print(j)
    }
}
# Internet search returns correct values for these five districts; fix manually:
master.df[48, "IRN"] <- 46383
master.df[152, "IRN"] <- 48512
master.df[178, "IRN"] <- 45914
master.df[316, "IRN"] <- 48553
master.df[564, "IRN"] <- 49148


### Fix 38 incorrect IRN values (discovered previously through anti-join) manually so that the join to the sf.base df works:
replace_values <- c("04360" = "04361", "04361" = "04360", "04975" = "10004", "04741" = "10002", "04356" = "10017", 
                    "04372" = "10030", "04960" = "10026", "04530" = "10020", "04534" = "10008", "04784" = "10013",
                    "04398" = "10000", "04980" = "10005", "04417" = "10025", "04544" = "10007", "05036" = "10022",
                    "04432" = "10019", "04692" = "10010", "13930" = "00094", "04442" = "10012", "05057" = "10033",
                    "04972" = "10021", "04658" = "10016", "04462" = "10015", "04462" = "10027", "04559" = "10028",
                    "04789" = "10014", "04473" = "10006", "04477" = "10009", "04478" = "10003", "09139" = "04926",
                    "04929" = "04900", "04562" = "10024", "04746" = "10001", "04496" = "10023", "04497" = "10029",
                    "05050" = "10018", "04501" = "10011", "04512" = "10032")                    
master.df$GEOID <- ifelse(master.df$ID %in% names(replace_values), replace_values[master.df$ID], master.df$ID)
master.df$GEOID <- paste("39", master.df$GEOID, sep = "")
master.df <- rename(master.df, District_name = "Cupp_District")

### All districts have enrollment numbers for '22, but five do not for '17
summary(master.df$Enrollment22)
summary(master.df$Enrollment17)

### Columbus City SD enrollment dropped from 72 K to 44 K
### Cleveland Municipal SD enrollment dropped from 55 K to 33 K
master.df[which.max(master.df$Enrollment17), c("District_name", "Enrollment22") ]
temp.df <- master.df[-113, ]
summary(temp.df$Enrollment22)
summary(temp.df$Enrollment17)
temp.df[which.max(temp.df$Enrollment17), c("District_name", "Enrollment22")]


### Initialize indicator variables
master.df$cen.17.ind = rep(NA, nrow(master.df))
master.df$cen.22.ind = rep(NA, nrow(master.df))
master.df$miss.ind = rep(NA, nrow(master.df))
# If either response variable is censored, set censored indicator to TRUE for that year, 
# then also for an either/or indicator:
for (i in 1:nrow(master.df)) {
    if ((!is.na(master.df$Homeschoolers17[i])) & (master.df$Homeschoolers17[i] == "<10")) {
        master.df$cen.17.ind[i] = TRUE
    } else {
        master.df$cen.17.ind[i] = FALSE
    } 
}
for (i in 1:nrow(master.df)) {
    if ((!is.na(master.df$Homeschoolers22[i])) & (master.df$Homeschoolers22[i] == "<10")) {
        master.df$cen.22.ind[i] = TRUE
    } else {
        master.df$cen.22.ind[i] = FALSE
    }
}
for (i in 1:nrow(master.df)) {
    if (master.df$cen.17.ind[i] | master.df$cen.22.ind[i]) {
        master.df$cen.ind[i] = TRUE
    } else {
        master.df$cen.ind[i] = FALSE
    }
}
# If either response variable is missing, set missing indicator to TRUE:
for (i in 1:nrow(master.df)) {
    if (is.na(master.df$Homeschoolers17[i]) | is.na(master.df$Homeschoolers22[i])) {
        master.df$miss.ind[i] = TRUE
    } else {
        master.df$miss.ind[i] = FALSE
    }
}


### Of the 606 districts, 467 have both response values, 102 are censored in one or both years, 
# 34 are missing in one or both years, and 3 are censored in one year and missing in the other
summary(master.df$cen.ind)
summary(master.df$miss.ind)
sum((master.df$cen.ind) & (master.df$miss.ind))
sum((!master.df$cen.ind) & (!master.df$miss.ind))
print(master.df[(master.df$cen.ind) & (master.df$miss.ind), 
                c("District_name", "Enrollment17", "Enrollment22", "Homeschoolers17", "Homeschoolers22")])

### Find out which districts have NA enrollments in data frame
# For now, impute 2022 enrollment for these 5 districts:
print(master.df[which(is.na(master.df$Enrollment17)), c("District_name", "Enrollment22", "Homeschoolers17", "Homeschoolers22")])
master.df[is.na(master.df$Enrollment17), ]$Enrollment17 <- master.df[is.na(master.df$Enrollment17), ]$Enrollment22


### Change homeschool counts to integer values (censored will show as NA)
# and create the ultimate response vector: Delta_Ratio
master.df <- master.df %>% 
    mutate(Home17 = as.integer(Homeschoolers17), Home22 = as.integer(Homeschoolers22)) %>%
    mutate(Delta_Ratio = Home22 / (Home22 + Enrollment22 + Est_private) - Home17 / (Home17 + Enrollment17 + Est_private))
summary(master.df$Delta_Ratio)

### Create bounds for censored imputation
# First confirm that all censored data has legitimate enrollment values
print(master.df[master.df$cen.ind, c("Enrollment17", "Enrollment22")], n = 105)
# Create columns of bounds for imputation of censored data:
master.df$max.Delta <- rep(NA, nrow(master.df))
master.df$min.Delta <- rep(NA, nrow(master.df))
for (i in 1:nrow(master.df)) {
    if (master.df$cen.17.ind[i] & !master.df$cen.22.ind[i]) {
        master.df$max.Delta[i] <- master.df$Home22[i] / (master.df$Home22[i] + master.df$Enrollment22[i] + master.df$Est_private[i]) 
        master.df$min.Delta[i] <- master.df$Home22[i] / (master.df$Home22[i] + master.df$Enrollment22[i] + master.df$Est_private[i]) - 9 / (9 + master.df$Enrollment17[i] + master.df$Est_private[i])
    } else if (!master.df$cen.17.ind[i] & master.df$cen.22.ind[i]) {
        master.df$max.Delta[i] <- 9 / (9 + master.df$Enrollment22[i] + master.df$Est_private[i]) - master.df$Home17[i] / (master.df$Home17[i] + master.df$Enrollment17[i] + master.df$Est_private[i])
        master.df$min.Delta[i] <-  - master.df$Home17[i] / (master.df$Home17[i] + master.df$Enrollment17[i] + master.df$Est_private[i])
    } else {
        master.df$max.Delta[i] <- 9 / (9 + master.df$Enrollment22[i] + master.df$Est_private[i]) 
        master.df$min.Delta[i] <- -9 / (9 + master.df$Enrollment17[i] + master.df$Est_private[i])
    }
}
 
# Handle three unique districts that were missing in '17 and censored in '22
for (j in c(172, 396, 492)) {
    master.df$max.Delta[j] <- 9 / (9 + master.df$Enrollment22[j] + master.df$Est_private[j]) 
    master.df$min.Delta[j] <- -9 / (9 + master.df$Enrollment17[j] + master.df$Est_private[j])
}
print(master.df[master.df$cen.ind, c("Enrollment17", "Enrollment22", "min.Delta", "max.Delta")], n = 105)



### See which rows in sf file don't have data in master.df
Extra.rows.df <- anti_join(OH.acs.df, master.df, by = "GEOID")
print(Extra.rows.df[,1:2])
Extra.rows.2.df <- anti_join(master.df, OH.acs.df, by = "GEOID")
print(Extra.rows.2.df[, c(2, 3, 30)])
### Fix 3 incorrect GEOID values manually so that a join to the sf.base df will work:
rows_without_match <- which(!master.df$GEOID %in% OH.acs.df$GEOID)
rows_without_match
master.df[306, "GEOID"] <- "3900537"
master.df[316, "GEOID"] <- "3910031"
master.df[439, "GEOID"] <- "3910027"
Extra.rows.2.df <- anti_join(master.df, OH.acs.df, by = "GEOID")
print(Extra.rows.2.df[, c(2, 3, 30)])

### Add column to indicate which 5 rows in map base do not have any data at all
### This will be helpful when attaching new columns that are model outputs to the map
OH.acs.df$invalid <- !OH.acs.df$GEOID %in% master.df$GEOID
summary(OH.acs.df$invalid)

saveRDS(OH.acs.df, "OH.acs.df")
saveRDS(master.df, "master.df")





