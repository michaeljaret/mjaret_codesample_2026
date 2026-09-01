#for data cleaning and manipulation
library(tidyverse)
#for geo-spatial data and functions
library(sf)
#for working with dates
library(lubridate)
#for file accessibility
library(here)

#define time period
period <- seq(as.Date("2000-01-01"), as.Date("2025-12-01"), by = "month")

#download city polygons from CAL FIRE
st_layers(here("raw_data", "incorp23_1.gdb"))
calfire_cities <- read_sf(here("raw_data", "incorp23_1.gdb"), layer = "incorp23_1")

#download ZHVI data from Zillow
zhvi <- read_csv(here("raw_data", "City_zhvi_uc_sfrcondo_tier_0.33_0.67_sm_sa_month (1).csv"))
zhvi_clean <- zhvi %>% 
  #filter to California data and retain necessary columns
  filter(State == "CA", RegionType == "city") %>%
  select(-RegionID, -SizeRank, -RegionType, -StateName, -State, -Metro, -CountyName) %>%
  #pivot to a panel structure
  pivot_longer(cols = -RegionName, names_to = "month_year", values_to = "price") %>%
  mutate(month_year = floor_date(as.Date(month_year), "month")) %>%
  rename(CITY = RegionName) %>%
  complete(CITY, month_year = period) %>%
  filter(month_year %in% period) %>%
  arrange(CITY, month_year)
#identify city spelling mismatches
setdiff(calfire_cities$CITY, zhvi_clean$CITY)
setdiff(zhvi_clean$CITY, calfire_cities$CITY)
zhvi_panel <- zhvi_clean %>%
  #match cities
  mutate(CITY = case_when(CITY == "Angels Camp" ~ "Angels", 
                          CITY == "California City" ~ "California", 
                          CITY == "La Canada Flintridge" ~ "La Cañada Flintridge",
                          CITY == "Saint Helena" ~ "St. Helena",
                          CITY == "Lone Pine" ~ "McFarland", .default = CITY)) %>%
  filter(CITY %in% calfire_cities$CITY) %>%
  arrange(CITY, month_year) %>%
  mutate(price_log = log1p(price))
#compile final city list
cities <- unique(zhvi_panel$CITY)

#download CAL FIRE FRAP Historical Fire Perimeters data
st_layers(here("raw_data", "fire25_1.gdb"))
calfire_frap <- read_sf((here("raw_data", "fire25_1.gdb")), layer = "firep25_1")
#clean data
fires_clean <- calfire_frap %>%
  mutate(start_date = as.Date(ALARM_DATE), end_date = as.Date(CONT_DATE)) %>%
  #choose 'GlobalID' as a unique fire identifier and rename as 'FIRE' for simplicity
  rename(FIRE = GlobalID) %>%
  filter(
    #remove fires where dates are missing or start date is past end date
    !(is.na(start_date)|is.na(end_date)), 
    !start_date > end_date, 
    #remove fires not active between January 2000 and December 2025
    !(end_date < as.Date("2000-01-01")|start_date > as.Date("2025-12-31"))) %>%
  select(FIRE, GIS_ACRES, start_date, end_date, Shape) %>% 
  mutate(start_month = as.Date(paste0(year(start_date), "-", month(start_date), "-01")), 
         end_month = as.Date(paste0(year(end_date), "-", month(end_date), "-01"))) %>% 
  #transform curved lines into multi-polygons to enable computation
  st_cast("MULTIPOLYGON") %>%
  distinct()

#download block data and polygons data set
st_layers(here("raw_data", "tl_2025_06_tabblock20 (2)"))
blocks <- read_sf(here("raw_data", "tl_2025_06_tabblock20 (2)"), layer = "tl_2025_06_tabblock20")
blocks_clean <- blocks %>% select(GEOID20, POP20, geometry) %>% rename(block = GEOID20, block_pop = POP20) %>% st_transform(crs = st_crs(calfire_cities))
#find geometric centroid and extract its coordinates
blocks_clean <- blocks_clean %>% mutate(centroid = st_centroid(geometry), 
                                        centroid_x = st_coordinates(centroid)[,1],
                                        centroid_y = st_coordinates(centroid)[,2],
                                        block_area = st_area(geometry) %>% as.numeric())
#find intersecting blocks and cities, and the proportion of the block that intersects
block_city <- st_intersection(blocks_clean, calfire_cities)
block_city <- block_city %>% mutate(
  intersection = st_area(geometry) %>% as.numeric(),
  area_weight = intersection/block_area,
  est_pop = block_pop*area_weight)
#find PWC coordinates using weighted means; combine coordinates into new column
city_pwc <- block_city %>% st_drop_geometry() %>% group_by(CITY) %>%
  summarise(pwc_x = weighted.mean(centroid_x, w = est_pop),
            pwc_y = weighted.mean(centroid_y, w = est_pop),
            total_pop = sum(est_pop)) %>%
  st_as_sf(coords = c("pwc_x", "pwc_y"), crs = 3310)
#find geometric centroid coordinates for each city
city_centroid <- calfire_cities %>%
  select(CITY) %>%
  mutate(centroid = st_centroid(SHAPE),
         centroid_x = st_coordinates(centroid)[,1],
         centroid_y = st_coordinates(centroid)[,2]) %>%
  st_drop_geometry()
#compare geometric centroid and PWC coordinates to ensure no outliers
city_compare <- city_pwc %>% 
  mutate(pwc_x = st_coordinates(geometry)[,1],
         pwc_y = st_coordinates(geometry)[,2]) %>%
  left_join(city_centroid, by = "CITY") %>%
  st_drop_geometry() %>%
  mutate(diff_x = centroid_x - pwc_x,
         diff_y = centroid_y - pwc_y)
summary(city_compare %>% select(diff_x, diff_y))

#define exponential decay parameters for continuous exposure calculation
theta <- c('0004' = 0.0004, '0002' = 0.0002, '0001' = 0.0001, '00005' = 0.00005, '000025' = 0.000025)

#find the maximum distance among all parameters in which the exponential decay function returns '0.001' (remove city-fire matches past this distance for computational efficiency)
cutoff <- function(x) {(-log(0.001)/x)}
theta_cutoff <- max(sapply(theta, cutoff))

#find city-fire pairs within the cutoff
within_cutoff <- st_is_within_distance(city_pwc, fires_clean, dist = theta_cutoff)
#creates a list of city-fire pairs using indices
cityfire_pairs <- tibble(city_i = rep(seq_along(within_cutoff), lengths(within_cutoff)),
                         fire_i = unlist(within_cutoff))
#creates a list of city-fire pairs using names and distances
cityfire_pairs <- cityfire_pairs %>% 
  mutate(CITY = city_pwc$CITY[city_i],
         FIRE = fires_clean$FIRE[fire_i],
         distance = st_distance(city_pwc$geometry[city_i],
                                fires_clean$Shape[fire_i],
                                by_element = TRUE) %>% as.numeric()) %>%
  select(CITY, FIRE, distance)
#repeat fires for each month of activity
fires_clean <- fires_clean %>%
  rowwise() %>%
  mutate(month_year = list(seq(start_month, end_month, by = "month"))) %>%
  unnest(month_year) %>%
  ungroup()
#calculate magnitude and duration as part of the continuous exposure measure
fires_summary <- fires_clean %>%
  #for each month find: 1) days in month, and 2) the days the fire was active per month
  mutate(days_month = days_in_month(month_year),
         #determine whether the month or the fire ended first
         month_ceiling = as.Date(pmin(ceiling_date(month_year, unit = "month") - days(1), end_date)),
         #determine whether the month or the fire started first
         month_floor = as.Date(pmax(floor_date(month_year, unit = "month"), start_date)),
         #calculate the number of days between the two values
         days_active = as.numeric(month_ceiling - month_floor) + 1,
         #calculate the ratio of the days the fire was active in the month and the total number of days in the month
         duration = days_active/days_month,
         #calculate magnitude using a log transformation
         magnitude = log1p(GIS_ACRES)) %>% 
  select(FIRE, month_year, start_date, end_date, magnitude, duration, Shape)
#join fire characteristics (duration, magnitude) with city-fire pairs and expand for each active month
exposure <- inner_join(cityfire_pairs, fires_summary %>% st_drop_geometry(), by = "FIRE", relationship = "many-to-many") %>%
  filter(CITY %in% cities)

#construct continuous exposure measure
cont_exp_panel <- exposure %>%
  #calculate exposure by row
  rowwise() %>%
  mutate(cont = list(map_dbl(theta, ~ magnitude*duration*exp(-.x*distance)))) %>%
  ungroup() %>%
  unnest_wider(cont, names_sep = "_") %>%
  group_by(CITY, month_year) %>%
  #combine into city-months observations and sum
  summarise(across(starts_with("cont_"), sum)) %>%
  ungroup() %>%
  #employ inverse hyperbolic sine (asinh) transformation on columns
  mutate(across(starts_with("cont_"), asinh, .names = "asinh_{.col}")) %>%
  #fill city-months with no exposure with zeroes
  complete(CITY = cities, month_year = period, 
           fill = list(cont_0004 = 0, cont_0002 = 0, cont_0001 = 0, cont_00005 = 0, cont_000025 = 0,
                       asinh_cont_0004 = 0, asinh_cont_0002 = 0, asinh_cont_0001 = 0, asinh_cont_00005 = 0, asinh_cont_000025 = 0)
  ) %>%
  arrange(CITY, month_year)

#calculate binary binned exposure variable
binary_exp_panel <- exposure %>%
  filter(distance <= 100000) %>%
  mutate(binary_0_10 = ifelse(distance >= 0 & distance <= 10000, 1, 0),
         binary_10_50 = ifelse(distance > 10000 & distance <= 50000, 1, 0),
         binary_50_100 = ifelse(distance > 50000 & distance <= 100000, 1, 0)) %>%
  group_by(CITY, month_year) %>%
  #combine into city-month observations (summarize with 1 showing at least one fire was active in the bin in the month)
  summarise(across(starts_with("binary_"), ~ as.integer(sum(.x) > 0))) %>%
  ungroup() %>%
  select(CITY, month_year, binary_0_10, binary_10_50, binary_50_100) %>%
  #complete the data set for all cities and months
  complete(CITY, month_year = period, fill = list(binary_0_10 = 0, binary_10_50 = 0, binary_50_100 = 0)) %>%
  arrange(CITY, month_year)

#calculate days binned exposure variable
days_exp_panel <- exposure %>%
  filter(distance <= 100000) %>%
  select(-magnitude, -duration, -FIRE) %>%
  #for each month, create variables defining the beginning and end of each month
  mutate(start_month = floor_date(month_year, unit = "month"),
         end_month = ceiling_date(month_year, unit = "month") - days(1)) %>%
  rowwise() %>%
  #for each month, expand data set to show each day
  mutate(days_active = list(seq(start_month, end_month, by = "day"))) %>%
  unnest(days_active) %>%
  ungroup() %>%
  #define bins with 1 showing the fire was active that day, and 0 otherwise
  mutate(
    active = ifelse(days_active >= start_date & days_active <= end_date, 1, 0),
    days_0_10 = ifelse(distance >= 0 & distance <= 10000, active, 0),
    days_10_50 = ifelse(distance > 10000 & distance <= 50000, active, 0),
    days_50_100 = ifelse(distance > 50000 & distance <= 100000, active, 0)) %>%
  #group by day, with 1 showing at least one fire was active in that bin for the day
  group_by(CITY, days_active) %>%
  summarise(
    days_0_10 = as.integer(sum(days_0_10) > 0),
    days_10_50 = as.integer(sum(days_10_50) > 0),
    days_50_100 = as.integer(sum(days_50_100) > 0),
    .groups = "drop"
  ) %>% ungroup() %>%
  mutate(month_year = floor_date(days_active, "month")) %>%
  #sum to find the number of days with fire activity for that month
  group_by(CITY, month_year) %>%
  summarise(
    days_0_10 = sum(days_0_10),
    days_10_50 = sum(days_10_50),
    days_50_100 = sum(days_50_100),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  #complete the data set for all cities and months
  filter(CITY %in% cities, month_year %in% period) %>%
  complete(CITY = cities, month_year = period, fill = list(days_0_10 = 0, days_10_50 = 0, days_50_100 = 0)) %>%
  arrange(CITY, month_year)

#prepare cities for PRISM data extraction
city_pwc_ll <- city_pwc %>%
  #convert PWCs to longitude and latitude
  st_transform(crs = 4326) %>%
  arrange(CITY) %>%
  filter(CITY %in% cities) %>%
  mutate(
    #uniquely identify cities by shortening and adding row number (PRISM cuts names to 12 characters, potentially creating duplicated names)
    CITY_short = substr(paste0(row_number(), CITY), 1, 12),
    lat = st_coordinates(geometry)[,2],
    lon = st_coordinates(geometry)[,1]) %>%
  st_drop_geometry()
#create a key to match cities to the shortened version later
city_key <- city_pwc_ll %>% select(CITY, CITY_short)
#simplify data frame to fit PRISM requirements
city_pwc_ll <- city_pwc_ll %>% select(lat, lon, CITY_short)
#ensure there are 478 unique city names to avoid mix-ups later
length(unique(city_pwc_ll$CITY_short))
#remove header to support PRISM requirements
names(city_pwc_ll) <- NULL
#create .csv file to upload into PRISM
write.csv(city_pwc_ll, here("source_data", "prism_cities_ll.csv"), row.names = FALSE)
#extract data from https://prism.oregonstate.edu/explorer/bulk.php, which permits 15 years per download
weather_00_14 <- read.csv(here("raw_data", "PRISM_ppt_tmean_stable_800m_200001_201412.csv"))
weather_15_25 <- read.csv(here("raw_data", "PRISM_ppt_tmean_stable_800m_201501_202512.csv"))
#clean PRISM data (start with 2000 to 2014)
col_names <- weather_00_14 %>%
  #remove top 9 rows
  slice(-(1:9)) %>%
  #data is placed into single column, so reorganize first 7 rows into the column names
  slice(1:7) %>% 
  pull(PRISM.Time.Series.Data)
weather_00_14 <- weather_00_14 %>%
  slice(-(1:9)) %>%
  slice(8:n()) %>%
  #disperse data past first 7 rows into the 7 columns
  mutate(
    row = (row_number() - 1) %/% 7 + 1,
    col = (row_number() - 1) %% 7 + 1) %>%
  pivot_wider(names_from = col, values_from = PRISM.Time.Series.Data) %>%
  select(-row) %>%
  #rename columns with the names from the first 7 rows
  rename_with(~ col_names)
#repeat for second data set (2014-2025)
weather_14_25 <- weather_14_25 %>%
  slice(-(1:9)) %>%
  slice(8:n()) %>%
  mutate(
    row = (row_number() - 1) %/% 7 + 1,
    col = (row_number() - 1) %% 7 + 1) %>%
  pivot_wider(names_from = col, values_from = PRISM.Time.Series.Data) %>%
  select(-row) %>%
  rename_with(~ col_names)
#combine two weather data sets into one 
weather_panel <- bind_rows(weather_00_14, weather_15_25) %>%
  rename(CITY_short = Name, month_year = Date, ppt = 'ppt (inches)', temp = 'tmean (degrees F)') %>%
  #add back in usual city names
  left_join(city_key, by = "CITY_short") %>%
  mutate(month_year = as.Date(paste0(month_year, "-01")),
         ppt = as.numeric(ppt),
         temp = as.numeric(temp)) %>%
  select(CITY, month_year, ppt, temp) %>%
  filter(CITY %in% cities, month_year %in% period) %>%
  arrange(CITY, month_year)

#combine all panels into large panel for regressions
panel_all <- left_join(zhvi_panel, cont_exp_panel, by = c("CITY", "month_year")) %>% 
  left_join(binary_exp_panel, by = c("CITY", "month_year")) %>%
  left_join(days_exp_panel, by = c("CITY", "month_year")) %>%
  left_join(weather_panel, by = c("CITY", "month_year")) %>%
  left_join(calfire_cities %>% 
              #add county and year clustered standard errors specification (Cameron and Miller, 2015)
              select(CITY, COUNTY) %>% 
              st_drop_geometry(), by = "CITY") %>%
  mutate(year = year(month_year)) %>%
  #add latitude and longitude for Conley (1999) spatial error robustness check
  left_join(city_pwc %>% 
              st_transform(crs = 4326) %>%
              mutate(lat = st_coordinates(geometry)[,2],
                     lon = st_coordinates(geometry)[,1]) %>% 
              st_drop_geometry() %>% 
              select(CITY, lat, lon), by = "CITY")

saveRDS(panel_all, here("source_data", "panel_all.rds"))
