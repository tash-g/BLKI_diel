# Evidence for latitude-driven changes in circadian rhythm in a wide-ranging seabird 
Natasha Gillies, Sébastien Descamps, Frederick McKendrick, Samantha C. Patrick

## Overview
This repository contains scripts and data to recreate the main results and figures of this paper (currently in prep). 

## Scripts
A short description of each script is given below.

- **1_extract-trips-and-behaviour.R** Extracts colony visitation/foraging trips and at-sea behaviour (rest, flight, foraging) from salt water immersion data collected from geolocator devices. 
- **2_prepare_attendance_data.R** Converts colony visitation data above into a time series for ACF analysis.
- **3_statistical_analysis.R** Code and functions to carry out all statistical analyses in the manuscript. 
- **3b_statistic_analysis-randomised-GAM.R** Randomly samples attendance data and fits to a GAM examining the predictive power of sun elevation angle on activity levels. Includes code to run the randomisation on a high-throughput system such as _Condor_.
- **BLKI_functions.R** Contains custom functions to run each of the above scripts. Must be run for the scripts to work effectively.

## Data inputs

These data are used in the above scripts. Note that all Rings/BirdIDs have been recoded and so cannot be linked to existing datasets. Please contact the authors if you would like to make use of these datasets, as we may be able to offer additional information, data, or advice. 

- **BLKI_attendance_2013-2022.Rdata** This is the main dataset for analysis. Each row represents an hour of a day for each individual bird. Data columns are as follows:
  - _datetime_: Date and hour of record
  - _Ring_: Factor encoding unique ID of bird
  - _colony_: Factor encoding the colony of provenance for the bird
  - _col_lat_: Numeric latitude of colony
  - _col_lon_: Numeric longitude of colony
  - _att_time_: Number of minutes of hour bird spent at colony
  - _year_: Factor encoding year of record
  - _num_visits_: Number of unique visits made by the bird
  - _time_ID_: Factor encoding time index of bird dataset
  - _att_poss_: Number of minutes of data recorded for the hourly record
  - _prop_att_: Numeric proportion of hour bird spent in attendance at colony
  - _n_flight_: Number of 10 minute flight bouts recorded
  - _n_rest_: Number of 10 minute rest bouts recorded
  - _prop_forage_: Proportion of hour bird spent foraging
  - _rest_flight_ratio_: Ratio of time spent in rest vs flight behaviour
  - _sun altitude_: Numeric sun altitude, in degrees. Sun altitude is the angular elevation of the center of the solar disk above the horizon
  - _temp_: Numeric temperature in degrees celsius
  - _hour_: Numeric value for hour of day; 24-hour format

- **BLKI_breeding-success.csv** Gives population-level annual breeding success for each colony. Data columns are as follows:
  - _year_: Factor encoding year of record
  - _n_nests_: Number of nests data collected from
  - _large_chicks_per_nest_: Numeric proportion of nests producing large chicks - this is our measure of breeding success
  - _colony_: Factor encoding the colony from which data were collected

- **BLKI_gls-daily-behaviour.RData** Dataset summarising daily behaviour for each bird. Data colums are as follows:
  - _ring_: Factor encoding unique ID of bird
  - _date_: Date of record
  - _n_recs24/day/night_: Number of 10 minute bouts per 24hour period/daylight hours/night hours for which data are available (indiviudal columns indicated by '/')
  - _prop_forage24/day/night_: Proportion of 24 hour period/daylight hours/night hours spent foraging (indiviudal columns indicated by '/'. 'nocol' in column heading indicates that colony attendance excluded from calculation)
  - _prop_rest24/day/night_: Proportion of 24 hour period/daylight hours/night hours spent resting (indiviudal columns indicated by '/'. 'nocol' in column heading indicates that colony attendance excluded from calculation)
  - _prop_flight24/day/night_: Proportion of 24 hour period/daylight hours/night hours spent in flight (indiviudal columns indicated by '/'. 'nocol' in column heading indicates that colony attendance excluded from calculation)
  - _col_att24/day/night_: Proportion of 24 hour period/daylight hours/night hours spent attending the colony (indiviudal columns indicated by '/')
  - _daylight.mins_: Number of minute of daylight per 24 hour period
  - _total_imm.day/night_: Total immersion recorded per daylight/night time hours (individual columns indicated by '/')
  - _colony_: Factor encoding the colony of provenance for the bird
  - _col_lat_: Numeric latitude of colony
  - _col_lon_: Numeric longitude of colony

- **BLKI_gls-hourly-behaviour.RData** Dataset summarising hourly daily behaviour for each bird. Data colums are as follows:
  - _colony_: Factor encoding the colony of provenance for the bird
  - _ring_: Factor encoding unique ID of bird
  - _date_: Date of record
  - _hour_: Numeric hour of record
  - _total_imm_: Numeric total immersion per hour
  - _n_flight_: Number of 10 minute flight bouts per hour
  - _n_rest_: Number of 10 minute rest bouts per hour
  - _n_forage_: Number of 10 minute foraging bouts per hour
  - _for_poss_: Number of 10 minute bouts for which data are available (maximum = 6).
  - _prop_forage_: Proportion of hour spent foraging
  - _rest_flight_ratio_: Ratio of time spent in rest vs flight behaviour

- **BLKI_gls-trips.RData** Condensed dataset summarising colony attendance or trip records for each bird; used to compare outputs of trips identified via GLS versus GPS data. Each row represents an individual bird record (trip or colony visit). Data columns are as follows:
  - _ring_: Factor encoding unique ID of bird
  - _type_: Whether record identified using geolocator (GLS) or GPS data
  - _start_: Datetime variable (POSIXct) indicating start time of record
  - _loc_: Whether the record is for a colony visit (colony) or trip to sea (trip)
  - _end_: Datetime variable (POSIXct) indicating end time of record
  - _duration_mins_: Duration of record in minutes
  - _colony_: Factor encoding the colony of provenance for the bird
  - _col_lat_: Numeric latitude of colony
  - _col_lon_: Numeric longitude of colony
 
- **BLKI_metadata.RData** Abridged dataset giving colony locations of each individual. Each row contains one individual bird. Data columns are as follows:
  - _ring_: Factor encoding unique ID of bird
  - _colony_: Factor encoding the colony of provenance for the bird
  - _col_lat_: Numeric latitude of colony
  - _col_lon_: Numeric longitude of colony
    
- **Grumant/NyAPyr/SEATRACK_GLS_processed.RData** Individual datasets separated by '/'. Datasets containing the raw immersion data for all birds, processed for consistent formatting and to calculate behavioural states. Each row indicates a 10 minute bout. Data columns are as follows:
  - _datetime_: Datetime variable (POSIXct) indicating start time of 10 minute bout
  - _immersion_: Numeric total salt water immersion recorded in the 10 minute bout. Geolocators sample every 3 seconds and log the total number of salt-water immersion events recorded every 10 minutes. Can therefore vary between 0 (no immersion) and 200 (complete immersion)
  - _ring_: Factor encoding unique ID of bird
  - _light_: Numeric light level recorded in lux
  - _behaviour_: Factor indicating most probable behavioural state based on behaviour-finding algorithm
