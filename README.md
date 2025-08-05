# Evidence for latitude-driven changes in diel rhythms in a wide-ranging seabird 
Natasha Gillies, Sébastien Descamps, Nicholas P. Huffeldt, Frederick McKendrick, Tycho Anker-Nilssen, Maria Bogdanova, Vegard Sandøy Bråthen, Olivier Chastel, Signe Christensen-Dalsgaard, Francis Daunt, Nina Dehnhard, Alexey Viktorovich Ezhov, Annette Fayet, Morten Frederiksen, Maria Gavrilo, April Hedd, Yann Kolbeinsson, Yuri Vladimirovich Krasnov, Aili Labansen, Jannie Linnebjerg, Svein-Håkon Lorentsen, Flemming R. Merkel, Børge Moe, Mark Newell, Stephen Newton, Bergur Olsen, Tone Reiertsen, Greg Robertson, Hallvard Strøm, Thorkell Lindberg Thorarinsson, Samantha C. Patrick

## Overview
This repository contains scripts and data to recreate the main results and figures of this paper (currently in submission). Analyses have been tested in R version 4.3.2 (2023-10-31 ucrt). Individual packages are listed under the 'Session Info' of each R Markdown file.

## Scripts
A short description of each script is given below. All files are provided as R markdown documents to aid interpretation and reproducibility. Note that the raw geolocator/salt-water immersion data are provided by the [SEATRACK](https://seapop.no/en/seatrack/) project and so are not shared here. Those interested in using these data should contact the authors/SEATRACK for further information.

- **1_data_processing.Rmd** Extracts colony visitation/foraging trips and at-sea behaviour (rest, flight, foraging) from saltwater immersion data collected via geolocators. The original data are not included (see above). Provided for transparency; not needed to reproduce manuscript results.
- **2_main_analysis.Rmd** Code and functions to carry out all statistical analyses in the manuscript. 
- **3_plot_actograms.Rmd** Plots actograms for all birds and the representative examples included in the manuscript.
- **BLKI_functions.R** Contains custom functions to run each of the above scripts. 
- **SI_validate_trips.Rmd** Compares trips identified via GPS with those identified from immersion data. For transparency only; not required for main analyses.

## Data inputs

These data are used in the above scripts. Note that all Rings/identities have been recoded and so cannot be linked to existing datasets. Please contact the authors if you would like to make use of these datasets, as we may be able to offer additional information, data, or advice. 

- **BLKI_breeding_dates.csv** Provides colony-specific breeding dates used to define behavioural analysis windows. Includes the estimated minimum lay date (min_lay) and the start (min_subset) and end (max_subset) dates for data subsetting.

- **BLKI_example_actograms.csv** Lists one representative individual per colony used to generate example actograms shown in the manuscript. Includes anonymised bird IDs (ring) matched to colony names.

- **BLKI_gls-hourly-behaviour.RData** Dataset summarising hourly daily behaviour for each bird. A .csv alternative is also provided. Data columns are as follows:
  - _colony_: Factor encoding the colony of origin for the bird
  - _year_: Year of record
  - _ring_: Factor encoding unique ID of bird (anonymised)
  - _col_lat_: Colony latitude
  - _col_lon_: Colony longitude
  - _datetime_: Date and time of record rounded to the nearest hour; POSIXct
  - _date_: Date of the record
  - _hour_: Numeric hour of record
  - _imm_total_: Numeric total immersion per hour
  - _activity_mins_: Total number of minutes in the activity state
  - _bout_dur.mins_: Total duration of the bout (record) in minutes
  - _sun_altitude_: Solar elevation angle (in degrees) at the recorded local datetime and location. Positive values indicate the sun is above the horizon; negative values indicate it is below the horizon (i.e., night)
  - _datetime_UTC_: Date and time of the record rounded to the nearest hour and converted to UTC; POSIXct
  - _hour_UTC_: Numeric hour of record converted to UTC

- **BLKI_metadata_anon.RData** Abridged dataset giving colony locations of each individual. A .csv alternative is also provided. Each row contains one individual bird. Data columns are as follows:
  - _ring_: Factor encoding unique ID of bird
  - _colony_: Factor encoding the colony of provenance for the bird
  - _col_lat_: Numeric latitude of colony
  - _col_lon_: Numeric longitude of colony