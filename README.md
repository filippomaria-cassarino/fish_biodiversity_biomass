# Fish diversity, biomass, and human pressures

This repository contains the code and data that I developed for 
my master's thesis project. In this current version, they are used to
model the response of fish biomass to fish diversity, fishing, and climate change 
in the Barents Sea.

## Data

Original data includes:

1) Fish abundance and biomass data from the Barents Sea Ecosystem Survey
(Eriksen et al., 2018). The data is stored in RData files named 
NOR-BTS_clean.RData (2004-2021) and BESS_data2022.RData (2022). The first file
was retrieved from the FISHGLOB_data dataset by Maureaud et al. (2024)
(https://github.com/fishglob/FishGlob_data), while the latter was shared
directly by Laurene Pecuchet since it has not been incorporated to FISHGLOB_data
yet.

2) Functional trait information (Beukhof et al., 2019) in a csv file named
trait_beukhof.csv

3) Oceanographic data from the E.U. Copernicus Marine Service, which is not 
available here because of its size. This data was compiled from
the Global Ocean Physics Reanalysis (https://doi.org/10.48670/moi-00021), where 
it was downloaded as mean monthly values at a 0.083° resolution, and from 
Global Ocean Colour (https://doi.org/10.48670/moi-00281) at a 4km resolution.
In both cases, the geographical range of the data to download was selected in 
order to cover all sampling locations of the biological data.

4) Fishing effort data from the Global Fishing Watch (Kroodsma et al., 2018),
which is not available here because of its excessive size, but can be 
obtained following the code section titled 4_fishing.R and the guide available 
here: https://github.com/GlobalFishingWatch/gfwr. 
