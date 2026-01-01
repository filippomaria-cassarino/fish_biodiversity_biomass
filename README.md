# Fish diversity, biomass, and human pressures

This repository contains the data and code that I developed for 
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

All intermediate and final data is available.

## Code

File 0 (0_functions.R) produces the custom functions used in the 
analysis.

Files 1-5 explore and prepare the data.

Files 6-7 model the relationships under investigation.

File 8 produces the main figures.

This section is currently under revision. Further details and references will
be made available once the code is finalized.

## References
Beukhof, E., Dencker, T. S., Palomares, M. L. D., & Maureaud, A. (2019). 
A trait collection of marine fish species from North Atlantic and
Northeast Pacific continental shelf seas [Data set]. PANGAEA.
https://doi.org/10.1594/PANGAEA.900866

Eriksen, E., Gjøsæter, H., Prozorkevich, D., Shamray, E., Dolgov, A., 
Skern-Mauritzen, M., Stiansen, J. E., Kovalev, Yu., & Sunnanå, K. (2018).
From single species surveys towards monitoring of the Barents Sea ecosystem. 
Progress in Oceanography, 166, 4–14. https://doi.org/10.1016/j.pocean.2017.09.007

Kroodsma, D. A., Mayorga, J., Hochberg, T., Miller, N. A., Boerder, K., Ferretti,
F., Wilson, A., Bergman, B., White, T. D., Block, B. A., Woods, P., Sullivan, B.,
Costello, C., & Worm, B. (2018). Tracking the global footprint of fisheries.
Science, 359(6378), 904–908. https://doi.org/10.1126/science.aao5646

Maureaud, A. A., Palacios-Abrantes, J., Kitchel, Z., Mannocci, L., Pinsky, M. L.,
Fredston, A., Beukhof, E., Forrest, D. L., Frelat, R., Palomares, M. L. D., Pecuchet,
L., Thorson, J. T., Van Denderen, P. D., & Mérigot, B. (2024).
FISHGLOB_data: An integrated dataset of fish biodiversity sampled with
scientific bottom-trawl surveys. Scientific Data, 11(1), 24.
https://doi.org/10.1038/s41597-023-02866-w
