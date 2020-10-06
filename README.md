## Tracking the global reduction of marine traffic during the COVID-19 pandemic

[![DOI](https://zenodo.org/badge/280381188.svg)](https://zenodo.org/badge/latestdoi/280381188)

This repository provides the R code that accompanies the article:

March D, Metcalfe K, Tintoré J, Godley BJ. Tracking the global reduction of marine traffic during the COVID-19 pandemic.


### Requirements
* R-studio with R >= 3.6.0


### Structure of this repostitory


Folder          |  Description    
--------------- | -------------------
analysis        | scripts used to analyse the data
data            | folder to store input and output data (see Data repository)
results         | figures and files generated during the analysis
scr             | scripts with custom fuctions that are sourced from different scripts



### Data repository: input data

Input data `data/input` contains all datasets required for the analysis. Because of size and license constrains to redistribute data, this folder is not included in the repository. I provide here the paths and sources of data.

Dataset (version) |  Path                         | Description    
----------------- | ----------------------------- | --------------------- 
GADM (v 3.6)      | `data/input/gadm`             | Database of Global Administrative Areas (https://gadm.org/)
EEZ (v 11)        | `data/input/marine_regions/World_EEZ_v11_20191118_gpkg`  | Exclusive Economic Zones (https://www.marineregions.org/)
S-AIS             | `data/input/exacEarth/`       | S-AIS datasets from ExactEarth (https://www.exactearth.com/)
FAO areas         | `data/input/FAO_AREAS/`       | FAO major fishing areas (http://www.fao.org/fishery/area/search/en)
Marine ecoregions | `data/input/MEOW/`            | Marine ecoregions of the world (https://www.worldwildlife.org/publications/marine-ecoregions-of-the-world-a-bioregionalization-of-coastal-and-shelf-areas)
High Seas         | `data/input/marine_regions/Intersect_EEZ_IHO_v4_2020` | Intersection between EEZ and IHO areas (https://www.marineregions.org/)
Income class      | `data/input/worldbank`        | Income class per country by World Bank (https://datahelpdesk.worldbank.org/knowledgebase/articles/906519-world-bank-country-and-lending-groups)



### License

Copyright (c) 2019 David March  
Licensed under the [MIT license](https://github.com/dmarch/abigoos/blob/master/LICENSE).



### Acknowledgements

We acknowledge support from the European Union's Horizon 2020 research and innovation programme under the Marie Skłodowska‐Curie grant agreement no 794938. We thank Marine Traffic and Exact Earth for their support with AIS data. We thank all the data providers for making their data available and making this study possible. 
