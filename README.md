# covid19-ais
Tracking the global reduction of marine traffic during the COVID-19 pandemic



# Structure


Folder          |  Description    
--------------- | -------------------
analysis        | scripts used to analyse the data
data            | folder to store input and output data (see Data repository)
results         | figures and files generated during the analysis
scr             | scripts with custom fuctions that are sourced from different scripts
data            | folder with input and output data (not included the repository)




# Input data repository

Input data `data/input` contains all datasets required for the analysis. Because of size and license constrains to redistribute data, I provide here the paths and sources of data

Dataset (version) |  Path                         | Description    
----------------- | ----------------------------- | --------------------- 
GADM (v 3.6)      | `data/input/gadm`             | Database of Global Administrative Areas (https://gadm.org/)
EEZ (v 11)        | `data/input/marine_regions/World_EEZ_v11_20191118_gpkg`  | Exclusive Economic Zones (https://www.marineregions.org/)
S-AIS             | `data/input/exacEarth/`       | S-AIS datasets from ExactEarth (https://www.exactearth.com/)
FAO areas         | `data/input/FAO_AREAS/`       | FAO major fishing areas (http://www.fao.org/fishery/area/search/en)
Marine ecoregions | `data/input/MEOW/`            | Marine ecoregions of the world (https://www.worldwildlife.org/publications/marine-ecoregions-of-the-world-a-bioregionalization-of-coastal-and-shelf-areas)
High Seas         | `data/input/marine_regions/Intersect_EEZ_IHO_v4_2020` | Intersection between EEZ and IHO areas (https://www.marineregions.org/)
Income class      | `data/input/worldbank`        | Income class per country by World Bank (https://datahelpdesk.worldbank.org/knowledgebase/articles/906519-world-bank-country-and-lending-groups)
