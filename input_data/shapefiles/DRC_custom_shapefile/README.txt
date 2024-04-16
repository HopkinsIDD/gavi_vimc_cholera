This directory includes files used to produce the custom vector used for the DRC Case study, which includes health zone locations
from Humanitarian Data Exchange instead of districts (administrative level 2) sourced from GADM.

COD_adm2.sf.rds: the shapefile used in the VIMC Core Model

DRC_targeting_raw.csv: a file with the targeting information for the DRC Case study (not the custom targeting table used in the pipeline) 

osm_rdc_sante_zones_211212.gpkg: a geopackage with DRC health zones from Humanitarian Data exchange

metadata-osm_rdc_sante_zones_211212-gpkg.csv: a file with matadata for drc_admin2_gadm.gpkg

process_custom_shapefile.R: an R Script that takes as input osm_rdc_sante_zones_211212.gpkg and COD_adm2.sf.rds and creates an sf (stored as an .rds file)
with the health zones following the format used for COD_adm2.sf.rds. This sf is used in the new pipeline for the DRC Case study.

custom_shapefile.rds: the output of process_custom_shapefile.R

