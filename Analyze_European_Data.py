#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script used to determines how much cells have been calculated on the domain and
creates a GeoJson containing the cells that have been calculated.

@author: Jérémy Bernard, Météo-France
#date: 2024-04-09
"""

import os
import geopandas as gpd
import pandas as pd
import random
import sys
import matplotlib.pylab as plt

##################################################################################
#########       DEFINE WHAT THE CODE IS GOING TO DO         ######################
##################################################################################

# Need to recalc, replot or reload GIS file
recalc = True
replot = True
reload = False

if not reload:
    ################################################################################
    ######### LOAD THE QGIS ENVIRONMENT TO USE QGIS FUNCTIONS ######################
    ################################################################################
    from qgis.core import QgsApplication

    # # Starts the qgis application without the graphical user interface
    gui_flag = False
    app = QgsApplication([], gui_flag)
    app.initQgis()

    sys.path.append('/usr/share/qgis/python/plugins/')

    import processing

    from processing.core.Processing import Processing
    Processing.initialize()


#############################################################################
#########       DECLARE SOME USEFUL PARAMETERS         ######################
#############################################################################
# LCZ definition by zone type
lcz_type = {"Urban": [1,2,3,4,5,6,7,8,9,105],
            "Rural": [101,102,103,104,106,107]}

# Names of the region containing european countries
country_regions = ['Western Europe', 'Northern Europe', 'Southern Europe',
                   'Eastern Europe']

# Size of one cell in °²
area_cell = 0.2 * 0.2 

# Number of cells by country
ncells_per_country = 25

# Minimum number of 1 km² pixels of a certain zone type to consider a cell
min_nb_pix = 5

# Extension name
extension_geo = ".fgb"
extension_txt = ".csv"

# File names
path_to_indic_data = "/cnrm/ville/DATA/OSM_RAW/osm_{0}_{1}_{2}_{3}/rsu_indicators"+extension_geo
path_to_lcz_data = "/cnrm/ville/DATA/OSM_RAW/osm_{0}_{1}_{2}_{3}/rsu_lcz"+extension_geo
file_name_indic = "rsu_indicators"  # Base name of the file containing RSU indicators
file_name_cells_calculated = "/cnrm/ville/USERS/bernardj/Data/GeoClimateEurope/LCZ/cells_calculated.geojson"
file_name_countries = "/cnrm/ville/USERS/Data/world-administrative-boundaries.geojson"
file_smod = {"Urban":"/cnrm/ville/USERS/bernardj/Data/GeoClimateEurope/GHS_SMOD/Only_Urban30.tif",
             "Rural":"/cnrm/ville/USERS/bernardj/Data/GeoClimateEurope/GHS_SMOD/Only_Rural11-13.tif"}
file_indic_out = {"Rural": "/cnrm/ville/USERS/bernardj/Data/GeoClimateEurope/Indics/indics_rural{0}",
                  "Urban": "/cnrm/ville/USERS/bernardj/Data/GeoClimateEurope/Indics/indics_urban{0}"}
file_boxplot = "/cnrm/ville/USERS/bernardj/Figures/GeoClimateEurope/DEODE_Europe/{0}_{1}.png"

###############################################################
#########       START THE SCRIPT         ######################
###############################################################
# Create a dictionary with an entry for each layer of indic (urban and rural)
vlayer_indic = {}
layer = {}

if recalc or replot:
    # Load file containing already calculated cells
    gdf_cells = gpd.read_file(file_name_cells_calculated)
    nb_rows = gdf_cells.index.size
    
    # Load the file containing World geometries and filter the European countries
    gdf_countries = gpd.read_file(file_name_countries)
    gdf_eu_countries = gdf_countries[gdf_countries.region.isin(country_regions)]
    country_ids = gdf_eu_countries.index

if recalc:
    # For each zone type (urban or rural)
    for zt in file_smod.keys():
        # Calculate cells that intersect Urban centre grid cells (GHS SMOD value 30)
        path_nb_zt_cells = processing.run("native:zonalstatisticsfb", 
                                          {'INPUT':file_name_cells_calculated,
                                           'INPUT_RASTER':file_smod[zt],
                                           'RASTER_BAND':1,'COLUMN_PREFIX':'_','STATISTICS':[0],
                                           'OUTPUT':"/tmp/tempo_file.fgb"})['OUTPUT']
        
        # Filter cells having a minimum amount of a given zone type cell
        gdf_cells_zt = gpd.read_file(path_nb_zt_cells)
        gdf_cells_zt = gdf_cells_zt[gdf_cells_zt["_count"] >= min_nb_pix]
        
        # Create a spatial index for the cell geometries
        if gdf_cells_zt.sindex is None:
            cells_zt_sindex = gdf_cells_zt.sindex
        else:
            cells_zt_sindex = gdf_cells_zt.sindex
                
        new_line = {}
        # Sample cells in each European country
        for j, c_id in enumerate(country_ids):
            print("Country {0}: {1}".format(j, c_id))
            country_polygon = gdf_eu_countries.loc[c_id, "geometry"]
            # Intersect cells with country geometry using the spatial index
            possible_matches_index = list(cells_zt_sindex.intersection(country_polygon.bounds))
            # Filter cells data based on the intersection
            country_cells = gdf_cells_zt.iloc[possible_matches_index]
            
    
            
            # Randomly choose cells within the country
            if ncells_per_country <= country_cells.index.size:
                ncells_country_j = ncells_per_country
            else:
                ncells_country_j = country_cells.index.size
            country_cells_sample = country_cells.sample(ncells_country_j)
            for i in country_cells_sample.index:
                print(i)
                # When a coordinate is an integer, the file is saved without the decimal --> need special treatment
                coord = pd.Series(country_cells.loc[i, "geometry"].bounds).to_dict()
                for k, nb in coord.items():
                    if nb % 1 == 0.0:
                        coord[k] = int(coord[k])
                
                # Select in the cell only the RSU that are of the zone type of interest
                gdf_rsu_lcz = gpd.read_file(path_to_lcz_data.format(coord[1], coord[0], coord[3], coord[2])).set_index("ID_RSU")
                rsu_lcz_zt = gdf_rsu_lcz[gdf_rsu_lcz["LCZ_PRIMARY"].isin(lcz_type[zt])].index
                
                # Calculate indicators only for RSU being the right zone type
                gdf_rsu_indic = gpd.read_file(path_to_indic_data.format(coord[1], coord[0], coord[3], coord[2])).to_crs(3857)
                gdf_rsu_indic = gdf_rsu_indic.reindex(rsu_lcz_zt)
                new_line[j * ncells_per_country + i] = \
                    pd.Series({"Country": gdf_eu_countries.loc[c_id,"name"],
                               "Undefined fraction": (gdf_rsu_indic.UNDEFINED_FRACTION*gdf_rsu_indic.area).sum()/(gdf_rsu_indic.area.sum()),
                               "Urban fraction": ((gdf_rsu_indic.BUILDING_FRACTION+gdf_rsu_indic.IMPERVIOUS_FRACTION+gdf_rsu_indic.ROAD_FRACTION\
                                                   +gdf_rsu_indic.HIGH_VEGETATION_BUILDING_FRACTION+gdf_rsu_indic.HIGH_VEGETATION_IMPERVIOUS_FRACTION+gdf_rsu_indic.HIGH_VEGETATION_ROAD_FRACTION)\
                                                  *gdf_rsu_indic.area).sum()/(gdf_rsu_indic.area.sum()),
                               "Individual housing area fraction": (gdf_rsu_indic.AREA_FRACTION_INDIVIDUAL_HOUSING*gdf_rsu_indic.area*\
                                                                 (gdf_rsu_indic.HIGH_VEGETATION_BUILDING_FRACTION+gdf_rsu_indic.BUILDING_FRACTION)).sum()/
                                                                   ((gdf_rsu_indic.area*(gdf_rsu_indic.HIGH_VEGETATION_BUILDING_FRACTION+gdf_rsu_indic.BUILDING_FRACTION)).sum()),
                               "Collective housing area fraction": (gdf_rsu_indic.AREA_FRACTION_COLLECTIVE_HOUSING*gdf_rsu_indic.area*\
                                                                 (gdf_rsu_indic.HIGH_VEGETATION_BUILDING_FRACTION+gdf_rsu_indic.BUILDING_FRACTION)).sum()/
                                                                   ((gdf_rsu_indic.area*(gdf_rsu_indic.HIGH_VEGETATION_BUILDING_FRACTION+gdf_rsu_indic.BUILDING_FRACTION)).sum()),
                               "Commercial area fraction": (gdf_rsu_indic.AREA_FRACTION_COMMERCIAL*gdf_rsu_indic.area*\
                                                                 (gdf_rsu_indic.HIGH_VEGETATION_BUILDING_FRACTION+gdf_rsu_indic.BUILDING_FRACTION)).sum()/
                                                                   ((gdf_rsu_indic.area*(gdf_rsu_indic.HIGH_VEGETATION_BUILDING_FRACTION+gdf_rsu_indic.BUILDING_FRACTION)).sum()),
                               "Tertiary area fraction": (gdf_rsu_indic.AREA_FRACTION_TERTIARY*gdf_rsu_indic.area*\
                                                                 (gdf_rsu_indic.HIGH_VEGETATION_BUILDING_FRACTION+gdf_rsu_indic.BUILDING_FRACTION)).sum()/
                                                                   ((gdf_rsu_indic.area*(gdf_rsu_indic.HIGH_VEGETATION_BUILDING_FRACTION+gdf_rsu_indic.BUILDING_FRACTION)).sum()),
                               "Education area fraction": (gdf_rsu_indic.AREA_FRACTION_EDUCATION*gdf_rsu_indic.area*\
                                                                 (gdf_rsu_indic.HIGH_VEGETATION_BUILDING_FRACTION+gdf_rsu_indic.BUILDING_FRACTION)).sum()/
                                                                   ((gdf_rsu_indic.area*(gdf_rsu_indic.HIGH_VEGETATION_BUILDING_FRACTION+gdf_rsu_indic.BUILDING_FRACTION)).sum()),
                               "Light industrial area fraction": (gdf_rsu_indic.AREA_FRACTION_LIGHT_INDUSTRIAL*gdf_rsu_indic.area*\
                                                                 (gdf_rsu_indic.HIGH_VEGETATION_BUILDING_FRACTION+gdf_rsu_indic.BUILDING_FRACTION)).sum()/
                                                                   ((gdf_rsu_indic.area*(gdf_rsu_indic.HIGH_VEGETATION_BUILDING_FRACTION+gdf_rsu_indic.BUILDING_FRACTION)).sum()),
                               "Heavy industrial area fraction": (gdf_rsu_indic.AREA_FRACTION_HEAVY_INDUSTRIAL*gdf_rsu_indic.area*\
                                                                 (gdf_rsu_indic.HIGH_VEGETATION_BUILDING_FRACTION+gdf_rsu_indic.BUILDING_FRACTION)).sum()/
                                                                   ((gdf_rsu_indic.area*(gdf_rsu_indic.HIGH_VEGETATION_BUILDING_FRACTION+gdf_rsu_indic.BUILDING_FRACTION)).sum()),
                               "Non heated area fraction": (gdf_rsu_indic.AREA_FRACTION_NON_HEATED*gdf_rsu_indic.area*\
                                                                 (gdf_rsu_indic.HIGH_VEGETATION_BUILDING_FRACTION+gdf_rsu_indic.BUILDING_FRACTION)).sum()/
                                                                   ((gdf_rsu_indic.area*(gdf_rsu_indic.HIGH_VEGETATION_BUILDING_FRACTION+gdf_rsu_indic.BUILDING_FRACTION)).sum()),
                               "Undefined area fraction": (gdf_rsu_indic.AREA_FRACTION_UNDEFINED*gdf_rsu_indic.area*\
                                                                 (gdf_rsu_indic.HIGH_VEGETATION_BUILDING_FRACTION+gdf_rsu_indic.BUILDING_FRACTION)).sum()/
                                                                   ((gdf_rsu_indic.area*(gdf_rsu_indic.HIGH_VEGETATION_BUILDING_FRACTION+gdf_rsu_indic.BUILDING_FRACTION)).sum()),   
                               "Weighted average building height": (gdf_rsu_indic.AVG_HEIGHT_ROOF_AREA_WEIGHTED*gdf_rsu_indic.area*\
                                                                 (gdf_rsu_indic.HIGH_VEGETATION_BUILDING_FRACTION+gdf_rsu_indic.BUILDING_FRACTION)).sum()/
                                                                   ((gdf_rsu_indic.area*(gdf_rsu_indic.HIGH_VEGETATION_BUILDING_FRACTION+gdf_rsu_indic.BUILDING_FRACTION)).sum()),                                                                  
                               "SVF": (gdf_rsu_indic.GROUND_SKY_VIEW_FACTOR*gdf_rsu_indic.area).sum()/(gdf_rsu_indic.area.sum()),
                               "High vegetation fraction": ((gdf_rsu_indic.HIGH_VEGETATION_FRACTION+gdf_rsu_indic.HIGH_VEGETATION_LOW_VEGETATION_FRACTION\
                                                             +gdf_rsu_indic.HIGH_VEGETATION_BUILDING_FRACTION+gdf_rsu_indic.HIGH_VEGETATION_IMPERVIOUS_FRACTION\
                                                                 +gdf_rsu_indic.HIGH_VEGETATION_ROAD_FRACTION+gdf_rsu_indic.HIGH_VEGETATION_WATER_FRACTION\
                                                                     +gdf_rsu_indic.HIGH_VEGETATION_RAIL_FRACTION)*gdf_rsu_indic.area).sum()/(gdf_rsu_indic.area.sum()),
                               "Low vegetation fraction": ((gdf_rsu_indic.LOW_VEGETATION_FRACTION+gdf_rsu_indic.HIGH_VEGETATION_LOW_VEGETATION_FRACTION)\
                                                           *gdf_rsu_indic.area).sum()/(gdf_rsu_indic.area.sum())})
    
        # Save the values in a csv file
        df = pd.concat(new_line, axis = 1).transpose()
        df.to_csv(file_indic_out[zt].format(extension_txt))
        
        # Save the median value and create a geographic file
        gdf_indic = gpd.GeoDataFrame(pd.concat([df.groupby("Country").median(), 
                                                pd.Series(df.groupby("Country")["SVF"].count(), name = "Count"),
                                                gdf_eu_countries.set_index("name").geometry],
                                               axis = 1).reset_index())
        indicators = df.columns.difference(pd.Index(["Country"])).tolist()
        gdf_indic[indicators] = gdf_indic[indicators].astype(float)
        gdf_indic.to_file(file_indic_out[zt].format(extension_geo))
    
if replot:
    for zt in file_smod.keys():
        df = pd.read_csv(file_indic_out[zt].format(extension_txt), header = 0, index_col = 0)
        # Get the list of indicators
        indicators = df.columns.difference(pd.Index(["Country"])).tolist()
        # Remove the non member states
        member_states = gdf_eu_countries[gdf_eu_countries.status == "Member State"].name.values
        df = df[df.Country.isin(member_states)]
        for ind in indicators:
            df2 = pd.concat([df[df.Country.isin([i])][ind].rename(i) for i in df.Country.unique()], axis = 1)
            # Sort the boxplots by decreasing order
            idx_sorted = df2.median().sort_values(ascending = False).index
            # Create a boxplot
            fig, ax = plt.subplots(figsize = (30, 20))
            df2.reindex(idx_sorted, axis = 1).boxplot(rot = 30, ax = ax)
            ax.tick_params(axis = "x", labelsize = 10)
            plt.setp(ax.xaxis.get_majorticklabels(), ha='right')
            ax.set_ylabel(ind, fontsize = 20)
            fig.savefig(file_boxplot.format(zt, ind))
        plt.close("all")
        
        
if reload:
    for zt in file_smod.keys():
        # Load the file directly into QGIS    
        vlayer_indic[zt] = QgsVectorLayer(file_indic_out[zt].format(extension_geo),
                                          "{0} indicators".format(zt),
                                          "ogr")
        vlayer_indic[zt] = iface.addVectorLayer(file_indic_out[zt].format(extension_geo),
                                        "{0} indicators".format(zt),
                                        "ogr")                                  
        """if not vlayer_indic[zt].isValid():
            print(f"Layer indicators failed to load!")
        else:
            QgsProject.instance().addMapLayer(vlayer_indic[zt])
            layer[zt] = iface.activeLayer()"""
"""

# Continue sampling within the calculated data as long as the amount of cell within each country is not met

    

# Define the size (in °) of the zone used for the calculation
x_res = 0.2
y_res = 0.2

# Get the centroid of the current map extent
current_map_extent = iface.mapCanvas().extent()
xcenter = (current_map_extent.xMaximum()+current_map_extent.xMinimum())/2
if xcenter > 0:
    xcenter = int(xcenter/x_res)*x_res
else:
    xcenter = int(xcenter/x_res)*x_res-x_res
ycenter = (current_map_extent.yMaximum()+current_map_extent.yMinimum())/2
if ycenter > 0:
    ycenter = int(ycenter/y_res)*y_res
else:
    ycenter = int(ycenter/y_res)*y_res-y_res

coord = {}
coord[0] = round(ycenter, 1)
coord[1] = round(xcenter, 1)

coord[2] = round(coord[0]+y_res, 1)
coord[3] = round(coord[1]+x_res, 1)

for i, nb in coord.items():
   if nb % 1 == 0:
       coord[i] = int(coord[i])

cells_to_map = [list(coord.values())]


url_lcz_style = "https://ncloud.zaclys.com/index.php/s/LSpwAqMG2MFNYFr/download/rsu_lcz_primary.qml"
filename_style = "/home/bernardj/Data/LCZ/rsu_lcz_primary.qml"
if not os.path.exists(filename_style):
    urlretrieve(url_lcz_style, filename_style)"""