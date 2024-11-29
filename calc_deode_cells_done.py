# -*- coding: utf-8 -*-
"""
Script used to determines how much cells have been calculated on the domain and
crates a GeoJson containing the cells that have been calculated.

@author: Jérémy Bernard, Météo-France
#date: 2023-09-15
"""

import glob
import numpy as np
import os
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon

# Input and output paths used
path_lcz = "/mnt/lfs/d40/ville/DATA/OSM_RAW/osm_*"
path_cells_done = "/cnrm/ville/USERS/bernardj/Data/GeoClimateEurope/LCZ/cells_calculated.fgb"

# Informations about the GeoClimate input data
epsg = 4326
grid_size = 0.2

# Get the list of folders (each one containing the coordinate of the cell)
list_folders = glob.glob(path_lcz)

# Convert the folder name into coordinates and calculate min and max lat and long of all cells
df_coord = pd.DataFrame({"lat": np.array([i.split(os.sep)[-1].split("_")[1] for i in list_folders]).astype(float),
                         "long": np.array([i.split(os.sep)[-1].split("_")[2] for i in list_folders]).astype(float)})
min_lat = df_coord["lat"].min()
max_lat = df_coord["lat"].max() + grid_size
min_long = df_coord["long"].min()
max_long = df_coord["long"].max() + grid_size

# Calculate and print the number of current calculated cells, expected number and % of cells calculated
nb_expected_cells = (max_lat - min_lat) * (max_long - min_long) / grid_size ** 2
nb_current_cells = len(list_folders)
print(f""""There is currently {int(nb_current_cells)} cells calculated against \
      {int(nb_expected_cells)} normally ({round(100*nb_current_cells/nb_expected_cells)}%)""")
  
# Create a file containing the geometry of all calculated cells
gdf_all = gpd.GeoSeries([Polygon([(df_coord.loc[i,"long"], df_coord.loc[i,"lat"]),
                                  (df_coord.loc[i,"long"] + grid_size, df_coord.loc[i,"lat"]), 
                                  (df_coord.loc[i,"long"] + grid_size, df_coord.loc[i,"lat"] + grid_size), 
                                  (df_coord.loc[i,"long"], df_coord.loc[i,"lat"] + grid_size),
                                  (df_coord.loc[i,"long"], df_coord.loc[i,"lat"])]) for i in df_coord.index]).set_crs(epsg)
gdf_all.to_file(path_cells_done)

vlayer = iface.addVectorLayer(path_cells_done, "DEODE European map", "ogr")