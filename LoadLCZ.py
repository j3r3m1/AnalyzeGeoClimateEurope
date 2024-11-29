"""
Script used to load LCZ or indicator maps from the European map

@author: Jérémy Bernard, Météo-France
#date: 2023-09-15
"""

from urllib.request import urlretrieve
import os

# Define which file you want to load
rsu_indic_bool = True

file_name_indic = "rsu_indicators"
file_name_lcz = "rsu_lcz"

# Extension name
extension = ".fgb"

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

path_to_data = "/cnrm/ville/DATA/OSM_RAW/osm_"
url_lcz_style = "https://ncloud.zaclys.com/index.php/s/LSpwAqMG2MFNYFr/download/rsu_lcz_primary.qml"
filename_style = "/cnrm/ville/USERS/bernardj/Data/GeoClimateEurope/LCZ/rsu_lcz_primary.qml"
if not os.path.exists(filename_style):
    urlretrieve(url_lcz_style, filename_style)

# Remove old layer
while len(QgsProject.instance().mapLayers().values()) > 2:
    QgsProject.instance().removeMapLayer(list(QgsProject.instance().mapLayers().keys())[0])

# Load vector layers
for bb in cells_to_map:
    print(bb)
    if rsu_indic_bool:
        filename_path_indic = path_to_data + str(bb)[1:-1].replace(",","_").replace(" ", "") + os.sep + file_name_indic + extension
        vlayer_indic = QgsVectorLayer(filename_path_indic, str(bb)[1:-1] + "_rsu_indic", "ogr")
        if not vlayer_indic.isValid():
            print(f"Layer {str(bb)[1:-1]} rsu indicators failed to load!")
        else:
            QgsProject.instance().addMapLayer(vlayer_indic)
            layer = iface.activeLayer()
    filename_path_lcz = path_to_data + str(bb)[1:-1].replace(",","_").replace(" ", "") + os.sep + file_name_lcz + extension
    vlayer_lcz = QgsVectorLayer(filename_path_lcz, str(bb)[1:-1] + "_rsu_lcz", "ogr")
    if not vlayer_lcz.isValid():
        print(f"Layer {str(bb)[1:-1]} LCZ failed to load!")
    else:
        QgsProject.instance().addMapLayer(vlayer_lcz)
        layer = iface.activeLayer()
        # Load LCZ style
        layer.loadNamedStyle(filename_style)
        layer.triggerRepaint()
        # Then zoom on the new tile
        extent = layer.extent()
        canvas = iface.mapCanvas()
        canvas.setExtent(extent)
        canvas.refresh()