# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 16:36:58 2023

@author: Administrador
"""

import pysheds
from pysheds.grid import Grid
import rasterio
import fiona
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np

Dir = 'D:/Curso/DEMD.tif' #Ingresa directorio de raster
Dir_guardar = r'D:/Curso/drenajes.shp' #Ingresa el directorio donde se guardará el shape

"Valores promedio para DEM de drone / resolución menor a 20 cm/px"
Area = 2000000 #Con este valor genera drenajes en un area de 800*550 m
Densidad = 15000 #Con este valor se pueden ver rios principales y tributarios en el area. A menor valor más densidad

with rasterio.open(Dir) as src:
    dem = src.read(1)
    
# Crea un objeto "grid" a partir de tu raster
grid = Grid.from_raster(Dir, data_dem='dem')
dem=grid.read_raster(Dir, data_dem='dem')

dirmap = (64,  128,  1,   2,    4,   8,    16,  32)

xm = xmax=src.bounds.right
ym = ymax=src.bounds.top

x, y = xm, ym

pit_filled_dem = grid.fill_pits(dem)

# Fill depressions in DEM
flooded_dem = grid.fill_depressions(pit_filled_dem)
    
# Resolve flats in DEM
inflated_dem = grid.resolve_flats(flooded_dem)

# Genera la dirección de flujo para cada celda en el grid
fdir = grid.flowdir(inflated_dem, dirmap=dirmap)

acc = grid.accumulation(fdir, dirmap=dirmap)

x_snap, y_snap = grid.snap_to_mask(acc > Area, (x, y))

branches = grid.extract_river_network(fdir, acc > Densidad, dirmap=dirmap)

sns.set_palette('flare')
fig, ax = plt.subplots(figsize=(8.5,6.5))

plt.xlim(grid.bbox[0], grid.bbox[2])
plt.ylim(grid.bbox[1], grid.bbox[3])
ax.set_aspect('equal')

for branch in branches['features']:
    line = np.asarray(branch['geometry']['coordinates'])
    plt.plot(line[:, 0], line[:, 1])
    
_ = plt.title('Drenajes', size=14)

# Define the Shapefile schema
schema = {
    'geometry': 'LineString',
    'properties': {
        'id': 'int'
    },
}

# Open the Shapefile for writing
with fiona.open(Dir_guardar, 'w', 'ESRI Shapefile', schema) as dst:
    for i, branch in enumerate(branches['features']):
        line = branch['geometry']['coordinates']
        feature = {
            'geometry': {
                'type': 'LineString',
                'coordinates': line
            },
            'properties': {
                'id': i
            },
        }
        dst.write(feature)

#%%

# # Si deseas plotear la acumulación de flujo puedes activar esta celda

fig, ax = plt.subplots(figsize=(32,24), dpi=300)
fig.patch.set_alpha(0)
plt.grid('on', zorder=0)
im = ax.imshow(acc, extent=grid.extent, zorder=2,
                cmap='cubehelix',
                norm=colors.LogNorm(1, acc.max()),
                interpolation='bilinear')
plt.colorbar(im, ax=ax, label='Upstream Cells')
plt.title('Flow Accumulation', size=14)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.tight_layout()

#%%

# calculando zonas de maxima inundación, falta trabajar
import numpy as np
from osgeo import gdal
from osgeo import ogr
from shapely.geometry import mapping

# Abrir el archivo raster
ds = gdal.Open("D:/Curso/DEMD.tif")

# Leer los datos del raster
dem = ds.ReadAsArray()

# Calcular la pendiente
slope = np.gradient(dem)

# Calcular la dirección de la pendiente
angle = np.arctan2(*slope)

# Calcular el umbral
threshold = 0.1 * np.max(np.abs(slope))

# Inicializar una lista para guardar los polígonos
polygons = []

# Recorrer cada píxel del raster
for x in range(dem.shape[1]):
    for y in range(dem.shape[0]):
        # Verificar si la pendiente supera el umbral
        if np.abs(slope[0][y, x]) > threshold or np.abs(slope[1][y, x]) > threshold:
            # Guardar el polígono en la lista
            polygon = [(x, y), (x + 5, y), (x + 5, y + 5), (x, y + 5), (x, y)]
            polygons.append(polygon)

# Crear un shapefile para guardar los polígonos
driver = ogr.GetDriverByName("ESRI Shapefile")
ds = driver.CreateDataSource("flood_zones.shp")
layer = ds.CreateLayer("flood_zones", geom_type=ogr.wkbPolygon)

# Crear una nueva columna para almacenar los polígonos
field_defn = ogr.FieldDefn("ID", ogr.OFTInteger)
layer.CreateField(field_defn)

# Añadir los polígonos al shapefile
for i, polygon in enumerate(polygons):
    # Crear una nueva feature
    feature = ogr.Feature(layer.GetLayerDefn())
    
    # Convertir el polígono a un objeto de tipo Polygon
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for point in polygon:
        ring.AddPoint(point[0], point[1])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    
    # Asignar el polígono a la feature
    feature.SetGeometry(poly)
    
    # Asignar un valor a la columna "ID"
    feature.SetField("ID", i)
    
    # Añadir la feature al shapefile
    layer.CreateFeature(feature)

schema = {
    'geometry': 'Polygon',
    'properties': {
        'id': 'int'
    },
}

with fiona.open('flood_zone.shp', 'w', 'ESRI Shapefile', schema) as c:
    c.write({
        'geometry': mapping(polygon),
        'properties': {
            'id': 1
        },
    })
    
# Cerrar el shapefile
ds = None




