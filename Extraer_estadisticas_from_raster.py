# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 12:03:56 2023

@author: Administrador
"""
Dir_shape = r'D:/Productos_PGA/DistritosC.shp'
Dir_raster1 = 'E:\Tesis_Briant_Ultima_versión\SM_Pachatusan_s2.tif'
Dir_raster2 = 'E:\Tesis_Briant_Ultima_versión\SM_Tambomachay_s.tif'
Dir_raster3 = 'E:\Tesis_Briant_Ultima_versión\SM_Qoricocha_s.tif'
rasters = [Dir_raster1, Dir_raster2, Dir_raster3]

import rasterio
from rasterio.features import geometry_mask
from shapely.geometry import shape
import fiona
import numpy as np
import matplotlib.pyplot as plt
import os

def graficar_valores(poly, rasters, feature):
    fig = plt.figure(figsize=(8,7))
    name = feature['properties']['NM_DIST']
    
    for raster in rasters:
        transform = raster.transform
        img = raster.read(1)
        out_shape = (raster.height, raster.width)

        mask = geometry_mask([poly], transform=transform, out_shape=out_shape, invert=True)

        img_mask = img[mask]
        
        if np.isnan(img_mask).all():
            print("No hay valores para graficar en esta máscara "+name)
            return
        if isinstance(raster.nodata, float):
            if (img_mask < raster.nodata).any() or (img_mask > raster.nodata).any():
                print("Los valores están fuera del rango del raster para "+name)
                return
                
        with rasterio.open(Dir_raster1) as raster1:
            mask1 = geometry_mask([poly], transform=raster1.transform, out_shape=(raster1.height, raster1.width), invert=True)
            img_mask1 = raster1.read(1)[mask1]
            min_val1 = img_mask1.min()
            mean_val1 = img_mask1.mean()
            max_val1 = img_mask1.max()

        with rasterio.open(Dir_raster2) as raster2:
            mask2 = geometry_mask([poly], transform=raster2.transform, out_shape=(raster2.height, raster2.width), invert=True)
            img_mask2 = raster2.read(1)[mask2]
            min_val2 = img_mask2.min()
            mean_val2 = img_mask2.mean()
            max_val2 = img_mask2.max()

        with rasterio.open(Dir_raster3) as raster3:
            mask3 = geometry_mask([poly], transform=raster3.transform, out_shape=(raster3.height, raster3.width), invert=True)
            img_mask3 = raster3.read(1)[mask3]
            min_val3 = img_mask3.min()
            mean_val3 = img_mask3.mean()
            max_val3 = img_mask3.max()
            

        plt.bar(x=[1, 2, 3], height=[min_val1, mean_val1, max_val1], width=0.2, color='r', label='raster1')
        plt.bar(x=[1.2, 2.2, 3.2], height=[min_val2, mean_val2, max_val2], width=0.2, color='g', label='raster2')
        plt.bar(x=[1.4, 2.4, 3.4], height=[min_val3, mean_val3, max_val3], width=0.2, color='b', label='raster3')
        plt.xticks([1.2, 2.2, 3.2], ['Mínimo', 'Media', 'Máximo'])
        for i, v in enumerate([round(min_val1, 2), round(mean_val1,2), round(max_val1,2)]):
            plt.text(i+0.932, v + 0.01, str(v), color='gray', fontweight='bold', size=10)
        for i, v in enumerate([round(min_val2, 2), round(mean_val2,2), round(max_val2,2)]):
            plt.text(i+1.12, v + 0.01, str(v), color='gray', fontweight='bold', size=10)  
        for i, v in enumerate([round(min_val3, 2), round(mean_val3,2), round(max_val3,2)]):
            plt.text(i+1.32, v + 0.01, str(v), color='gray', fontweight='bold', size=10)   
            
    plt.ylabel('PGA (g)', font='Gabriola', fontsize=20)
    plt.xlabel('Valores', font='Gabriola', fontsize=20)
    plt.title(f'Valores PGA para ({name})', fontsize=28, font='Gabriola', fontweight='bold')
    plt.axhline(y = 0.45, color='red', ls='-.') #Norma tecnica max
    plt.axhline(y = 0.35, color='orange', ls='-.') #Norma tecnica med
    plt.axhline(y = 0.25, color='yellow', ls='-.') #Norma tecnica min
    plt.axhline(y = 0.98, color='gray', ls='-') #Norma tecnica min
    plt.legend(['PGA Alto - N.T.E030','PGA Medio - N.T.E030', 'PGA Bajo - N.T.E030', '0.98 (g)-Gravedad Estándar', 'Fuente sísmica:F. Pachatusan', 'Fuente sísmica:F. Tambomachay', 'Fuente sísmica:F. Qoricocha'], loc='center', bbox_to_anchor=(0.82, 0.122), prop={'size': 8})
    # plt.legend(['PGA alto', 'Falla Tambomachay', 'Falla Qoricocha'], loc='center', bbox_to_anchor=(0.15, 0.93), prop={'size': 8.5})
   
    if not os.path.exists("Figuras"):
       os.mkdir("Figuras")
    fig.savefig(f"Figuras/{name}.png")
    plt.show()

with fiona.open(Dir_shape) as src:
    rasters = [rasterio.open(Dir_raster1), rasterio.open(Dir_raster2), rasterio.open(Dir_raster3)]
    for feature in src:
        poly = shape(feature['geometry'])
        graficar_valores(poly, rasters, feature)