# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 21:12:55 2019
@author: INGEMMET - GRUPO DE NEOTECTÓNICA
"""
import numpy as np
import fiona
import pyproj
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
from scipy import stats
from shapely.geometry import Point, LineString, Polygon, mapping
from shapely.ops import linemerge, unary_union, polygonize
from descartes import PolygonPatch
from AS08_const_coef import *
from osgeo import gdal
from osgeo import ogr
from osgeo import gdalconst
import multiprocessing as mp
import rasterio
import os.path as path
#%%
def plot_coords(ax, ob,color='#999999'):
    x, y = ob.xy
    ax.plot(x, y, 'o', color=color, zorder=1)
            
def plot_line(ax, ob, color='black',linewidth=2):
    x, y = ob.xy
    ax.plot(x, y, alpha=1.0, linewidth=linewidth, solid_capstyle='round', zorder=2, color=color)
def distanceWGS84(puntoA,puntoB): #Calcula distancias en el sistema WGS84 - "Geográficas"
    long1,lat1 = puntoA# Y,X (-15.996333373904548,-76.81000486854916)
    long2,lat2 = puntoB# Y,X (-13.701694863680594,-13.701694863680594)
    geod = pyproj.Geod(ellps="WGS84")
    angle1,angle2,distance = geod.inv(long1, lat1, long2, lat2)
    return distance #is in meters

def inPuntos(Linea,nsamples): #modificado para tener los segmentos en partes iguales
    puntos = []
    cellsize = Linea.length/nsamples
    ptosInLine = np.arange(0,Linea.length,cellsize) ####
    for punto in ptosInLine:
        x,y = Linea.interpolate(punto).coords[0]
        puntos.append((x,y))
    x,y = Linea.interpolate(Linea.length).coords[0]
    puntos.append((x,y))
    
    puntos2 = []
    ptosInLine2 = np.arange(cellsize/2,Linea.length,cellsize) ####
    for punto in ptosInLine2:
        x,y = Linea.interpolate(punto).coords[0]
        puntos2.append((x,y))
    
    return puntos,puntos2

def cIntervalo(lBasePerfil,muestreo):
    Dm = 0
    for i in range(0,len(lBasePerfil.coords)-1):
        Dm = Dm + distanceWGS84(lBasePerfil.coords[i],lBasePerfil.coords[i+1])
#    Dm = distanceWGS84(lBasePerfil.coords[0],lBasePerfil.coords[1]) #Longitud en metros
    Das = lBasePerfil.length #Longitud en Arcos de Segundo
    return muestreo*(Das / Dm)

def cut_polygon_by_line(polygon, line):
    merged = linemerge([polygon.boundary, line])
    borders = unary_union(merged)
    polygons = polygonize(borders)
    return list(polygons)
def getboxesPos(segmento,l):
    box= segmento.buffer(l,cap_style=2)
    boxes = cut_polygon_by_line(box, segmento)
    box1 = boxes[0]
    b1minx, b1miny, b1maxx, b1maxy = box1.bounds
    box2 = boxes[1]
    b2minx, b2miny, b2maxx, b2maxy = box2.bounds

    boxN = None
    boxS = None
    boxE = None
    boxW = None
    if b1miny != b2miny:
        if b1miny > b2miny:
            boxN = box1
            boxS = box2
        else:
            boxN = box2
            boxS = box1
    else:
        if b1minx > b2minx:
            boxE = box1
            boxW = box2
        else:
            boxE = box2
            boxW = box1
    print(boxN)
    print(boxS)
    print(boxE)
    print(boxW)
    return boxN,boxS,boxE,boxW
#%%
def getboxZones(segmento,l,dirBz=None):
#def setboxZones(segmento,dirBz=None,Bz=None,tipo=None,Mw=None,Hypo=None):# (Linea; Norte('N') o Sur('S'); Bz<90; 'NM','IN','SL')
    # Pendiente en el plano de la recta definida por el segmento
    try:
        m_seg, b_seg, r_value, p_value, std_err = stats.linregress(segmento.xy[0],segmento.xy[1])
        #m_seg,b_seg= np.polyfit(segmento.xy[0],segmento.xy[1],1)
        print(m_seg)
    except:
        pass
    print(segmento)
    boxNorte,boxSur,boxEste,boxOeste=getboxesPos(segmento,l)
    
    if dirBz == 'S':
        if np.isnan(m_seg):
            surfaceArea=boxEste# Caso 3
            dipDir='E'
        elif m_seg >= 0:
            surfaceArea=boxSur # Caso 4 y 5
            dipDir='SE'
            if m_seg == 0:
                dipDir='S'
        elif m_seg < 0:
            surfaceArea=boxSur # Caso 6
            dipDir='SO'
    elif dirBz == 'N':
        if np.isnan(m_seg):
            surfaceArea=boxOeste # Caso 7
            dipDir='O'
        elif m_seg >= 0:
            surfaceArea=boxNorte # Caso 8 y 1
            dipDir='NO'
            if m_seg == 0:
                dipDir='N'
        elif m_seg < 0:
            surfaceArea=boxNorte # Caso 2
            dipDir='NE'
        print(surfaceArea)
        print(dipDir)
        print(m_seg)
        print(b_seg)
        return surfaceArea,dipDir,m_seg,b_seg

def getGeometryFault(shape,
                     nroTrazo_shp=0,
                     UTM=True,
                     sam_profile=100,
                     dirBz=None, # Norte o Sur
                     Bz=None, # en grados
                     tipo=None, # NM, IN, *SL StrikeSlip falta implementar
                     Mw=None, # Magnitud en Mw
                     Hypo=None, # Hipocentro en Km
                     Ztor=None, # Profundidad del topa de la ruptura en Km
                     AS=False): 
    
    #eDEM,gteDEM = openRaster(rasterDEM)
    shp = fiona.open(shape)
    noTrazo= nroTrazo_shp
    if type(shp[0]['geometry']['coordinates'][0]) == tuple:
        lTrazoFalla=LineString(shp[noTrazo]['geometry']['coordinates'])
    elif type(shp[0]['geometry']['coordinates'][0]) == list:
        ltf = []
        for llin in shp[0]['geometry']['coordinates']:
            ltf = ltf + llin
        lTrazoFalla=LineString(ltf)
    else:
        print('El shapefile de la Falla esta mal!!!')
    #ccc = shp[noTrazo]['geometry']['coordinates']
    crs = shp.crs
    shp.close()
    
#### EL ANCHO INCLINADO(buzamiento) DEL RECTANGULO 
    W = 10.0**(CONS[tipo][Wa]+(CONS[tipo][Wb]*Mw))
#### EL ANCHO HORIZONTAL(en superficie) DEL RECTANGULO
    l=(W*np.cos(np.deg2rad(Bz)))*1000.0 #l es el ancho del rectangulo 
#### EL HIPOCENTRO
    if Hypo == None:
        Hypo = CONS[tipo][Ha]+CONS[tipo][Hb]*Mw
    else:
        Hypo = Hypo
#### PROFUNDIDAD DEL TOPE DE LA RUPTURA
    if Ztor == None:
        Ztor = max([(Hypo - 0.6 *W * np.sin(np.deg2rad(Bz))),0])
    else:
        Ztor = Ztor
    
    
    if UTM == False:
        cellsize = cIntervalo(lTrazoFalla,sam_profile)
    elif UTM == True:
        cellsize = sam_profile
    puntos,puntosh=inPuntos(lTrazoFalla,cellsize)
    lmTrazoFalla = LineString(puntos)
    lmTrazoFallah = LineString(puntosh)
    nsegmentos = len(lmTrazoFalla.coords)-1
    Segmentos=[]
    AreasRuptura = []
    DipDirs = []
    Ms_seg = []
    Bs_seg = []            
    for segmento in range(0,nsegmentos,1):
        if segmento <= (nsegmentos-1):
            lsegmento = LineString([lmTrazoFalla.coords[segmento],lmTrazoFalla.coords[segmento+1]])
            Segmentos.append(lsegmento)
            try:
                lsegmentoh = LineString([lmTrazoFallah.coords[segmento],lmTrazoFallah.coords[segmento+1]])
                Segmentos.append(lsegmentoh)
            except:
                pass
            # print(l)
            surfaceArea,dipDir,m_seg,b_seg = getboxZones(lsegmento,l,dirBz=dirBz)
            try:
                surfaceAreah,dipDirh,m_segh,b_segh = getboxZones(lsegmentoh,l,dirBz=dirBz)
            except:
                pass
            
            AreasRuptura.append(surfaceArea)
            try:
                AreasRuptura.append(surfaceAreah)
            except:
                pass
            
            DipDirs.append(dipDir)
            try:
                DipDirs.append(dipDirh)
            except:
                pass
            
            Ms_seg.append(m_seg)
            try:
                Ms_seg.append(m_segh)
            except:
                pass
            
            Bs_seg.append(b_seg)
            try:
                Bs_seg.append(b_segh)
            except:
                pass
            
    Ftipo  = tipo
    Fnor = 1*(tipo=='NM')
    Finv = 1*(tipo=='IN')
    Dip    = np.deg2rad(Bz)
    M      = Mw
    Zhipo  = Hypo
    Ztop   = Ztor
    W      = W
    FRV    = Finv
    FN     = Fnor
    AS     = 1*AS
    
    Falla={}
    Falla['tipo']=tipo
    Falla['dip']=Dip
    Falla['mw']=M
    
    Falla['zhipo']=Zhipo
    Falla['ztop']=Ztop
    Falla['w']=W
    
    Falla['oTrazo']=lTrazoFalla #linestring => Trazo Original
    Falla['sTrazo']=lmTrazoFalla #linestring => Trazo Simplificado
    Falla['segs']={}
    for segmento in range(0,len(Segmentos),1):
        Falla['segs'][segmento]={}
        Falla['segs'][segmento]['sfalla']=Segmentos[segmento] #linestring
        Falla['segs'][segmento]['aRupt']=AreasRuptura[segmento] #polygon
        Falla['segs'][segmento]['dipDir']=DipDirs[segmento] #string
        Falla['segs'][segmento]['m']=Ms_seg[segmento] #double
        Falla['segs'][segmento]['b']=Bs_seg[segmento] #double
    
    Falla['FRV']=FRV
    Falla['FN']=FN
#    Falla['AS']=AS
    Falla['crs']=shp.crs
    
#    return lTrazoFalla,lmTrazoFalla,Segmentos,AreasRuptura,DipDirs,Dip,Ztop,W,Ms_seg,Bs_seg,shp.crs
    return Falla
          #TrazoFalla,mTrazoFalla,Segmentos,AreasRuptura,DipDirs,Dip,Ztop,W,Ms_seg,Bs_seg,shpCRS
          
def getPointZona(lon,lat,Falla,segmento):
#def getPointZona(lon,lat,segmento,areaRuptura,dipDir,dip,ztop,w,m_seg,b_seg):
    # Retorna la zona y los parametros del punto con respecto a la falla
    # Se construye el punto
    Punto = Point(lon,lat)
    #segmento = Falla['segs'][segmento]['sfalla']
    areaRuptura= Falla['segs'][segmento]['aRupt']
    dipDir= Falla['segs'][segmento]['dipDir']
    dip= Falla['dip']
    ztop= Falla['ztop']
    w= Falla['w']
    m_seg= Falla['segs'][segmento]['m']
    b_seg= Falla['segs'][segmento]['b']
    
#### Se define el nombre y el orden de los vertices del area de Ruptura
    if m_seg == 0 or np.isnan(m_seg):
        minx, miny, maxx, maxy = areaRuptura.bounds
        a = Point(minx,maxy)
        b = Point(maxx,maxy)
        c = Point(maxx,miny)
        d = Point(minx,miny)
    else:
        df_box = pd.DataFrame()
        df_box['lon'],df_box['lat'] = areaRuptura.boundary.xy
        df_box= df_box.drop([4],axis=0)

        sortLon= df_box.sort_values('lon').reset_index(drop=True)
        minX = sortLon.loc[ 0 , : ]
        maxX = sortLon.loc[ 3 , : ]
        
        sortLat= df_box.sort_values('lat').reset_index(drop=True)
        minY = sortLat.loc[ 0 , : ]
        maxY = sortLat.loc[ 3 , : ]
        
        a = Point(minX.lon,minX.lat)            
        b = Point(maxY.lon,maxY.lat)
        c = Point(maxX.lon,maxX.lat)
        d = Point(minY.lon,minY.lat)
#### Construyendo los lados del rectangulo que define la superficie de ruptura
    ab = LineString([a, b])
    bc = LineString([b, c])
    cd = LineString([c, d])
    da = LineString([d, a])
    abcd = Polygon([[p.x, p.y] for p in [a,b,c,d]])
#### Cálculo de distancias del punto a la geometría de la falla
    dPuntoBox = Punto.distance(areaRuptura) #Distancia más corta entre el punto y el area de ruptura (vertice o segmento)
    dPa = Punto.distance(a) #Distancia al vertice a
    dPb = Punto.distance(b) #Distancia al vertice b
    dPc = Punto.distance(c) #Distancia al vertice c
    dPd = Punto.distance(d) #Distancia al vertice d
    dPab = Punto.distance(ab) #Distancia más corta entre el punto y el segmento ab
    dPbc = Punto.distance(bc) #Distancia más corta entre el punto y el segmento bc
    dPcd = Punto.distance(cd) #Distancia más corta entre el punto y el segmento cd
    dPda = Punto.distance(da) #Distancia más corta entre el punto y el segmento da
    
    # Se determina el vertice más cercano al punto
    df_disPP = pd.DataFrame()
    df_disPP['nS'],df_disPP['dP'] = ['dPa','dPb','dPc','dPd'],[dPa,dPb,dPc,dPd]
    sortdPP = df_disPP.sort_values('dP').reset_index(drop=True)
    minDisPP = sortdPP.loc[ 0 , : ]
    # Se determina el segmento más cercano al punto
    df_disPL = pd.DataFrame()
    df_disPL['nS'],df_disPL['dP'] = ['dPab','dPbc','dPcd','dPda'],[dPab,dPbc,dPcd,dPda]
    sortdPL = df_disPL.sort_values('dP').reset_index(drop=True)
    minDisPL = sortdPL.loc[ 0 , : ]
    
#### Se determina en que zona, respecto a la falla, se encuentra el punto  ####   
    if Punto.intersects(areaRuptura)==True:
        zz = 'E'
    else:       
        if dPuntoBox == minDisPP['dP']:
            if 'dPa' == minDisPP['nS']:
                zz = 'A'
            elif 'dPb' == minDisPP['nS']:
                zz = 'C'
            elif 'dPc' == minDisPP['nS']:
                zz = 'I'
            elif 'dPd' == minDisPP['nS']:
                zz = 'G'
        elif minDisPP['dP'] > dPuntoBox:
            if 'dPab' == minDisPL['nS'] :
                zz = 'B'
            elif 'dPbc' == minDisPL['nS'] :
                zz = 'F'
            elif 'dPcd' == minDisPL['nS']:
                zz = 'H'
            elif 'dPda' == minDisPL['nS']:
                zz = 'D'
    zona = ZONA[dipDir][zz]
    
#### Cálculo del ángulo alpha  ####     
    alpha = getAlpha(dipDir,zona,Punto,a,b,c,d,m_seg)
#### Cálculo distancia dRjb #### esta en m, para km dividir entre 1000
    if zona == 5:
        Rjb = 0
    else:
        Rjb = dPuntoBox/1000
#### Cálculo distancia dRx ####
    # y=m_seg*x+b_seg => ECUACION DE LA RECTA DEFINIDA POR EL SEGMENTO DE FALLA, /1000 para KM
    if zona == 1 or zona == 4 or zona == 7:
        fz = -1
    else:
        fz = 1
    Rx = abs((m_seg*Punto.x - Punto.y + b_seg)/(np.sqrt((m_seg**2.000) + 1.000)))*fz/1000
    
    #dF=(self.para_m*x-y+self.para_b)/(math.sqrt((self.para_m**2.000)+1.000))
    
#### Cálculo distancia dRy ####        
    if abs(alpha) == np.pi/2:
        Ry = 0
    elif abs(alpha) == 0 or abs(alpha) == np.pi:
        Ry = Rjb
    else:
        Ry = abs(Rx*(1/np.tan(alpha)))        
#### Cálculo distancia drRup' #### 
    if Rx < (ztop * np.tan(dip)):
        rRup_p = np.sqrt((Rx**2)+(ztop**2))
    elif ztop*np.tan(dip)<=Rx<= (ztop*np.tan(dip) + w*1/np.cos(dip)):
        rRup_p = Rx*np.sin(dip)+ztop*np.cos(dip)
    elif Rx > (ztop*np.tan(dip) + w*1/np.cos(dip)):
        rRup_p = np.sqrt((Rx-w*np.cos(dip))**2 + (ztop+w*np.sin(dip))**2)        
#### Cálculo distancia drRup #### 
    rRup = np.sqrt((rRup_p**2)+(Ry**2))
    return zona,alpha,Rjb,Rx,Ry,rRup_p,rRup,a,b,c,d,minDisPP,minDisPL,dip,w,ztop
    #return zona,alpha,Rjb,Rx,Ry,rRup_p,rRup,a,b,c,d
    
def getAlpha(dipDir,zona,punto,a,b,c,d,m_seg):#### Cálculo del ángulo alpha  ####
    
    angle_md = get_angle(m_seg,get_m(punto,d))
    angle_mc_p = np.pi - get_angle(m_seg,get_m(punto,c))
    angle_ma = get_angle(m_seg,get_m(punto,a))
    angle_md_p = np.pi - get_angle(m_seg,get_m(punto,d))
    angle_mb = get_angle(m_seg,get_m(punto,b))
    angle_ma_p = np.pi - get_angle(m_seg,get_m(punto,a))
    angle_mc = get_angle(m_seg,get_m(punto,c))
    angle_mb_p = np.pi - get_angle(m_seg,get_m(punto,b))
    
    mAlpha={}
    mAlpha['4N']=mAlpha['4NO']=mAlpha['4NE']=mAlpha['4E']=mAlpha['4SE']=mAlpha['4S']=mAlpha['4SO']=mAlpha['4O']=-np.pi/2 # -90
    mAlpha['5N']=mAlpha['5NO']=mAlpha['5NE']=mAlpha['5E']=mAlpha['5SE']=mAlpha['5S']=mAlpha['5SO']=mAlpha['5O']=np.pi/2 # 90
    mAlpha['6N']=mAlpha['6NO']=mAlpha['6NE']=mAlpha['6E']=mAlpha['6SE']=mAlpha['6S']=mAlpha['6SO']=mAlpha['6O']=np.pi/2 # 90
    
    mAlpha['1N']=mAlpha['1NO']= -angle_md # -(get_angle(m_seg,get_m(punto,d)))  
    mAlpha['2N']=mAlpha['2NO']= angle_md # get_angle(m_seg,get_m(punto,d))  
    mAlpha['3N']=mAlpha['3NO']= angle_md # get_angle(m_seg,get_m(punto,d))  
    
    mAlpha['7N']=mAlpha['7NO']= -angle_mc_p # -(np.pi - get_angle(m_seg,get_m(punto,c)))  
    mAlpha['8N']=mAlpha['8NO']= angle_mc_p # np.pi - get_angle(m_seg,get_m(punto,c))  
    mAlpha['9N']=mAlpha['9NO']= angle_mc_p # np.pi - get_angle(m_seg,get_m(punto,c))  
 
    mAlpha['1NE']=mAlpha['1E']= -angle_ma # -(get_angle(m_seg,get_m(punto,a)))  
    mAlpha['2NE']=mAlpha['2E']= angle_ma # get_angle(m_seg,get_m(punto,a))  
    mAlpha['3NE']=mAlpha['3E']= angle_ma # get_angle(m_seg,get_m(punto,a))  
    
    mAlpha['7NE']=mAlpha['7E']= -angle_md_p # -(np.pi - get_angle(m_seg,get_m(punto,d)))  
    mAlpha['8NE']=mAlpha['8E']= angle_md_p # np.pi - get_angle(m_seg,get_m(punto,d))  
    mAlpha['9NE']=mAlpha['9E']= angle_md_p # np.pi - get_angle(m_seg,get_m(punto,d))  
    
    mAlpha['1SE']=mAlpha['1S']= -angle_mb # -(get_angle(m_seg,get_m(punto,b)))  
    mAlpha['2SE']=mAlpha['2S']= angle_mb # get_angle(m_seg,get_m(punto,b))  
    mAlpha['3SE']=mAlpha['3S']= angle_mb # get_angle(m_seg,get_m(punto,b))  
    
    mAlpha['7SE']=mAlpha['7S']= -angle_ma_p # -(np.pi - get_angle(m_seg,get_m(punto,a)))  
    mAlpha['8SE']=mAlpha['8S']= angle_ma_p # np.pi - get_angle(m_seg,get_m(punto,a))  
    mAlpha['9SE']=mAlpha['9S']= angle_ma_p # np.pi - get_angle(m_seg,get_m(punto,a))  
    
    mAlpha['1SO']=mAlpha['1O']= -angle_mc # -(get_angle(m_seg,get_m(punto,c)))  
    mAlpha['2SO']=mAlpha['2O']= angle_mc # get_angle(m_seg,get_m(punto,c))  
    mAlpha['3SO']=mAlpha['3O']= angle_mc # get_angle(m_seg,get_m(punto,c))  
    
    mAlpha['7SO']=mAlpha['7O']= -angle_mb_p # -(np.pi - get_angle(m_seg,get_m(punto,b))) 
    mAlpha['8SO']=mAlpha['8O']= angle_mb_p # np.pi - get_angle(m_seg,get_m(punto,b))  
    mAlpha['9SO']=mAlpha['9O']= angle_mb_p # np.pi - get_angle(m_seg,get_m(punto,b))  
    return mAlpha[str(zona)+dipDir]

def get_m(punto,extremo):  #algo pasa midiendo los angulos
    l = LineString([punto, extremo])
    m, b, r, p, std_err = stats.linregress(l.xy[0],l.xy[1])
    return m
def get_angle(m1,m2):
    return np.arctan(abs((m1-m2)/(1+m1*m2)))

#%%
class modeloAS08:
    def __init__(self,Falla,tiempoPeriodo=0.02,AS=False):
    #### PARAMETROS INICIALES OBTENIDOS A PARTIR DE LA FALLA ############################
        self.Falla  = Falla
        
        self.dip    = Falla['dip']
        self.M      = Falla['mw']
        self.Zhipo  = Falla['zhipo']
        self.Ztop   = Falla['ztop']
        self.W      = Falla['w']
        self.FRV    = Falla['FRV']
        self.FN     = Falla['FN']
        self.AS     = 1*AS
        self.Ftipo  = Falla['tipo']
        
        ################################
        if '%.3f'%(tiempoPeriodo) in Tcoeficientes:
            self.T = '%.3f'%(tiempoPeriodo)
        else:
            print("tiempoPeriodo no esta en la lista de Periodos")
            raise ValueError
        ################################
    def __call__(self,longitud,latitud,Vs30,segmento):
#    def __call__(self,Falla,segmento,tiempoPeriodo,AS=False):
    #### PARAMETROS PROPIOS #############################
        self.lat    = latitud
        self.lon    = longitud
        self.Vs30   = Vs30
        self.segmento = segmento
    #### PARAMETROS A PARTIR DE LA FALLA ############################
#        self.Rx,self.Rjb,self.Alpha,self.Rrup,self.Ry= self.Falla.Zona(self.lon,self.lat)
        self.zona, self.Alpha, self.Rjb, self.Rx, self.Ry, self.rRup_p, self.Rrup, self.Ara, self.Arb, self.Arc, self.Ard,self.minDisPP,self.minDisPL,self.dip,self.W,self.Ztop = getPointZona(self.lon, self.lat, self.Falla, self.segmento)
    #### PARAMETROS CALCULADOS A PARTIR DE LOS ANTERIORES
        self.z1     = self.calculo_Z1(self.Vs30)
        self.z25    = self.calculo_Z25(self.z1)
        # Calculamos si el Punto esta en el hanging wall o en el foot wall
        if self.Rx >= 0 and self.dip !=np.pi/2: 
            self.Fhw = 1
        else:
            self.Fhw = 0
        if self.dip == np.pi/2:
            self.Fhw = 0
        #### PARA F1 ####################################
        self.R = np.sqrt((self.Rrup**2.0)+(c4**2.0))
        #### PARA F5 ####################################
        self.PGA = self.calculo_INTENSIDAD()
        return self.lon,self.lat,self.PGA            
    #Calculamos el Z1 que es la velocidad a la Vs = 1.0Km
    def calculo_Z1(self,Vs30):
        if Vs30 < 180.0:
            return np.exp(6.745)
        elif 180.0 <= Vs30 <= 500.0:
            return np.exp(6.745 - 1.35*np.log(Vs30/180.0))
        elif 500.0 < Vs30:
            return np.exp(5.394 - 4.48*np.log(Vs30/500.0))
            
    #Calculamos el Z2.5 que es la velocidad a la Vs = 2.5Km
    def calculo_Z25(self,Z1):
        return 519.0 + 3.595*Z1

# EMPEZAMOS EL CALCULO DE CADA UNO DE LOS TERMINOS DEL MODELO DE ATENUACIÓN:
#### F1(M,Rrup) - Base del Modelo ###########################
    def F1(self,T=None): # el parámetro T es un Texto
        if T == None:
            T=self.T
            
        a1 = Tcoeficientes[T]['a1']
        a2 = Tcoeficientes[T]['a2']
        a3 = Tcoeficientes[T]['a3']
        a4 = Tcoeficientes[T]['a4']
        a5 = Tcoeficientes[T]['a5']
        a8 = Tcoeficientes[T]['a8']
        
        if self.M <= c1:
            return a1+a4*(self.M-c1)+a8*((8.5 - self.M)**2.0)+(a2+a3*(self.M-c1))*np.log(self.R)
        elif self.M > c1:
            return a1+a5*(self.M-c1)+a8*((8.5 - self.M)**2.0)+(a2+a3*(self.M-c1))*np.log(self.R)
            
#### F5(PGA1100,VS30*) - Modelo de Respuesta del Lugar #################################################################
    # Hallando V1
    def calculo_V1(self,T): # el parámetro T es un Número
        if T == -2.0:
            return 862.0 # m/s
        else:
            if T <= 0.50:
                return 1500.0 # m/s
            elif 0.50 < T <= 1.0:
                return np.exp(8.0 - 0.795*np.log(T/0.21)) # m/s
            elif 1.0 < T < 2.0:
                return np.exp(6.76 - 0.297*np.log(T)) # m/s
            elif T >= 2.0:
                return 700.0 # m/s
            
    # Hallando Vs30*    
    def calculo_Vs30star(self,Vs30,T): # el parámetro T es un Número
        V1 = self.calculo_V1(T)
        if Vs30 < V1:
            return Vs30
        elif Vs30 >= V1:
            return V1
            
    # La Funcion F5 propiamente definida        
    def F5(self,PGA_1100,Vs30=None,T=None): # el parámetro T es un Texto
        if Vs30 == None:
            Vs30 = self.Vs30
        if T == None:
            T = self.T
            
        Vs30s= self.calculo_Vs30star(Vs30,float(T))
        a10  = Tcoeficientes[T]['a10']
        b    = Tcoeficientes[T]['b']
        Vlin = Tcoeficientes[T]['Vlin']
        if Vs30 < Vlin:
            return a10*np.log(Vs30s/Vlin) - b*np.log(PGA_1100+c) + b*np.log(PGA_1100+c*((Vs30s/Vlin)**n))
        elif Vs30 >= Vlin:
            return (a10+b*n)*np.log(Vs30s/Vlin)

#### F4(Rjb,Rrup,dip,Ztop,M) - Modelo del Bloque que se mueve ##########################################################
    # La funcion T1
    def calculo_T1(self):
        if self.Rjb < 30.0: #en Kilometros
            return 1.0 - (self.Rjb/30.0)
        elif self.Rjb >= 30.0:
            return 0.0
    # La funcion T2
    def calculo_T2(self):
        wCos = self.W*np.cos(self.dip)
        if self.Rx <= wCos:
            return 0.5 + (self.Rx/(2*wCos))
        if self.Rx > wCos or self.dip == np.pi/2:
            return 1.0
    # La funcion T3
    def calculo_T3(self):
        if self.Rx >= self.Ztop:
            return 1.0
        elif self.Rx < self.Ztop:
            if self.Ztop !=0:
                return (self.Rx/self.Ztop)
            else:
                return 0 # Corregido por que no existe la división entre 0
    # La funcion T4
    def calculo_T4(self):
        if self.M <= 6.0:
            return 0.0
        elif 6.0 < self.M < 7.0:
            return self.M - 6.0
        elif self.M >= 7.0:
            return 1.0
    # La funcion T5
    def calculo_T5(self): #corregido con AS08_NGA errata 2009
        idip = np.rad2deg(self.dip )
        if idip >= 30.0:
            return 1.0 - ((idip - 30.0)/60.0)
        elif idip < 30.0:
            return 1.0
    # La Funcion F4 propiamente definida
    def F4(self,T=None): # el parámetro T es un Texto
        if T == None:
            T = self.T
            
        a14=Tcoeficientes[T]['a14']
        return a14*self.calculo_T1()*self.calculo_T2()*self.calculo_T3()*self.calculo_T4()*self.calculo_T5()
        
#### F6(Ztop) - Modelo de la Ruptura (depth - top) #####################################################################
    def F6(self,T = None): # el parámetro T es un Texto
        if T == None:
            T = self.T
            
        a16=Tcoeficientes[T]['a16']
        if self.Ztop < 10.0: #en Kilometros
            return (a16*self.Ztop)/10.0
        elif self.Ztop >= 10.0:
            return a16

#### F8(Rrup,M) - Modelo a gran distancia ##############################################################################
    # La Función T6    
    def calculo_T6(self):
        if self.M < 5.5:
            return 1.0
        elif 5.5 <= self.M <= 6.5:
            return 0.5*(6.5-self.M)+0.5
        elif self.M > 6.5:
            return 0.5
    # La Función F8 propiamente definida
    def F8(self,T = None): # el parámetro T es un Texto
        if T == None:
            T = self.T
            
        a18=Tcoeficientes[T]['a18']
        if self.Rrup < 100.0:
            return 0.0
        elif self.Rrup >= 100.0:
            return a18*(self.Rrup - 100.0)*self.calculo_T6()*self.M
            
#### F10(Z1.0,Vs30) - Modelo del comportamiento del suelo a profundidad ################################################
    # Función Ln(Z10(Vs30))
    def calculo_Z1s(self,Vs30):
        if Vs30 < 180.0: #en metros/segundo
            return np.exp(6.745)
        elif 180.0 <= Vs30 <= 500.0:
            return np.exp(6.745 - 1.35 * np.log(Vs30/180.0))
        elif Vs30 > 500.0:
            return np.exp(5.394 - 4.48 * np.log(Vs30/500.0))
            
    # Función e2        
    def calculo_e2(self,Vs30,T): # el parámetro T es un Número
        if T < 0.35 or Vs30 > 1000.0:
            return 0.0
        elif 0.35 <= T <= 2.0:
            return (-0.25)*np.log(Vs30/1000.0)*np.log(T/0.35)
        elif T > 2.0:
            return (-0.25)*np.log(Vs30/1000.0)*np.log(2/0.35)
    # Función para hallar a21        
    def calculo_a21(self,Vs30,T): # PGV se computan en T=1segundo ###### el parámetro T es un Texto
        a10  = Tcoeficientes[T]['a10']
        b    = Tcoeficientes[T]['b']
        T = float(T)
        V1 = self.calculo_V1(T)
        Vs30s= self.calculo_Vs30star(Vs30,T)
        nume=(a10+b*n)*np.log(Vs30s/min(V1,1000.0))
        deno=np.log((self.calculo_Z1(Vs30)+c2)/(self.calculo_Z1s(Vs30)+c2))
        e2 = self.calculo_e2(Vs30,T)
        if Vs30 >= 1000.0:
            return 0.0
        elif (nume + e2*deno) < 0.0:
            return  (-1*nume)/deno
        else:
            return e2
    # Función para hallar a22
    def calculo_a22(self,T): # PGV se computan en T=1segundo ##### el parámetro T es un Número
        if T < 2.0:
            return 0.0
        elif T >= 2.0:
            return 0.0625*(T-2)
    # Función F10 propiamente definida
    def F10(self,Vs30=None,T=None): # el parámetro T es un Texto
        if Vs30 == None:
            Vs30 = self.Vs30
        if T == None:
            T = self.T
            
        Z1  = self.calculo_Z1(Vs30)
        Z1s = self.calculo_Z1s(Vs30)
        p_1erTermino = self.calculo_a21(Vs30,T)*np.log((Z1+c2)/Z1s+c2)
        if Z1 >= 200.0:
            return p_1erTermino + self.calculo_a22(float(T))*np.log(Z1/200.0)
        elif Z1 < 200.0:
            return p_1erTermino + 0.0
        
#### Calculamos el PGA1100 Que Pasaremos al F5 #####################################################################
    def calculo_PGA1100(self):
        PGA1100Roca = 0.0
        Vs30Roca = 1100.0
        T = '%.3f'%(-1.0)
        return np.exp(self.F1(T) + self.TRES(T) + self.F5(PGA1100Roca,Vs30Roca,T)+self.Fhw*self.F4(T)+self.F6(T)+self.F8(T)+self.F10(Vs30Roca,T))

#Otras Funciones
    def loglineal(self,x1,x2,y1,y2,x):
        # Interpolacion Lineal
        k = (y2-y1)/(x2-x1)
        C = y1-k*x1
        y = k*x+C
        return y
#Terminos Intermedios
    def TRES(self,T=None):
        if T == None:
            T = self.T
        return Tcoeficientes[T]['a12']*self.FRV + Tcoeficientes[T]['a13']*self.FN + Tcoeficientes[T]['a15']*self.AS
        
#CALCULO DE LA INTENSIDAD
    def calculo_INTENSIDAD(self):
        PGA1100 = self.calculo_PGA1100()
        Td = 10**(-1.25 + 0.3 * self.M)
        Td = min(Td,10)
        if float(self.T) <= Td:
            LnSA = self.F1() + self.TRES() + self.F5(PGA1100) + self.Fhw*self.F4() + self.F6() + self.F8() + self.F10()
        elif float(self.T) > Td:
            periodos = np.array([ -1.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0, -2.0])
            lim_menor = (periodos<Td).nonzero()[0] 
            lim_mayor = (periodos>Td).nonzero()[0]
            Tdmenor = max(periodos[lim_menor])
            Tdmayor = min(periodos[lim_mayor])
            Tdm = '%.3f'%(Tdmenor)
            TdM = '%.3f'%(Tdmayor)
            
            LnSAmenor = self.F1(Tdm) + self.TRES(Tdm) + self.F5(PGA1100,1100,Tdm) + self.Fhw*self.F4(Tdm) + self.F6(Tdm) + self.F8(Tdm) + self.F10(1100,Tdm)
            LnSAmayor = self.F1(TdM) + self.TRES(TdM) + self.F5(PGA1100,1100,TdM) + self.Fhw*self.F4(TdM) + self.F6(TdM) + self.F8(TdM) + self.F10(1100,TdM)
            
            LnSaTd = self.loglineal(np.log(Tdmenor),np.log(Tdmayor),LnSAmenor,LnSAmayor,np.log(Td))
            #Acorde a la formula 22 de Abrahamson & Silva 2008
            LnSA = LnSaTd*((Td/float(self.T))**2) + self.F5(PGA1100) + self.F10()
        return np.exp(LnSA)
    
def calc(segmto,lon,lat,vs30,Falla):
        as08_Modelo = modeloAS08(Falla)
#        pga_Pto = as08_Modelo(lon,lat,vs30)
        return as08_Modelo(lon,lat,vs30,segmento=segmto)
    
def worker(row):
    return calcPGA_pto(row[1],row[2],row[3],row[4])    
       
def calcPGA_pto(lon,lat,vs30,Falla):
    if vs30 > 0:
        as08_Modelo = modeloAS08(Falla)
        PGAs_Pto = [as08_Modelo(lon,lat,vs30,segmto) for segmto in Falla['segs']]
        return max(PGAs_Pto)
    elif vs30 == 0:
        return lon,lat,0
    else:
        return lon,lat,-999999.0
    
def Vs30_SHP_to_Raster(shpVs30,rasterVs30,cellsize):
    gpd_VS30 = gpd.read_file(shpVs30,encoding = 'utf-8')
    VS30_extension = gpd_VS30.unary_union.envelope
    x_min = min(VS30_extension.boundary.xy[0])
    y_max = max(VS30_extension.boundary.xy[1])
    x_max = max(VS30_extension.boundary.xy[0])
    y_min = min(VS30_extension.boundary.xy[1])
    pixelSizeX=cellsize
    pixelSizeY=cellsize
    x_res = int((x_max - x_min) / pixelSizeX)
    y_res = int((y_max - y_min) / pixelSizeY)
    Vs30_raster = rasterVs30
    mb_v = ogr.Open(shpVs30)
    mb_l = mb_v.GetLayer()
    target_ds = gdal.GetDriverByName('GTiff').Create(Vs30_raster, x_res, y_res, 1, gdal.GDT_Float32)
    target_ds.SetGeoTransform((x_min, pixelSizeX, 0, y_max, 0, -pixelSizeY))

    band = target_ds.GetRasterBand(1)
    NoData_value = -999999
    band.SetNoDataValue(NoData_value)
    band.FlushCache()
    gdal.RasterizeLayer(target_ds, [1], mb_l, options=["ATTRIBUTE=VS30"])
    target_ds = None
    return VS30_extension.boundary

def rasterizar(shpVs30,rasterVs30,cellsize):
    if path.exists(rasterVs30):
        rVs30 = rasterio.open(rasterVs30)
        elt=Point(rVs30.bounds.left,rVs30.bounds.top)
        ert=Point(rVs30.bounds.right,rVs30.bounds.top)
        elb=Point(rVs30.bounds.left,rVs30.bounds.bottom)
        erb=Point(rVs30.bounds.right,rVs30.bounds.bottom)
        Extension = LineString([elt,ert,erb,elb,elt])
    else:
        Extension = Vs30_SHP_to_Raster(shpVs30,rasterVs30,cellsize)
        rVs30 = rasterio.open(rasterVs30)
    return rVs30,Extension    
   
