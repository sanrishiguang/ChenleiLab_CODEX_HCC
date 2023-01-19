# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 14:46:23 2022

@author: dell
"""
#Python 3.9.12 
import numpy as np  #1.22.4
import matplotlib.pyplot as plt  #3.5.1
import seaborn as sns  #0.11.2
from shapely.geometry import MultiPoint, Point, Polygon  #1.8.2
from scipy.spatial import Voronoi #1.7.3
import shapely.geometry #1.8.2
import pandas as pd  #1.4.2
from sklearn.preprocessing import LabelEncoder  #1.0.2


def voronoi_finite_polygons_2d(vor, radius=None):

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1] # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)


path_to_data = r'D:\Qiu Xinyao\CODEX\LS\TP voronoi\All_CN_TP_voronoi.csv'
cells2 = pd.read_csv(path_to_data)
a=cells2['All_CN'].unique()
b=range(0,26)
celltype_anno=['otherCNs', 
               'CN2_others_LII',  ##inner
               'CN34_others_LII',  ##inner
               'CN5_others', 
               'CN10_others',
               'CN0_Hepar1_LII', ##inner
               'CN16_Hepar1', 
               'CN38_Hepar1_HII',
               'CN13_Hepar1_LII', 
               'CN33_c_Myc_LII',  ##inner
               'CN3_c_Myc', 
               'CN17_p_mTOR', ##tumor
               'CN32_p_mTOR', ##stromal  
               'CN8_CD107a', ##inner
               'CN31_CD107a',   
               'CN36_p53', ##inner
               'CN4_p53',       
               'CN23_p_S6', ##inner
               'CN11_p_S6',        
               'CN6_E_Cad_LII',   ##inner
               'CN25_E_Cad',    
               'CN12_Glypican3', 
               'CN24_Twist1',  
               'CN37', 
               'CN26', 
               'CN22']

celltype_map=dict(zip(celltype_anno,b))
color_code=["#787b7d",
            "#FA339A",
            "#ff6699",
            "#ff99cc",
            "#cc6699",
            "#7E8EE3",
            "#99ccff",
            "#669999",
            "#D083FF",
            "#8358CF",
            "#cc99ff",
            "#FBF88B",
            "#ffffcc",
            "#AE2519",
            "#cc6666",
            "#EE6600",
            "#ff9966",
            "#f1dd60",
            "#d8cb93",
            "#62FFFF",
            "#ccffff",
            "#27A7FA",
            "#D92B00",
            "#74F882",
            "#9C9B99",
            "#cc99cc"]
            

label_code=['A','B','C','D','E','F',
            'G', 'H','I','J','K','L','M','N','O','P',
            'Q', 'R','S','T','V','W','X','Y','Z',
            'n']

voronoi_hue = 'All_CN'
X = 'XMin'
Y = 'YMin'
edge_color = 'facecolor'
line_width = .1
alpha = 1
figsize = (30,30)

tps=cells2['TP'].unique()
for tp in tps:
    TP = cells2[cells2['TP']==tp]
    class_=TP['Class'].unique()
    for m in class_:
        spot = cells2[cells2['Class']==m]
        spot['All_CN'] = spot['All_CN'].map(celltype_map) 
        labels  = [label_code[i] for i in spot[voronoi_hue]]
        colors  = [color_code[i] for i in spot[voronoi_hue]]
        points=spot[[X,Y]].values
        size_max=4000 ##voronoi_kwargs={'size_max':4000}
        points[:,1] = max(points[:,1])-points[:,1]
        vor = Voronoi(points)
        label=[]
        regions, vertices = voronoi_finite_polygons_2d(vor)
        pts = MultiPoint([Point(i) for i in points])
        mask = pts.convex_hull
        new_vertices = []
        if type(alpha)!=list:
            alpha = [alpha]*len(points)
        areas = []
        plt.figure(figsize=figsize)
        for i,(region,alph) in enumerate(zip(regions,alpha)):
            polygon = vertices[region]
            shape = list(polygon.shape)
            shape[0] += 1
            p = Polygon(np.append(polygon, polygon[0]).reshape(*shape)).intersection(mask)
            areas+=[p.area]
            if p.area <size_max:
                poly = np.array(list(zip(p.boundary.coords.xy[0][:-1], p.boundary.coords.xy[1][:-1])))
                new_vertices.append(poly)
                if edge_color == 'facecolor':
                    plt.fill(*zip(*poly), alpha=alph,edgecolor= 'black',linewidth = line_width , facecolor = colors[i])
                    plt.text(poly.mean(0)[0],poly.mean(0)[1],s=labels[i],fontsize=8)
                else:
                    plt.fill(*zip(*poly), alpha=alph,edgecolor=  edge_color,linewidth = line_width, facecolor = colors[i])
        mu=0
        for color in color_code:
            plt.scatter([],[],c=color,s=30,label=label_code[mu]+':'+celltype_anno[mu])
            mu+=1
        font = {'family': 'Arial', 'weight': 'normal', 'size': 20}
        #plt.legend(prop=font,loc="upper right")
        plt.title(m, fontsize = 100)
        plt.axis('off')
        plt.savefig('D:\\Qiu Xinyao\\CODEX\\LS\\TP voronoi\\single\\hiPPI_{}_{}.tif'.format(str(m), str(tp)))##需要师姐修改保存的路径
        plt.show() 




