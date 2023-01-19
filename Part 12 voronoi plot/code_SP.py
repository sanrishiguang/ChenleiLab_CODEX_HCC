# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 14:39:30 2022

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


path_to_data = r'D:\Qiu Xinyao\CODEX\LS\SP voronoi\All_CN_TP_voronoi.csv'
cells2 = pd.read_csv(path_to_data)

a=cells2['All_CN'].unique()
b=range(0,16)

celltype_anno=['otherCNs', 
               'CN39_ImmuneMix_B',
               'CN15_ImmMix_B', 
               'CN18_Fibro',
               'CN7_Fibro_HII',
               'CN14_ImmMix_Tumor', 
               'CN20_ImmMix_Mac',
               'CN9_CD44',
               'CN28_CD44_ImmMix',
               'CN35_Ki67_LII', 
               'CN30_Ki67', 
               'CN21_ImmMix_Lym', 
               'CN27_Endo', 
               'CN19_Biliary',     
               'CN29_Hif1a',
               'CN1_HLA_DR']

celltype_map=dict(zip(celltype_anno,b))

color_code=["#787b7d",
            "#8358CF",
            "#cc99ff",
            "#74F882",
            "#ccffcc",
            "#FA339A",
            "#62FFFF",
            "#ffcccc",
            "#f57f7f",
            "#406E30",
            "#83b372",
            "#EE6600",
            "#F6D026",
            "#AE2519",
            "#73904D",
            "#149684"]
label_code=['A','B','C','D','E','F',
            'G', 'H','K','N','P',
            'Q', 'R','S','T',
            'x']
voronoi_hue = 'All_CN'
X = 'XMin'
Y = 'YMin'
edge_color = 'facecolor'
line_width = .1
alpha = 1
figsize = (30,30)

class_=cells2['Class'].unique()
tps=cells2['SP'].unique()

for tp in tps:
    TP = cells2[cells2['SP']==tp]
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
        plt.title(m, fontsize = 100)
        plt.axis('off')
        plt.savefig('...\\hiPPI_{}_{}.tif'.format(str(m), str(tp)))
        plt.show() 


