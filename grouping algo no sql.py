# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 14:06:01 2023

@author: benk
"""

from scipy.spatial.distance import pdist, squareform
from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import copy


class Point():
    def __init__(self, x, y):
        self.x, self.y = x, y
        self.points = [(self.x, self.y)]
        self.center = self.x, self.y
        self.x_points = [x]
        self.y_points = [y]
        
    def __add__(self, other_point):
        self.points += other_point.points
        self.update()
        
    def update(self):
        self.x_points = [point[0] for point in self.points]
        self.y_points = [point[1] for point in self.points]
        self.x  = sum(self.x_points)/len(self.x_points)
        self.y = sum(self.y_points)/len(self.y_points)
        self.center = self.x, self.y
        
    def check_max_distance(self, other_point = None):
        if other_point is not None:
            points = self.points + other_point.points
        else:
            points = self.points
        dist_matrix = pdist(points)
        max_dist = np.max(dist_matrix)
        return max_dist
    
    def radius(self):
        points = self.points + [(self.center)]
        dist_matrix =  squareform(pdist(points))
        radius = np.max(dist_matrix[:, -1])
        return radius
    
    
def main(points, dist_limit):
    
    print('main called with points list of length ', len(points))
    point_centers = [point.center for point in points]
    all_dist = squareform(pdist(point_centers))
    all_dist[np.tril_indices(all_dist.shape[0])] = np.inf


    points_to_use = list(combinations(range(len(points)), 2))
    dist_by_idx = [(val[0], val[1], all_dist[val[0], val[1]]) for val in points_to_use]
    dist_by_idx = sorted(dist_by_idx, key = lambda val: val[2])
    
    for par in dist_by_idx:
        
        point1 = points[par[0]]
        point2 = points[par[1]]
        
        max_distance = point1.check_max_distance(point2)
        if max_distance <= dist_limit:
            point1 + point2
            points.remove(point2)
            return main(points, dist_limit)
            
            
    return points



points = np.random.randint(0,500, (500,2))
# points = list(zip(df.x, df.y))
dist_limit = 50
points = [Point(par[0], par[1]) for par in points]
og_points = copy.deepcopy(points)
groups = main(points, dist_limit)


x = [point.x for point in og_points]
y = [point.y for point in og_points]
centers_x = [point.center[0] for point in groups]
centers_y = [point.center[1] for point in groups]

fig, ax = plt.subplots()
ax.scatter(x, y, c='red')
ax.scatter(centers_x, centers_y, c='blue')
print('groups', len(groups))
print('points', len(x))
dups = []
for group in groups:
    if len(group.x_points)>1:
        dups.append(group)
        ax.add_patch(plt.Circle( (group.x, group.y ),
                                          group.radius() ,
                                          fill = False ))
plt.show()

