import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

import pandas as pd

import glob

import os

from vedo import Text3D
import vedo
import numpy as np
from scipy.interpolate import splprep, splev

from brainrender import Scene , actors, settings, actor, camera, cameras
from brainrender.actors import Neuron, Points, Point, Line


# # in brainrender lf = z, ap = x, dv = y


# This document uses master_detailed_comment.csv for neuron groups and soma locations. There are several
# adjustments to be made to get slightly different figures.
# 1. Centroids and border points are optional add-ons applied per-plot via generate_soma_plot()'s
#    include_border / include_centroids arguments (see bottom of file).
# 2. Soma size and transparencty can be adjusted at lines 103 and 104, respectively.



settings.CHECK_FOR_UPDATES = False
settings.SHADER_STYLE = "flat"
settings.SHADER_STYLE = "ambient"
settings.SHOW_AXES = False
settings.DEFAULT_CAMERA = "top"


# function extracted from Brainrender source code to plot the ID of the soma
def add_label(
    point_actor,
    cord,
    size=50,
    color=None
):
   
    cord[0], cord[2] = cord[2], cord[0]
    cord = cord * np.array([-1, -1, 1])
    
    new_actors = []
    color = point_actor.color
    label = point_actor.name
    txt = Text3D(
            label, cord, s=size, c="black", depth=0.1
            )
    new_actors.append(txt.rotate_x(180).rotate_y(180))

    return new_actors



path = str(Path(__file__).parent)
print(path)
screenshots_dir = f"{path}/plots/soma_plots"
os.makedirs(screenshots_dir, exist_ok=True)
df = pd.read_csv(path + "/preprocessed_data/master_detailed_comment.csv", usecols=["neuron_ID", "proj", "x", "y", "z"])
df.loc[df["neuron_ID"]=="1114_10","proj"] = "rsp_orb"
df.loc[df["neuron_ID"]=="1849931_7","proj"] = "rsp_orb"



# RSP and ORB projecting neurons 
rsp_orb = df.loc[df["proj"]=="rsp_orb", "neuron_ID"].to_list()
print(rsp_orb)
# RSP projecting neurons 
rsp = df.loc[df["proj"]=="rsp", "neuron_ID"].to_list()
# ORB projecting neurons
orb = df.loc[df["proj"]=="orb", "neuron_ID"].to_list()





###### CENTROIDS
def calucalte_centroid(group_list):   
    centroid = (
        df[df["neuron_ID"].isin(group_list)]
        .select_dtypes(include="number")
        .sum(axis=0)
        .div(len(group_list))
    )
    return centroid.to_numpy()


## for somata points the alpha is set to 0.6 and the radius to 20 
centroid_rsp_orb = calucalte_centroid(rsp_orb)
centroid_rsp = calucalte_centroid(rsp)
centroid_orb = calucalte_centroid(orb)


# add centroids
def add_centroids(scene):
    centroids_plot1 = actors.Point(centroid_orb, color="#3b83a3", radius = 70, alpha=1)
    centroids_plot2 = actors.Point(centroid_rsp, color="#db2787", radius = 70, alpha=1)
    centroids_plot3 = actors.Point(centroid_rsp_orb, color="#e3bf02", radius = 70, alpha=1)
    scene.add(centroids_plot1, centroids_plot2, centroids_plot3)
    return [centroids_plot1, centroids_plot2, centroids_plot3]



##### CUSTOM BORDER POINTS

def add_borders(scene):
    border_actors = []

    # for custom lateral boarded
    x = np.array([9975, 9650,9025, 7875, 7200])
    y = np.array([2225,3000, 3900, 4575, 5500])
    z = np.array([9125, 9825,10325, 10300, 10175])

    # Fit parametric spline
    tck, u = splprep([x, y, z], s=0)

    # Evaluate smooth curve
    u_fine = np.linspace(0, 1, 20)
    x_s, y_s, z_s = splev(u_fine, tck)

    for xi, yi, zi in zip(x_s, y_s, z_s):
        point_coord = [xi, yi, zi]
        border_point = Point(point_coord, color="dark grey", radius=20) #to add lateral border!!!!
        scene.add(border_point)
        border_actors.append(border_point)

    # for custom medial border
    x = np.array([9786.942980, 9600.088825, 9542.595239, 9542.595239, 9528.221842, 9528.221842, 9506.661747, 9456.354860, 9341.367688])
    y = np.array([2450.0, 2925.0, 3225.0, 3500.0, 3825.0, 4300.0, 4650.0, 4875.0, 5250.0])
    z = np.array([8729.504374, 8585.770409, 8657.637392, 8722.317676, 8779.811262, 8837.304848, 8822.931452, 8808.558055, 8751.064469])

    # Fit parametric spline
    tck, u = splprep([x, y, z], s=0)

    # Evaluate smooth curve
    u_fine = np.linspace(0, 1, 20)
    x_s, y_s, z_s = splev(u_fine, tck)

    for xi, yi, zi in zip(x_s, y_s, z_s):
        point_coord = [xi, yi, zi]
        border_point = Point(point_coord, color="dark grey", radius=20) #to add medial border!!!
        scene.add(border_point)
        border_actors.append(border_point)

    return border_actors





side_angle_view = {'pos': (28112.846874489107, 3501.493190957214, -16800.163794753542), 'focal_point': (5434.172000237177, 4186.207857314123, -7353.979126679797), 'viewup': (-0.017424366656200126, -0.9993795800862043, 0.030607945918536997), 'distance': 24576.849533854038, 'clipping_range': (11375.652346347175, 43146.76954489577)}

zoom_in_on_soma = {
    "pos": (8699.956534515386, -9539.068724921435, -9184.347059174846),
    "viewup": (-1.0, 0.0, 0.0),
    "clipping_range": (5783.754263567886, 23037.430492592688),
    'focal_point': (8726.844246108687, 4317.75295241172, -9088.53876886988)
}

three_quarters_camera = {
    "pos": (-20169, -7298, 14832),
    "viewup": (0, -1, 0),
    "clipping_range": (16955, 58963),
}




def generate_soma_plot(filename, camera, include_border=False, include_centroids=False, zoom=1.5, scale=10):
    # set the scene with allen mouse atlas
    scene = Scene(title="" ,  inset=False, atlas_name='allen_mouse_25um', check_latest=False,screenshots_folder= screenshots_dir)
    if include_centroids:
        radius = 20 
        alpha = 0.6
    else:
        radius = 30 
        alpha = 1 
    number=0

    # # go over each of the neuron soma and extract the coordinates
    for index, soma in df.iterrows():

        x_soma = soma["x"]
        y_soma = soma["y"]
        z_soma = soma["z"]
        
        coordinates = [x_soma, y_soma, z_soma]

        soma_id = soma["neuron_ID"]
        #print(soma_id[-5:])
        #print(coordinates)



        # colour soma with colours based on the group the neuron belongs to
        if soma_id in rsp_orb:
            colour = "#FFD600"
            soma_plot = actors.Point(coordinates, color=colour, radius=radius, name = soma_id[-5:], alpha = alpha)
            soma_number = add_label(point_actor = soma_plot, cord = coordinates)
            scene.add(soma_plot)
            #print(soma_id)
            number+=1
            
            
            #scene.add(*soma_number)

        
        elif soma_id in orb:
            colour = "#49A2C9"
            soma_plot = actors.Point(coordinates, color=colour, radius=radius, name = soma_id[-5:], alpha = alpha)
            soma_number = add_label(point_actor = soma_plot, cord = coordinates)
            scene.add(soma_plot)
            #print(soma_id)
            number+=1
            #scene.add(*soma_number)
        
        elif soma_id in rsp:
            colour = "#EA4FA2"
            soma_plot = actors.Point(coordinates, color=colour, radius=radius, name = soma_id[-5:], alpha = alpha)
            soma_number = add_label(point_actor = soma_plot, cord = coordinates)
            scene.add(soma_plot)
            #print(soma_id)
            number+=1
            #scene.add(*soma_number)
        else:
            #print(f"not in the list: {soma_id}")
            continue


    if include_border:
        add_borders(scene)
    if include_centroids:
        add_centroids(scene)

    scene.screenshot(name=filename, zoom=zoom, scale=scale, camera=camera)

   


generate_soma_plot("soma_plot_from_back_with_border.png", side_angle_view, include_border=True)
generate_soma_plot("soma_plot_zoomed_without_border.png", zoom_in_on_soma)
generate_soma_plot("soma_plot_zoomed_with_centroids_without_border.png", zoom_in_on_soma, include_centroids=True)

#scene.export("soma_full_brain.html")




