import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

import pandas as pd

import glob

import os

from vedo import Text3D
import vedo

from brainrender import Scene , actors, settings, actor, camera, cameras
from brainrender.actors import Neuron, Points, Point

settings.CHECK_FOR_UPDATES = False
settings.SHADER_STYLE = "ambient"
settings.SHOW_AXES = False
settings.DEFAULT_CAMERA = "top"
settings.BACKGROUND_COLOR = "white"

path = os.path.dirname(os.path.abspath(__file__))


scene = Scene(title="", atlas_name='allen_mouse_25um', check_latest=False, screenshots_folder=path)


soma_file = pd.read_csv(f"{path}/input_data/merfish_67gene_plus_extra_genes.csv")



# # go over each of the neuron soma and extract the coordinates
for index, soma in soma_file.iterrows():

    x_soma = soma["x_ccf"]
    y_soma = soma["y_ccf"]
    z_soma = soma["z_ccf"]
    
    coordinates = [x_soma, y_soma, z_soma]

    soma_id = soma["supertype"]
    #print(soma_id[-5:])
    print(coordinates)


    # colour soma with colours based on the group the neuron belongs to
    if soma_id == "0007 L5/6 IT TPE-ENT Glut_1":
        colour = "orange" #"#CBBCA1" 
        soma_plot = actors.Point(coordinates, color=colour, radius=20, name = soma_id[-5:], alpha = 1)
        scene.add(soma_plot)

       
    elif soma_id == "0008 L5/6 IT TPE-ENT Glut_2":
        colour = "#F1D68C"
        soma_plot = actors.Point(coordinates, color=colour, radius=20, name = soma_id[-5:], alpha = 1)
        scene.add(soma_plot)
       
    
    elif soma_id == "0010 L5/6 IT TPE-ENT Glut_4":
        colour = "#73969A"
        soma_plot = actors.Point(coordinates, color=colour, radius=20, name = soma_id[-5:], alpha = 1)
        scene.add(soma_plot)
       
    elif soma_id == "0129 NP PPP Glut_1":
        colour = "#E9B9A9"
        soma_plot = actors.Point(coordinates, color=colour, radius=20, name = soma_id[-5:], alpha = 1)
        scene.add(soma_plot)
        
        continue

scene.add_brain_region("ENT",alpha=0.4, color = "#5977D3")

plane = scene.atlas.get_plane(plane = "horizontal", pos = (0,4000,0), norm = (0,1,0))
plane2 = scene.atlas.get_plane(plane = "horizontal", pos = (0,5000,0), norm = (0,-1,0))
# plane = scene.atlas.get_plane(plane = "horizontal", pos = (0,3100,0), norm = (0,1,0))
# plane2 = scene.atlas.get_plane(plane = "horizontal", pos = (0,4000,0), norm = (0,-1,0))

cam_angled = dict(
    pos=(31808.5, -18205.9, 9684.36),
    focal_point=(8230.24, 4523.06, -5674.10),
    viewup=(-0.538962, -0.777720, -0.323530),
    roll=-178.678,
    distance=36172.1,
    clipping_range=(23289.7, 54878.8),
)
cam_zoomed_in = dict(
    pos=(21810.3, -6837.64, 4177.79),
    focal_point=(6401.61, 4805.27, -2449.16),
    viewup=(-0.395295, -0.789951, -0.468743),
    roll=164.072,
    distance=20418.2,
    clipping_range=(6824.34, 39891.0),
)
scene.slice(plane)
scene.slice(plane2)

# scene.render(camera=cam_zoomed_in)
scene.screenshot(name = "subtypes_4000_5000.png", zoom = 1.5, scale=10, camera = cam_zoomed_in)
#scene.screenshot(name = "subtypes_31000_4000.png", zoom = 1.5, scale=10, camera = cam_zoomed_in)


