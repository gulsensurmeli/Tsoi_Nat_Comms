### PLOT ALL NEURONS (AXONS ONLY) THAT WERE USED IN ANALYSES

from pathlib import Path

import pandas as pd
import numpy as np

import os


from brainrender import Scene , actors, settings, actor, camera, cameras
from brainrender.actors import Neuron, Points, Point


settings.SHADER_STYLE = "ambient"
settings.SHOW_AXES = False



# IF ONLY SPECIFIC NEURONS ARE WANTED FOR PLOTTING:
# 1. Uncomment line 51 and 52
# 2. In line 51 and add neuron IDs within brackets e.g.: wanted_neurons = ["1114_11","1135_10"]
# 3. Run the code; lines 85 and 105 determine where the images will be saved

# to remove small brain at the bottom (the inset) got to -> class Scene and change instet=True to False

path = Path(__file__).parent
screenshots_dir = f"{path}/plots"
os.makedirs(screenshots_dir, exist_ok=True)


top_side_camera = cam = dict(
    pos=(3390.71, -31489.4, -5386.91),
    focal_point=(7215.95, 4630.65, -5650.02),
    viewup=(-0.0147359, 8.84362e-3, 0.999852),
    roll=85.1943,
    distance=36323.0,
    clipping_range=(33317.6, 39751.0),
)

sagittal_camera2 = {
    "pos": (9782, 1795, -40999),
    "viewup": (0, -1, 0),
    "clipping_range": (23256, 51031),
}



df = pd.read_csv(f"{path}/master_detailed_comment.csv", usecols=["neuron_ID", "proj", "x", "y", "z"])
# wanted_neurons = ["1056_16"]
# df = df.loc[df["neuron_ID"].isin(wanted_neurons)]


# RSP and ORB projecting neurons 
rsp_orb = df.loc[df["proj"]=="rsp_orb", "neuron_ID"].to_list()
# RSP projecting neurons 
rsp = df.loc[df["proj"]=="rsp", "neuron_ID"].to_list()

# ORB projecting neurons
orb = df.loc[df["proj"]=="orb", "neuron_ID"].to_list()
print(orb)



df.loc[df["neuron_ID"]=="1114_10","proj"] = "rsp_orb"
df.loc[df["neuron_ID"]=="1849931_7","proj"] = "rsp_orb"


groupings = {
    "orb": orb,
    "rsp_orb": rsp_orb,
    "rsp": rsp
}


def neuron(group):
    if group == "orb":
         colour = "#49A2C9"
    elif group == "rsp":
         colour = "#EA4FA2"
    elif group == "rsp_orb":
         colour = "#FFB300"
    for i in groupings[group]:
        i
        full_id_list = i.split("_")
        mouseID = full_id_list[0]
        neuronNR = full_id_list[1]
        neuron_file_path = f"{path}/preprocessed_data/{mouseID}/{mouseID}_{neuronNR}_neurites.xlsx"
        
        neuron_file = pd.read_excel(neuron_file_path)
        neuron_file = neuron_file[neuron_file['Ntype'].isin([2,5,6])]
        #neuron_file = neuron_file[~neuron_file['Abbreviation'].astype(str).str.startswith("ENT")] 
        print(neuron_file_path)

        cortical_axon_coords = []
        coords = neuron_file[['x', 'y', 'z']].to_numpy() * 25
        
        if not os.path.exists(f"{screenshots_dir}/{i}_top.png"):
            scene = Scene(title="", inset=False,atlas_name='allen_mouse_25um', check_latest=False, screenshots_folder= screenshots_dir)
            scene.add(Points(coords, radius = 20, colors=colour, alpha=1)) # to add the points
            scene.screenshot(name = i + "_top", zoom = 1.5, scale=2, camera = top_side_camera)
            
        if not os.path.exists(f"{screenshots_dir}/{i}_side.png"):
            scene = Scene(title="", inset=False,  atlas_name='allen_mouse_25um', check_latest=False, screenshots_folder= screenshots_dir)
            scene.add(Points(coords, radius = 20, colors=colour, alpha=1))
            scene.screenshot(name = i + "_side", zoom = 1.5, scale=2, camera = sagittal_camera2)
        
        #scene.render()
         


neuron("rsp")
neuron("orb")
neuron("rsp_orb")


