from importlib.resources import path

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

import pandas as pd

import glob

import os

from vedo import *

from brainrender import Scene , actors, settings, actor, camera, cameras
from brainrender.actors import Neuron, Points, Point, PointsDensity



# areas of interest
isocortical_areas = [
    'Isocortex'
]


select_isocortical = ['ILA','PL', 'ORB', 'AUD', 'TEa', 'ECT', 'PERI', 'VISC', 
    'AI', 'MO', 'SS', 'ACA', 'PTLp','VIS',  'RSP'] # GU and frontal pole excluded


subcortical_areas = ['OLF', 'CTXsp', 
        'CLA', 'EP', 'BLA', 'LA', 'BMA', 'PA', 
    'STR',     
    'PAL', 
    'HPF']

remove_ec = ['ENT'] # we don't want to include ec in the plotting becuase of how dense the local axons are and it will bias the visual


# path to this script's folder, where structure_tree_safe.csv lives
path = str(Path(__file__).parent)

# load the structure tree (same file/columns as used in preprocessing_of_files.ipynb)
structure_tree = pd.read_csv(f"{path}/structure_tree_safe.csv")


# creates a list of id paths for the areas of interest
def get_id_paths(areas, structure_tree):
    return [structure_tree.loc[structure_tree["acronym"] == area, "structure_id_path"].values[0]
            for area in areas]

# makes the list into regex pattern
def make_pattern(id_paths):
    return '|'.join(map(str, id_paths))


isocortical_area_id_paths = get_id_paths(isocortical_areas, structure_tree)
iso_pattern = make_pattern(isocortical_area_id_paths)

select_isocortical_area_id_paths = get_id_paths(select_isocortical, structure_tree)
select_isocortical_pattern = make_pattern(select_isocortical_area_id_paths)

subcortical_area_id_paths = get_id_paths(subcortical_areas, structure_tree)
sub_pattern = make_pattern(subcortical_area_id_paths)

ec_id = get_id_paths(remove_ec, structure_tree)





# in brainrender scene.py file set inset=False to remove little brain icon at the bottom of the scene
settings.CHECK_FOR_UPDATES = False
settings.SHOW_AXES = False
settings.SHADER_STYLE = "ambient"




print(path)

screenshots_dir = f"{path}/plots/density_maps"
os.makedirs(screenshots_dir, exist_ok=True)

df = pd.read_csv(path + "/preprocessed_data/master_detailed_comment.csv", usecols=["neuron_ID", "proj", "x", "y", "z"])


# RSP and ORB projecting neurons 
rsp_orb = df.loc[df["proj"]=="rsp_orb", "neuron_ID"].to_list()
print(rsp_orb)
# RSP projecting neurons 
rsp = df.loc[df["proj"]=="rsp", "neuron_ID"].to_list()
# ORB projecting neurons
orb = df.loc[df["proj"]=="orb", "neuron_ID"].to_list()


# grouping based on endpoinbt count and axon length
groupings = {
    "orb": orb,
    "rsp_orb": rsp_orb,
    "rsp": rsp
}




# get all the neurite files with all neuron puncta coordinates
# everything is in parent folder, code is in folder Code and all brain data is in folders starting ENT*

list_of_files_full = glob.glob(str(path) + '/preprocessed_data/*/*neurites.xlsx') # all the forlders must start with ENT and then the brain ID


# need to make 12 different plots with these variables: sagittal2/top_side, endpoint/all axon, isocortical/all areas

def density_map (camera_angle, axon_part, areas, group):
    """
    The function to automaticaly generate 12 different plots for each condition

    camera_angle can be "sagittal2" or "top_side"

    axon_part can be "full_axon"  or "endpoints"

    areas can be "subcortical" or "isocortex"

    group can be "orb", "rsp", "rsp_orb"
    """
    settings.DEFAULT_CAMERA = camera_angle
    scene = Scene(title="", atlas_name='allen_mouse_25um',inset=False, check_latest=False, screenshots_folder= screenshots_dir)

    if axon_part == "full_axon":
        axon_ntype = [2,5,6]
    elif axon_part== "endpoints":
        axon_ntype = [6]
  
    cortical_axon_coords =  [] # empty list to fill with coordinates
    for neuron in groupings[group]:
        full_id_list = neuron.split("_")
        mouseID = full_id_list[0]
        neuron_file_path = path + "/preprocessed_data/" + mouseID + "/" + neuron + "_neurites.xlsx"
        axon_points = pd.read_excel(neuron_file_path, usecols=['Ntype', 'structure_id_path', 'x', 'y', 'z'])
        print(neuron_file_path)
        axon_points_filtered = axon_points[axon_points['Ntype'].isin(axon_ntype)]
        axon_points_filtered = axon_points_filtered[axon_points_filtered['z']> 228.5]
        
        
        if areas == "subcortical":
            axon_points_filtered = axon_points_filtered[axon_points_filtered['structure_id_path'].astype(str).str.contains(sub_pattern)] # to filter for subcortical areas
            axon_points_filtered = axon_points_filtered[~axon_points_filtered['structure_id_path'].astype(str).str.startswith(ec_id[0])] # to remove all points that are in the entorhinal area, which is the only isocortical area in the subcortical list
        elif areas == "isocortex":
            axon_points_filtered = axon_points_filtered[axon_points_filtered['structure_id_path'].astype(str).str.contains(iso_pattern)] # to filter for isocortical areas
        
        coords = axon_points_filtered[['x', 'y', 'z']].to_numpy() * 25
        cortical_axon_coords.append(coords)
        
        
    cortical_axon_coords_np = np.concatenate(cortical_axon_coords, axis=0) # to create a numpy array which is accepted to Points 
    #scene.add(Points(cortical_axon_coords_np, radius = 20)) # to add the points
    scene.add(PointsDensity(cortical_axon_coords_np, colors = "gnuplot", radius=669))

    name = group + "(" + areas + ", " + axon_part + ", " + camera_angle + ")"
    scene.screenshot(name = name, zoom = 1.5, scale=2) # change name accordingly to what is being plotted
    #scene.render()




    
density_map("sagittal2","full_axon","subcortical","orb") 
density_map("sagittal2","full_axon","subcortical","rsp")
density_map("sagittal2","full_axon","subcortical","rsp_orb")

# density_map("top_side","full_axon","subcortical","orb") 
# density_map("top_side","full_axon","subcortical","rsp")
# density_map("top_side","full_axon","subcortical","rsp_orb")

# density_map("sagittal2","endpoints","subcortical","orb")
# density_map("sagittal2","endpoints","subcortical","rsp")
# density_map("sagittal2","endpoints","subcortical","rsp_orb")

# density_map("top_side","endpoints","subcortical","orb") 
# density_map("top_side","endpoints","subcortical","rsp")
# density_map("top_side","endpoints","subcortical","rsp_orb")


density_map("sagittal2","full_axon","isocortex","orb") 
density_map("sagittal2","full_axon","isocortex","rsp")
density_map("sagittal2","full_axon","isocortex","rsp_orb")

# density_map("top_side","full_axon","isocortex","orb") 
# density_map("top_side","full_axon","isocortex","rsp")
# density_map("top_side","full_axon","isocortex","rsp_orb")

# density_map("sagittal2","endpoints","isocortex","orb")
# density_map("sagittal2","endpoints","isocortex","rsp")
# density_map("sagittal2","endpoints","isocortex","rsp_orb")

# density_map("top_side","endpoints","isocortex","orb") 
# density_map("top_side","endpoints","isocortex","rsp")
# density_map("top_side","endpoints","isocortex","rsp_orb")

    


