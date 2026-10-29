import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

import pandas as pd

import glob

import os

from vedo import Text3D
import vedo
import anndata as ad
from brainrender import Scene , actors, settings, actor, camera, cameras
from brainrender.actors import Neuron, Points, Point
import matplotlib
from matplotlib.colors import to_hex
from tqdm import tqdm
import shutil



settings.CHECK_FOR_UPDATES = False
settings.SHADER_STYLE = "ambient"
settings.SHOW_AXES = False
settings.DEFAULT_CAMERA = "top"
settings.BACKGROUND_COLOR = "white"

cam_zoomed_in = dict(
    pos=(21810.3, -6837.64, 4177.79),
    focal_point=(6401.61, 4805.27, -2449.16),
    viewup=(-0.395295, -0.789951, -0.468743),
    roll=164.072,
    distance=20418.2,
    clipping_range=(6824.34, 39891.0),
)

path = os.path.dirname(os.path.abspath(__file__))



# read the MERFISH gene file to know cells expression of each gene
soma_file = pd.read_csv(f"{path}/input_data/merfish_67gene_plus_extra_genes.csv")
# adata = ad.read_h5ad(f"{path}/input_data/merfish_zhuang_ec_supertypes_of_interest.h5ad")
# soma_file = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
# soma_file = soma_file.join(adata.obs).reset_index()
# soma_file.to_csv("merfish_gene_expression_coords.csv")

# read the gene file with the gene names to link short name with ID
categories = pd.read_csv(f"{path}/output_data/merfish_final_gene_categories.csv")
gene_names = pd.read_csv(f"{path}/output_data/bulk_gene_classes.csv", usecols=["gene_id", "gene_name"])
genes_of_interest_ensambleID = gene_names.loc[gene_names["gene_name"].isin(['Pou3f3','Prss12','Gpr88','Boc',"Fibcd1","Cyp26b1"]),"gene_id"].to_list()

remaining = genes_of_interest_ensambleID

file_list = glob.glob(str(path)+'/*.png')



for column in tqdm(remaining, desc="Processing genes"):
# # go over each of the neuron soma and extract the coordinates
    if "ENS" in column and not os.path.exists(f"{path}/{column}_4000_6000_.png"):

        # set the scene with allen mouse atlas
        scene = Scene(title="", atlas_name='allen_mouse_25um', check_latest=False, screenshots_folder=path)
        soma_file_filtered = soma_file[[column, "x_ccf", "y_ccf", "z_ccf"]]
        cmap_norm = matplotlib.colors.Normalize(vmin=0, vmax=soma_file_filtered[column].max())
        print(column)
        for index, soma in soma_file_filtered.iterrows():
            x_soma = soma["x_ccf"]
            y_soma = soma["y_ccf"]
            z_soma = soma["z_ccf"]
            cmap = matplotlib.colormaps["Reds"]


            coordinates = [x_soma, y_soma, z_soma]

            colour = to_hex(cmap(cmap_norm(soma[column])), keep_alpha=False)

            soma_plot = actors.Point(coordinates, color=colour, radius=20,  alpha = 1)
            scene.add(soma_plot)
        scene.add_brain_region("ENT",alpha=0.4, color = "#5977D3")
        plane = scene.atlas.get_plane(plane = "horizontal", pos = (0,4000,0), norm = (0,1,0))
        plane2 = scene.atlas.get_plane(plane = "horizontal", pos = (0,6000,0), norm = (0,-1,0))
        scene.slice(plane)
        scene.slice(plane2)
        file_name = f"{column}_4000_6000_.png"
        scene.screenshot(name = file_name, zoom = 1.5, scale=10, camera = cam_zoomed_in)
        scene.close()
        scene.plotter.close()
        del scene.actors
        del scene
        import gc
        gc.collect()
    if "ENS" in column and not os.path.exists(f"{path}/{column}_3100_4000_.png"):

        # set the scene with allen mouse atlas
        scene = Scene(title="", atlas_name='allen_mouse_25um', check_latest=False, screenshots_folder=path)
        soma_file_filtered = soma_file[[column, "x_ccf", "y_ccf", "z_ccf"]]
        cmap_norm = matplotlib.colors.Normalize(vmin=0, vmax=soma_file_filtered[column].max())
        print(column)
        for index, soma in soma_file_filtered.iterrows():
            x_soma = soma["x_ccf"]
            y_soma = soma["y_ccf"]
            z_soma = soma["z_ccf"]
            cmap = matplotlib.colormaps["Reds"]


            coordinates = [x_soma, y_soma, z_soma]

            colour = to_hex(cmap(cmap_norm(soma[column])), keep_alpha=False)

            soma_plot = actors.Point(coordinates, color=colour, radius=20,  alpha = 1)
            scene.add(soma_plot)
        scene.add_brain_region("ENT",alpha=0.4, color = "#5977D3")
        plane = scene.atlas.get_plane(plane = "horizontal", pos = (0,3100,0), norm = (0,1,0))
        plane2 = scene.atlas.get_plane(plane = "horizontal", pos = (0,4000,0), norm = (0,-1,0))
        scene.slice(plane)
        scene.slice(plane2)
        file_name = f"{column}_3100_4000_.png"
        scene.screenshot(name = file_name, zoom = 1.5, scale=10, camera = cam_zoomed_in)
        scene.close()
        scene.plotter.close()
        del scene.actors
        del scene
        import gc
        gc.collect()
    if "ENS" in column and not os.path.exists(f"{path}/{column}_full.png"):

        # set the scene with allen mouse atlas
        scene = Scene(title="", atlas_name='allen_mouse_25um', check_latest=False, screenshots_folder=path)
        soma_file_filtered = soma_file[[column, "x_ccf", "y_ccf", "z_ccf"]]
        cmap_norm = matplotlib.colors.Normalize(vmin=0, vmax=soma_file_filtered[column].max())
        print(column)
        for index, soma in soma_file_filtered.iterrows():
            x_soma = soma["x_ccf"]
            y_soma = soma["y_ccf"]
            z_soma = soma["z_ccf"]
            cmap = matplotlib.colormaps["Reds"]


            coordinates = [x_soma, y_soma, z_soma]

            colour = to_hex(cmap(cmap_norm(soma[column])), keep_alpha=False)

            soma_plot = actors.Point(coordinates, color=colour, radius=20,  alpha = 1)
            scene.add(soma_plot)
        scene.add_brain_region("ENT",alpha=0.4, color = "#5977D3")


        file_name = f"{column}_full.png"
        scene.screenshot(name = file_name, zoom = 1.5, scale=10, camera = cam_zoomed_in)
        scene.close()
        scene.plotter.close()
        del scene.actors
        del scene
        import gc
        gc.collect()
    else:
        continue



gene_expression_stats_nb_categories = categories.copy()
gene_expression_stats_nb_categories = gene_expression_stats_nb_categories.loc[gene_expression_stats_nb_categories["category"]!="none"]
gene_categories = (
    gene_expression_stats_nb_categories.groupby("category")["gene_id"]
    .apply(list)
    .to_dict()
)
for category in tqdm(gene_categories.keys()):
    dest_dir = Path(path) / category
    dest_dir.mkdir(exist_ok=True)

    for gene in gene_categories[category]:
        matches = glob.glob(f"{path}/{gene}_4000_6000_.png")
        if matches:
            shutil.move(matches[0], dest_dir)
        matches = glob.glob(f"{path}/{gene}_3100_4000_.png")
        if matches:
            shutil.move(matches[0], dest_dir)
        matches = glob.glob(f"{path}/{gene}_full.png")
        if matches:
            shutil.move(matches[0], dest_dir)
        else:
            print(f"Image file not found for {gene}.")
