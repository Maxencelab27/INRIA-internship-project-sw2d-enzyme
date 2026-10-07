#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mesh_cost_plot.py — trace la fonction cout (misfit) par cellule sur le maillage
pour un couple de fichiers de résultats spécifique (comparaison directe instantanée).

Usage :
  python3 utils/python/plot/mesh_plot/mesh_cost_plot.py <mesh.2dm> <ref_file.txt> <fwd_file.txt> [--log] [--save cost_per_cell.csv]

Exemple :
  python3 utils/python/plot/mesh_plot/mesh_cost_plot.py \
      sp_gardon_travail/primary_files/Gardon_sp.2dm \
      sp_gardon_travail/reference/ResultHydro_p0000d03h55m00s000.txt \
      sp_gardon_copy_2/output/ResultHydro_p0000d03h55m00s000.txt
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

def read_2dm(file_path):
    nodes = {}
    elements = []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            if parts[0] == "ND":
                nodes[int(parts[1])] = (float(parts[2]), float(parts[3]))
            elif parts[0] == "E3T":
                _, _, n1, n2, n3, *_ = parts
                elements.append([int(n1), int(n2), int(n3)])
            elif parts[0] == "E4Q":
                _, _, n1, n2, n3, n4, *_ = parts
                elements.append([int(n1), int(n2), int(n3), int(n4)])
    return nodes, elements


def read_frame_hq(file_path, ncells):
    """Lit h, qx, qy (colonnes 3,4,5) indexes par numero de ligne (= cell_id)."""
    h  = np.zeros(ncells); qx = np.zeros(ncells); qy = np.zeros(ncells)
    with open(file_path, 'r') as f:
        row = 0
        header_seen = False
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            if 'x' in parts and 'y' in parts:
                header_seen = True
                continue
            if not header_seen:
                continue
            if row >= ncells:
                break
            try:
                h[row]  = float(parts[3])
                qx[row] = float(parts[4])
                qy[row] = float(parts[5])
            except (ValueError, IndexError):
                pass
            row += 1
    return h, qx, qy


def panel(ax, nodes, elements, values, title, cmap='jet', log=False,
          vmin=None, vmax=None):
    polys, vals = [], []
    for i, elem in enumerate(elements):
        if all(n in nodes for n in elem):
            polys.append([nodes[n] for n in elem])
            vals.append(values[i] if i < len(values) else 0.0)
    vals = np.array(vals)
    plot_vals = np.log10(vals + 1e-12) if log else vals
    pc = PolyCollection(polys, cmap=cmap, edgecolors='black', linewidths=0.1)
    pc.set_array(plot_vals)
    if vmin is not None or vmax is not None:
        pc.set_clim(vmin, vmax)
    ax.add_collection(pc)
    cb = plt.colorbar(pc, ax=ax)
    cb.set_label(("log10 " if log else "") + title)
    ax.set_aspect('equal', 'box')
    ax.set_title(title)
    ax.set_xlabel("X (m)"); ax.set_ylabel("Y (m)")
    ax.autoscale()
    return pc


if __name__ == "__main__":
    args = sys.argv[1:]
    
    log = '--log' in args
    save_path = None
    if '--save' in args:
        try:
            save_path = args[args.index('--save') + 1]
        except IndexError:
            pass
            
    positional_args = []
    skip = False
    for arg in args:
        if skip:
            skip = False
            continue
        if arg == '--log':
            continue
        if arg == '--save':
            skip = True
            continue
        positional_args.append(arg)

    if len(positional_args) < 3:
        print(__doc__)
        sys.exit(1)

    mesh_path = positional_args[0]
    ref_file = positional_args[1]
    fwd_file = positional_args[2]

    print(f"Maillage : {mesh_path}")
    nodes, elements = read_2dm(mesh_path)
    ncells = len(elements)
    print(f"--> {len(nodes)} noeuds, {ncells} elements")

    print(f"Comparaison directe :\n  Ref = {ref_file}\n  Sim = {fwd_file}")
    
    if not os.path.exists(ref_file):
        print(f"Erreur : Le fichier de référence n'existe pas : {ref_file}")
        sys.exit(1)
    if not os.path.exists(fwd_file):
        print(f"Erreur : Le fichier simulé n'existe pas : {fwd_file}")
        sys.exit(1)

    h_o, qx_o, qy_o = read_frame_hq(ref_file, ncells)
    h_s, qx_s, qy_s = read_frame_hq(fwd_file, ncells)

    cost = (h_s - h_o)

    print(f"--> J total = {cost.sum():.6e}   max/cellule = {cost.max():.6e} (cellule {int(cost.argmax())})")

    if save_path:
        with open(save_path, 'w') as f:
            f.write("cell_id,cost\n")
            for i in range(ncells):
                f.write(f"{i},{cost[i]:.10e}\n")
        print(f"--> cout par cellule ecrit dans {save_path}")

    fig, ax = plt.subplots(figsize=(11, 8))
    base_ref = os.path.basename(ref_file)
    panel(ax, nodes, elements, cost, f"Cout par cellule (misfit) - {base_ref}", log=log)

    plt.tight_layout()
    out_png = "cost_mesh_single.png"
    plt.show()
