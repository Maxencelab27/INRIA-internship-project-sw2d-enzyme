import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.collections import PolyCollection
import numpy as np
import sys
import csv
import os
import re

def read_2dm(file_path):
    nodes = {}
    elements = []

    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue

            if parts[0] == "ND":
                node_id = int(parts[1])
                x = float(parts[2])
                y = float(parts[3])
                nodes[node_id] = (x, y)

            elif parts[0] == "E3T":
                _, elem_id, n1, n2, n3, *_ = parts
                elements.append([int(n1), int(n2), int(n3)])

            elif parts[0] == "E4Q":
                _, elem_id, n1, n2, n3, n4, *_ = parts
                elements.append([int(n1), int(n2), int(n3), int(n4)])

    return nodes, elements


def _read_sp_map(file_path, n_elements=None):
    """
    Lit une carte de porosite au format SP_cell_porosity_map.txt (Unif/Distrib, 1 colonne) :
        Unif <0|1>
        ==== Default Param
        <nom>           (ex: Phi, SpPhiW)
        <valeur defaut>
        ==== Distrib    (uniquement si Unif 0)
        <v0> <v1> ...   (1 valeur par ligne, dans l'ordre des elements)
    Retourne {index_element: valeur}. Pour Unif 1, remplit n_elements cellules.
    """
    values = {}
    with open(file_path, 'r', encoding='utf-8') as f:
        raw = f.readlines()
    toks0 = raw[0].split() if raw else []
    unif = int(toks0[1]) if len(toks0) >= 2 else 0
    if unif == 1:
        try:
            v = float(raw[3].strip())
        except (IndexError, ValueError):
            v = 0.0
        if n_elements:
            for i in range(n_elements):
                values[i] = v
    else:
        started = False
        idx = 0
        for line in raw:
            s = line.strip()
            if not started:
                if 'distrib' in s.lower():
                    started = True
                continue
            if s == '':
                continue
            try:
                values[idx] = float(s)
                idx += 1
            except ValueError:
                continue
    return values


def read_data(file_path, param_type, n_elements=None):
    """
    Lit le fichier de données associé et retourne {index_cellule: valeur}.
    phi/grad/grad_signe : CSV "cell_id,phi[,gradient]" OU format SP Unif/Distrib (auto-detecte).
    h/q/z               : fichier ResultHydro (TXT).
    """
    values = {}

    if param_type in ['phi', 'grad', 'grad_signe', 'cost']:
        with open(file_path, 'r', encoding='utf-8') as f:
            first = ''
            for line in f:
                if line.strip():
                    first = line.strip()
                    break
        if first.lower().startswith('unif'):
            return _read_sp_map(file_path, n_elements)

        with open(file_path, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            header = next(reader, None)
            
            val_idx = 1
            if header:
                if param_type in ['grad', 'grad_signe'] and 'gradient' in header:
                    val_idx = header.index('gradient')
                elif param_type == 'phi' and 'phi' in header:
                    val_idx = header.index('phi')
                elif param_type in ['grad', 'grad_signe']:
                    val_idx = 2
            
            for row in reader:
                if not row or len(row) <= val_idx:
                    continue
                try:
                    cell_id = int(row[0])
                    val = float(row[val_idx])
                    values[cell_id] = val
                except ValueError:
                    continue
    else:
        if param_type in ['z', 'zb']:
            target_col = 'zb'
            default_idx = 2
        elif param_type == 'h':
            target_col = 'water_depth'
            default_idx = 3
        else:
            target_col = 'unit_discharge_norm'
            default_idx = 6
        
        with open(file_path, 'r', encoding='utf-8') as f:
            header_line = None
            col_idx = None
            row_counter = 0
            
            for line in f:
                parts = line.strip().split()
                if not parts:
                    continue
                
                if 'x' in parts and 'y' in parts:
                    header_line = parts
                    if target_col in header_line:
                        col_idx = header_line.index(target_col)
                    else:
                        col_idx = default_idx
                    continue
                
                if header_line is None:
                    continue
                
                try:
                    if col_idx is not None and col_idx < len(parts):
                        val = float(parts[col_idx])
                        values[row_counter] = val
                        row_counter += 1
                except ValueError:
                    continue
                    
    return values


def plot_mesh(nodes, elements, values, param_type, file_path=None,
              save_path=None, show=True, dpi=130):
    fig, ax = plt.subplots(figsize=(11, 8))

    vmin, vmax = None, None
    norm = None

    if param_type == 'zb':
        param_type = 'z'

    if param_type == 'h':
        cmap = 'Blues'
        label_str = "Hauteur d'eau (h)"
        vmin_h = 1e-4
        vmax_h = 0.2
    elif param_type == 'q':
        cmap = 'Purples'
        label_str = "Débit unitaire (q)"
    elif param_type == 'z':
        cmap = 'terrain'
        label_str = "Topographie (z)"
    elif param_type == 'phi':
        cmap = 'gray'
        label_str = "Porosité (phi)"
        vmin, vmax = 0, 1.0
    elif param_type == 'grad':
        cmap = 'jet'
        label_str = "Gradient"
    elif param_type == 'grad_signe':
        cmap = mcolors.ListedColormap(['red', 'white', 'green'])
        norm = mcolors.BoundaryNorm([-1.5, -0.5, 0.5, 1.5], cmap.N)
        label_str = "Signe du Gradient"
    elif param_type == 'cost':
        cmap = 'YlOrRd'
        label_str = "Poids du cout (w)"
    else:
        cmap = 'viridis'
        label_str = param_type

    polygons = []
    poly_values = []

    for i, elem in enumerate(elements):
        if all(n in nodes for n in elem):
            pts = [nodes[n] for n in elem]
            polygons.append(pts)
            
            val = values.get(i, 0.0)
            poly_values.append(val)

    if param_type == 'grad_signe':
        mapped_values = []
        for v in poly_values:
            if v > 1e-12:
                mapped_values.append(1)
            elif v < -1e-12:
                mapped_values.append(-1)
            else:
                mapped_values.append(0)
        poly_values = mapped_values

    if param_type == 'h':
        arr = np.asarray(poly_values, dtype=float)
        h_dry = 1e-4
        arr = np.where(arr < h_dry, np.nan, arr)
        poly_values = arr
        finite = arr[np.isfinite(arr)]
        if finite.size:
            env_hmax = float(os.environ.get('SW2D_HMAX', 0) or 0)
            vmax_h = env_hmax if env_hmax > 0 else float(np.percentile(finite, 99))
            vmin_h = max(h_dry, float(finite.min()))
            if vmax_h > vmin_h:
                norm = mcolors.PowerNorm(gamma=0.4, vmin=vmin_h, vmax=vmax_h)
        cmap = plt.get_cmap('Blues').copy()
        cmap.set_bad(color='#e8e8e8')

    if polygons:
        pc = PolyCollection(polygons, cmap=cmap, norm=norm, edgecolors='black', linewidths=0.2)
        pc.set_array(np.array(poly_values))
        
        if vmin is not None and vmax is not None:
            pc.set_clim(vmin, vmax)
            
        ax.add_collection(pc)

        if param_type == 'grad_signe':
            cbar = fig.colorbar(pc, ax=ax, ticks=[-1, 0, 1])
            cbar.ax.set_yticklabels(['Négatif', 'Nul', 'Positif'])
        else:
            cbar = fig.colorbar(pc, ax=ax)
        cbar.set_label(label_str, fontsize=12)
    else:
        print("Erreur : Aucun élément valide n'a pu être construit.")

    title_suffix = ""
    if file_path:
        filename = os.path.basename(file_path)
        match = re.search(r'_p(\d{4})d(\d{2})h(\d{2})m(\d{2})s(\d{3})', filename)
        if match:
            days, hours, minutes, seconds, ms = match.groups()
            d_int = int(days)
            if d_int > 0:
                title_suffix = f" - t = {d_int}d {hours}h {minutes}m {seconds}s"
            else:
                title_suffix = f" - t = {hours}h {minutes}m {seconds}s"

    ax.set_aspect('equal', 'box')
    ax.set_title(f"Champ de - {label_str}{title_suffix}", fontsize=14, fontweight='bold')
    ax.set_xlabel("X (m)", fontsize=11)
    ax.set_ylabel("Y (m)", fontsize=11)
    
    ax.autoscale()
    plt.grid(True, linestyle="--", alpha=0.5)
    
    if save_path:
        out_path = os.path.abspath(save_path)
        out_dir = os.path.dirname(out_path)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        fig.savefig(out_path, dpi=dpi, bbox_inches='tight')
        print(f"--> Figure sauvegardée : {out_path}")

    if show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(
        description="Trace un champ (phi, h, q, z, gradient, cost) sur un maillage .2dm.",
        epilog="Exemples :\n"
               "  python3 mesh_all_phi.py Gardon_sp.2dm phi optimized_porosity.txt\n"
               "  python3 mesh_all_phi.py Gardon_sp.2dm z ResultHydro.txt\n"
               "  python3 mesh_all_phi.py Gardon_sp.2dm porosity.txt --out /tmp/phi.png --no-show\n"
               "  python3 mesh_all_phi.py Gardon_sp.2dm cost input/w_h.csv --save",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mesh", help="fichier maillage .2dm")
    ap.add_argument("param", help="phi | grad | grad_signe | cost | h | q | z | zb")
    ap.add_argument("data", help="fichier de données associé")
    ap.add_argument("--out", default=None,
                    help="chemin PNG de sortie (sauvegarde directe, sans prompt)")
    ap.add_argument("--save", action="store_true",
                    help="sauvegarde dans /home/chapumax/Documents/résultats avec un nom auto")
    ap.add_argument("--no-show", action="store_true",
                    help="ne pas ouvrir la fenêtre (mode batch)")
    ap.add_argument("--dpi", type=int, default=130)
    args = ap.parse_args()

    file_2dm = args.mesh
    param_type = args.param.lower()
    data_file = args.data

    if param_type not in ['phi', 'grad', 'grad_signe', 'cost', 'h', 'q', 'z', 'zb']:
        print("Erreur : le paramètre doit être 'phi', 'grad', 'grad_signe', 'cost', 'h', 'q', 'z' ou 'zb'.")
        sys.exit(1)

    save_path = None
    if args.out:
        save_path = args.out
    elif args.save:
        base_name = os.path.splitext(os.path.basename(data_file))[0]
        out_param_name = 'z' if param_type in ['z', 'zb'] else param_type
        save_path = os.path.join("/home/chapumax/Documents/résultats",
                                 f"{base_name}_{out_param_name}.png")

    print(f"Lecture du maillage : {file_2dm} ...")
    nodes, elements = read_2dm(file_2dm)
    print(f"--> {len(nodes)} nœuds et {len(elements)} éléments chargés.")

    print(f"Lecture des données ({param_type}) : {data_file} ...")
    values = read_data(data_file, param_type, len(elements))
    print(f"--> {len(values)} valeurs lues.")

    print("Affichage du maillage coloré...")
    plot_mesh(nodes, elements, values, param_type, data_file,
              save_path=save_path, show=not args.no_show, dpi=args.dpi)
