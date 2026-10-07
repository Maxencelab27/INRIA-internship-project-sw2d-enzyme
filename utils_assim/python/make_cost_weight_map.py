#!/usr/bin/env python3
"""
make_cost_weight_map.py — cartes de ponderation du cout (data-driven, levier 3).

Construit, par cellule, des poids egaux au CARRE de la difference entre DEUX cartes
ResultHydro PRISES AU MEME INSTANT : l'EBAUCHE (output/, φ=1) et la VERITE (reference/, φ=0).

    w_h[i] = ( h_ref(i)  - h_ebauche(i)  )^2
    w_q[i] = ( qx_ref - qx_eb )^2 + ( qy_ref - qy_eb )^2

Forts la ou verite et ebauche divergent a cet instant = zone influencee par le barrage
(retenue + ouvrage + aval) -> le cout "regarde" le barrage sans presupposer sa geometrie.

Sortie : <out-dir>/w_h.csv et <out-dir>/w_q.csv (format cell_id,phi, relisible par
h_weight / q_weight via load_phi_csv).

Usage :
    python3 make_cost_weight_map.py <ebauche.txt> <reference.txt> [options]
      <ebauche.txt>   : ex. sp_gardon/output/ResultHydro_p0000d01h30m00s000.txt   (φ=1)
      <reference.txt> : ex. sp_gardon/reference/ResultHydro_p0000d01h30m00s000.txt (φ=0)
    options :
      --out-dir DIR : ou ecrire w_h.csv / w_q.csv (defaut : <dossier_simu>/input,
                      deduit du chemin de la reference ; sinon dossier courant)
      --norm mean|max|none : normalisation (defaut mean -> moyenne=1, garde l'echelle de J)
      --floor F     : poids plancher ajoute partout (defaut 0.0)
"""
import sys, os, glob, argparse
import numpy as np

def read_2dm(file_path):
    """Lit un maillage .2dm -> (nodes{id:(x,y)}, elements[[n1,n2,n3(,n4)]])."""
    nodes, elements = {}, []
    with open(file_path) as f:
        for line in f:
            p = line.split()
            if not p: continue
            if p[0] == "ND":
                nodes[int(p[1])] = (float(p[2]), float(p[3]))
            elif p[0] == "E3T":
                elements.append([int(p[2]), int(p[3]), int(p[4])])
            elif p[0] == "E4Q":
                elements.append([int(p[2]), int(p[3]), int(p[4]), int(p[5])])
    return nodes, elements

def plot_weight_on_mesh(w, label):
    """Trace la carte arc-en-ciel (jet) des poids w (par cellule) sur le maillage
    auto-detecte dans ./primary_files/*.2dm (dossier courant)."""
    meshes = glob.glob(os.path.join("primary_files", "*.2dm"))
    if not meshes:
        print("  [plot] aucun maillage ./primary_files/*.2dm trouve."); return
    mesh = meshes[0]
    print(f"  [plot] maillage : {mesh}")
    try:
        import matplotlib.pyplot as plt
        from matplotlib.collections import PolyCollection
    except ImportError:
        print("  [plot] matplotlib indisponible."); return
    nodes, elements = read_2dm(mesh)
    polys, vals = [], []
    for i, elem in enumerate(elements):
        if all(n in nodes for n in elem):
            polys.append([nodes[n] for n in elem])
            vals.append(w[i] if i < len(w) else 0.0)
    fig, ax = plt.subplots(figsize=(11, 8))
    pc = PolyCollection(polys, cmap="jet", edgecolors="black", linewidths=0.2)
    pc.set_array(np.array(vals))
    ax.add_collection(pc)
    fig.colorbar(pc, ax=ax).set_label(label, fontsize=12)
    ax.set_aspect("equal", "box"); ax.autoscale()
    ax.set_title(f"Carte de cout - {label}", fontsize=14, fontweight="bold")
    ax.set_xlabel("X (m)"); ax.set_ylabel("Y (m)")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.show()

def read_hydro(path):
    """Retourne (h[], qx[], qy[]) en ordre de cellule depuis un ResultHydro_*.txt."""
    h, qx, qy = [], [], []
    ih = iqx = iqy = None
    with open(path) as f:
        hdr = None
        for line in f:
            p = line.split()
            if not p: continue
            if 'x' in p and 'y' in p:
                hdr = p
                ih  = hdr.index('water_depth')     if 'water_depth'     in hdr else 3
                iqx = hdr.index('unit_discharge_x') if 'unit_discharge_x' in hdr else 4
                iqy = hdr.index('unit_discharge_y') if 'unit_discharge_y' in hdr else 5
                continue
            if hdr is None: continue
            try:
                h.append(float(p[ih])); qx.append(float(p[iqx])); qy.append(float(p[iqy]))
            except (ValueError, IndexError):
                pass
    return h, qx, qy

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ebauche", help="ResultHydro de l'ebauche (φ=1, output/)")
    ap.add_argument("reference", help="ResultHydro de la verite (φ=0, reference/), MEME instant")
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--norm", choices=["mean", "max", "none"], default="mean")
    ap.add_argument("--floor", type=float, default=0.0)
    a = ap.parse_args()

    for f in (a.ebauche, a.reference):
        if not os.path.isfile(f):
            print(f"ERREUR : fichier introuvable {f}"); sys.exit(1)

    out_dir = a.out_dir if a.out_dir is not None else os.path.join(os.getcwd(), "input")
    os.makedirs(out_dir, exist_ok=True)

    he, qxe, qye = read_hydro(a.ebauche)
    hr, qxr, qyr = read_hydro(a.reference)
    m = min(len(he), len(hr))
    if m == 0:
        print("ERREUR : aucune donnee lue (verifier le format ResultHydro)."); sys.exit(1)
    if len(he) != len(hr):
        print(f"ATTENTION : nb cellules different (ebauche={len(he)}, ref={len(hr)}) -> on prend {m}.")

    wh = [np.exp(hr[i] - he[i]) for i in range(m)]
    wq = [((qxr[i]-qxe[i])**2 + (qyr[i]-qye[i])**2) for i in range(m)]

    def normalize(w):
        w = [x + a.floor for x in w]
        if a.norm == "none": return w
        ref_val = (sum(w)/len(w)) if a.norm == "mean" else max(w)
        if ref_val <= 0: return [1.0]*len(w)
        return [x/ref_val for x in w]


    def write_csv(path, w, label):
        with open(path, "w") as f:
            f.write("cell_id,phi\n")
            for i, x in enumerate(w):
                f.write(f"{i},{x:.10g}\n")
        nz = sum(1 for x in w if x > 1e-9)
        print(f"  {label} -> {path}  ({len(w)} cellules, {nz} non nulles, max={max(w):.3g}, moy={sum(w)/len(w):.3g})")

    print(f"Instant : ebauche={os.path.basename(a.ebauche)}  vs  reference={os.path.basename(a.reference)}")
    print("Cartes de poids (carre de la difference a cet instant) :")
    write_csv(os.path.join(out_dir, "w_h.csv"), wh, "w_h")
    write_csv(os.path.join(out_dir, "w_q.csv"), wq, "w_q")
    print("\nDans input.solv :")
    print("  h_weight   input/w_h.csv")
    print("  q_weight   input/w_q.csv   ")

    ans = input("\nTracer la carte arc-en-ciel des couts sur le maillage ? [h/q/n] : ").strip().lower()
    if ans in ("h", "y", "yes", "o", "oui"):
        plot_weight_on_mesh(wh, "w_h (cout pondere h)")
    elif ans == "q":
        plot_weight_on_mesh(wq, "w_q (cout pondere q)")

if __name__ == "__main__":
    main()
