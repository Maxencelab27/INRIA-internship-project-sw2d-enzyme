#!/usr/bin/env python3
"""Trace la convergence d'une minimisation de porosite a partir de descent_log.csv.

Figure a EXACTEMENT 3 panneaux (1x3), en fonction de l'iteration :
  (1) J            : fonction cout J_total (echelle log)
  (2) sse_phi_ref  : ecart quadratique a la carte VERITE
  (3) rel_init_pct : pourcentage de descente relatif a J_init (col 5 du CSV,
                     = le 'rel = X %' affiche dans le terminal)

Usage :
  python3 plot_descent.py [descent/descent_log.csv] [-o sortie.png]
"""
import argparse
import csv
import os
import sys
import numpy as np

import matplotlib.pyplot as plt


def load_log(path):
    """Retourne dict {colonne: [valeurs float]} ; ignore les cellules non numeriques."""
    cols = {}
    with open(path) as f:
        reader = csv.reader(f)
        header = next(reader)
        keys = [h.strip() for h in header]
        for h in keys:
            cols[h] = []
        for row in reader:
            if not row or not row[0].strip():
                continue
            for k, v in zip(keys, row):
                try:
                    cols[k].append(float(v))
                except ValueError:
                    cols[k].append(float("nan"))
    return cols


def has(cols, name):
    return name in cols and len(cols[name]) > 0


def main():
    ap = argparse.ArgumentParser(description="Plot de convergence descent_log.csv (J | sse_phi_ref | rel_init_pct)")
    ap.add_argument("csv", nargs="?", default="descent/descent_log.csv",
                    help="chemin du descent_log.csv (defaut: descent/descent_log.csv)")
    ap.add_argument("-o", "--out", default=None,
                    help="image de sortie (defaut: descent_convergence.png a cote du csv)")
    args = ap.parse_args()

    if not os.path.isfile(args.csv):
        sys.exit(f"[Err] fichier introuvable : {args.csv}")

    cols = load_log(args.csv)
    it = cols.get("iter") or list(range(1, len(next(iter(cols.values()))) + 1))
    out = args.out or os.path.join(os.path.dirname(os.path.abspath(args.csv)),
                                   "descent_convergence.png")

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    fig.suptitle(f"Convergence minimisation porosite - {args.csv}", fontsize=13, fontweight='bold')

    bbox_style = dict(boxstyle="round,pad=0.4", fc="#fff9c4", ec="#fbc02d", alpha=0.9, lw=1.5)
    arrow_style = dict(arrowstyle="->", color="black", lw=1.2)

    ax = axes[0]
    if has(cols, "J_total"):
        j_vals = cols["J_total"]
        ax.semilogy(it, j_vals, "o-", color="tab:blue", ms=4, label="J", alpha=0.7)
        
        valid_idx = [i for i, v in enumerate(j_vals) if not np.isnan(v)]
        if valid_idx:
            best_i = min(valid_idx, key=lambda i: j_vals[i])
            best_it = it[best_i]
            best_j = j_vals[best_i]
            
            ax.plot(best_it, best_j, ".", color="#fbc02d", ms=12, mec="black", zorder=5, label="Meilleur")
            ax.annotate(f"Min: {best_j:.4e}\n(it {int(best_it)})",
                        xy=(best_it, best_j), xytext=(20, 20), textcoords="offset points",
                        arrowprops=arrow_style, bbox=bbox_style, fontweight="bold", fontsize=9)

        jpos = [v for v in j_vals if v > 0]
        if jpos:
            ax.set_ylim(bottom=min(jpos) * 0.5, top=j_vals[0] * 2.0)
            
    ax.set_title("Fonction cout J")
    ax.set_xlabel("iteree")
    ax.set_ylabel("J (log)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()

    ax = axes[1]
    if has(cols, "sse_phi_ref"):
        sse = cols["sse_phi_ref"]
        ax.plot(it, sse, "o-", color="tab:red", ms=4, alpha=0.7, label="SSE")
        
        valid_idx = [i for i, v in enumerate(sse) if not np.isnan(v)]
        if valid_idx:
            best_i = min(valid_idx, key=lambda i: sse[i])
            best_it = it[best_i]
            best_sse = sse[best_i]
            
            ax.plot(best_it, best_sse, ".", color="#fbc02d", ms=12, mec="black", zorder=5, label="Meilleur")
            ax.annotate(f"Min: {best_sse:.4e}\n(it {int(best_it)})",
                        xy=(best_it, best_sse), xytext=(20, 20), textcoords="offset points",
                        arrowprops=arrow_style, bbox=bbox_style, fontweight="bold", fontsize=9)
            
        ax.set_ylim(bottom=max(0.0, min(sse) * 0.95), top=sse[0] * 1.5)
        ax.set_ylabel("sse_phi_ref")
    else:
        ax.text(0.5, 0.5, "sse_phi_ref absent", ha="center", va="center")
        
    ax.set_title("SSE porosite vs reference")
    ax.set_xlabel("iteree")
    ax.grid(True, alpha=0.3)
    ax.legend()

    ax = axes[2]
    if has(cols, "rel_init_pct"):
        rel = cols["rel_init_pct"]
        ax.plot(it, rel, "o-", color="tab:green", ms=4, alpha=0.7, label="rel %")
        ax.axhline(0.0, color="0.5", lw=0.8)
        
        valid_idx = [i for i, v in enumerate(rel) if not np.isnan(v)]
        if valid_idx:
            best_i = min(valid_idx, key=lambda i: rel[i])
            best_it = it[best_i]
            best_rel = rel[best_i]
            
            ax.plot(best_it, best_rel, ".", color="#fbc02d", ms=12, mec="black", zorder=5, label="Meilleur")
            ax.annotate(f"Max descente:\n{best_rel:.2f}%\n(it {int(best_it)})",
                        xy=(best_it, best_rel), xytext=(20, -35), textcoords="offset points",
                        arrowprops=arrow_style, bbox=bbox_style, fontweight="bold", fontsize=9)
            
        ax.set_ylabel("rel_init_pct (%)")
    else:
        ax.text(0.5, 0.5, "rel_init_pct absent", ha="center", va="center")
        
    ax.set_title("% de descente vs J_init")
    ax.set_xlabel("iteree")
    ax.grid(True, alpha=0.3)
    ax.legend()

    fig.tight_layout(rect=[0, 0, 1, 0.93])

    plt.show()

    try:
        rep = input(f"Sauvegarder la figure en PNG ? [{out}] (o/N) ").strip().lower()
    except EOFError:
        rep = ""
    if rep in ("o", "oui", "y", "yes"):
        fig.savefig(out, dpi=130)
        print(f"[OK] figure ecrite : {out}")
    else:
        print("[skip] figure non sauvegardee.")


if __name__ == "__main__":
    main()
