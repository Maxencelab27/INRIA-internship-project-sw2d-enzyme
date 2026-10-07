#!/usr/bin/env python3
"""Courbe d'identifiabilite temporelle d'une structure (barrage, ville, ...) :
innovation |reference - ebauche| agregee sur les mailles de la structure, en fonction
du temps. Sert a REGLER LA FENETRE d'assimilation sur le pic d'identifiabilite :
la ou l'innovation est forte, le parametre (porosite) est observable ; ailleurs il ne
l'est pas (signal nul avant l'arrivee de l'onde, etc.).

Idee (data-driven, meme logique que make_weight_maps.py) :
    innov_h(t) = agg_{i dans structure} | h_ref(i,t) - h_eb(i,t) |
    innov_q(t) = agg_{i dans structure} | q_ref(i,t) - q_eb(i,t) |   (norme (dqx,dqy))
La "structure" est definie par la carte verite (cell_id,phi) : mailles phi < --dam-thr.
Sans --dam, l'innovation est agregee sur TOUTES les mailles (identifiabilite globale).

La fenetre observable retournee = intervalle des t ou innov_h >= --window-frac * pic.

Prerequis (deux dossiers de ResultHydro aux MEMES instants) :
  reference/ : sorties verite            (ex.  sw2d_enzyme . gen-ref-e)
  output/    : sorties ebauche phi=ones   (      sw2d_enzyme . forward)
Format ResultHydro : entete + 1 ligne/maille (x, y, zb, h, qx, qy, |q|).
Temps encode dans le nom ...dHHhMMmSSsMMM : d=heures, h=minutes, m=secondes, s=ms.

Exemples :
  # barrage Gardon, fenetre >=50% du pic, trace h ET q :
  python3 identifiability_plot.py --ref-dir sp_gardon_travail/reference --eb-dir sp_gardon_travail/output

  # Ne tracer QUE la hauteur d'eau h :
  python3 identifiability_plot.py --ref-dir sp_gardon_travail/reference --eb-dir sp_gardon_travail/output h

  # Ne tracer QUE le débit q :
  python3 identifiability_plot.py --ref-dir sp_gardon_travail/reference --eb-dir sp_gardon_travail/output q
"""
import argparse, glob, math, os, re, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TIME_RE = re.compile(r"d(\d+)h(\d+)m(\d+)s(\d+)")


def time_from_name(path):
    """ResultHydro_..dHHhMMmSSsMMM -> secondes (d=heures, h=minutes, m=secondes, s=ms)."""
    m = TIME_RE.search(os.path.basename(path))
    if not m:
        return None
    hh, mm, ss, ms = (int(x) for x in m.groups())
    return hh * 3600 + mm * 60 + ss + ms / 1000.0


def index_dir(d):
    out = {}
    for p in glob.glob(os.path.join(d, "ResultHydro_*.txt")):
        t = time_from_name(p)
        if t is not None:
            out[round(t, 3)] = p
    return out


def load_hydro(path):
    """Lit un ResultHydro -> (h[], qx[], qy[]) par maille (ordre du fichier = cell_id)."""
    h, qx, qy = [], [], []
    with open(path) as f:
        next(f)
        for line in f:
            c = line.split()
            if len(c) < 6:
                continue
            h.append(float(c[3])); qx.append(float(c[4])); qy.append(float(c[5]))
    return h, qx, qy


def read_structure(path, thr):
    """carte verite (cell_id,phi) -> set des cell_id de la structure (phi < thr)."""
    cells = set()
    with open(path) as f:
        for line in f:
            c = line.replace(",", " ").split()
            if len(c) < 2:
                continue
            try:
                i, v = int(c[0]), float(c[1])
            except ValueError:
                continue
            if v < thr:
                cells.add(i)
    return cells


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref-dir", default="reference", help="ResultHydro de reference (verite)")
    ap.add_argument("--eb-dir", default="output", help="ResultHydro d'ebauche (forward phi=ones)")
    ap.add_argument("--dam", default=None,
                    help="carte verite (cell_id,phi) : agrege l'innovation sur les mailles phi<--dam-thr")
    ap.add_argument("--dam-thr", type=float, default=0.5, help="seuil phi de la structure (defaut 0.5)")
    ap.add_argument("--agg", choices=["mean", "max", "sum"], default="mean",
                    help="agregation spatiale de l'innovation sur la structure (defaut mean)")
    ap.add_argument("--window-frac", type=float, default=0.5,
                    help="fraction du pic definissant la fenetre observable (defaut 0.5)")
    ap.add_argument("--mark-peak", type=float, default=None, help="instant a marquer (ex. pic hydrogramme)")
    ap.add_argument("--mark-old", type=float, default=None, help="instant a marquer (ex. ancien resultat)")
    ap.add_argument("--title", default="Identifiabilite temporelle")
    ap.add_argument("--out", default="identifiability.png")
    ap.add_argument("--dpi", type=int, default=140)
    ap.add_argument("--csv", default=None, help="ecrit aussi la serie temporelle en CSV")
    
    ap.add_argument("variable", nargs="?", choices=["h", "q"], default=None,
                    help="Variable à tracer : 'h' (uniquement la hauteur d'eau), 'q' (uniquement le débit) ou les deux si non spécifié")
    
    a = ap.parse_args()

    ref_idx, eb_idx = index_dir(a.ref_dir), index_dir(a.eb_dir)
    if not ref_idx:
        sys.exit(f"[Err] aucun ResultHydro dans {a.ref_dir}")
    if not eb_idx:
        sys.exit(f"[Err] aucun ResultHydro dans {a.eb_dir} (genere l'ebauche : sw2d_enzyme . forward)")
    times = sorted(t for t in ref_idx if t in eb_idx)
    if not times:
        sys.exit("[Err] aucun instant commun entre reference et ebauche.")

    struct = read_structure(a.dam, a.dam_thr) if a.dam else None
    if struct is not None:
        print(f"[info] structure = {len(struct)} mailles (phi<{a.dam_thr}) depuis {a.dam}")
    print(f"[info] {len(times)} instants communs : {times[0]} .. {times[-1]} s (agg={a.agg})")

    def aggregate(vals, sel):
        xs = [vals[i] for i in sel] if sel is not None else vals
        if not xs:
            return 0.0
        if a.agg == "max":
            return max(xs)
        if a.agg == "sum":
            return sum(xs)
        return sum(xs) / len(xs)

    th, tq = [], []
    for t in times:
        hr, qxr, qyr = load_hydro(ref_idx[t])
        he, qxe, qye = load_hydro(eb_idx[t])
        n = min(len(hr), len(he))
        sel = [i for i in struct if i < n] if struct is not None else None
        ih = [abs(hr[i] - he[i]) for i in range(n)]
        iq = [math.hypot(qxr[i] - qxe[i], qyr[i] - qye[i]) for i in range(n)]
        th.append(aggregate(ih, sel)); tq.append(aggregate(iq, sel))

    plot_h = a.variable in [None, "h"]
    plot_q = a.variable in [None, "q"]

    active_series = tq if a.variable == "q" else th
    active_label = "q" if a.variable == "q" else "h"

    peak = max(active_series) if active_series else 0.0
    thr = a.window_frac * peak
    obs = [times[k] for k in range(len(times)) if active_series[k] >= thr]
    win = (obs[0], obs[-1]) if obs else None
    if win:
        print(f"[fenetre] innovation {active_label} >= {a.window_frac:.0%} du pic sur [{win[0]:.0f}, {win[1]:.0f}] s "
              f"(pic={peak:.3g} a t={times[active_series.index(peak)]:.0f} s)")

    if a.csv:
        with open(a.csv, "w") as f:
            f.write("t_s,innov_h,innov_q\n")
            for t, xh, xq in zip(times, th, tq):
                f.write(f"{t:.3f},{xh:.6g},{xq:.6g}\n")
        print("[csv]", a.csv)

    fig, ax = plt.subplots(figsize=(13, 6.2), constrained_layout=True)
    lbl = "au barrage" if a.dam else "global"
    
    if plot_h:
        ax.plot(times, th, "-o", color="#c0392b", ms=4, lw=2,
                label=r"innovation $h$ %s  $|h_{verite}-h_{ebauche,\phi=1}|$" % lbl)
    if plot_q:
        ax.plot(times, tq, "-^", color="#27ae60", ms=4, lw=1.6, label=r"innovation $|q|$ %s" % lbl)
        
    if win:
        ax.axvspan(win[0], win[1], color="#f6c85f", alpha=0.30,
                   label=f"fenetre observable (sur {active_label}) >={a.window_frac:.0%} pic  [{win[0]:.0f},{win[1]:.0f}]s")
                   
    ax.axvline(times[-1], color="#c0392b", ls="--", lw=1.5, label=f"maximum a t={times[-1]:.0f}s (= t_max)")
    
    if a.mark_peak is not None:
        ax.axvline(a.mark_peak, color="#1f3a93", ls=":", lw=2, label=f"pic hydrogramme amont ({a.mark_peak:.0f}s)")
        
    if a.mark_old is not None:
        pct = 100 * (active_series[min(range(len(times)), key=lambda k: abs(times[k] - a.mark_old))] / peak) if peak else 0
        ax.axvline(a.mark_old, color="0.35", ls=":", lw=1.3, label=f"ancien resultat t={a.mark_old:.0f} ({pct:.0f}% du pic {active_label})")
        
    ax.set_xlabel("temps (s)")
    
    ylabel_suffix = " au barrage" if a.dam else ""
    if a.variable == "h":
        ax.set_ylabel(f"innovation h{ylabel_suffix} (m)")
    elif a.variable == "q":
        ax.set_ylabel(f"innovation q{ylabel_suffix} (m²/s)")
    else:
        ax.set_ylabel(f"innovation{ylabel_suffix} (unités mixtes)")
        
    ax.set_title(a.title, fontsize=13, fontweight="bold")
    ax.grid(True, alpha=0.3); ax.legend(loc="upper left", fontsize=9)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    fig.savefig(a.out, dpi=a.dpi, bbox_inches="tight")
    print("SAVED", a.out)


if __name__ == "__main__":
    main()
