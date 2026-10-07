#!/usr/bin/env python3
"""usage : python3 utils/python/plot/compare_h_phimaps.py h_opt.txt h_ref.txt --mesh mesh.2dm
        --label1 "optimisee" --label2 "reference" --title "t=4200s" --out fig.png
Compare 2 champs de hauteur d'eau h (meme maillage, meme instant, seule la carte de
porosite change), coloriage par CELLULE (polygones du .2dm). 4 panneaux :

    [ h(map1) | h(map2) | dh absolu (m) | dh relatif (%) ]

- h : palette Blues, plage [0, hmax] (hmax = 99.5e percentile, option --hmax).
- dh absolu  = h1 - h2  (m),  RdBu_r symetrique +/- amax (99.5e pct de |dh|).
- dh relatif = (h1 - h2)/h2 * 100  (%),  masque les cellules seches (h2 < --hwet),
  echelle symetrique +/- --relcap  ->  montre que l'erreur ne vaut qu'un faible % de la
  profondeur locale la ou h est grand (evite l'effet "bleu fonce = grosse erreur").

  python3 compare_h_phimaps.py h_opt.txt h_ref.txt --mesh Gardon_sp.2dm \
      --label1 "phi optimisee" --label2 "phi reference" --title "Gardon t=4200s" --out fig.png
"""
import argparse, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

def read_2dm(path):
    """Lit un maillage SMS .2dm -> (nodes {id:(x,y)}, elements [[n1,n2,n3(,n4)], ...])."""
    nodes={}; elems=[]
    with open(path) as f:
        for l in f:
            p=l.split()
            if not p: continue
            if p[0]=="ND":
                nodes[int(p[1])]=(float(p[2]), float(p[3]))
            elif p[0]=="E3T":
                elems.append([int(p[2]), int(p[3]), int(p[4])])
            elif p[0]=="E4Q":
                elems.append([int(p[2]), int(p[3]), int(p[4]), int(p[5])])
    return nodes, elems

def read_h(path):
    """h par cellule (ordre du fichier ResultHydro = ordre des elements du .2dm)."""
    H=[]
    with open(path) as f:
        hdr=None
        for l in f:
            p=l.split()
            if not p: continue
            if 'x' in p and 'y' in p: hdr=p; ih=hdr.index('water_depth'); continue
            if hdr is None: continue
            try: H.append(float(p[ih]))
            except ValueError: pass
    return np.array(H)

ap=argparse.ArgumentParser()
ap.add_argument("h1", help="champ h de la carte 1 (ex: porosite optimisee)")
ap.add_argument("h2", help="champ h de la carte 2 (ex: reference / verite)")
ap.add_argument("--mesh",required=True,help="maillage .2dm (fournit les polygones de cellules)")
ap.add_argument("--label1",default="carte 1"); ap.add_argument("--label2",default="carte 2")
ap.add_argument("--short1",default="h1",help="nom court de la carte 1 dans les panneaux d'erreur")
ap.add_argument("--short2",default="h2",help="nom court de la carte 2 dans les panneaux d'erreur")
ap.add_argument("--title",default=""); ap.add_argument("--out",default="compare_h.png")
ap.add_argument("--dpi",type=int,default=140)
ap.add_argument("--hmax",type=float,default=None,help="borne haute de l'echelle h (defaut: 99.5e pct)")
ap.add_argument("--relcap",type=float,default=30.0,help="borne (%) de l'echelle d'erreur relative")
ap.add_argument("--hwet",type=float,default=0.05,help="seuil h_ref (m) sous lequel l'erreur relative est masquee")
ap.add_argument("--no-rel",dest="no_rel",action="store_true",
                help="ne PAS tracer le panneau d'erreur relative (3 panneaux : h1 | h2 | difference)")
ap.add_argument("--dhmax",type=float,default=None,
                help="borne (m) de l'echelle de difference (defaut: 99.5e pct de |dh|). "
                     "Ex: --dhmax 5 -> meme plage de couleurs que la hauteur d'eau")
a=ap.parse_args()

nodes,elems=read_2dm(a.mesh)
h1=read_h(a.h1); h2=read_h(a.h2)

polys=[[nodes[n] for n in e if n in nodes] for e in elems]
n=min(len(polys),len(h1),len(h2))
polys=polys[:n]; h1=h1[:n]; h2=h2[:n]

d      = h1 - h2
wet1=int((h1>1e-3).sum()); wet2=int((h2>1e-3).sum())
nrm2=float((d**2).sum())

wetmask = h2 > a.hwet
drel = np.full_like(d, np.nan)
drel[wetmask] = d[wetmask] / h2[wetmask] * 100.0
mape = float(np.mean(np.abs(drel[wetmask]))) if wetmask.any() else 0.0

hmax = a.hmax if a.hmax is not None else float(np.percentile(np.concatenate([h1,h2]),99.5))
amax = float(np.percentile(np.abs(d),99.5)) or (float(np.abs(d).max()) or 1e-9)
relcap = a.relcap

s1, s2 = a.short1, a.short2
dmax = a.dhmax if a.dhmax is not None else amax

if a.no_rel:
    fig,axg=plt.subplots(1,3,figsize=(21,6.2),constrained_layout=True)
else:
    fig,axg=plt.subplots(2,2,figsize=(14,11.5),constrained_layout=True)
ax=axg.ravel()
fig.suptitle(f"Comparaison de la hauteur d'eau — {a.title}",fontsize=15,fontweight='bold')
def draw(a_,vals,cmap,vmin,vmax):
    pc=PolyCollection(polys,array=np.asarray(vals,dtype=float),cmap=cmap,
                      edgecolors='face',linewidths=0.0)
    pc.set_clim(vmin,vmax); a_.add_collection(pc); a_.autoscale_view(); a_.set_aspect('equal')
    return pc

p0=draw(ax[0],np.clip(h1,0,hmax),'Blues',0,hmax)
ax[0].set_title(a.label1,fontsize=12,fontweight='bold')
fig.colorbar(p0,ax=ax[0],shrink=.85,label="h (m)")
p1=draw(ax[1],np.clip(h2,0,hmax),'Blues',0,hmax)
ax[1].set_title(a.label2,fontsize=12,fontweight='bold')
fig.colorbar(p1,ax=ax[1],shrink=.85,label="h (m)")
p2=draw(ax[2],d,'RdBu_r',-dmax,dmax)
ax[2].set_title(f"Difference  ({s1} $-$ {s2})   "
                f"(max|dh|={np.abs(d).max():.2f} m,  echelle $\\pm${dmax:g} m)",
                fontsize=12,fontweight='bold')
fig.colorbar(p2,ax=ax[2],shrink=.85,label="dh (m)",extend='both')
if not a.no_rel:
    cmap_rel=plt.cm.RdBu_r.copy(); cmap_rel.set_bad('white')
    p3=draw(ax[3],drel,cmap_rel,-relcap,relcap)
    ax[3].set_title(f"erreur relative  $({s1}-{s2})/{s2}$   "
                    f"(MAPE={mape:.1f}%,  masque $h_{{ref}}$<{a.hwet:g} m)",fontsize=12,fontweight='bold')
    fig.colorbar(p3,ax=ax[3],shrink=.85,label="erreur relative (%)",extend='both')

for a_ in ax: a_.set_xlabel("x"); a_.set_ylabel("y")
fig.savefig(a.out,dpi=a.dpi,bbox_inches='tight')
print(f"SAVED {a.out} | wet1={wet1} wet2={wet2} | ||dh||^2={nrm2:.4g} | "
      f"hmax={hmax:.3g} dmax={dmax:.3g}")
