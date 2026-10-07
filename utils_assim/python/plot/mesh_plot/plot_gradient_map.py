#!/usr/bin/env python3
"""Trace la carte du gradient dJ/dphi (mode gradient par defaut -> gradient_porosity.csv)
sur le maillage .2dm, coloriage par CELLULE. Montre l'espace nul : le gradient est
non nul seulement sur les cellules OBSERVABLES (front/axes d'ecoulement), ~0 ailleurs.
  python3 plot_gradient_map.py gradient_porosity.csv --mesh ToceAPC.2dm --out grad.png
"""
import argparse, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

def read_2dm(path):
    nodes={}; elems=[]
    with open(path) as f:
        for l in f:
            p=l.split()
            if not p: continue
            if p[0]=="ND": nodes[int(p[1])]=(float(p[2]),float(p[3]))
            elif p[0]=="E3T": elems.append([int(p[2]),int(p[3]),int(p[4])])
            elif p[0]=="E4Q": elems.append([int(p[2]),int(p[3]),int(p[4]),int(p[5])])
    return nodes, elems

ap=argparse.ArgumentParser()
ap.add_argument("csv"); ap.add_argument("--mesh",required=True)
ap.add_argument("--title",default="Gradient dJ/dphi"); ap.add_argument("--out",default="gradient_map.png")
ap.add_argument("--dpi",type=int,default=150)
a=ap.parse_args()

nodes,elems=read_2dm(a.mesh)
g=[]
with open(a.csv) as f:
    hdr=None
    for l in f:
        p=l.strip().split(',')
        if hdr is None: hdr=p; ig=hdr.index('gradient'); continue
        try: g.append(float(p[ig]))
        except (ValueError,IndexError): pass
g=np.array(g)
polys=[[nodes[n] for n in e if n in nodes] for e in elems]
m=min(len(polys),len(g)); polys=polys[:m]; g=g[:m]

nz=int((np.abs(g)>1e-3).sum())
amax=float(np.percentile(np.abs(g),99.8)) or (np.abs(g).max() or 1e-9)

fig,ax=plt.subplots(1,2,figsize=(16,6.5),constrained_layout=True)
pc=PolyCollection(polys,array=g,cmap='RdBu_r',edgecolors='face',linewidths=0.0)
pc.set_clim(-amax,amax); ax[0].add_collection(pc); ax[0].autoscale_view(); ax[0].set_aspect('equal')
ax[0].set_title(f"{a.title}  (signe)\ncellules |g|>1e-3 : {nz} / {m}",fontsize=12,fontweight='bold')
fig.colorbar(pc,ax=ax[0],shrink=.85,label="dJ/dphi")
ag=np.abs(g); ag_plot=np.log10(ag+1e-6)
pc2=PolyCollection(polys,array=ag_plot,cmap='inferno',edgecolors='face',linewidths=0.0)
pc2.set_clim(-6,np.log10(amax+1e-6)); ax[1].add_collection(pc2); ax[1].autoscale_view(); ax[1].set_aspect('equal')
ax[1].set_title("log10 |dJ/dphi|  (espace nul = zones sombres)\n"
                f"mean|g|={ag.mean():.2e}  max|g|={ag.max():.2e}",fontsize=12,fontweight='bold')
fig.colorbar(pc2,ax=ax[1],shrink=.85,label="log10|dJ/dphi|")
for a_ in ax: a_.set_xlabel("x"); a_.set_ylabel("y")
fig.suptitle("Sensibilite de J a la porosite : seules les cellules OBSERVABLES ont un gradient",
             fontsize=13,fontweight='bold')
fig.savefig(a.out,dpi=a.dpi,bbox_inches='tight'); print("SAVED",a.out,"| cellules actives (|g|>1e-3):",nz,"/",m)
