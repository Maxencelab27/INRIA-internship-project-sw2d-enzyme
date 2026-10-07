# Couche d'assimilation / AD Enzyme — découplée du SW2D natif

Ce dossier `assim/` contient **toute** la couche de différentiation automatique
(Enzyme) et d'assimilation de données. Le reste de l'arbre (`../src`, `../cmake`,
`../CMakeLists.txt`) est le **SW2D natif strictement pristine** (commit pré-Enzyme
`4e0199f`) : aucun fichier natif — pas même un header — n'est modifié.

## Principe du découplage

La classe `swSpE` **dérive** de la classe native `swSp` :

```cpp
class swSpE : public swSp { ... };   // assim/swSpE.h
```

C'est possible sans toucher au natif parce que :

- `swSp` n'a **aucune section `private:`** (uniquement `public:` puis `protected:`)
  et son constructeur est public ;
- `swParent` (classe de base) est **entièrement `public:`**.

Une classe dérivée voit donc **tout** l'état natif (`water_depth`,
`unit_discharge_x/y`, `F1_int`, `edge_bnd_factor_*`, `z_riemann_int`,
`bc_val_interpolated`, …) par héritage, et peut appeler les méthodes natives
(`computeMusclSlopes`, `specificForce`, `lambdaMinus/Plus`, `phi`).

Toutes les régularisations Enzyme-safe (vitesse désingularisée `desingVelocity`,
gardes cellule sèche, `min` lissé C¹, recalcul des facteurs `1/phi` dans la tape)
vivent **ici**, jamais dans le natif. Le solveur natif de référence reste identique
à l'amont (les références d'assimilation sont générées par le forward `_E` via le
mode *gen-ref-e*, donc l'identité bit-à-bit classique≡`_E` n'est pas requise).

## Sens de couplage autorisé

- **`_E` → natif** (une méthode `swSpE` appelle du code natif non modifié) : **autorisé**,
  non contaminant.
- **natif ← `_E`** (éditer un fichier natif pour l'AD) : **proscrit**. Rien ici ne le fait.

## Instanciation

Le factory natif `swParent::create` fabrique un `swSp`. La couche assim
(`main_enzyme.cpp`) instancie directement un `swSpE` puis appelle `init()` hérité —
le factory natif n'est pas touché.

## Points de vigilance

- La macro AD `#define phi(i) ...` (accès porosité via `phi_ptr`) doit rester
  **confinée aux unités de compilation `assim/`** et ne jamais fuir dans un `.cpp`
  natif.

## Plan de migration (par étapes)

- **Étape 1 (faite)** : copie pristine + structure `assim/` + squelette `swSpE`
  (méthodes `specificForce_E`, `hypComputeBalance_E`) qui **compile** contre le natif.
- **Étape 2 (à venir)** : migration de toutes les fonctions `_E` en méthodes
  `swSpE::` (hyp/flux/HLLC/BC/sources/checkpoint/minimize) + membres de données AD.
- **Étape 3 (à venir)** : `main_enzyme.cpp` + `CMakeLists.txt` autonome
  (natif compilé en bibliothèque statique, cible assim liée dessus, flags Enzyme) +
  build complet avec la toolchain clang+Enzyme.

## Validation étape 1

```
g++ -std=c++17 -c swSpE.cpp     # OK : la classe dérivée voit l'état protected natif
```
