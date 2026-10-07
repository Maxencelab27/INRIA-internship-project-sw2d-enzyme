# Journal de migration — découplage strict SW2D natif / couche AD (Enzyme)

Ce journal documente, étape par étape, la mise en œuvre du **découplage strict** :
un SW2D natif **pristine** (`../src`, commit pré-Enzyme `4e0199f`) + une couche
d'assimilation/AD entièrement contenue dans `assim/`, **sans aucune modification
d'un fichier natif** (headers compris).

---

## Décisions d'architecture (validées)

- **Classe dérivée** : `class swSpE : public swSp`. Possible sans toucher au natif
  car `swSp` n'a aucune section `private:` (uniquement `public:`/`protected:`) et
  `swParent` est tout `public:`. Une classe dérivée voit donc tout l'état natif par
  héritage.
- **CMake** : `assim/` autonome ; le natif est compilé en bibliothèque, la cible
  assim est liée dessus. Zéro fichier natif touché.
- **Natif strictement pristine** : les gardes « cellule sèche » et régularisations
  vivent uniquement côté `_E`.
- **Livraison par étapes**, validée à chaque palier.

---

## Étape 1 — TERMINÉE

- Copie pristine `sw2d_natif/` → `sw2d_enzyme_release/` (sans `.git`, `test`).
  Vérifié : `grep -rl enzyme src/` = 0.
- Création de `assim/` + squelette `swSpE.{h,cpp}` avec deux méthodes `_E`
  représentatives (`specificForce_E`, `hypComputeBalance_E`).
- **Validation** : `g++ -std=c++17 -c swSpE.cpp` compile — la classe dérivée
  accède bien aux ~12 tableaux **protected** du natif par héritage.

---

## Étape 2 — EN COURS

Objectif : migrer TOUTES les fonctions `_E` en méthodes `swSpE::` (noms de
fonctions inchangés ; noms de fichiers .cpp libres), reloger les ~500 lignes de
déclarations AD depuis l'ancien `swSp.h` vers `swSpE.h`.

### Méthode de reconstruction de `swSpE.h`

Le diff natif↔courant de `swSp.h` mêle natif et AD (blocs `c`), et certains
membres natifs relocalisés réapparaissent « ajoutés » (commentaire retiré). Une
diff ligne-à-ligne est donc dangereuse (risque de **shadowing silencieux** d'un
membre natif). Méthode retenue, en deux temps :

1. **Méthodes AD** : déclarées à la main dans `swSpE.h` (signatures non ambiguës,
   relevées dans l'ancien `swSp.h`). Direction sûre : une omission = erreur de
   compilation (corrigible) ; jamais un shadow.
2. **Membres de données AD** : extraits par différence de lignes exacte
   (courant − natif), puis **classés par présence du nom dans le `swSp.h` natif** :
   - NATIFS (exclus de `swSpE`) : `phi`, `water_depth`, `water_depth_old`, …
   - AD (mis dans `swSpE`) : `dz_dx`, `dz_dy`, `z_cell`, `swept_data`, `phi_ptr`,
     `obs_times`, `h_init`, `qx_init`, `qy_init`, `water_depth_obs`, `bc_flat`,
     `dt_seq`, `obs_trig`, `segs_record_`, `assim_window_tiles_`, tous les
     `alpha_*`, `lbfgsb_*`, `adam_*`, `phi_norm_*`, poids `h/q_weight_*`, etc.

### Règles de transformation des .cpp

- `swSp::foo_E(...)` → `swSpE::foo_E(...)`.
- `swParent::computeBC_E` / `interpolateTimeSeries_E` → méthodes `swSpE` (héritent
  de `swParent`, membres tout `public:`).
- **`ckpt_segment_E(swSp* solver, …)` → `ckpt_segment_E(swSpE* solver, …)`** :
  RÈGLE C++ CRUCIALE — une méthode (même `static`) de `swSpE` ne peut accéder aux
  membres *protected* du natif qu'à travers un `swSpE*`, pas un `swSp*`. Les
  fonctions statiques différenciées par Enzyme doivent donc prendre `swSpE*`.
- Includes relatifs reciblés vers `../src/...` ; nouvel en-tête
  `swSpDesingVelocity.h` inclus par le seul `_E`.

### Types AD à reloger

`struct SegPhysData` et `enum class PhiTransform` (+ helpers inline `phi_correct`,
`m_from_phi`) sont définis en tête de l'ancien `swSp.h` (types AD, absents du
natif). Ils vont dans `swSpE.h` (ou un petit en-tête AD inclus par `swSpE.h`).

### Point de vigilance : macro `#define phi(i)`

Macro AD (`phi_ptr ? phi_ptr[i] : phi[i]`) confinée aux unités de compilation
`assim/`. N'expanse que la syntaxe d'appel `phi(i)`, jamais `phi[i]` ni la
déclaration `std::vector<double> phi` — donc coexiste avec le membre natif `phi`.

---

## ⚠ BLOCAGE IDENTIFIÉ — couplage de CONFIGURATION (à trancher)

Le découplage des **fonctions physiques** est propre (sous-classe). En revanche,
la **configuration** de l'assimilation est tissée dans DEUX fichiers natifs :

1. `common/swSolverParam.{h,cpp}` : ~40 champs AD + setters/getters + parsing
   (diff de **977 lignes**) — `minimize_algo`, `alpha_reg`, `alpha_spatial`,
   `adam_step_size`, `lbfgsb_*`, `t_state_assimilation`, `dt_map_ad`,
   `phi_transform`, `windowed_assimilation`, `n_assimilation_windows`, `obs_type`,
   poids `h/q_weight`, `zonation`, chemins de cartes de référence/zone, …
2. `models/sp/swSp.cpp::importAttributesFromSolverParam` : lit ces champs
   (`solver_param->getMinimizeAlgo()`, `getAlphaReg()`, `getTStateAssimilation()`,
   `getDtMapAd()`, `getPhiTransform()`, …) dans les membres de `swSp`.

En natif **pristine**, ces champs et getters n'existent pas → la couche AD ne peut
pas lire sa configuration comme aujourd'hui. Les fonctions physiques se
décollent ; la config, non. **Décision requise avant de continuer** (voir la
question posée à l'utilisateur). Options envisagées :

- (A) Lecteur de config AD dans `assim/` qui lit les clés AD depuis le même
  `input.solv` et peuple les membres de `swSpE` (zéro modif native ; ~40 clés à
  parser).
- (B) Fichier de config compagnon séparé (`assim.solv`) parsé uniquement par
  `assim/` (séparation physique maximale ; `input.solv` natif intact).
- (C) Extension **additive** de `swSolverParam` natif (petite modif native, contre
  l'objectif strict).

### DÉCISION (retenue) — COPIE FILTRÉE

**Un seul `input.solv`, zéro modif native, réutilisation MAXIMALE du natif.**
Obstacle confirmé : le parseur natif pristine fait `logFatal` sur toute clé
inconnue (le `else` terminal de la chaîne `if/else if (parameter_name == …)`).

Solution retenue (choix utilisateur, motif : réutiliser la fonction de lecture
native `readParameters()`) :

> Un lecteur **NOUVEAU** dans `assim/` (`assimConfig.{h,cpp}`) :
> 1. ouvre le **seul** `input.solv` rédigé par l'utilisateur ;
> 2. extrait les valeurs des clés AD (liste figée d'~40) → les stocke pour `swSpE` ;
> 3. écrit une **copie FILTRÉE** d'`input.solv` (clés AD retirées) ;
> 4. le **parseur natif `swSolverParam::readParameters()` lit cette copie filtrée
>    telle quelle** → aucune clé inconnue → pas de `logFatal`.
>
> `swSolverParam` reste 100% pristine et fait tout le parsing des clés natives.
> `assimConfig` ne parse QUE les ~40 clés AD (pas de duplication de la table native).

Mécanisme fichier (à finaliser en étape 3, options non destructives) :
- (i) dossier temporaire `input/` en miroir (symlinks des autres fichiers d'entrée
  + `input.solv` filtré), `output/` renvoyé vers le dossier réel ; ou
- (ii) filtrage en place avec sauvegarde `.assimbak` puis restauration juste après
  le `readParameters()`.
Recommandation : (i) — non destructif, aucune écriture dans le dossier utilisateur.

### Sous-paliers réalisés

- **2a — TERMINÉ** : `swSpE.h` (types AD `SegPhysData`/`PhiTransform`, helpers inline
  `phi_correct`/`dphi_dm`/`m_from_phi`/`get_phi`/`interpolateTimeSeries_E`, ~100 membres
  AD classés, déclarations de toutes les méthodes `_E`/AD, macro `phi(i)` confinée) +
  `swSpDesingVelocity.h` (copie standalone). Compile-check header : OK.
- **2b — TERMINÉ** : migration des 5 fichiers PHYSIQUES en `swSpE::`, tous compilés
  contre le natif pristine :
  - `swSpE_core.cpp`     (ex-`swSp_E.cpp` : `specificForce_E`, `allocate_as_adjoint`,
    `update_phi_norm_bounds`, `resolve_cost_weight`, `storeReferenceMap`)
  - `swSpE_hyp.cpp`      (ex-`swSpHyp_E.cpp` : LRStates/Fluxes/HLLC/BoundaryFlux/
    TimeStep/EdgeFactors/Balance `_E`)
  - `swSpE_muscl.cpp`    (ex-`swSpHypMusclEvr_E.cpp`)
  - `swSpE_friction.cpp` (ex-`swSpSourceFriction_E.cpp`)
  - `swSpE_bc.cpp`       (ex-`swParent_E.cpp` : `computeBC_E`)

  Technique d'includes retenue : compiler chaque fichier avec
  `-I<dossier_natif_d'origine>` (`.../models/sp`, et `.../models` pour le BC) → tous
  les includes relatifs natifs (`../common/…`, `../../../common/…`, `../swParent.h`)
  se résolvent vers le natif sans réécriture. Seules transformations sed : `swSp.h`→
  `swSpE.h`, `swParent.h`→`swSpE.h` (fichier BC), `swSp::`→`swSpE::`,
  `swParent::`→`swSpE::`, `swSp*`→`swSpE*`.

  Corrections header déclenchées par compile-check : ajout de `storeReferenceMap`
  (AD, absente du natif) et de l'inline `interpolateTimeSeries_E` (ex-`swParent`).

### Sous-palier 2c — TERMINÉ

Migration des 5 fichiers PIPELINE, tous compilés contre le natif pristine :
- `swSpE_checkpoint.cpp` (ex-`porosity_checkpoint.cpp` : `build_checkpoint_records`,
  `ckpt_segment_E`, `run_forward_E`, `generate_references_E`, `get_spatial_gradient_*`,
  `replay_trajectory_cost`, `run_taylor_test_*` ; `__enzyme_autodiff` sur `ckpt_segment_E`)
- `swSpE_ad.cpp`        (ex-`porosity_AD.cpp` : `execute_E`, FD checks)
- `swSpE_minimize.cpp`  (ex-`porosity_minimize.cpp` : `run_minimization`, PSQN/LBFGSpp/optim)
- `swSpE_gauges.cpp`    (ex-`porosity_gauges.cpp` : `locate_gauges`, `compute_cost_gauges_ptr`)
- `swSpE_spAD.cpp`      (ex-`sp_AD.cpp` : wrapper `__enzyme_autodiff` de démonstration)

Intrinsèques Enzyme (`__enzyme_autodiff`, `enzyme_dup`, `enzyme_const`, `enzyme_out`)
déclarés LOCALEMENT dans chaque fichier AD -> aucune dépendance au `param.h` natif.

Règle sed enrichie : ajout de `new swSp` -> `new swSpE` (allocation du solveur adjoint/
shadow). Corrections header via compile-check : ajout des membres `dt_map_ad`,
`t_state_assimilation`, `assimilation_window_duration` ; des méthodes
`run_pure_simulation` (x2), `test_enzyme_integration`, `compare_versions` ; include
`swBoundaryCondition.h` (type complet). **Aucun symbole `swSpE::` dupliqué** entre
fichiers (vérifié) -> pas de collision au link.

Bilan : **10 fichiers AD** (`swSpE_*.cpp`) + `swSpE.h` + `swSpDesingVelocity.h`,
tous compilent en C++17 contre le natif pristine. Noms de FONCTIONS inchangés.

---

## Architecture de BUILD (dédiée, natif intact)

Décision : cmake **DÉDIÉ** dans `assim/`, natif (`cmake/BuildTypes.cmake`,
`cmake/MainSetup.cmake`, `CMakeLists.txt`) **NON modifié**.

```
sw2d_enzyme_release/
├── src/ cmake/ CMakeLists.txt      <- natif PRISTINE (build de reference, compilo normal)
└── assim/
    ├── CMakeLists.txt              <- build AD autonome
    └── cmake/EnzymeFlags.cmake     <- flags Enzyme, gardes par ENABLE_ENZYME
```

Le mot-clé **`-DENABLE_ENZYME=ON`** est LE commutateur : sans lui, aucun flag Enzyme
n'est appliqué. `assim/CMakeLists.txt` :
- `option(ENABLE_ENZYME "..." OFF)` ;
- liste les sources NATIVES nécessaires (référencées `../src/...`, **non modifiées**)
  + les sources AD (`swSpE_*.cpp`, `main_enzyme.cpp`, `assimConfig.cpp`) ;
- include dirs : natif + `optim/include` + `LBFGSpp/include` + `psqn` + eigen3 ; link
  `liboptim.so` ;
- `if(ENABLE_ENZYME)` -> applique `${ENZYME_COMPILE_FLAGS}` (`-flto -O2 -fno-vectorize
  -fno-slp-vectorize -fno-unroll-loops`) et `${ENZYME_LINK_FLAGS}` (`-flto
  -fuse-ld=lld-18 -fpass-plugin=LLDEnzyme-18.so -Wl,--load-pass-plugin=...`) à la SEULE
  cible AD. Sinon : compilo/flags standard.

POINT DE CORRECTION IMPORTANT vs l'idée initiale « natif = lib statique » : les
fonctions natives appelées DANS la tape (`computeMusclSlopes`, `lambdaMinus/Plus`,
`specificForce`) doivent être vues par Enzyme en bitcode LTO. La cible AD **compile
donc ces sources natives elle-même** (en entrée, non modifiées) avec les flags gardés,
plutôt que de lier une lib pré-compilée en objets natifs. Le natif reste pristine
(fichiers non touchés) ; il est juste compilé une 2e fois, avec les flags AD, dans
l'unité LTO de la cible AD. L'exécutable de référence natif, lui, est bâti par le
CMake natif avec le compilo normal — les flags Enzyme ne le touchent JAMAIS.

## Étape 3 — TERMINÉE (build complet réussi)

- `assimConfig.{h,cpp}` : lit les 32 clés AD d'`input.solv`, écrit une copie FILTRÉE
  (clés AD retirées, commentaires inline nettoyés) dans un dossier temporaire
  `<case>/.assim_tmp/input/` (miroir par symlink du reste d'`input/`), lue TELLE
  QUELLE par `swSolverParam::readParameters()` natif. `apply_to(swSpE*)` pousse les
  valeurs AD dans les membres (le `importAttributesFromSolverParam` natif est `final`
  -> non surchargeable, d'où l'application après coup).
- `main_enzyme.cpp` : instancie `swSpE` directement (`new swSpE()+setIOPointer`,
  réplique de `swParent::create`), câble `assimConfig`, appelle `resolve_cost_weight`
  côté AD (l'init native pristine ne le fait pas). Accès à `phi` (protected natif) via
  l'accesseur public `swSpE::phi_data()`.
- `assim/cmake/EnzymeFlags.cmake` + `assim/CMakeLists.txt` : build autonome, natif
  compilé DANS la cible AD (bitcode LTO), flags Enzyme gardés par `-DENABLE_ENZYME=ON`.

### Résultat du build

```
cmake -S assim -B build_assim -DENABLE_ENZYME=ON -DCMAKE_CXX_COMPILER=clang++-18
cmake --build build_assim -j4
```
-> `build_assim/sw2d_enzyme` produit et EXÉCUTABLE (`-v` OK). Link Enzyme (plugin
LLDEnzyme-18) réussi.

Vérifications finales :
- **Mot-clé** : sans `-DENABLE_ENZYME=ON`, la cible AD n'est PAS construite (aucun flag
  Enzyme). Avec, flags LTO+plugin appliqués à la seule cible AD.
- **Natif pristine** : `diff -rq src/ sw2d_natif/src/` = **0 différence** ; 0 fichier
  natif contenant `enzyme`/`_E`.

### BUG de migration trouvé & corrigé — `computeBC_E` (inflow)

Validation runtime sur le cas Toce : le forward `_E` du release produisait un
domaine **entièrement sec** (0 cellule mouillée) alors que le binaire original
inonde (2820 cellules, vol 0→1.35). Même code physique `_E` (diff sémantique = 0).

CAUSE : il existe **deux** `computeBC_E`. `swParent::computeBC_E` (migré depuis
`swParent_E.cpp`) itère `bc_idx < bc_type_nb`, or **`bc_type_nb` vaut toujours 0**
(jamais initialisé) -> il ne remplit JAMAIS `bc_val_interpolated` -> l'hydrogramme
d'entrée (BC `q`) n'est jamais injecté -> domaine sec. La bonne version est le
`swSp::computeBC_E` **inline de swSp.h**, qui remplit `bc_val_interpolated` **indexé
par arête** (`edge_bnd_nb`), tel que le lit `hypComputeFluxes_E`. `swSp::computeBC_E`
masque celui de `swParent` -> le binaire original utilise la version par arête.

FIX : `computeBC_E` réécrit en inline (indexé par arête) dans `swSpE.h` ; suppression
de `swSpE_bc.cpp` (retiré des `AD_SOURCES`) ; ajout de l'include
`swBoundaryConditionTimeSeries.h` (son include dans `swSp.h` était un ajout AD absent
du pristine). Après fix, le forward `_E` du release reproduit l'original **à
l'identique** (vol=0.424/1.335/1.351/1.242 à t=10/15/20/25).

LEÇON : lors d'une redéfinition masquant une méthode de base (`swSp` masque
`swParent`), migrer la version EFFECTIVEMENT appelée (la plus dérivée), pas
l'homonyme de la classe de base.

### VALIDATION RUNTIME (Toce, reconstruction de porosité) — FAITE, −75.8 %

Cas isolé (aucune écriture dans le cas de travail ; `reference/` en symlink
lecture seule). Ébauche φ=1 partout, L-BFGS-B pavé (48 fenêtres).

| iter | RELEASE découplé | −75bak (réf) | run live utilisateur |
|------|------------------|--------------|----------------------|
| 3    | −6.41 %          | −6.41 %      | ✓ |
| 4    | −17.32 %         | −17.32 %     | ✓ |
| 5    | −9.79 %          | −9.79 %      | ✓ |
| 6    | −49.32 %         | −49.32 %     | ✓ |
| 7    | −64.59 %         | −64.6 %      | ✓ |
| 8    | −73.61 %         | −73.6 %      | ✓ |
| 9    | **−75.80 %**     | **−75.8 %**  | **−75.807 %** |

- **Le release découplé reproduit la descente profonde −75.80 % au chiffre près**,
  trajectoire itér.-par-itér. identique au run de référence ET au binaire original
  (build_enzyme). Forward `_E` par ailleurs bit-à-bit (vol=0.424/1.335/1.351/1.242).
- Natif pristine après tous les runs : `diff src/ sw2d_natif/src/` = **0**.

#### ⚠ Piège de configuration (leçon)

**Il n'y a AUCUNE régression de source.** La bonne config d'assimilation utilise
`h_weight 1` (pondération uniforme du coût). Le fichier `descent_-75bak/input.solv`
a été **édité ~4 h APRÈS** le run −75.8 % (dates : `descent_log.csv`=10:18,
`input.solv`=14:24) pour y activer la **carte de poids** `h_weight input/w_h_city.csv` :
cette carte **casse la descente** (−16.5 % au lieu de −75.8 %). Toujours reproduire
depuis la config RÉELLEMENT active (`sp_toce_travail/input/input.solv`), pas depuis un
`input.solv` de dossier de résultats potentiellement post-édité (vérifier les mtimes).

### Reste (optionnel)

- Lancer `sw2d_enzyme <cas> forward` / `gen-ref-e` / `taylor` sur un cas (ex. Gardon)
  et vérifier le test de Taylor (gradient adjoint == différences finies).
- Vendorer `optim/`, `LBFGSpp/`, `psqn/` dans le release (actuellement pointés vers le
  dépôt de dev via les caches `OPTIM_DIR`/`LBFGSPP_DIR`/`PSQN_DIR`).

---

Getters AD requis : `getMinimizeAlgo`, `getAlphaReg`, `getAlphaSpatial`,
`getAlphaL1`, `getAlphaPhiRef`, `getAdamStepSize`, `getAdamSafe`, `getLbfgsbM/…`,
`getTStateAssimilation`, `getDtMapAd`, `getPhiTransform`, `getWindowedAssimilation`,
`getAssimilationWindowDuration/Steps`, `getNAssimilationWindows`, `getZonation`,
`getHWeight`, `getQWeight`, `getObsType`, chemins ref/zone, … → tous fournis
désormais par `assimConfig` (côté assim), plus par `swSolverParam`.
