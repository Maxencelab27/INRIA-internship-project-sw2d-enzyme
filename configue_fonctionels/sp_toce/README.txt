CAS FONCTIONNEL — sp_toce (comparaison SP / SP, assimilation de porosité)
=========================================================================

But
----
Cas de reconstruction de porosité par assimilation variationnelle (Enzyme),
destiné à la comparaison SP/SP : les observations de reference/ ont été
générées par le MÊME modèle SP (single porosity) que celui optimisé, sur le
même maillage (twin experiment). On part d'une ébauche de porosité uniforme
(phi = 1 partout) et on reconstruit la carte de porosité qui explique les
observations (la "ville" / structure du cas Toce).

Contenu
-------
  input/       maillage (geometry.sw2d), cartes (friction, head loss, CI),
               hydrogramme amont, sondes, et input.solv (la configuration).
  reference/   observations SP de référence (ResultHydro_p*.txt), auto-cohérentes
               (mêmes équations SP que le forward assimilé).
  descent_log_reference_-75.8pct.csv   trace de descente attendue (référence).

Lancer
------
  # binaire découplé (depuis la racine du release, build avec -DENABLE_ENZYME=ON) :
  ./build_assim/sw2d_enzyme  configue_fonctionels/sp_toce  minimize

  Modes utiles : forward | gen-ref-e | taylor | minimize

Résultat attendu
----------------
  Ébauche phi=1 -> J_init ≈ 3.104. La minimisation L-BFGS-B (pavée, 48 fenêtres)
  descend à J ≈ 0.751, soit ~ -75.8 % du coût, en ~9 itérations (cf. le
  descent_log de référence fourni).

Point de configuration IMPORTANT
--------------------------------
  input.solv utilise  h_weight 1  (pondération UNIFORME du coût). C'est ce qui
  donne la descente -75.8 %. NE PAS activer la carte de poids w_h_city.csv
  (h_weight input/w_h_city.csv) : elle bloque la descente vers ~ -16 %.

(Documentation générale du code AD/Enzyme et du découplage : voir ../../assim/.)
