PBEsol

## H 

  * H.in: Uses one p-projector to improve transferability. Can be used for GW
    Using v_p as local part produces vloc(q) with high Fourier components at large q
    Phonons are noisy.

## He 

 *  He.in: Uses one p-projector to improve transferability. Can be used for GW.
    Using v_p as local part produces vloc(q) with high Fourier components at large q
    phonons are noisy.

## Li

  * Li-s.in: Default version with 1s in valence.

  * Li-s-high.in: Similar to Be-s but with smaller core radius.

## Be

  * Be-s.in: Default version with 1s in valence.

  * Be-s-high.in: Similar to Li-s but with smaller core radius.

## B

  * B.in: Default version with (2s, 2p) in valence.

  * B-gw.in: Version with improved logarithmic derivatives, recommended for GW.

## C No PHGAMMA

  * C.in: Default version with (2s, 2p) in valence
  * C-gw.in: Version with improved logarithmic derivatives, recommended for GW.

## N   

  * N.in: Default version with (2s, 2p) in valence
  Phonons at Gamma present some noise

## O NO PHGAMMA

  * O.in: Default version with (2s, 2p) in valence and one d-projector for unbound 3d.
  * O-high.in: Similar to O.in, smaller core radii for high-accuracy applications.

## F:
    TODO: PHONONS are missing

Ne: Pseudo without nlcc since modeling 1s core density without spoiling convergence is not trivial
    Ne.in: Default version with (2s, 2p) in valence.
    Ne-high.in: Similar to Ne.in but with smaller core radii.
-->

## Na

  * Na-sp.in: (2s, 2p) semicore in valence. No model core charge (ae core too localized).
    Phonons woasr are noisy
 
## Mg

  * Mg.in: has ghost at +80 eV, it's the fast low-accuracy version without semicore.
  * Mg-sp.in: (2s, 2p) semicore in valence. No model core charge (ae core too localized)

## Al 

  * Al.in: Oscillations in from factors, Requires GW version.
  Has 2 projectors in the d-channel to improve transferability.

## Si  REQUIRES validation (see PBE)

  * Si.in: Requires GW version
  Has 2 projectors in the d-channel to improve transferability.

## P 
  * P.in: requires GW version 
  Has 2 projectors in the d-channel to improve transferability.
  $$ the model core charge should be improved
  $$ the algorithm detects dispersionless states but inspection of the BS does not show any

<!--
[new pseudo, TO BE TESTED Complete PHGAMMA]
Cl:** (vloc(d), slow convergence)
  $$ lets use new, some phonons are missing though
-->

## Ar
  * Ar.in: (3s, 3p) in valence. Has 2 projectors in the d-channel to improve transferability.
    Phonons woasr are noisy

## K
  * K-sp.in: (3s, 3p) semicore states in valence. 
  Has 2 projectors in the d-channel to improve transferability.

## C
  * Ca-sp.in: (3s, 3p) semicore states in valence. 
  Has 2 projectors in the d-channel to improve transferability.

## Sc
  * Sc-sp.in: (3s, 3p) semicore states in valence. 
    Phonons woasr are noisy

## Ti
  * Ti-sp.in: (3s, 3p) semicore states in valence. 

## V
  * V-sp.in: (3s, 3p) semicore states in valence. 

## Cr the mc could be improved
  * Cr-sp.in: Default version with (3s, 3p) semicore states in valence. 
    TODO: PHGAMMA is missing RUNNING

  * Cr-sp-high.in: Similar to Cr-sp.in but with smaller core radii (and larger cutoff)
    for high-accuracy calculations.

## Mn

  * Mn-sp.in: Default version with (3s, 3p) semicore states in valence. 

  * Mn-sp-high.in: Similar to Mn-sp.in but with smaller core radii (and larger cutoff)
    for high-accuracy calculations.
    TODO: Phgamma RUNNING

## Fe model core charge too hard, try the one from Fe-sp-high.

  * Fe-sp.in: Default version with (3s, 3p) semicore states in valence.
    Uses `ncon 3`, the standard is 4

  * Fe-sp-high.in: Smaller core radii for high accuracy calculations.
    Uses `ncon 3`, the standard is 4

## Co mc could be improved

  * Co-sp.in: Default version with (3s, 3p) semicore states in valence. deltafactor ~= 1.2.
  * Co-sp-high.in: Smaller core radii for high accuracy calculations. deltafactor ~= 0.7

## Ni
  * Ni-sp.in: Default version with (3s, 3p) semicore states in valence. deltafactor ~= 1.2
  * Ni-sp-high.in: Smaller core radii for high accuracy calculations. deltafactor ~= 1.2

## Cu
  * Cu-sp.in: Default version with (3s, 3p) semicore states in valence. deltafactor ~= 0.3
  * Cu-sp-high.in: Smaller core radii for high accuracy calculations. deltafactor ~= 0.3

Zn-sp: 
    I've added the GW tag (ok but not "perfect", ask Michiel if he has specialized version)
    Ask about mc params
    $$ Michiel will look at this. the total energy convergence look a bit suspicious...
    ** phonons are noisy.

## Ga
  * Ga-d.in: 3d in valence, default for GS applications.

## Ge
  * Ge-d: Require GW version 

<!--
[DONE, RUNNING] phonons are missing
As:
    As.in:** No convergence (v(d)?) Try As-new with modcore from As
    As-d:* Require GW version 
    As-spd-high has good logders but it's hard)

-->
## Se 
  * Se-d.in: Require GW version 
  TODO: PHGAMMA RUNNING

## Br needs GW version
   TODO: PHGAMMA RUNNING

## Kr needs GW version
   phonons are noisy but it could be due to the particular system

## Rb  
    Rb-sp.in: use two d-projectors for unbound d
    use dvloc0=2.5 to improve scattering properties at high energy

Sr:
-->

# Y
  * Y-sp.in: 4s-4p semicore in valence
    Phonons woasr are noisy


# Zr

  * Zr-sp.in: 4s-4p semicore in valence.
  Phonons woasr are noisy

# Nb
  * Nb-sp.in: 4s-4p semicore in valence.

# Mo
  * Mo-sp.in: 4s-4p semicore in valence.

# Tc
  * Tc-sp.in: 4s-4p semicore in valence.

# Ru
  * Ru-sp.in: 4s-4p semicore in valence.

# Rh
  * Rh-sp.in: 4s-4p semicore in valence.

<!--
Pd: Previous attempt to generate Pd without semicore lead to ghost state at +1 eV  

Ag: It seems difficult to get good logders without semicore

Cd:
    Cd.in:* Has ghost at +73 Ha $$ to discuss
-->

In-d: ghost at +30eV (In-spd-high is much better for GW)

Sn-d: (ghost at +60, use Sn-spd-high for GW)

Sb-d: ghost at +20 eV, Use Sb-spd-high for GW

Te-d:** ghost at +77 eV  (Should try Te-spd-high for GW)

I: needs GW version

Xe: Require GW version
    Xe.in: 5s-5p in valence. Two projectors for f to improve trasferability.

<!--
## Cs

  * Cs-sp.in: (5s, 5p, 6s) in valence, 4f in core.
    Use two f-projectors to improve transferability.
    Require GW version.
    PHGAMMA is missing

## Ba

  * Ba-sp-new.in: seems to have a ghost around 6 eV
    Ba-sp-new adds two projectors for f to improve deltafactor and GBRV tests.

## Hf

  Needs a GW version fsp is probably unavoidable but log ders are not optimal.
  the p version is the only one that does not have an explicit f projector, did we never try?

  * Hf-p.in:
    Use dvloc0 8.0

  * Hf-sp.in:

  * Hf-fsp.in:

## Ta

  * Ta-sp.in:
  The missing f, only at 1 H below the fermi level is probably too much to classify the sp as GW

  * Ta-fsp.in:

## W

  * W-sp.in:

  * W-fsp.in:

## Re:

  * Re-sp.in: needs GW version

## Os

  * Os-sp.in: needs GW version

## Ir

  * Ir-sp: removed the GW tag, 4f to close

## Pt
  fsp softer model core charge? the two highest energy points seem off...
  add an f projector to sp?
  $ removed the GW tag from sp, its high in energy and the fsp shows a serious overlap
  Pt-sp: Regenerated with modcore 3
  (Pt-p gives much better deltafactor but has ghosts at +25 eV,
  Pt-ex-sp does not improve, new df with pt-sp-icmod3 is 4.8)
  Ask michiel about MC

  * Pt-sp.in:

  * Pt-fsp.in:

## Au
  $ removed the GW tag from sp, its high in energy and the fsp shows a serious overlap
  $ we actually did the GW calcualtion on gold both SigmaX and SigmaC are seriously affected by unfreezing
  $ the 4f electrons ~10% total effect. since they are les deep in the Hf-Pt series I don't trust the GW tag on any
  $ of them on non f potentials

  * Au-sp.in

## Hg: ??

  * Hg.in: (5d, 6s) in valence.
    Not recommended for GW. Ghost at +66 eV 
    PHGAMMA is missing (problems to converge)

  * Hg-sp.in: (5s, 5p, 5d, 6s) in valence.
    Recommended version for GW.

## Tl:
  
  * Tl-d.in: Default version for GS applications with 5d in valence. 
    Not recommended for GW. Ghost at ~ +60 eV.

  * Tl-spd-high.in: Include full (6s, 6p, 6d) shell in valence.
    Recommended version for GW.
    4f electrons in core may still be to high in energy (4.2H) for GW.

## Pb:

  * Pb-d.in: Default version for GS applications with 5d in valence. 
    Not recommended for GW. Ghost at +70 eV.

  * Pb-spd-high.in: Include full (6s, 6p, 6d) shell in valence.
    Recommended version for GW.
    the 4f may still be to high in energy (4.8H) for GW

## Bi:

  * Bi-d.in: Default version for GS applications with 5d in valence. 
    Not recommended for GW. Ghost at +70 eV.
    TODO: PHGAMMA Missing

  * Bi-spd-high.in: Include full (6s, 6p, 6d) shell in valence.
    Recommended version for GW. $ 4f at -5.6H

## Po

  * Po-d.in: Default version for GS applications with 5d in valence. 
    Not recommended for GW. Ghost at +58 eV

  * Po-spd-high.in: Include full (6s, 6p, 6d) shell in valence.
    Recommended version for GW. $ 4f at 6.4H

## Rn: 

  * Rn.in: Default version for GS applications with (6s, 6p) in valence. 
    Not recommended for GW. d is only 1.7H below ef

  * Rn-d.in: (5d, 6s, 6p) in valence for high accuracy calculations. 
    Not recommended for GW. 


## TODO
    * test for the 6th row the addition of an explicit f projector
    * Pd --> Cd (look at Pd/Ag/Cd carefully)
    
    Cu-sp high: TO BE INVESTIGATED:
    Test that there's no fr pseudo with local part taken from v(L)

    W-sp, Re-sp (f channel could be improved)

    In-d is very sensitive, changing XC (GGa-->LDA) gives a ghost very close to 0.
    The first s projectors presents oscillations I don't see in the other pseudos.
-->


 pseudo_dojo/pseudos/ONCVPSP-PBEsol-PDv0.3/As/As-d.djrepo        |    137 +
 pseudo_dojo/pseudos/ONCVPSP-PBEsol-PDv0.3/F/F-new.djrepo        | 158536 ++++++++++++++++++++++++++++++++++++++++++++++
 pseudo_dojo/pseudos/ONCVPSP-PBEsol-PDv0.3/Lu/Lu-fsp-soft.djrepo |     84 +
 pseudo_dojo/pseudos/ONCVPSP-PBEsol-PDv0.3/Lu/Lu-sp.djrepo       |     63 +
