## GGA-PBE table (oncvpsp NC pseudos)

Otherwise explicitly noted, all pseudos have been generated using the neutral ground-state
configuration as reference, non-linear core correction in included except for
pseudos in which the all-electron (AE) core is too localized or pseudos having 1s in valence.
Pseudos with the GW tag have good logarithmic derivatives up to

## H

  * H.in: Uses one p-projector to improve transferability. Can be used for GW.
    Using v_p as local part produces vloc(q) with high Fourier components at large q

## He

  * He.in: Uses one p-projector to improve transferability. Can be used for GW.
    Using v_p as local part produces vloc(q) with high Fourier components at large q

## Li

  * Li-s.in: Default version with 1s in valence. Can be used for GW.

  * Li-s-high.in: Similar to Li-s but with smaller core radii. Can be used for GW

## Be

  * Be-s.in: Default version with 1s in valence. Can be used for GW

  * Be-s-high.in: Similar to Be-s but with smaller core radii. Can be used for GW.

## B:

  * B.in: Default version with (2s, 2p) in valence.

  * B-gw.in: Version with improved logarithmic derivatives, based on B.in.
    Recommended for GW.

## C

  * C.in: Default version with (2s, 2p) in valence

  * C-gw.in: Version with improved logarithmic derivatives, based on C.in.
    Recommended for GW.

## N

  * N.in: Default version with (2s, 2p) in valence. Can be used for GW.

## O

  * O.in: Default version with (2s, 2p) in valence and one d-projector for unbound 3d.
    Can be used for GW.

  * O-high.in: Similar to O.in but with smaller core radii for high-accuracy applications.
    Can be used for GW.

<!--
[Generated new pseudo, TO BE TESTED]
F:***   
    F.psp8 can be improved (model core a bit too hard, phonons don't converge)
    PHONONS are not smooth
-->

## Ne

These pseudos are without nlcc since modeling 1s core density 
without spoiling convergence is not trivial

  * Ne.in: Default version with (2s, 2p) in valence.

  * Ne-high.in: Similar to Ne.in but with smaller core radii. Requires GW version.

## Na

  * Na-sp.in: (2s, 2p) semicore in valence. No model core charge (AE core too localized).
    Can be used for GW.

## Mg

  * Mg.in: Low-cutoff low-accuracy version without semicore.
    Has ghost at +80 eV, not recommended for GW.

  * Mg-sp.in: (2s, 2p) semicore states in valence. No model core charge (AE core too localized).
    Can be used for GW.

## Al

  * Al.in: Has 2 d-projectors to improve transferability.
    Oscillations in form factors. Requires GW version.

## Si

  * Si.in: Has 2 d-projectors to improve transferability.
    Requires GW version

## P

  * P.in: Has 2 d-projectors to improve transferability.
    Requires GW version
    $$ the model core charge should be improved
    $$ the algorithm detects dispersionless states but inspection of the BS does not show any

##Cl
    Cl-new
    Has 2 d-projectors to improve transferability.

## Ar

  * Ar.in: (3s, 3p) in valence. Has 2 d-projectors to improve transferability.
    Can be used for GW.

## K

  * K-sp.in: (3s, 3p) semicore states in valence. Has 2 d-projectors to improve transferability.
    Can be used for GW.

## Ca

  * Ca-sp.in: (3s, 3p) semicore states in valence. Has 2 d-projectors to improve transferability.
    Can be used for GW.

## Sc

  * Sc-sp.in: (3s, 3p) semicore states in valence. Can be used for GW.

## Ti

  * Ti-sp.in: (3s, 3p) semicore states in valence. Can be used for GW.

## V

  * V-sp.in: (3s, 3p) semicore states in valence. Can be used for GW.

## Cr

  the mc could be improved

  * Cr-sp.in: Default version with (3s, 3p) semicore states in valence. 
    Can be used for GW.

  * Cr-sp-high.in: Similar to Cr-sp.in but with smaller core radii (and larger cutoff) 
    for high-accuracy calculations. Can be used for GW.

## Mn

  * Mn-sp.in: Default version with (3s, 3p) semicore states in valence.
    Can be used for GW.

  * Mn-sp-high.in: Similar to Mn-sp.in but with smaller core radii (and larger cutoff)
    for high-accuracy calculations. Can be used for GW.

## Fe

  model core charge too hard, try the one from Fe-sp-high.

  * Fe-sp.in: Default version with (3s, 3p) semicore states in valence.
    Can be used for GW. Uses `ncon 3`, the standard is 4.

  * Fe-sp-high.in: Smaller core radii for high accuracy calculations
    Can be used for GW. Uses `ncon 3`, the standard is 4.

## Co:

  mc could be improved

  * Co-sp.in: Default version with (3s, 3p) semicore states in valence.
    Can be used for GW.

  * Co-sp-high.in: Smaller core radii for high accuracy calculations.
    Can be used for GW.

## Ni

  * Ni-sp.in: Default version with (3s, 3p) semicore states in valence.

  * Ni-sp-high.in: Smaller core radii for high accuracy calculations.

## Cu

  * Cu-sp.in: Default version with (3s, 3p) semicore states in valence.
    Can be used for GW.

  * Cu-sp-high.in: Smaller core radii for high accuracy calculations.
    Can be used for GW.

## Zn

  * Zn.in

  * Zn-sp.in
    I've added the GW tag (ok but not "perfect", ask Michiel if he has specialized version)
    Ask about mc params
    $$ Michiel will look at this. the total energy convergence look a bit suspicious...

## Ga

  * Ga-low.in: Low-accuracy version with 3d in core, (4s,4p) in valence.
    2 d-projectors for unbound states.

  * Ga-d.in: 3d in valence, default for GS applications. NOT recommended for GW.

  * Ga-spd-high.in: Include full (3s, 3p, 3d) shell in valence.
    Recommended version for GW and high accuracy calculations but it's hard

## Ge

  * Ge-low.in: Low-accuracy version with 3d in core, (4s,4p) in valence.
    2 d-projectors for unbound states.

  * Ge-d: 3d in valence, default for GS applications. NOT recommended for GW.

  * Ge-spd-high.in: Include full (3s, 3p, 3d) shell in valence. Recommended version for GW
    and high accuracy calculations but it's hard

## As
  Trying As-new-mc1.in
  [DONE, RUNNING] phonons are missing

  * As.in:** No convergence (v(d)?) 
    Try As-new with modcore from As

  * As-d: Require GW version

  * As-spd-high has good logders but it's hard

## Se

  * Se.in: Not recommended for GW

  * Se-d.in: Include 3d in valence. Not recommended for GW

  * Se-spd.in: Se-spd-high has good logders but it's hard)
    TODO: Remove model core-charge too localized

## Br

  * Br.in: Not recommended for GW

  * Br-d.in: Include 3d in valence. Not recommended for GW

  * Br-spd.in: Include full (3s, 3p, 3d) shell in valence.
    Recommended version for GW and high accuracy calculations but it's hard
    TOO HARD (model core)

## Kr

  * Kr.in: needs GW version

## Rb  
    Rb-sp.in: use two d-projectors for unbound d
    use dvloc0 2.5 to improve scattering properties at high energy
    Can be used for GW.

## Sr

  * Sr-sp.in: use two d-projectors for d (4d is bound in GGA-PBE)
    Require GW version.

## Y

  * Y-sp.in: 4s-4p semicore in valence. Requires GW version?

## Zr

  * Zr-sp.in: 4s-4p semicore in valence. Can be used for GW.

## Nb

  * Nb-sp.in: 4s-4p semicore in valence. Can be used for GW

## Mo

  * Mo-sp.in: 4s-4p semicore in valence. Can be used for GW

## Tc

  * Tc-sp.in: 4s-4p semicore in valence. Can be used for GW

## Ru

  * Ru-sp.in: 4s-4p semicore in valence. Can be used for GW

## Rh

  * Rh-sp.in: 4s-4p semicore in valence.

## Pd

  * Pd-sp.in: 4s-4p semicore in valence. 
    Require GW version
    Previous attempt to generate Pd without semicore lead to ghost state at +1 eV  

## Ag

  * Ag-sp.in. 
    It seems difficult to get good logders without semicore

## Cd

  * Cd.in: (4d, 5s) in valence. Has ghost at +73 Ha

  * Cd-sp.in. (4s, 4p, 4d, 5s) in valence
    Can be used for GW.

## In

  * In-d: Default version for GS application with 4d in valence.
    GHOST_In-20eV. Not recommended for GW.

  * In-spd-high.in: Include full (4s, 4p, 4d) shell in valence.
    Recommended version for GW and high accuracy calculations.

## Sn

  * Sn.in: (5s,5p) in valence + 2 d-projectors for unbound state to improve transferability
    Not recommended for GW.

  * Sn-d.in: Include 4d in valence. Default version for GS applications. 
    Not recommended for GW. Ghost at +60 eV

  * Sn-spd-high.in: Include full (4s, 4p, 4d) shell in valence.
    Recommended version for GW and high accuracy calculations.

## Sb

  * Sb.in: (5s,5p) in valence + 2 d-projectors for unbound state to improve transferability
    Not recommended for GW.

  * Sb-d.in: Include 4d in valence. Default version for GS applications. 
    Not recommended for GW. Ghost at +20 eV,

  * Sb-spd-high.in: Include full (4s, 4p, 4d) shell in valence.
    Recommended version for GW and high accuracy calculations.

## Te

  * Te.in: (5s,5p) in valence + 2 d-projectors for unbound state to improve transferability
    Not recommended for GW.

  * Te-d.in: Include 4d in valence. 
    Not recommended for GW. Ghost at +77 eV.

  * Te-spd-high.in: TODO Complete tests

## I

  * I.in: Not recommended for GW

  * I-d.in: Include 4d in valence. Not recommended for GW

  * I-spd-high: Include full (4s, 4p, 4d) shell in valence.
    Recommended version for GW and high accuracy calculations.

## Xe

  * Xe.in: 5s-5p in valence. Two projectors for f to improve trasferability.
    Require GW version

## Cs

  * Cs-sp.in: (5s, 5p, 6s) in valence, 4f in core.
    Use two f-projectors to improve transferability.
    Require GW version.

## Ba

  * Ba-sp-new.in: seems to have a ghost around 6 eV
    Ba-sp-new adds two projectors for f to improve deltafactor and GBRV tests.
    Uses excited configuration

## Hf

  Needs a GW version fsp is probably unavoidable but log ders are not optimal.
  the p version is the only one that does not have an explicit f projector, did we never try?

  * Hf-p.in:
    Use dvloc0 8.0

  * Hf-sp.in:
    4f frozen in core. Has 2 f-projectors to improve transferability.
    Recommended for GS applications.

  * Hf-fsp.in:

## Ta 
  Needs GW version

  * Ta-sp.in:
    4f frozen in core. Has 2 f-projectors to improve transferability.
    Recommended for GS applications.
    The missing f, only at ~1 H below the fermi level is probably too much to classify the sp as GW

  * Ta-fsp.in:

## W
  Needs GW version

  * W-sp.in:
    4f frozen in core. Has 2 f-projectors to improve transferability.
    Recommended for GS applications.
    The missing f, only at ~1 H below the fermi level is probably too much to classify the sp as GW

  * W-fsp.in:

## Re:

  * Re-sp.in
    4f frozen in core. Has 2 f-projectors to improve transferability.
    Recommended for GS applications.
    The missing f, only at ~1.5 H below the fermi level is probably too much to classify the sp as GW

## Os

  * Os-sp.in: needs GW version
    4f frozen in core. Has 2 f-projectors to improve transferability.
    Recommended for GS applications.

## Ir

  * Ir-sp: removed the GW tag, 4f to close
    4f frozen in core. Has 2 f-projectors to improve transferability.
    Recommended for GS applications.

## Pt
  fsp softer model core charge? the two highest energy points seem off...
  add an f projector to sp?
  $ removed the GW tag from sp, its high in energy and the fsp shows a serious overlap
  (Pt-p gives much better deltafactor but has ghosts at +25 eV,
   Pt-ex-sp does not improve, new df with pt-sp-icmod3 is 4.8)

  * Pt-sp.in:

  * Pt-fsp.in:

## Au
  $ removed the GW tag from sp, its high in energy and the fsp shows a serious overlap
  $ we actually did the GW calcualtion on gold both SigmaX and SigmaC are seriously affected by unfreezing
  $ the 4f electrons ~10% total effect. since they are les deep in the Hf-Pt series I don't trust the GW tag on any
  $ of them on non f potentials

  * Au-sp.in
    4f frozen in core. Has 1 f-projectors to improve transferability.
    Recommended for GS applications.

## Hg: ??

  * Hg.in: (5d, 6s) in valence.
    Not recommended for GW. Ghost at +66 eV 

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

<!--
GGA_X_PBE_SOL = 116 GGA_C_PBE_SOL = 133
LDA_X = 1 LDA_C_PW = 12
7: dict(x=xcf.LDA_X, c=xcf.LDA_C_PW),          # PW 001012
-->

 pseudo_dojo/pseudos/ONCVPSP-PBE-PDv0.3/F/F-new.djrepo        |      5 -
La/La-sp.in:# lloc lpopt rc5 dvloc0
La/La-sp.in-4 5 1.6 2.0
