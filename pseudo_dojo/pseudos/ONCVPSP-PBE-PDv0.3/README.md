## GGA-PBE table (oncvpsp NC pseudos)

Unless explicitly noted otherwise, all pseudos have been generated using the neutral ground-state
configuration as reference, non-linear core correction is always included except for
pseudos in which the all-electron (AE) core is too localized (mentioned in the text) or pseudos with 1s in valence.
Pseudos with the GW tag have good logarithmic derivatives up to 6 H.

## H

  * H.in: Default version with one p-projector to improve transferability. Can be used for GW.
    Using v_p as local part produces vloc(q) with high Fourier components at large q.

## He

  * He.in: Default version with one p-projector to improve transferability. Can be used for GW.
    Using v_p as local part produces vloc(q) with high Fourier components at large q.

## Li

  * Li-s.in: Default version with 1s in valence. Can be used for GW.

  * Li-s-high.in: Similar to Li-s but with smaller core radii. Can be used for GW

## Be

  * Be-s.in: Default version with 1s in valence. Can be used for GW

  * Be-s-high.in: Similar to Be-s but with smaller core radii. Can be used for GW.

## B:

  * B.in: Default version with (2s, 2p) in valence. Can be used for GW

## C

  * C.in: Default version with (2s, 2p) in valence

## N

  * N.in: Default version with (2s, 2p) in valence. Can be used for GW.

## O

  * O.in: Default version with (2s, 2p) in valence and one d-projector for unbound 3d.
    Can be used for GW.

  * O-high.in: Similar to O.in but with smaller core radii for high-accuracy applications.
    Can be used for GW.

## F: 

  * F.in: Default version with (2s, 2p) in valence and one d-projector for unbound 3d.
    Needs GW version

## Ne

  These pseudos are without nlcc since modeling 1s core density 
  without spoiling convergence is not trivial

  * Ne.in: Default version with (2s, 2p) in valence.
    Needs GW version

  * Ne-high.in: Similar to Ne.in but with smaller core radii.
    Needs GW version

## Na

  * Na-sp.in: Default version with (2s, 2p) semicore in valence. No model core charge (AE core too localized).
    Can be used for GW.

## Mg

  * Mg.in: Low-cutoff low-accuracy version without semicore.
    Has ghost at +80 eV, not recommended for GW.

  * Mg-sp.in: Default version with (2s, 2p) semicore states in valence. No model core charge (AE core too localized).
    Recommended for GS. Can be used for GW.

## Al

  * Al.in: Default version with 2 d-projectors to improve transferability.
    Needs GW version
    TODO: Oscillations in form factors. 

## Si

  * Si.in: Default version with 2 d-projectors to improve transferability.
    Needs GW version

## P

  * P.in: Default version with 2 d-projectors to improve transferability.
    TODO: the model core charge and GBRV should be improved
    Needs GW version

## S

  * S.in: Default version with 2 d-projectors to improve transferability.
    Needs GW version

## Cl

  * Cl.in: Default version with 2 d-projectors to improve transferability.
    Needs GW version

## Ar

  * Ar.in: Default version with (3s, 3p) in valence. Has 2 d-projectors to improve transferability.
    Can be used for GW.

## K

  * K-sp.in: Default version with (3s, 3p) semicore states in valence. Has 2 d-projectors to improve transferability.
    Can be used for GW.

## Ca

  * Ca-sp.in: Default version with (3s, 3p) semicore states in valence. Has 2 d-projectors to improve transferability.
    Can be used for GW.

## Sc

  * Sc-sp.in: Default version with (3s, 3p) semicore states in valence. 
    Can be used for GW.

## Ti

  * Ti-sp.in: Default version with (3s, 3p) semicore states in valence. 
    Can be used for GW.

## V

  * V-sp.in: Default version with (3s, 3p) semicore states in valence. 
    Can be used for GW.

## Cr

  * Cr-sp.in: Default version with (3s, 3p) semicore states in valence. 
    Can be used for GW.

  * Cr-sp-high.in: Similar to Cr-sp.in but with smaller core radii (and larger cutoff) 
    for high-accuracy calculations e.g. magnetic systems.
    Can be used for GW.

## Mn

  * Mn-sp.in: Default version with (3s, 3p) semicore states in valence.
    Can be used for GW.

  * Mn-sp-high.in: Similar to Mn-sp.in but with smaller core radii (and larger cutoff)
    for high-accuracy calculations e.g. in magnetic systems.
    Can be used for GW.

## Fe

  * Fe-sp.in: Default version with (3s, 3p) semicore states in valence.
    Can be used for GW. Uses `ncon 3`, the standard is 4.

  * Fe-sp-high.in: Smaller core radii for high accuracy calculations e.g. in magnetic systems.
    Can be used for GW. Uses `ncon 3`, the standard is 4.

## Co

  AE core very localized --> Had to find a compromise for model core charge 

  * Co-sp.in: Default version with (3s, 3p) semicore states in valence.
    Can be used for GW.

  * Co-sp-high.in: Smaller core radii for high accuracy calculations e.g. in magnetic systems.
    Can be used for GW.

## Ni
  Needs GW version

  * Ni-sp.in: Default version with (3s, 3p) semicore states in valence.

  * Ni-sp-high.in: Smaller core radii for high accuracy calculations.

## Cu

  * Cu-sp.in: Default version with (3s, 3p) semicore states in valence.
    Can be used for GW.

  * Cu-sp-high.in: Smaller core radii for high accuracy calculations.
    Can be used for GW.

## Zn

  * Zn.in:  Low-accuracy version.

  * Zn-sp.in: Default version. Has the GW tag but mind ZnO.

## Ga

  * Ga-low.in: Low-accuracy version with 3d in core, (4s,4p) in valence.
    2 d-projectors for unbound states.

  * Ga-d.in: Default version with 3d in valence, for GS applications. NOT recommended for GW.

  * Ga-spd-high.in: Include full (3s, 3p, 3d) shell in valence.
    Recommended version for GW and high accuracy calculations but it's hard
    TODO: model core-charge too localized

## Ge

  * Ge-low.in: Low-accuracy version with 3d in core, (4s,4p) in valence.
    2 d-projectors for unbound states.

  * Ge-d.in: Default version with 3d in valence, for GS applications. NOT recommended for GW.

  * Ge-spd-high.in: Include full (3s, 3p, 3d) shell in valence. Recommended version for GW
    and high accuracy calculations but it's hard
    TODO: model core-charge too localized

## As
  Needs GW version

  * As.in:
    Version with d frozen.

  * As-d.in: Default version with d electrons in valence. Recommended for GS applications.


  * As-spd-high has good logders but it's hard and model core is peaked

## Se

  * Se.in: Not recommended for GW

  * Se-d.in: Default version with 3d in valence. Not recommended for GW

  * Se-spd-high.in: Se-spd-high has good logders but it's hard)
    TODO: model core-charge too localized

## Br

  * Br.in: Default version for GS applications. Not recommended for GW

  * Br-d.in: Include 3d in valence. Not recommended for GW

  * Br-spd-high.in: Include full (3s, 3p, 3d) shell in valence.
    Recommended version for GW and high accuracy calculations but it's hard

    TODO: model core-charge too localized (hints are not enough to converge GBRV
    compounds with normal)

## Kr
  Needs GW version

  * Kr.in: Default version, needs GW version

## Rb  

  * Rb-sp.in: Default version with two d-projectors for unbound d
    uses dvloc0 2.5 to improve scattering properties at high energy
    Can be used for GW.

## Sr
  Needs GW version

  * Sr-sp.in: Default version with two d-projectors for d (4d is bound in GGA-PBE)
    Requires GW version.

## Y
  Needs GW version

  * Y-sp.in: Default version with 4s-4p semicore in valence.

## Zr

  * Zr-sp.in: Default version with 4s-4p semicore in valence. 
    Can be used for GW.

## Nb

  * Nb-sp.in: Default version with 4s-4p semicore in valence. 
    Can be used for GW

## Mo

  * Mo-sp.in: Default version with 4s-4p semicore in valence. 
    Can be used for GW

## Tc

  * Tc-sp.in: Default version with 4s-4p semicore in valence. 
    Can be used for GW

## Ru

  * Ru-sp.in: Default version with 4s-4p semicore in valence. 
    Can be used for GW

## Rh
  Needs GW version

  * Rh-sp.in: Default version with 4s-4p semicore in valence.

## Pd
  Needs GW version

  * Pd-sp.in: Default version with 4s-4p semicore in valence. 
    Previous attempt to generate Pd without semicore lead to ghost state at +1 eV  

## Ag
  Needs GW version

  * Ag-sp.in: Default version. 
    It seems difficult to get good logders without sp semicore

## Cd

  * Cd.in: version with (4d, 5s) in valence. Has ghost at +73 Ha

  * Cd-sp.in: (4s, 4p, 4d, 5s) in valence
    Default version. Can be used for GW.

## In
  PROBLEMATIC?

  * In-d.in: Default version for GS application with 4d in valence.
    Has GHOST_In-20eV. Not recommended for GW.

  * In-spd-high.in: Include full (4s, 4p, 4d) shell in valence.
    Recommended version for GW and high accuracy calculations.

## Sn

  * Sn.in: (5s,5p) in valence + 2 d-projectors for unbound state to improve transferability
    Not recommended for GW.

  * Sn-d.in: Default version with 4d in valence. Default version for GS applications. 
    Not recommended for GW. Ghost at +60 eV

  * Sn-spd-high.in: Include full (4s, 4p, 4d) shell in valence.
    Recommended version for GW and high accuracy calculations.

## Sb
  PROBLEMATIC? 

  * Sb.in: (5s,5p) in valence + 2 d-projectors for unbound state to improve transferability
    Not recommended for GW.

  * Sb-d.in: Default version with 4d in valence. Default version for GS applications. 
    Not recommended for GW. Ghost at +20 eV,

  * Sb-spd-high.in: Include full (4s, 4p, 4d) shell in valence.
    Recommended version for GW and high accuracy calculations.

## Te

  * Te.in: Version (5s,5p) in valence + 2 d-projectors for unbound state to improve transferability
    Not recommended for GW.

  * Te-d.in: Include 4d in valence. Default version.
    Not recommended for GW. Ghost at +77 eV.

  * Te-spd-high.in: Include full (4s, 4p, 4d) shell in valence.
    Recommended version for GW and high accuracy calculations.
    TODO Complete tests

## I

  * I.in: Default version for GS applications. Not recommended for GW

  * I-d.in: Include 4d in valence. Not recommended for GW

  * I-spd-high.in: Include full (4s, 4p, 4d) shell in valence.
    Recommended version for GW and high accuracy calculations.
    TODO: Complete tests

## Xe
  Needs GW version

  * Xe.in: Default version with 5s-5p in valence. Two projectors for f to improve trasferability.

## Cs
  Needs GW version

  * Cs-sp.in: Default version with (5s, 5p, 6s) in valence, 4f in core.
    Use two f-projectors to improve transferability.
    Require GW version.

## Ba
  Needs GW version

  * Ba-sp.in:Default version with 5s, 5p in valence. Use ground-state as reference.
    TODO: Ba-sp with f?

  * Ba-sp-exc.in: 
    Ba-sp-new adds two projectors for f to improve deltafactor and GBRV tests.
    Uses excited state (5d1, 6s1) as reference configuration
    seems to have a ghost around 6 eV

## Hf
  Needs a GW version fsp is probably unavoidable but log ders are not optimal.

  * Hf-sp.in:
    Default version with 4f frozen in core. Has 2 f-projectors to improve transferability.
    Recommended for GS applications.

  * Hf-fsp.in:

## Ta 
  Needs GW version

  * Ta-sp.in: Default version with 4f frozen in core. Has 2 f-projectors to improve transferability.
    The missing f, only at ~1 H below the fermi level is probably too much to classify the sp as GW

  * Ta-fsp.in:

## W
  Needs GW version

  * W-sp.in: Default version with  4f frozen in core. Has 2 f-projectors to improve transferability.
    The missing f, only at ~1 H below the fermi level is probably too much to classify the sp as GW

  * W-fsp.in:

## Re:
  Needs GW version
 
  * Re-sp.in: Default version with 4f frozen in core. Has 2 f-projectors to improve transferability.
    The missing f, only at ~1.5 H below the fermi level is probably too much to classify the sp as GW

## Os
  Needs GW version

  * Os-sp.in:  Default version with 4f frozen in core. Has 2 f-projectors to improve transferability.

## Ir
  Needs GW version

  * Ir-sp.in: Default version with 4f frozen in core. Has 2 f-projectors to improve transferability.
    Recommended for GS applications.
    removed the GW tag, 4f to close

## Pt

  Needs GW version

  * Pt-sp.in: Default version with 2 projectors for unbound f states
    Recommended for GS applications.

  * Pt-fsp.in:
    Experimental version!

## Au

  Needs GW version

  * Au-sp.in: Default version with  5s, 5p, 5d, 6s in valence,  4f frozen in core. Has 1 f-projector for unbound 5f to improve transferability.
    Recommended for GS applications.
    We actually did GW calculations for gold and both SigmaX and SigmaC are seriously affected by unfreezing
    the 4f electrons ~10% total effect. 

## Hg:

  * Hg.in: (5d, 6s) in valence.
    Not recommended for GW. Ghost at +66 eV 

  * Hg-sp.in: Default version with (5s, 5p, 5d, 6s) in valence. Add f projector
    Recommended for GS and GW.

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
  Needs GW version

  * Rn.in: Version wtih (6s, 6p) in valence and 5d frozen
    Not recommended for GW. d is only 1.7H below ef

  * Rn-d.in: Default version for GS applications with (5d, 6s, 6p) in valence 
    Not recommended for GW. 

## TODO
    * test for the 6th row the addition of an explicit f projector
    * Pd --> Cd (look at Pd/Ag/Cd carefully)
    
    Cu-sp high: TO BE INVESTIGATED:
    Test that there's no fr pseudo with local part taken from v(L)

    W-sp, Re-sp (f channel could be improved)

    In-d is very sensitive, changing XC (GGa-->LDA) gives a ghost very close to 0.
    The first s projectors presents oscillations I don't see in the other pseudos.

    add hints for 
    Te-spd-high
    Pt-spd


<!--
GGA_X_PBE_SOL = 116 GGA_C_PBE_SOL = 133
LDA_X = 1 LDA_C_PW = 12
7: dict(x=xcf.LDA_X, c=xcf.LDA_C_PW),          # PW 001012
-->
 pseudo_dojo/pseudos/ONCVPSP-PBE-PDv0.3/F/F-new.djrepo        |      5 -
La/La-sp.in:# lloc lpopt rc5 dvloc0
La/La-sp.in-4 5 1.6 2.0

In-d, Sb-d
