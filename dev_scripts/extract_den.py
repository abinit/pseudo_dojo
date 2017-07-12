#!/usr/bin/env python
from __future__ import division, print_function, unicode_literals

import sys
import os
import numpy as np

def psp8_get_densities(path):
    """
    Extract densities from a psp8 files.
    """
    from pymatgen.io.abinit.pseudos import _dict_from_lines

    # !(1) title (character) line
    # !(2) znucl,zion,pspdat
    # !(3) pspcod,pspxc,lmax,lloc,mmax,r2well  (r2well not used)
    # !(4) rchrg,fchrg,qchrg  (fchrg /=0 if core charge, qchrg not used)
    # !(5) nproj(0:lmax)  (several projectors allowed for each l)
    # !Then, for ll=0,lmax :
    # !if(nproj(ll)>0)
    # !1/<u1|vbkb1>, 1/<u2|vbkb2>, ...
    # !for  irad=1,mmax  : irad, r(irad), vbkb1(irad,ll), vbkb2(irad,ll), ...
    # !elseif ll=lloc
    # !for  irad=1,mmax  : irad, r(irad), vloc(irad)
    # !end if
    # !
    # !If(lloc>lmax)
    # !for  irad=1,mmax  : irad, r(irad), vloc(irad)
    # !end if
    # !
    # !vbkb are Bloechl-Kleinman-Bylander projectors,(vpsp(r,ll)-vloc(r))*u(r,ll),
    # !unnormalized
    # !Note that an arbitrary local potential is allowed.  Set lloc>lmax, and
    # !provide projectors for all ll<=lmax
    # !
    # !Finally, if fchrg>0,
    # !for  irad=1,mmax  : irad, r(irad), xccc(irad),
    # !xccc'(irac), xccc''(irad), xccc'''(irad), xccc''''(irad)
    # !
    # !Model core charge for nonlinear core xc correction, and 4 derivatives

    with open(path, "rt") as fh:
        lines = [fh.readline() for _ in range(6)]
        header = _dict_from_lines(lines[1:3], [3, 6])

        # Number of points on the linear mesh.
        mmax = int(header["mmax"])

        # Non-linear core correction parameters.
        rchrg, fchrg, qchrg = [float(t) for t in lines[3].split()[:3]]

        # Read Number of projectors(l) and extension switch
        nproj = [int(t) for t in lines[4].split()[:5]]
        assert len(nproj) == 5
        tokens = lines[5].split()
        extension_switch = int(tokens[0])

        # Old format, Densities are not available.
        if len(tokens) == 1:
            print("psp8 file does not contain density records")
            return {}

        den_flag = int(tokens[1])
        if den_flag != 1:
            raise ValueError("Expecting den_flag 1 but got %s" % den_flag)

        # Read SOC projectors
        has_soc = extension_switch in (2, 3)
        if has_soc:
            line = fh.readline()
            # Start at l=1
            nproj_soc = [int(t) for t in line.split()[:4]]
            nproj_soc.insert(0, 0)
            print("nproj_soc", nproj_soc)
            raise NotImplementedError("SOC not tested")

        lmax = int(header["lmax"])
        lloc = int(header["lloc"])
        nso = 1 if not has_soc else 2

        # Will now proceed at the reading of pots and projectors
        # rad(:)=radial grid r(i)
        # vpspll(:,1),...,vpspll(:,lnmax)=nonlocal projectors
        # vloc(:)=local potential

        #for nn in range(nso):
        # Skip projectors (scalar relativistic, always present).
        for l, npl in enumerate(nproj):
            if npl == 0 and l != lloc: continue

            line = fh.readline() # l, ekb[:npl]
            l_file = int(line.split()[0])
            if l != l_file:
                print("For l=%s, npl=%s" % (l, npl), "wrong line", line)
                raise RuntimeError("l != l_file (%s != %s)" % (l, l_file))

            for ir in range(mmax):
                fh.readline()

        # Skip local potential.
        if lloc == 4:
            lloc_file = int(fh.readline())
            assert lloc_file == lloc
            for ir in range(mmax):
                fh.readline()

        # Skip model core charge function and derivatives, if present.
        if fchrg > 1e-15:
            for ir in range(mmax):
                fh.readline()

        # Read pseudo valence charge in real space on the linear mesh.
        # [i, r, PS_val, AE_val, AE_core]
        mesh, psval, aeval, aecore= [np.empty(mmax) for _ in range(4)]
        for ir in range(mmax):
            l = fh.readline()
            #print("denline", l)
            findx, rad, v1, v2, v3 = l.split()
            assert ir + 1 == int(findx)
            mesh[ir] = float(rad)
            psval[ir] = float(v1)
            aeval[ir] = float(v2)
            aecore[ir] = float(v3)

        #return {
        #    "psval": RadialDensity(mesh, psval),
        #    "aeval": RadialDensity(mesh, aeval),
        #    "aecore": RadialDensity(mesh, aecore),
        #}


def main():
    path = sys.argv[1]
    if os.path.isfile(path):
	d = psp8_get_densities(path)
	print(d)
    else:
        for dirpath, dirnames, filenames in os.walk(path):
            for f in filenames:
                if not f.endswith(".psp8"): continue
                path = os.path.join(dirpath, f)
                psp8_get_densities(path)


if __name__ == "__main__":
    sys.exit(main())
