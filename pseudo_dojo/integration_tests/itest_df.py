"""Integration tests for deltafactor calculations with Abinit and AbiPy."""
import pytest
import pseudo_dojo.data as pdj_data
import abipy.abilab as abilab

from pseudo_dojo.dojo.works import DeltaFactory, GbrvFactory


def itest_deltafactor_gga_pawxml(fwp, tvars):
    """
    Testing the flow used for the computation of the deltafactor with PAW and GGA XC.
    """
    # Path of the pseudopotential to test.
    pseudo = pdj_data.pseudo("Si.GGA_PBE-JTH-paw.xml").as_tmpfile(tmpdir=fwp.workdir)
    assert pseudo is not None
    assert pseudo.has_dojo_report
    assert not pseudo.dojo_report.exceptions

    flow = abilab.Flow(workdir=fwp.workdir, manager=fwp.manager)

    # Build the workflow for the computation of the deltafactor.
    # The workflow will produce a pdf file with the equation of state
    # and a file deltafactor.txt with the final results in the
    # outdir directory DELTAFACTOR/work_0/outdir.
    kppa = 20  # this value is for testing purpose (6570 is the correct one)
    ecut = 2
    pawecutdg = ecut * 2 if pseudo.ispaw else None

    assert not pseudo.dojo_report.has_trial("deltafactor", ecut=ecut)
    work = DeltaFactory(pseudo.xc).work_for_pseudo(pseudo,
                                          kppa=kppa, ecut=ecut, pawecutdg=pawecutdg,
                                          include_soc=False,
                                          paral_kgb=tvars.paral_kgb)

    # Register the workflow.
    flow.register_work(work)
    flow.build_and_pickle_dump(abivalidate=True)

    for task in flow[0]:
        task.start_and_wait()

    flow.check_status(show=True)
    assert flow.all_ok
    assert all(work.finalized for work in flow)
    results = flow[0].get_results()

    #20.453 ang^3 88.545 GPa 4.31 20.8658081501 336.680999051 GPa -35.681897152
    #delta Equation of State: deltafactor_polyfit
    #Minimum volume = 20.87 Ang^3
    #modulus = 2.10 eV/Ang^3 = 336.68 GPa, b1 = -35.68
    #Deltafactor = 15.681 meV
    assert pseudo.dojo_report.has_trial("deltafactor", ecut=ecut)
    assert not pseudo.dojo_report.exceptions


def itest_deltafactor_gga_ncsoc(fwp, tvars):
    """
    Testing the flow used for the computation of the deltafactor with GGA and NC+SOC.
    """
    # return
    # Path of the pseudopotential to test.
    pseudo = pdj_data.pseudo("Pb-d-3_r.psp8").as_tmpfile(tmpdir=fwp.workdir)
    assert pseudo.has_dojo_report
    assert pseudo.supports_soc
    assert not pseudo.dojo_report.exceptions

    flow = abilab.Flow(workdir=fwp.workdir, manager=fwp.manager)

    # Build the workflow for the computation of the deltafactor.
    # The workflow will produce a pdf file with the equation of state
    # and a file deltafactor.txt with the final results in the
    # outdir directory DELTAFACTOR/work_0/outdir.
    kppa = 40  # this value is for testing purpose (6570 is the correct one)
    ecut = 8
    pawecutdg = ecut * 2 if pseudo.ispaw else None

    assert pseudo.dojo_report.has_trial("deltafactor", ecut=12)
    assert not pseudo.dojo_report.has_trial("deltafactor_soc", ecut=ecut)

    work = DeltaFactory(pseudo.xc).work_for_pseudo(pseudo,
                                          kppa=kppa, ecut=ecut, pawecutdg=pawecutdg,
                                          include_soc=True,
                                          paral_kgb=tvars.paral_kgb)

    # Register the workflow.
    flow.register_work(work)
    flow.build_and_pickle_dump(abivalidate=True)

    for task in flow[0]:
        task.start_and_wait()

    flow.check_status(show=True)
    assert flow.all_ok
    assert all(work.finalized for work in flow)
    results = flow[0].get_results()

    assert pseudo.dojo_report.has_trial("deltafactor_soc", ecut=ecut)
    assert not pseudo.dojo_report.exceptions
