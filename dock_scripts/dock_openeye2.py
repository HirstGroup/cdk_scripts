import sys

from openeye import oechem
from openeye import oedocking


def main(argv=[__name__]):
    dockOpts = oedocking.OEDockOptions()
    opts = oechem.OERefInputAppOptions(dockOpts, "DockMolecules", oechem.OEFileStringType_Mol3D,
                                       oechem.OEFileStringType_Mol3D, oechem.OEFileStringType_DU, "-receptor")
    if oechem.OEConfigureOpts(opts, argv, False) == oechem.OEOptsConfigureStatus_Help:
        return 0
    dockOpts.UpdateValues(opts)

    ifs = oechem.oemolistream()
    if not ifs.open(opts.GetInFile()):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % opts.GetInFile())

    rfs = oechem.oeifstream()
    if not rfs.open(opts.GetRefFile()):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % opts.GetRefFile())

    ofs = oechem.oemolostream()
    if not ofs.open(opts.GetOutFile()):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % opts.GetOutFile())

    du = oechem.OEDesignUnit()
    if not oechem.OEReadDesignUnit(rfs, du):
        oechem.OEThrow.Fatal("Failed to read design unit")
    if not du.HasReceptor():
        oechem.OEThrow.Fatal("Design unit %s does not contain a receptor" % du.GetTitle())

    dock = oedocking.OEDock(dockOpts)
    dock.Initialize(du)

    for mcmol in ifs.GetOEMols():
        print("docking", mcmol.GetTitle())
        dockedMol = oechem.OEGraphMol()
        retCode = dock.DockMultiConformerMolecule(dockedMol, mcmol)
        if (retCode != oedocking.OEDockingReturnCode_Success):
            oechem.OEThrow.Fatal("Docking Failed with error code " + oedocking.OEDockingReturnCodeGetName(retCode))

        sdtag = oedocking.OEDockMethodGetName(dockOpts.GetScoreMethod())
        oedocking.OESetSDScore(dockedMol, dock, sdtag)
        dock.AnnotatePose(dockedMol)
        oechem.OEWriteMolecule(ofs, dockedMol)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))