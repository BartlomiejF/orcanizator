from typing import List
from rdkit import Chem
import pandas as pd
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
import subprocess


ORCA_EXECUTABLE = "orca"


def orca_generator(smiles: str, nb: int) -> None:
    """Generate ORCA input file from SMILES. ORCA arguments are hardcoded into
    the function.

    Arguments:
        smiles {str} -- SMILES structure of the molecule.
        nb {int} -- number. This will be the name of the input file.
    """

    molecule = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(molecule)

    AllChem.EmbedMolecule(mol)

    output = Chem.rdmolfiles.MolToXYZBlock(mol)

    charge = Chem.GetFormalCharge(mol)

    spin_multiplicity = Descriptors.NumRadicalElectrons(mol)*2*0.5 + 1

    newline = "\n"
    file_output = f"""! DFT B3LYP/G RIJCOSX def2-SVP Opt FREQ BOHRS xyzfile

    %scf
    MaxIter 500
    end

    %geom
    MaxIter 500
    end

    %elprop Dipole true
    Polar 1
    end

    * xyz {charge} {spin_multiplicity}
    {newline.join(output.split(newline)[2:-1])}
    *"""

    with open(f"{nb}.inp", "w") as f:
        f.write(file_output)


def orca_reader(smiles: str, nb: int) -> dict:
    """Read following properties of ORCA output file:
    1. HOMO energy
    2. LUMO energy
    3. Dipole moment
    4. Isotropic polarizability
    5. Total thermal energy, U
    6. Total enthalpy, H
    7. Entropy term, S
    8. Gibbs free energy, G
    9. Rotational constant, A
    10. Rotational constant, B
    11. Rotational constant, C

    Arguments:
        smiles {str} -- SMILES structure of the molecule.
        nb {int} -- number. This is the name of the output file to be read.

    Returns:
        {dict} -- a dictionary of read properties.
    """

    with open(f"{nb}.out", encoding="utf-16") as f:
        output = f.readlines()

    props = {"smiles": smiles}

    parse = False
    orbitals_start = False
    for line in output:
        if "ORBITAL ENERGIES".encode("utf-16").decode("utf-16") in line:
            parse = True
            continue

        if "E(Eh)".encode("utf-16").decode("utf-16") in line:
            orbitals_start = True
            continue

        if parse & orbitals_start:
            splitted = line.split()

            if float(splitted[1]) > 0:
                lumo_Eh = splitted[2]
                lumo_eV = splitted[3]
                homo_Eh = splitted[2]
                homo_eV = splitted[3]
            else:
                lumo_Eh = splitted[2]
                lumo_eV = splitted[3]
                parse = False
                orbitals_start = False

            props["LUMO (Eh)"] = lumo_Eh
            props["LUMO (eV)"] = lumo_eV
            props["HOMO (Eh)"] = homo_Eh
            props["HOMO (eV)"] = homo_eV

        if "DIPOLE MOMENT".encode("utf-16").decode("utf-16") in line:
            parse = True

        if parse & ("Magnitude (Debye)".encode("utf-16").decode("utf-16") in line):
            moment = line.split()[3]
            props["Dipole Moment"] = moment

        if parse & ("Isotropic polarizability".encode("utf-16").decode("utf-16") in line):
            polarizability = line.split()[3]
            props["polarizability"] = polarizability
            parse = False

        if "INNER ENERGY".encode("utf-16").decode("utf-16") in line:
            parse = True

        if parse & ("Total thermal energy".encode("utf-16").decode("utf-16") in line):
            U = line.split()[3]
            props["U"] = U

        if parse & ("Total Enthalpy".encode("utf-16").decode("utf-16") in line):
            H = line.split()[3]
            props["H"] = H

        if parse & ("Final entropy term".encode("utf-16").decode("utf-16") in line):
            S = line.split()[4]
            props["S"] = S

        if parse & ("Final Gibbs".encode("utf-16").decode("utf-16") in line):
            G = line.split()[5]
            props["G"] = G

        if parse & ("Rotational constants in cm-1:".encode("utf-16").decode("utf-16") in line):
            A, B, C = line.split()[4:7]
            props["A"] = A
            props["B"] = B
            props["C"] = C
    return props


def orcanize(smiles: str, nb: int) -> dict:
    """This function does the following steps:
    1. Generate ORCA input file based on molecule SMILES by calling
    orca_generator function.
    2. Call ORCA job with generated input file.
    3. Read calculated properties from generated output file and return
    them as dictionary. This is done by calling orca_reader function.

    SMILES of successfully processed molecules are saved into done.log file.
    The done.log also contains the number which is input and output files name.
    If there was an error in calculations, then the molecule's SMILES is
    saved into undone.log file.

    Arguments:
        smiles {str} -- [SMILES molecule structure]
        nb {int} -- [Number. This is ORCA input and output file name]

    Returns:
        {dict} -- [Dictionary of properties]
    """

    #  generate input file
    orca_generator(smiles, nb)

    #  run ORCA job
    subprocess.run([ORCA_EXECUTABLE, f"{smiles}.inp", ">", f"{smiles}.out"])

    #  try to read the calculated data; if there is an error(mostly because
    # of calculation error, then write number and smiles to undone.log file)
    try:
        data = orca_reader(smiles, nb)
        with open("done.log", "a") as f:
            f.write(f"{nb} {smiles}\n")
        return data
    except:
        with open("undone.log", "a") as f:
            f.write(f"{nb} {smiles}\n")


def orcanize_many(
                molecules: List[str],
                save: bool = True
                ) -> pd.DataFrame:
    """Run orcanize function with a collection of SMILES.

    Arguments:
        molecules {Union[str]} -- collection of SMILES strings.
        save {bool} -- save the pandas dataframe to csv(True, default) or
                    not(False); by convention it is saved everytime ORCA
                    successfully finishes calculation

    Returns:
        {pd.DataFrame} -- pandas dataframe with read properties
    """
    outdata = pd.DataFrame()
    for nb, smiles in enumerate(molecules):
        print(nb,smiles)
        dat = orcanize(smiles, nb)
        outdata = outdata.append(pd.DataFrame(dat, index=[0]), ignore_index=True)
        if save:
            outdata.to_csv("orkanized.csv")
    return outdata
