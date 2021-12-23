import os
import glob
import shutil
import fnmatch
from rdkit import Chem


parent_dir = os.getcwd()


class pKaCalculator:
    def __init__(
        self,
        path,
        charge,
        spin,
        dry_run: bool = False,
        conformers: bool = True,
        tautomers: bool = True,
    ):

        self.path = os.path.dirname(path) + "/"
        self.molecule = os.path.basename(path).strip(".xyz")
        self.moleculename = os.path.basename(path).strip(".xyz")
        self.workdir = self.path + "../molecules/" + self.molecule + "/"

        self.charge = charge
        self.spin = spin

        self.dry_run = dry_run
        self.conformers = conformers
        self.tautomers = tautomers

        self.radicalcation_energy: float
        self.neutralradical_energy: float
        self.neutralsinglet_energy: float

        self.deprotomer: str

        self.pKa: float

    def cleanup(self):

        blacklist = [
            "coord",
            "coord.original",
            "fixpositions",
            "struc.xyz",
            "tautomerize_*.xyz",
            "deprotonate_*.xyz",
            "wbo",
            ".UHF",
            ".CHRG",
            ".history.xyz",
            ".xtboptok",
            "charges",
            "cregen_0.tmp",
            "crest.energies",
            "cre_members",
            "g98.out",
            "gfnff_charges",
            "gfnff_topo",
            "hessian",
            "vibspectrum",
            "xtbopt.log",
            "xtbrestart",
            "*.inp",
            "*.cpcm",
            "*.densities",
            "*.gbw",
            "*.ges",
            "*.smd.out",
            "*property.txt",
        ]

        for root, dirs, files in os.walk("."):
            for file in files:
                for item in blacklist:
                    if fnmatch.fnmatch(file, item):
                        os.remove(f"{root}/{file}")

    def tautomer_search(self, charge, spin):

        workdir = self.workdir + "tautomer_search/"

        if not os.path.exists(workdir):
            os.makedirs(workdir)
        os.chdir(workdir)

        if self.dry_run is False:
            os.system(
                f"crest {self.path}{self.molecule}.xyz --gfn2 --chrg {charge} --uhf {spin} --alpb water --mquick --fstrict --tautomerize > {self.molecule}_taut.out 2>> {self.molecule}_taut.out"
            )

        # Looking for cyclization
        os.system(
            f"crest --testtopo {self.path}{self.molecule}.xyz > start_topo.out 2>> start_topo.out"
        )
        os.system(
            f"crest --testtopo tautomers.xyz > tautomers_topo.out 2>> tautomers_topo.out"
        )
        with open("start_topo.out", "r") as out:
            start_rings = 0
            for line in out:
                if "Total number of rings in the system" in line:
                    start_rings = int(line.split()[-1])
                    break
        with open("tautomers_topo.out", "r") as out:
            tautomer_rings = 0
            for line in out:
                if "Total number of rings in the system" in line:
                    tautomer_rings = int(line.split()[-1])
                    break
                if "No (valid) input file!" in line:
                    print(
                        f"INFO: No topomers found for {self.molecule}. Ignoring tautomer search."
                    )
                    break

        if start_rings < tautomer_rings:
            print(
                f"WARNING: Cyclization spotted for {self.molecule}. Ignoring tautomer search.\n"
            )

        else:
            shutil.copyfile("tautomers.xyz", "../" + self.molecule + "_taut.xyz")

            self.path = self.workdir
            self.molecule = self.molecule + "_taut"

        self.cleanup()

        os.chdir(parent_dir)

        return

    def optimization(self, type, charge, spin):

        workdir = self.workdir + type + "/"

        if not os.path.exists(workdir):
            os.makedirs(workdir)
        os.chdir(workdir)

        if self.dry_run is False:
            os.system(
                f"xtb {self.path}{self.molecule}.xyz --gfn2 --chrg {charge} --uhf {spin} --alpb water --ohess  > {self.molecule}_opt.out 2>> {self.molecule}_opt.out"
            )
            if self.conformers is True:
                os.system(
                    f"crest xtbopt.xyz --gfn2 --chrg {charge} --uhf {spin} --alpb water --mquick > {self.molecule}_conf.out 2>> {self.molecule}_conf.out"
                )
                os.system(
                    f"xtb crest_best.xyz --gfn2 --chrg {charge} --uhf {spin} --alpb water --ohess  > {self.molecule}_opt.out 2>> {self.molecule}_opt.out"
                )

        mol_file = [f for f in os.listdir(".") if f.endswith(".mol")][-1]

        end_mol = Chem.MolFromMolFile(
            mol_file, sanitize=False, removeHs=False, strictParsing=False
        )
        end_smiles = Chem.MolToSmiles(end_mol)

        with open(self.molecule + "_opt.out", "r") as out:
            for line in out:
                if "TOTAL FREE ENERGY" in line:
                    energy = float(line.split()[-3])

        if "." in end_smiles:
            print(
                f"ERROR: {self.molecule} (charge {charge} uhf {spin}) has undergone dissociation.\n"
            )

        self.cleanup()

        os.chdir(parent_dir)

        return energy

    def single_point_B97(self, type, charge, spin):

        workdir = self.workdir + type + "/"

        if not os.path.exists(workdir):
            os.makedirs(workdir)
        os.chdir(workdir)

        with open(f"{self.molecule}.inp", "w") as inp:
            inp.write(
                f"""! RIJCOSX B97-D3 D3BJ
! def2-TZVP def2/J

%CPCM
    SMD True
    SMDsolvent "water"
end

* xyzfile {charge} {spin+1} xtbopt.xyz
"""
            )

        if self.dry_run is False:
            os.system(f"$ORCADIR/orca {self.molecule}.inp > {self.molecule}_orca.out")

        with open(self.molecule + "_opt.out", "r") as out:
            for line in out:
                if "G(RRHO) contrib." in line:
                    vib = float(line.split()[-3])

        with open(self.molecule + "_orca.out", "r") as out:
            for line in out:
                if "FINAL SINGLE POINT ENERGY" in line:
                    scf = float(line.split()[-1])

        energy = scf + vib

        self.cleanup()

        os.chdir(parent_dir)

        return energy

    def deprotonate(self, charge, spin):

        workdir = self.workdir + "deprotonation/"

        if not os.path.exists(workdir):
            os.makedirs(workdir)
        os.chdir(workdir)

        if self.dry_run is False:
            os.system(
                f"crest {self.workdir}neutral_singlet/xtbopt.xyz --gfn2 --chrg {charge} --uhf {spin} --alpb water -deprotonate > {self.molecule}_deprot.out 2>> {self.molecule}_deprot.out"
            )

        self.evaluate_deprotomers()

        for deprotomer in glob.glob("deprotomer*.xyz"):
            shutil.copyfile(deprotomer, "../" + deprotomer)

        self.path = self.workdir

        self.cleanup()

        os.chdir(parent_dir)

        return

    def evaluate_deprotomers(self):

        with open("deprotonated.xyz", "r") as f:
            molsize = int(f.readline().strip())
        with open("deprotonated.xyz", "r") as f:
            numlines = int(sum(1 for line in f))

        j = 0
        deprotomer = 1
        with open("deprotonated.xyz", "r") as f:
            while j < numlines:
                i = 0
                with open(f"deprotomer_{deprotomer}.xyz", "w") as out:
                    while i < molsize + 2:
                        out.write(f.readline())
                        i += 1
                        j += 1
                deprotomer += 1

    def try_calculate_pka(self, return_self: bool = True):

        # tautomer search
        if self.tautomers is True:
            self.tautomer_search(self.charge, self.spin)

        # radical cation optimization
        self.radicalcation_energy = self.optimization(
            "radical_cation", self.charge + 1, self.spin + 1,
        )

        # neutral singlet optimization
        self.neutralsinglet_energy = self.optimization(
            "neutral_singlet", self.charge, self.spin,
        )

        # deprotonating neutral molecule
        self.deprotonate(self.charge, self.spin)

        # neutral radicals optimization and iterative pKa calculation
        self.neutralradical_energies = []

        pka_correction = 164.22  # kcal/mol
        self.pKas = []

        deprotomers = sorted(glob.glob(f"{self.workdir}deprotomer*.xyz"))
        best_deprotomer_energy = 0

        for deprotomer in deprotomers:

            self.molecule = os.path.basename(deprotomer).strip(".xyz")

            neutralradical_energy = self.optimization(
                f"neutral_radical_{self.molecule[-1]}", self.charge, self.spin + 1,
            )

            pKa = -(
                (self.radicalcation_energy - neutralradical_energy) * 627.5
                + 270.29
                - pka_correction
            ) / (2.303 * 1.98720425864083 / 1000 * 298.15)

            if len(deprotomers) > 1:

                new_deprotomer_energy = self.single_point_B97(
                    f"neutral_radical_{self.molecule[-1]}", self.charge, self.spin + 1,
                )

                if new_deprotomer_energy < best_deprotomer_energy:

                    best_deprotomer_energy = new_deprotomer_energy

                    self.deprotomer = self.molecule
                    self.neutralradical_energy = neutralradical_energy
                    self.pKa = pKa

            else:
                self.deprotomer = self.molecule
                self.neutralradical_energy = neutralradical_energy
                self.pKa = pKa

            if pKa > 14:
                break

        if return_self is True:
            return self.moleculename, self.pKa, self
        else:
            return self.moleculename, self.pKa

    def calculate_pka(self, return_self: bool = True):
        try:
            return self.try_calculate_pka(return_self)

        except Exception as error:
            print(f"ERROR: {self.moleculename}! {error}")
            return None
