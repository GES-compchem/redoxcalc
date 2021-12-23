import pickle


def calculate_potentials(pickle_file):

    data = pickle.load(open(pickle_file, "rb"))

    redox_potentials = []

    for entry in data:

        if entry is None:
            continue

        moleculename = entry[2].moleculename
        pKa = entry[2].pKa
        deprotomer = entry[2].deprotomer

        radicalcation_energy = entry[2].radicalcation_energy
        neutralradical_energy = entry[2].neutralradical_energy
        neutralsinglet_energy = entry[2].neutralsinglet_energy

        reference_potential = 4.28  # V
        pH = 0  # pH units
        self_interaction_shift = 4.846  # V
        PCET_shift = 7.12  # V

        for i in range(15):
            if pKa < pH:
                redox_potentials.append(
                    [
                        moleculename,
                        deprotomer,
                        pKa,
                        pH,
                        -(
                            (neutralsinglet_energy - neutralradical_energy) * 627.5
                            + 270.29
                        )
                        * (4184 / 96485)
                        - reference_potential
                        - 0.059 * pH
                        - self_interaction_shift
                        + PCET_shift,
                    ]
                )

            else:
                redox_potentials.append(
                    [
                        moleculename,
                        deprotomer,
                        pKa,
                        pH,
                        -((neutralsinglet_energy - radicalcation_energy) * 627.5)
                        * (4184 / 96485)
                        - reference_potential
                        - self_interaction_shift,
                    ]
                )
            pH += 1

    return redox_potentials

