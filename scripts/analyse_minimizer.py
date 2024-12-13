import argparse
import matplotlib.pyplot as plt


ALPHABET = "ACTG"


def main():
    parser = argparse.ArgumentParser(description="Detect length of extended minimizers")

    parser.add_argument(
        "input_file",  # name of the argument
        type=str,  # type of the argument
        help="input file containing reads with the minimizer in it",
    )

    parser.add_argument(
        "minimizer",  # name of the argument
        type=str,  # type of the argument
        help="minimizer",
    )

    args = parser.parse_args()  # Parse the arguments

    input_file = args.input_file
    minimizer = args.minimizer

    match_reads = []
    print("Collecting...", end="")
    with open(input_file, "r") as fichier:
        for ligne in fichier:
            ligne = ligne.strip()
            if minimizer not in ligne:
                # print("error")
                # exit()
                pass
            match_reads.append(ligne)
    print(" Done.")
    nb_match = []
    nb_match.append(len(match_reads))

    while True:
        potential_extended_minimizers = []
        counts = []
        for letter_start in ALPHABET:
            for letter_end in ALPHABET:
                potential_extended_minimizers.append(
                    letter_start + minimizer + letter_end
                )
                counts.append(0)

        index_to_remove = []
        for index, match_read in enumerate(match_reads):
            found = False
            for i, candidate in enumerate(potential_extended_minimizers):
                if candidate in match_read:
                    counts[i] += 1
                    found = True
            if not found:
                index_to_remove.append(index)

        # remove elements in a list
        # we iterate in reverse, put the leement to remove at the end, then remove it
        # (fatser to remove en element at the end of a list in python)
        for index in index_to_remove[::-1]:
            if index != len(match_reads):
                match_reads[index], match_reads[-1] = (
                    match_reads[-1],
                    match_reads[index],
                )
                match_reads.pop(-1)

        max_nb_match = max(counts)
        best_extension_id = counts.index(max_nb_match)
        best_extension = potential_extended_minimizers[best_extension_id]

        nb_match.append(max_nb_match)

        minimizer = best_extension
        # print(best_extension)
        print(max_nb_match)
        # print()
        if max_nb_match == 0:
            break
    # print(nb_match)
    minimizer_size = [len(args.minimizer) + 2 * x for x in range(len(nb_match))]
    plt.plot(minimizer_size, nb_match)
    plt.title("number of match when extending a highly present minimizer")
    plt.xlabel("size of the minimizer")
    plt.ylabel("number of match")
    plt.show()


if __name__ == "__main__":
    main()
