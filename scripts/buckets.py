import matplotlib.pyplot as plt

import yaml
import argparse


def plot_minimizers_t1(data, colors, les_x, k, m, error_rate):
    for i, x in enumerate(m):
        plt.plot(
            les_x,
            data[error_rate][k][x]["minimizers"],
            "o--",
            color=colors[i],
            label=f"number of minimizers (m = {x})",
        )


def plot_entries_t1(data, colors, les_x, k, m, error_rate):
    for i, x in enumerate(m):
        plt.plot(
            les_x,
            data[error_rate][k][x]["entries"],
            "x--",
            color=colors[i],
            label=f"total entries (m = {x})",
        )


def main():

    parser = argparse.ArgumentParser(
        description="Show minimizer repartition in KFC buckets"
    )

    parser.add_argument(
        "input_file",  # name of the argument
        type=str,  # type of the argument
        help="input file containing the minimizer repartition",
    )

    args = parser.parse_args()  # Parse the arguments

    input_file = args.input_file

    with open(input_file, "r") as fichier:
        data = yaml.safe_load(fichier)

        error_rate = 0.1
        k = 99
        m = [0, 10, 30]
        colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]

        les_x = range(len(data[error_rate][k][m[0]]["minimizers"]))
        max_minimizers = max((max(data[error_rate][k][x]["minimizers"]) for x in m))
        max_entries = max((max(data[error_rate][k][x]["entries"]) for x in m))

    plot_entries_t1(data, colors, les_x, k, m, error_rate)
    plot_minimizers_t1(data, colors, les_x, k, m, error_rate)

    plt.xlim([0, len(les_x) - 1 + 0.5])
    plt.ylim([0, max_entries + 0.5])
    plt.xlabel("bucket ID")
    plt.ylabel("number of occurences")

    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
