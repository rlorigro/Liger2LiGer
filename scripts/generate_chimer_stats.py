from collections import defaultdict
from modules.IterativeHistogram import IterativeHistogram
from matplotlib import pyplot
import matplotlib
import argparse
import numpy
import os

matplotlib.use("Agg")


class ChimerStats:
    def __init__(self, name, n50, n_chimers, n_non_chimers):
        self.name = name
        self.n50 = n50
        self.n_chimers = n_chimers
        self.n_non_chimers = n_non_chimers

    def __str__(self):
        s = self.name + ',' + str(int(self.n50)) + ',' + str(self.n_chimers) + ',' + str(self.n_non_chimers)

        return s

    def get_header(self):
        s = "name" + ',' + "n50" + ',' + "n_chimers" + ',' + "n_non_chimers"

        return s


def generate_chimer_stats(paths):

    path_pairs = defaultdict(lambda: [None,None])

    results = list()

    for path in paths:
        name = '_'.join(os.path.basename(path).split('.')[:-2])

        print(name, "non_chimer_lengths" in path, "chimer_lengths" in path)

        if "non_chimer_lengths" in path:
            if path_pairs[name][0] is None:
                path_pairs[name][0] = path
            else:
                exit("ERROR: path with name prefix already encountered: " + path_pairs[name][0])
        elif "chimer_lengths" in path:
            if path_pairs[name][1] is None:
                path_pairs[name][1] = path
            else:
                exit("ERROR: path with name prefix already encountered: " + path_pairs[name][1])
        else:
            exit("ERROR: path lacks identifying suffix: " + path)

    for name,pair in path_pairs.items():
        r = ChimerStats(name, None, 0, 0)

        name = name.split("_VS_")[0]
        print(name)

        non_chimer_histogram = IterativeHistogram(start=0, stop=500_000, n_bins=2000)
        chimer_histogram = IterativeHistogram(start=0, stop=500_000, n_bins=2000)

        with open(pair[0], 'r') as file:
            for line in file:
                l = int(line.strip())

                non_chimer_histogram.update(l)

                r.n_non_chimers += 1

        with open(pair[1], 'r') as file:
            for line in file:
                l = int(line.strip())

                chimer_histogram.update(l)

                r.n_chimers += 1

        # -------

        # find N50

        all_frequencies = chimer_histogram.get_histogram() + non_chimer_histogram.get_histogram()
        bins = chimer_histogram.get_bin_centers()
        total_bp = sum([bins[i]*all_frequencies[i] for i in range(len(all_frequencies))])

        cumulative_sum = 0
        for i in range(len(all_frequencies)):
            cumulative_sum += bins[i]*all_frequencies[i]

            if cumulative_sum > total_bp/2:
                print("N50\t" + str(bins[i]))
                r.n50 = bins[i]
                break

        results.append(r)

        # -------

        f = pyplot.figure()
        axes = pyplot.axes()

        chimer_frequencies = chimer_histogram.get_histogram()
        x = chimer_histogram.get_bin_centers()
        axes.bar(x=x, height=chimer_frequencies, width=x[1]-x[0], color="C1")

        non_chimer_frequencies = non_chimer_histogram.get_histogram()
        x = non_chimer_histogram.get_bin_centers()
        axes.bar(x=x, height=non_chimer_frequencies, bottom=chimer_frequencies, width=x[1]-x[0], color="C0")

        axes.set_title(name + " raw")
        axes.set_xlim([0,200_000])
        axes.set_ylim([0,numpy.mean(sorted(non_chimer_frequencies, reverse=True)[:3])])

        output_path = os.path.join(os.path.dirname(pair[0]),name + "_chimer_distribution_raw.png")

        pyplot.legend(["chimers","non_chimers"])

        f.set_size_inches([10,8])
        axes.set_xlabel("Read Length (bp)")
        axes.set_ylabel("Frequency (#)")

        pyplot.savefig(output_path, dpi=200)
        pyplot.close()

        # ----------------

        f = pyplot.figure()
        axes = pyplot.axes()

        x = chimer_histogram.get_bin_centers()
        y = 100*(numpy.array(chimer_frequencies) + 0.000001)/(numpy.array(non_chimer_frequencies) + numpy.array(chimer_frequencies) + 0.001)
        axes.bar(x=x, height=y, width=x[1]-x[0], color="C1")

        axes.set_title(name + " percentages")
        axes.set_xlim([0,200_000])
        axes.set_ylim([0,100])

        output_path = os.path.join(os.path.dirname(pair[0]),name + "_chimer_distribution_percent.png")

        pyplot.legend(["chimers","non_chimers"])

        f.set_size_inches([10,8])
        axes.set_xlabel("Read Length (bp)")
        axes.set_ylabel("Frequency (%)")

        pyplot.savefig(output_path, dpi=200)
        pyplot.close()

        # ----------------

        f = pyplot.figure()
        axes = pyplot.axes()

        chimer_frequencies = chimer_histogram.get_histogram()
        x = chimer_histogram.get_bin_centers()
        chimer_frequencies *= x
        axes.bar(x=x, height=chimer_frequencies, width=x[1]-x[0], color="C1")

        non_chimer_frequencies = non_chimer_histogram.get_histogram()
        x = non_chimer_histogram.get_bin_centers()
        non_chimer_frequencies *= x
        axes.bar(x=x, height=non_chimer_frequencies, bottom=chimer_frequencies, width=x[1]-x[0], color="C0")

        axes.set_title(name + " coverage")
        axes.set_xlim([0,200_000])

        output_path = os.path.join(os.path.dirname(pair[0]),name + "_chimer_distribution_coverage.png")

        pyplot.legend(["chimers","non_chimers"])

        f.set_size_inches([10,8])
        axes.set_xlabel("Read Length (bp)")
        axes.set_ylabel("Frequency (coverage in bp)")

        pyplot.savefig(output_path, dpi=200)
        pyplot.close()

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input","-i",
        type=str,
        required=True,
        help="Comma separated paths of length txt files containing lengths of chimers and non-chimers \
              (one length [in bp] per line). Can run any number of pairs of txt files as long as they have the suffix \
              non_chimer_lengths.txt or chimer_lengths.txt and pairs share a prefix."
        )

    args = parser.parse_args()

    generate_chimer_stats(paths=args.input)
