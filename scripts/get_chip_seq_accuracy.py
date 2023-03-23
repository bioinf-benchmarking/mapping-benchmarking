import numpy as np
import bionumpy as bnp

truth = bnp.open(snakemake.input.truth).read().sort_by("start")
peaks = bnp.open(snakemake.input.peaks).read()

print("%d peaks found" % len(peaks))

# pick the best scoring peaks (as many as in truth)
peaks = peaks[np.argsort(-peaks.score)][0:len(truth)].sort_by("start")

print(truth)
print(peaks)

n_found = 0
for peak in peaks:
    for truth_peak in truth:
        truth_summit = truth_peak.start + (truth_peak.stop-truth_peak.start)/2
        peak_summit = peak.start + peak.summit
        if abs(truth_summit-peak_summit) < 200:
            n_found += 1
            break

accuracy = n_found / len(truth)
print("accuracy", accuracy)
with open(snakemake.output.result, "w") as f:
    f.write(str(accuracy) + "\n")
