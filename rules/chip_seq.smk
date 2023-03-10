

rule macs2_callpeaks:
    input:
        treatment = "{path}/mapped.bam"
    output:
        multiext("{path}/peaks",
            "_peaks.xls",### required
            ### optional output files
            "_peaks.narrowPeak",
            "_summits.bed"
        )
    params:
        "-f BAM -g hs --nomodel -q 0.2"
    wrapper:
        "v1.23.5-27-g1638662a/bio/macs2/callpeak"



rule get_accuracy:
    input:
        peaks=f"data/{reference_genome}/{chip_seq}/peaks_peaks.narrowPeak",
        truth=f"data/{reference_genome}/{chip_seq.until('n_peaks')}/simulated_peaks.bed"
    output:
        peaks=f"data/{reference_genome}/{chip_seq}/accuracy.txt",
    run:
        import numpy as np
        import bionumpy as bnp

        truth = bnp.open(input.truth).read()
        peaks = bnp.open(input.peaks).read()

        # pick the best scoring peaks (as many as in truth)
        peaks = peaks[np.argsort(peaks.score)][0:len(truth)]

        print(truth)
        print(peaks)

        n_found = 0
        for peak in peaks:
            for truth_peak in truth:
                truth_summit = truth_peak.start + (truth_peak.stop-truth_peak.start)/2
                peak_summit = peak.start + peak.summit
                if abs(truth_summit-peak_summit) < 100:
                    n_found += 1
                    break

        accuracy = n_found / len(truth)
        print("accuracy", accuracy)
        with open(output.peaks, "w") as f:
            f.write(str(accuracy) + "\n")
