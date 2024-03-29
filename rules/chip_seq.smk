from mapping_benchmarking.config import CalledPeaks, SimulatedChipSeqPeaks, PeakCallingAccuracy

rule macs2_callpeaks:
    input:
        treatment = "{path}.bam"
    output:
        multiext("{path}",
            "_peaks.xls",### required
            ### optional output files
            "_peaks.narrowPeak",
            "_summits.bed"
        )
    params:
        "-f BAM -g hs --nomodel -q 0.05"
    wrapper:
        "v1.23.5-27-g1638662a/bio/macs2/callpeak"



rule get_accuracy:
    input:
        #peaks=f"data/{reference_genome}/{chip_seq}/peaks_peaks.narrowPeak",
        peaks = CalledPeaks.path() + "_peaks.narrowPeak",
        truth = SimulatedChipSeqPeaks.path()
        #truth=f"data/{reference_genome}/{chip_seq.until('n_peaks')}/simulated_peaks.bed"
    output:
        #peaks=f"data/{reference_genome}/{chip_seq}/peak_calling_accuracy.txt",
        #peaks = CalledPeaks.path() + "/peak_calling_accuracy.txt"
        result = PeakCallingAccuracy.path()
    script: "../scripts/get_chip_seq_accuracy.py"
