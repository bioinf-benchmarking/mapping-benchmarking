# Snakemake pipeline for benchmarking read mappers

This is an attempt at an open Snakemake pipeline for benchmarking read mappers. The idea behind this repository is:

* It can be used to benchmark new read mappers (no need to spend time writing your own benchmarking scripts).
* It is difficult to test and present accuracy/performance across all combinations of paramaters (read length, error rate, error profiles, genome diversity, etc.) when presenting benchmark results in a paper. this pipeline, each person can instead run the benchmarks they find interesting when assessing read mappers.  
* With the help of Snakemake and Conda, anyone should be able to clone this repository and easily test read mappers using a configuration of their choice. In addition, we try to run benchmarks on a limited set of cases each night and present them here so that one can get an overview of how mappers perform on the most common cases.

## How to use
There are typically three ways to use this repository:

1) Look at the [latest benchmarking results](...).
2) Run benchmarks locally for a given set of parameters.
3) Contribute by adding read-mapper to this repository or changing configurations/settings

Below are more detailed guides for these:

### Latest benchmarking results
You will find the latest results [here](). If you want these to include other parameters/settings, feel free to edit the configuration files and make a pull-request (see guide further down).

### Run benchmarks locally
1. Install Snakemake and Conda if you havn't already.
2. Clone this repository: `git clone ...`
3. Check that  


### Contribute 

#### 1) Add a new read-mapper
Read-mappers that are not currently listed in config.yaml can be added. Follow these steps below. We assume you already have some experience with writing Snakemake rules, if not check out the Snakemake documentation first.

1. Add the mapper to the `method` list in `config/config.yaml`
2. Add the name of the mapper (this is what is shown in plots) to `method_names`
3. Create a `.smk` file in rules for the mapper (e.g. `bwa.smk`)
4. Implement two rules: `NAME_map_single_end` and `NAME_map_paired_end` (replace NAME with the mapper ID, e.g. `bwa`). The single-end rule should take `{data}/whole_genome_single_end/{config}/reads.fq.gz` as input. The paired-end rule should take this list of two fq-files as input: `["{data}/whole_genome_paired_end/{config}/reads" + n + ".fq.gz" for n in ("1", "2")]` 
5. The rules can also take an index as input if that is needed, you will then need to implement a rule for creating the index (see the `bwa_index` rule for reference).
6. As output, the rules should produce this bam-file: `{data}/whole_genome_single_end/{config}/NAME/{n_threads}/mapped.bam`
7. The rule should use the `n_threads` wildcard as a parameter to specify number of threads. 
8. Either use a Snakemake wrapper or specify a conda-environment for the rule. Nothing should be needed to be installed manually. Note that some Snakemake wrappers do not use the `n_threads` parameter correctly (e.g. the BWA MEM wrapper). If that is the case, you need to implement the run-command yourself and use a conda environment.

#### 2) Add plots/cases
All plots are specified in `config/plots.yaml`. Feel free to edit this file to add plots you believe are useful, and make a pull request. The file contains documentation that should explain how to edit the file, but feel free to reach out with questions if anything is unclear. 

#### 3) General help/contribution
This is a simple first attempt at a pipeline that tries to be flexible and allow benchmarking across what we think are the relevant parameters. However, we want this pipeline to be shaped by what the community believe is important. Feel free to open an issue to discuss things that can be changed or added, e.g. cases or benchmarks that are not currently supported.


## Developer guide

### Creating a plot

This pipeline follows the Snakemake principles, meaning that the user defines what the final result should be, and then the pipeline tries to run the necessary jobs for creating that output. For instance, you can ask for a plot where the x-axis is something, the y-axis is something and the pipeline will try to run what is needed to generate that plot. "Something" needs to be a valid *parameter* or *result_type*. For instance, the x-axis could be `method` (i.e. read mapper) and the y-axis can be `memory_usage` and pipeline will then run all methods and capture the memory usage and present that.

You can add a plot type by adding a configuration under `plot_types` in `config/plots.yaml`. Example:

```yaml
plot_types:
  accuracy_vs_read_length:
    type: "line"
    x: "read_length"
    y: "f1_score"
    color: "method"
    facet_col: "variant_filter"
    facet_row: "error_profile"
```

The above defines a plot type with the name `accuracy_vs_read_length`. We tell the pipeline to make a plot type where the x-axis is `read_length` and the y-axis is `f1_score`. The "color" is `method`, which means that we want one line for each available method (color is the term that **Plotly** uses). We want to repeat this plot for different "variant_filters" along the columns (specified by `facet_col`) and for different error profiles along the rows (specified by `facet_row`). Note that only `x` and `y` are mandatory, the rest can be ommited (in that case only a single plot is created).

The above **only** specifies a **plot_type**, not a plot. For instance, we say that the x-axis should have `method`, but we don't tell it *which* methods to include. Using this plot type, we can specify actual plots with data, under the `plots`-section in the same YAML-file:

```yaml
plots:
  my_plot:
    plot_type: accuracy_vs_read_length
```

We can now ask Snakemake to generate my_plot:

```bash
snakemake reports/presets/my_plot.png
```

Running the above will generate a plot with **default values** for all variables (method, read length etc). This is because no parameters were specified. The default values are specified in `config/config.yaml`. For instance, the default `parameter_group` for read length is `[75, 150, 300]`, specified with the `read_lengths_rough` name. The default individual is `hg002`. If you want to change any of the defaults, you can easily do that. Example: 

```yaml
plots:
  my_plot:
    plot_type: accuracy_vs_read_length
    parameters:
      min_mapq: 30
      individual: hg003
```

Now, if creating `my_plot` it will be run for hg003 and not hg002 which is the default. The minimum mapq for reads considered will be 30 (not 0 which is the default value).


### Add a new read-mapper
Read-mappers that are not currently listed in config.yaml can be added. Follow these steps below. We assume you already have some experience with writing Snakemake rules, if not check out the Snakemake documentation first.

1. Add the mapper to the `method` list in `config/config.yaml`
2. Add the name of the mapper (this is what is shown in plots) to `method_names`
3. Create a `.smk` file in rules for the mapper (e.g. `bwa.smk`)
4. Implement two rules: `NAME_map_single_end` and `NAME_map_paired_end` (replace NAME with the mapper ID, e.g. `bwa`). The single-end rule should take `{data}/whole_genome_single_end/{config}/reads.fq.gz` as input. The paired-end rule should take this list of two fq-files as input: `["{data}/whole_genome_paired_end/{config}/reads" + n + ".fq.gz" for n in ("1", "2")]` 
5. The rules can also take an index as input if that is needed, you will then need to implement a rule for creating the index (see the `bwa_index` rule for reference).
6. As output, the rules should produce this bam-file: `{data}/whole_genome_single_end/{config}/NAME/{n_threads}/mapped.bam`
7. The rule should use the `n_threads` wildcard as a parameter to specify number of threads. 
8. Either use a Snakemake wrapper or specify a conda-environment for the rule. Nothing should be needed to be installed manually. Note that some Snakemake wrappers do not use the `n_threads` parameter correctly (e.g. the BWA MEM wrapper). If that is the case, you need to implement the run-command yourself and use a conda environment.


### Add a new read type
Reads are currently simulated using Mason in the rule `simulate_reads_for_chromosome_and_haplotype` but any tool that is able to take a fasta file (reference) and produce simulated reads and truth positions of those reads will work.

The idea behind read simulation is currently that two haploid references are created and that reads are simulated independently on those. The truth positions are then "translated" to the real reference genome coordinates. This is a bit complicated, and not how read simulation is normally done, but this lets us simulate from a diploid genome and do benchmarking of variant calling and genotyping since we can have heterozygous variants. All the conversion of coordinates is taken care of as long as the read simulation rule outputs the correct files.

The pipeline does currently not support long read sequences, but this could be included by and named `whole_genome_long_reads` for the parameter called `read_type` with a rule similar to `simulate_reads_for_chromosome_and_haplotype`. 