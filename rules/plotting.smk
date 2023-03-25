from mapping_benchmarking.config import *
import itertools
from snakehelp.plotting import PlotType
import tabulate


result_to_classes_mapping = {
    "runtime": Runtime,
    "memory_usage": MemoryUsage,
    "mapping_recall": MappingRecall,
    "mapping_one_minus_precision": MappingOneMinusPrecision,
    "mapping_f1_score": MappingF1Score,
    "variant_calling_recall": VariantCallingRecall,
    "variant_calling_one_minus_precision": VariantCallingOneMinusPrecision,
    "variant_calling_f1_score": VariantCallingF1Score,
    "peak_calling_accuracy": PeakCallingAccuracy
}


def get_plot_type_parameters(plot_type, plot_type_object):
    """
    Returns a dict with parameters for the plot type.
    Uses default range values if a parameter is not specified in the config.
    """
    plot_type_config = config["plot_types"][plot_type]

    # parse those specified in yaml config
    if "parameters" in plot_type_config:
        parameters = plot_type_config["parameters"]
    else:
        parameters = {}

    # set defaults for those not specified
    for required_parameter in plot_type_object.parameter_types():
        if required_parameter not in parameters:
            parameters[required_parameter] = config["default_parameter_sets"][required_parameter]

    return parameters


def get_plot(plot_name):
    plot_config = config["plots"][plot_name]
    print("Plot config", str(plot_config))
    plot_type_config = config["plot_types"][plot_config["plot_type"]]

    parsed_config = {}
    for name, val in plot_type_config.items():
        print(name, val)
        if name != "parameters" and name != "layout":
            if val in result_to_classes_mapping:
                val = result_to_classes_mapping[val]
        parsed_config[name] = val

    # replace literal values with classes

    print("Parsed config")
    print(parsed_config)
    plot_type = PlotType.from_yaml_dict(parsed_config)

    out_base_name = f"plots/{plot_name}"

    parameters = get_plot_type_parameters(plot_config["plot_type"],plot_type)
    plot = plot_type.plot(out_base_name,**parameters)
    return plot


def get_plot_input_files(wildcards):
    plot_name = wildcards.plot_name
    plot = get_plot(plot_name)
    files = plot.file_names()
    return files


rule make_plot:
    input: get_plot_input_files
    output:
        plot="plots/{plot_name}.png",
        csv="plots/{plot_name}.csv",
        md="plots/{plot_name}.md",
    run:
        plot = get_plot(wildcards.plot_name)
        df = plot._parameter_combinations.get_results_dataframe()
        plot.plot()

        df.to_csv(output.csv, index=False)

        markdown_table = tabulate.tabulate(df,headers=df.columns,tablefmt="github")
        print(markdown_table)
        with open(output.md, "w") as f:
            f.write(markdown_table + "\n")


def get_report_input(wildcards):
    plots = config["test_plots"] if wildcards.type == "test" else config["nightly_plots"]
    return list(itertools.chain.from_iterable(zip(["plots/" + name + ".png" for name in plots],
               ["plots/" + name + ".md" for name in plots])))


rule report:
    input: get_report_input
    output:
        "reports/{type, test|main}.md"

    run:
        import time
        timestamp = time.strftime("%Y-%m-%d %H:%M")
        out = "# Report, created  " + timestamp + "\n\n"
        out += "Note: This pipeline is under development, and the results may be wrong or inaccurate.\n\n"

        # remove reports/
        files = ['../' + file for file in input]

        for image, table in zip(files[0::2], files[1::2]):
            print(image, table)
            name = image.split("/")[-1].split(".")[0]
            plot_config = config["plots"][name]
            description = ""
            if "description" in plot_config:
                description = plot_config["description"] + "\n\n"
            out += "![](" + image + ") \n\n " + description + " [Link to plot data](" + table + ") \n\n"

        with open(output[0],"w") as f:
            f.write(out)



rule dummy_report:
    output: "reports/dummy.md"
    shell:
        "echo 'test' > {output}"