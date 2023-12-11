import logging
logging.basicConfig(level=logging.INFO)
from mapping_benchmarking.config import Runtime, MemoryUsage, MappingRecall, MappingOneMinusPrecision, MappingF1Score, VariantCallingRecall, VariantCallingOneMinusPrecision, VariantCallingF1Score, PeakCallingAccuracy
import itertools
from snakehelp.plotting import PlotType
import tabulate
from snakehelp import ResultLike


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


def pretty_name(name):
    logging.info("Checking pretty name %s against %s" % (name, config["pretty_names"]))
    if name not in config["pretty_names"]:
        return name.capitalize()
    return config["pretty_names"][name]


def get_plot_type_parameters(plot_name, plot_type_object):
    """
    Returns a dict with parameters for the plot type.
    Uses default range values if a parameter is not specified in the config.
    """
    plot_type = config["plots"][plot_name]["plot_type"]
    plot_type_config = config["plot_types"][plot_type]
    plot_config = config["plots"][plot_name]

    # parse those specified in yaml config
    if "parameters" in plot_config:
        parameters = plot_config["parameters"]
    else:
        parameters = {}

    # set defaults for those not specified
    for required_parameter in plot_type_object.parameter_types():
        if required_parameter not in parameters:
            parameters[required_parameter] = config["default_parameter_sets"][required_parameter]

    return parameters


def _is_result_class(name):
    try:
        t = eval(name)
        if issubclass(t, ResultLike):
            return True
        return False
    except Exception as e:
        return False

def get_plot(plot_name):
    plot_config = config["plots"][plot_name]
    try:
        plot_type_config = config["plot_types"][plot_config["plot_type"]]
    except KeyError:
        print(f"Plot type {plot_config['plot_type']} not found in config {config['plot_types']}")
        raise

    parsed_config = {}
    result_types = {}
    for name, val in plot_type_config.items():
        if name != "parameters" and name != "layout":
            if val in result_to_classes_mapping:
                val = result_to_classes_mapping[val]
            elif _is_result_class(val):
                val = eval(val)
                result_types[name] = val
        parsed_config[name] = val

    # Limit ResultType union types if specified
    # (WHen a Result is unwrapped and there are union-types, we need to know
    # which type to choose for this plot
    if "type_limits" in plot_config:
        for name, result_type in result_types.items():
            for field_name, new_field in plot_config["type_limits"].items():
                if name == "type":
                    continue
                new_field = eval(new_field)
                #print("New field", new_field)
                #print("  Replacing parsed config %s, field name %s to %s on %s" % (name, field_name, new_field, result_type))
                new_type = result_type.replace_field(field_name,(field_name, new_field, None))
                #print(new_type)
                parsed_config[name] = new_type
        """
        for base_class, limits in plot_config["type_limits"].items():
            cls = eval(base_class)
            cls.clear_union_choices()
            for limit in limits:
                cls.limit_union_choice(limit)
                print("  Limiting result %s to %s" % (cls, limit))
        """

    # replace literal values with classes

    plot_type = PlotType.from_yaml_dict(parsed_config)

    out_base_name = f"plots/{plot_name}"

    parameters = get_plot_type_parameters(plot_name, plot_type)
    plot = plot_type.plot(out_base_name, **parameters)
    return plot, parameters



def get_plot_input_files(wildcards):
    plot_name = wildcards.plot_name
    plot, parameters = get_plot(plot_name)
    files = plot.file_names()
    return files


rule make_plot:
    input: get_plot_input_files
    output:
        plot="plots/{plot_name}.png",
        csv="plots/{plot_name}.csv",
        md="plots/{plot_name}.md",
    run:
        plot, parameters = get_plot(wildcards.plot_name)
        df = plot._parameter_combinations.get_results_dataframe(**parameters)
        plot.plot(pretty_names_func=pretty_name)
        #plot.plot()

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
            title = plot_config["title"] if "title" in plot_config else name
            out += "## " + title + "\n\n"

            description = ""
            if "description" in plot_config:
                description = plot_config["description"] + "\n\n"

            out += description
            out += "![](" + image + ") \n\n [Link to plot data](" + table + ") \n\n"

        with open(output[0],"w") as f:
            f.write(out)



rule dummy_report:
    output: "reports/dummy.md"
    shell:
        "echo 'test' > {output}"