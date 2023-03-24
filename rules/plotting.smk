import itertools
from snakehelp.plotting import PlotType, Plot
import plotly.express as px
import tabulate
from hierarchical_results.hierarchical_results import HierarchicalResults, ParameterCombinations
#hr = HierarchicalResults(parameters_wgs, config["result_types"], prefix="data/")
from mapping_benchmarking.config import *



def parse_plot_specification(plot_type):
    specification = {}
    for dimension in config["plotting_dimensions"]:
        if dimension in config["plot_types"][plot_type]:
            specification[dimension] = config["plot_types"][plot_type][dimension]
        else:
            specification[dimension] = None

    specification["labels"] = config["pretty_names"]
    return specification


result_to_classes_mapping = {
    "runtime": Runtime,
    "memory_usage": MemoryUsage,
    "recall": MappingRecall,
    "one_minus_precision": MappingOneMinusPrecision,
    "f1_score": MappingF1Score,
    "variant_calling_precision": None,
    "variant_calling_one_minus_precision": None,
    "variant_calling_recall": None,
    "variant_calling_f1score": None,
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

    print(parameters)
    return parameters


def get_plot(plot_name):
    plot_config = config["plots"][plot_name]
    plot_type_config = config["plot_types"][plot_config["plot_type"]]

    parsed_config = {}
    for name, val in plot_type_config.items():
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


def get_plot_input_files2(wildcards):
    plot_name = wildcards.plot_name
    plot = get_plot(plot_name)
    files = plot.file_names()
    return files


rule make_plot2:
    input: get_plot_input_files2
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

"""
rule make_plot:
    input: get_plot_input_files
    output:
        plot="reports/plots/{plot_type, \w+}/{path}/plot.png",
        plot_html="reports/plots/{plot_type, \w+}/{path}/plot.html",
        data="reports/plots/{plot_type, \w+}/{path}/plot.csv",
        table="reports/plots/{plot_type, \w+}/{path}/table.md"
    run
        def pretty_name(name):
            if name not in config["pretty_names"]:
                return name.capitalize()
            return config["pretty_names"][name]


        hr = get_hierarchical_results(wildcards)
        parameter_combinations, result_names = get_parameter_combinations_and_result_names(wildcards)
        df = hr.get_results_dataframe(parameter_combinations, result_names)
        df.to_csv(output.data, index=False)

        markdown_table = tabulate.tabulate(df, headers=df.columns, tablefmt="github")
        print(markdown_table)
        with open(output.table, "w") as f:
            f.write(markdown_table + "\n")

        if wildcards.plot_type not in config["plot_types"]:
            print("Invalid plot type ", wildcards.plot_type, " not specified in config")

        plot_config = config["plot_types"][wildcards.plot_type]
        plot_type = plot_config["type"]
        title = wildcards.plot_type.capitalize().replace("_", " ")
        if "title" in plot_config:
            title = plot_config["title"]

        specification = parse_plot_specification(wildcards.plot_type)
        markers = False
        if plot_type != "scatter" and "markers" in plot_config and plot_config["markers"]:
            specification["markers"] = True
            assert "labels" in plot_config, "When markers: True, you need to define labels in the plot config"
            specification["text"] = plot_config["labels"]

        assert plot_type in plotting_functions, "Plot type %s not supported"
        func = plotting_functions[plot_type]
        fig = func(df, **specification, template="simple_white", title=title)

        # prettier facet titles, names, etc
        fig.for_each_annotation(lambda a: a.update(text=pretty_name(a.text.split("=")[-1])))
        fig.for_each_trace(lambda t: t.update(name=pretty_name(t.name)))
        if "text" in specification:
            fig.update_traces(textposition="bottom right")

        #fig.update_annotations(font=dict(size=20))
        #fig.update_layout(font=dict(size=20))
        if "layout" in plot_config:
            fig.update_layout(**plot_config["layout"])

        #if plot_type != "bar":
        #    fig.update_traces(marker_size=15)
        fig.show()
        fig.write_image(output.plot)
        fig.write_html(output.plot_html)
"""


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
        files = ['/'.join(file.split("/")[1:]) for file in input]

        for image, table in zip(files[0::2], files[1::2]):
            name = image.split("/")[1].split(".")[0]
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