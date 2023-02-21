import itertools

import plotly.express as px
import tabulate
from hierarchical_results.hierarchical_results import HierarchicalResults, ParameterCombinations
hr = HierarchicalResults(config["parameter_types"],config["result_types"], prefix="data/")

plotting_functions = {
    "bar": px.bar,
    "line": px.line,
    "scatter": px.scatter
}


def get_parameter_from_config_path(parameter, path):
    print(parameter, path)
    assert parameter in config["parameter_types"]
    return path.split("/")[config["parameter_types"].index(parameter)]


def permute_files(files, parameter, values):
    print("Permuting files %s with parameter %s and values %s" % (files, parameter, values))
    new_files = []
    for file in files:
        splitted = file.split("/")
        for value in values:
            splitted[config["parameter_types"].index(parameter)] = str(value)
            new_files.append("/".join(splitted))
    return new_files


def get_parameter_combinations_and_result_names(wildcards):
    parameter_combinations = ParameterCombinations.from_path(hr.get_names(), wildcards.path)

    type = wildcards.plot_type

    dimensions = []
    for possible_dimension in config["plotting_dimensions"]:
        if possible_dimension in config["plot_types"][type]:
            dimensions.append(config["plot_types"][type][possible_dimension])

    result_names = []

    for dimension in dimensions:
        if dimension in config["parameter_types"]:
            parameter_group = get_parameter_from_config_path(dimension, wildcards.path)
            assert parameter_group in config["parameter_sets"], "Parameter group %s invalid" % parameter_group
            values = config["parameter_sets"][parameter_group]["values"]
            parameter_name = config["parameter_sets"][parameter_group]["parameter_type"]
            parameter_combinations.set(parameter_name, values)
        elif dimension in config["result_types"]:
            result_names.append(dimension)
        else:
            raise Exception("Could not link %s to a parameter type or result type" % dimension)

    return parameter_combinations, result_names


def get_plot_input_files(wildcards):
    parameter_combinations, result_names = get_parameter_combinations_and_result_names(wildcards)
    files = hr.get_result_file_names(parameter_combinations, result_names)
    files = [f for f in files]
    return files


def parse_plot_specification(plot_type):
    specification = {}
    for dimension in config["plotting_dimensions"]:
        if dimension in config["plot_types"][plot_type]:
            specification[dimension] = config["plot_types"][plot_type][dimension]
        else:
            specification[dimension] = None

    specification["labels"] = config["pretty_names"]
    return specification


rule make_plot:
    input: get_plot_input_files
    output:
        plot="reports/plots/{plot_type, \w+}/{path}/plot.png",
        plot_html="reports/plots/{plot_type, \w+}/{path}/plot.html",
        data="reports/plots/{plot_type, \w+}/{path}/plot.csv",
        table="reports/plots/{plot_type, \w+}/{path}/table.md"
    run:
        def pretty_name(name):
            if name not in config["pretty_names"]:
                return name.capitalize()
            return config["pretty_names"][name]


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

        assert plot_type in plotting_functions, "Plot type %s not supported"
        func = plotting_functions[plot_type]
        fig = func(df, **specification, template="simple_white", title=title)

        # prettier facet titles, names, etc
        fig.for_each_annotation(lambda a: a.update(text=pretty_name(a.text.split("=")[-1])))
        fig.for_each_trace(lambda t: t.update(name=pretty_name(t.name)))

        #fig.update_annotations(font=dict(size=20))
        #fig.update_layout(font=dict(size=20))
        if "layout" in plot_config:
            fig.update_layout(**plot_config["layout"])

        #if plot_type != "bar":
        #    fig.update_traces(marker_size=15)
        fig.show()
        fig.write_image(output.plot)
        fig.write_html(output.plot_html)


def get_plot_name(wildcards):
    name = wildcards.name

    if name in config["plots"]:
        #assert name in config["plots"], "Plot name %s not defined in plots.yaml" % name
        plot_config = config["plots"][name]
        assert plot_config["plot_type"] in config["plot_types"], "Plot specifies a plot type %s that is not in config.plot_types" % plot_config["plot_type"]
        plot_type_config = config["plot_types"][plot_config["plot_type"]]
    else:
        # allow direct use of a plot type without having it defined as a plot
        plot_config = config["plots"]["generic"]
        print(plot_config)
        plot_config["plot_type"] = name
        plot_type_config = config["plot_types"][name]

    # Parameters that can vary for this plot:
    variables = [plot_type_config[dimension] for dimension in config["plotting_dimensions"] if dimension in plot_type_config]
    print(variables)

    plot_path = []
    for parameter in config["parameter_types"]:
        if "parameters" in plot_config and parameter in plot_config["parameters"]:
            parameter = str(plot_config["parameters"][parameter])
        else:
            # not specified, use default value
            if parameter in variables:
                # use parameter_group default value
                parameter = config["default_parameter_sets"][parameter]
            else:
                # use default parameter
                parameter = config["default_parameter_values"][parameter]
        plot_path.append(parameter)

    base_name = "reports/plots/" + plot_config["plot_type"] + "/" + "/".join(plot_path) + "/"
    endings = ["plot.png", "table.md"]
    return [base_name + ending for ending in endings]


# Wrapper around the make_plot rule that uses default parameters
# so that less stuff needs to be specified
rule plot_from_name:
    input: get_plot_name
    output:
        png="reports/presets/{name}.png",
        table="reports/presets/{name}.md",
    shell:
        "cp {input[0]} {output[0]} && "
        "cp {input[1]} {output[1]}"


def get_report_input(wildcards):
    plots = config["test_plots"] if wildcards.type == "test" else config["nightly_plots"]
    return list(itertools.chain.from_iterable(zip(["reports/presets/" + name + ".png" for name in plots],
               ["reports/presets/" + name + ".md" for name in plots])))


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