from collections import defaultdict

import plotly.express as px
from hierarchical_results.hierarchical_results import HierarchicalResults, Parameters
hr = HierarchicalResults(config["parameter_types"],config["result_types"], prefix="data/")

plotting_functions = {
    "bar": px.bar,
    "line": px.line
}


# reports/plots/f1_score_bar/hg38/hg002/small/whole_genome_single_end/low_error/150/1000/all_methods/4/5/variants/plot.png


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
    parameter_combinations = Parameters.from_path(hr.get_names(), wildcards.path)

    type = wildcards.plot_type

    dimensions = []
    for possible_dimension in config["plotting_dimensions"]:
        if possible_dimension in config["plot_types"][type]:
            dimensions.append(config["plot_types"][type][possible_dimension])

    result_names = []

    for dimension in dimensions:
        if dimension in config["parameter_types"]:
            parameter_group = get_parameter_from_config_path(dimension, wildcards.path)
            assert parameter_group in config["parameter_groups"], "Parameter group %s invalid" % parameter_group
            values = config["parameter_groups"][parameter_group]["values"]
            parameter_name = config["parameter_groups"][parameter_group]["parameter_type"]
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
    return specification


rule make_plot:
    input: get_plot_input_files
    output:
        plot="reports/plots/{plot_type, \w+}/{path}/plot.png",
        plot_html="reports/plots/{plot_type, \w+}/{path}/plot.html",
        data="reports/plots/{plot_type, \w+}/{path}/plot.csv",
    run:
        parameter_combinations, result_names = get_parameter_combinations_and_result_names(wildcards)
        df = hr.get_results_dataframe(parameter_combinations, result_names)
        df.to_csv(output.data)
        print(df)

        if wildcards.plot_type not in config["plot_types"]:
            print("Invalid plot type ", wildcards.plot_type, " not specified in config")

        plot_config = config["plot_types"][wildcards.plot_type]
        plot_type = plot_config["type"]
        if plot_type == "scatter_and_line":
           specification = parse_plot_specification(wildcards.plot_type)
           fig = px.scatter(df, **specification)
           fig.add_traces(px.line(df, **specification).data)
        else:
            assert plot_type in plotting_functions, "Plot type %s not supported"
            func = plotting_functions[plot_type]
            fig = func(df, **parse_plot_specification(wildcards.plot_type))

        fig.update_layout(font=dict(size=20))
        fig.show()
        fig.write_image(output.plot)
        fig.write_html(output.plot_html)



def get_plot_name(wildcards):
    name = wildcards.name
    assert name in config["plots"], "Plot name %s not defined in plots.yaml" % name
    plot_config = config["plots"][name]
    assert plot_config["plot_type"] in config["plot_types"], "Plot specifies a plot type %s that is not in config.plot_types" % plot_config["plot_type"]
    plot_type_config = config["plot_types"][plot_config["plot_type"]]

    # Parameters that can vary for this plot:
    variables = [plot_type_config[dimension] for dimension in config["plotting_dimensions"] if dimension in plot_type_config]
    print(variables)

    plot_path = []
    for parameter in config["parameter_types"]:
        print(parameter)
        if parameter in plot_config["parameters"]:
            parameter = plot_config["parameters"][parameter]
        else:
            # not specified, use default value
            if parameter in variables:
                # use parameter_group default value
                parameter = config["default_parameter_groups"][parameter]
            else:
                # use default parameter
                parameter = config["default_parameter_values"][parameter]
        plot_path.append(parameter)

    file = "reports/plots/" + plot_config["plot_type"] + "/" + "/".join(plot_path) + "/plot.png"
    return file



# Wrapper around the make_plot rule that uses default parameters
# so that less stuff needs to be specified
rule make_plot_from_name:
    input: get_plot_name
    output:
        "reports/presets/{name}.png"
    shell:
        "cp {input} {output}"

