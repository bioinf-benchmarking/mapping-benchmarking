import plotly.express as px
from hierarchical_results.hierarchical_results import HierarchicalResults, Parameters
hr = HierarchicalResults(config["parameter_types"],config["result_types"], prefix="data/")


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


def get_plot_input_files(wildcards):
    type = wildcards.plot_type
    print("Plot type: ", type)
    x_axis = config["plot_types"][type]["x_axis"]
    y_axis = config["plot_types"][type]["y_axis"]
    print("x axis / yaxis", x_axis, y_axis)

    files = [wildcards.path + "/result_type.txt"]
    result_type = None

    for axis in [x_axis, y_axis]:
        if axis in config["parameter_types"]:
            print("AXIS: ", axis)
            parameter_group = get_parameter_from_config_path(axis, wildcards.path)
            print("Parameter group: ", parameter_group)
            assert parameter_group in config["parameter_groups"]
            values = config["parameter_groups"][parameter_group]["values"]
            print("Values: ", values)
            files = permute_files(files, config["parameter_groups"][parameter_group]["parameter_type"], values)
        elif axis in config["result_types"]:
            result_type = axis

    assert result_type is not None, "One of the axis for plot must be a result type"
    files = [f.replace("result_type", result_type) for f in files]
    files = ["data/" + f for f in files]
    print(files)
    return files


def get_parameter_combinations_and_result_names(wildcards):
    parameter_combinations = Parameters.from_path(hr.get_names(), wildcards.path)

    type = wildcards.plot_type
    x_axis = config["plot_types"][type]["x_axis"]
    y_axis = config["plot_types"][type]["y_axis"]

    result_names = []

    for axis in [x_axis, y_axis]:
        if axis in config["parameter_types"]:
            parameter_group = get_parameter_from_config_path(axis, wildcards.path)
            assert parameter_group in config["parameter_groups"]
            values = config["parameter_groups"][parameter_group]["values"]
            parameter_name = config["parameter_groups"][parameter_group]["parameter_type"]
            parameter_combinations.set(parameter_name, values)
        elif axis in config["result_types"]:
            result_names.append(axis)

    return parameter_combinations, result_names


def get_plot_input_files2(wildcards):
    parameter_combinations, result_names = get_parameter_combinations_and_result_names(wildcards)
    files = hr.get_result_file_names(parameter_combinations, result_names)
    files = [f for f in files]
    return files


def get_input_file_for_x_and_y_parameters(path, x, y):
    assert y in config["result_types"], "Y must be a result type, now %s" % y
    file_name = path + "/" + y + ".txt"
    if x in config["parameter_types"]:
        pass

def get_result_file(path, x_parameter, x_value, y_parameter):
    path = path.split("/")
    # replace x value
    path[config["parameter_types"].index(x_parameter)] = x_value
    return "data/" + "/".join(path) + "/"  + y_parameter + ".txt"


rule make_plot:
    input: get_plot_input_files2
    output:
        #"reports/plots/{plot_type}/{genome_build}/{individual}/{dataset_size}/{read_type}/" + \
        #"{error_profile}/{read_length}/{n_reads}/{method}/{n_threads}/{min_mapq}/{variant_filter}/plot.png"
        plot="reports/plots/{plot_type, \w+}/{path}/plot.png",
        plot_html="reports/plots/{plot_type, \w+}/{path}/plot.html",
        data="reports/plots/{plot_type, \w+}/{path}/plot.csv",
    run:
        parameter_combinations, result_names = get_parameter_combinations_and_result_names(wildcards)
        df = hr.get_results_dataframe(parameter_combinations, result_names)
        df.to_csv(output.data)
        print(df)

        plot_config = config["plot_types"][wildcards.plot_type]
        fig = px.bar(df, x=plot_config["x_axis"], y=plot_config["y_axis"])
        fig.write_image(output.plot)
        fig.write_html(output.plot_html)
