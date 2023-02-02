

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
    input: get_plot_input_files
    output:
        #"reports/plots/{plot_type}/{genome_build}/{individual}/{dataset_size}/{read_type}/" + \
        #"{error_profile}/{read_length}/{n_reads}/{method}/{n_threads}/{min_mapq}/{variant_filter}/plot.png"
        "reports/plots/{plot_type, \w+}/{path}/plot.png"
    run:
        plot_config = config["plot_types"][wildcards.plot_type]
        if plot_config["type"] == "bar":
            y = plot_config["y_axis"]
            x = plot_config["x_axis"]
            assert y in config["result_types"], "y value %s" % y

            # x can also potentially be a result_type
            assert x in config["parameter_types"], "x values %s" % x

            x_values = get_parameter_from_config_path(x, wildcards.path)
            assert x_values in config["parameter_groups"], "x values %s not in config"

            x_values = config["parameter_groups"][x_values]["values"]

            y_files = [get_result_file(wildcards.path, x, x_value, y) for x_value in x_values]
            y_values = [float(open(y_file).read().strip()) for y_file in y_files]

            import plotly.express as px
            fig = px.bar(x=x_values, y=y_values)
            fig.write_image(output[0])

        else:
            assert False, "not implemented"

