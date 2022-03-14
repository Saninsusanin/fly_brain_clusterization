from __future__ import print_function
from mpl_toolkits.mplot3d import Axes3D
from collections import defaultdict

import os
import numpy as np
import imageio
import pandas as pd
import matplotlib.pyplot as plt

path_to_data = '../data/'
names = ['zonetab-mel-animal-1.csv',
         'zonetab-mel-animal-2.csv',
         'zonetab-mel-animal-3.csv',
         'zonetab-sim-animal-2.csv',
         'zonetab-sim-animal-3.csv',
         'zonetab-sim-animal-4.csv']


def draw(data, labels, title, name, elev=90):
    plt.close('all')
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(elev=elev, azim=90)
    ax.scatter(data[:, 0],
               data[:, 1],
               data[:, 2],
               c=labels)
    ax.set_title(title)

#    # Hide grid lines
#    ax.grid(False)
#
#    # Hide axes ticks
#    ax.set_xticks([])
#    ax.set_yticks([])
#    ax.set_zticks([])
    plt.savefig(name)
    plt.show()


def get_labels_for_graph_partition(path_to_partition):
    cluster_id = 1
    data = np.array([])

    with open(path_to_partition, 'r') as partition:
        partition.readline()

        for cluster in partition.readlines():
            vertices = np.array(list(map(int, cluster.split())))
            labels = np.array([cluster_id] * len(vertices))
            cluster_id += 1
            tmp = np.concatenate(([vertices], [labels]), axis=0)

            if 0 in data.shape:
                data = tmp
            else:
                data = np.append(data, tmp, axis=1)

    return pd.DataFrame(np.transpose(data), columns=['vertex_ID', 'cluster_ID'])


def reorder_labels(data, vertex_IDs):
    labels = [0] * len(vertex_IDs)
    data_dict = pd.Series(data.cluster_ID.values, index=data.vertex_ID).to_dict()

    for index in range(len(vertex_IDs)):
        labels[index] = data_dict[vertex_IDs[index] - 1]

    return np.array(labels)


def write_table(name, target, labels):
    target_keys = np.unique(target)
    target_values = [i for i in range(len(target_keys))]
    target_map = {target_keys[i]: target_values[i] for i in range(len(target_keys))}

    labels_keys = np.unique(labels)
    labels_values = [i for i in range(len(labels_keys))]
    labels_map = {labels_keys[i]: labels_values[i] for i in range(len(labels_keys))}

    table = np.array([[0.] * len(target_keys)] * len(labels_keys))

    number_of_points = [0] * len(target_keys)

    for index in range(len(target)):
        target_index = target_map[target[index]]
        label_index = labels_map[labels[index]]
        number_of_points[target_index] += 1
        table[label_index][target_index] += 1

    for target_key in target_keys:
        target_index = target_map[target_key]

        for index in range(len(labels_keys)):
            table[index, target_index] /= number_of_points[target_index]

    table = np.transpose(table)
    data = pd.DataFrame(table, index=target_keys, columns=labels_keys)
    data.to_csv(name)


def csv_to_tex(input, output, label, caption):
    os.system(f'python3 tably.py {input} -o {output} -l {label} -c {caption}')


def test_genetic():

    for name in names:
        brain = pd.read_csv(path_to_data + name)

        for hemi_ID in ['l', 'r']:
            dir_name = name.split('.')[0] + '-' + hemi_ID
            os.system(f'mkdir -p {dir_name}')
            hemisphere = brain[brain.h == hemi_ID]
            hemisphere['id'] = [i for i in range(len(hemisphere['id']))]
            ides = hemisphere.id.values + 1
            target = hemisphere.zone.values
            hemisphere_points = hemisphere.iloc[:, [1, 2, 3]].to_numpy()
            hemisphere.to_csv(f'{dir_name}/data.csv', index=False)

            os.system(f'./graph_builder -i {dir_name}/data.csv -o {dir_name}/graph.in '
                      '-c graph_builder_config --distance minkowski --p 1')
            print('graph was built')

            os.system(f'./algorithm --input {dir_name}/graph.in --output {dir_name}/genetic_partition.out '
                      f'--multi --algo genetic --seed 42 --gen louvain --mut separating extracting merging '
                      f'--engine girvan-newman --config genetic_config')
            print('partition was generated')

            with open(f'{dir_name}/real_partition.out', 'w') as real_partition:
                real_partition.write('1\n')

                for label in np.unique(hemisphere.zone):
                    real_partition.write(' '.join(list(map(str,
                                                           hemisphere[hemisphere.zone == label].id.values))) + '\n')
            print('real partition was saved')
            os.system(f'./metric -r {dir_name}/real_partition.out -f {dir_name}/genetic_partition.out '
                      f'-o {dir_name}/genetic_metric_value')
            draw(hemisphere_points, target, 'real', dir_name + os.sep + 'real')
            genetic_labels = reorder_labels(get_labels_for_graph_partition(f'{dir_name}/genetic_partition.out'), ides)
            draw(hemisphere_points, genetic_labels,
                 'genetic', dir_name + os.sep + 'genetic')
            write_table(f'{dir_name}/genetic_vs_real.csv', target, genetic_labels)


def test_graph_builder():

    for name in names:
        brain = pd.read_csv(path_to_data + name)

        for hemi_ID in ['l', 'r']:
            dir_name = name.split('.')[0] + '-' + hemi_ID
            os.system(f'mkdir -p {dir_name}')
            hemisphere = brain[brain.h == hemi_ID]
            hemisphere['id'] = [i for i in range(len(hemisphere['id']))]
            ides = hemisphere.id.values + 1
            target = hemisphere.zone.values
            hemisphere_points = hemisphere.iloc[:, [1, 2, 3]].to_numpy()
            hemisphere.to_csv(f'{dir_name}/data.csv', index=False)
            weight_dict = {'uni': 'uni',
                           'linear': 'linear',
                           'quadratic': 'quadratic',
                           'qubiq': 'qubiq',
                           'exponential': 'exponential'}

            for weight_type, path_to_graph_config in weight_dict.items():
                os.system(f'./graph_builder -i {dir_name}/data.csv -o {dir_name + os.sep + weight_type}_graph.in '
                          f'-c {path_to_graph_config} --distance minkowski --p 1')
                print('graph was built')

                os.system(f'./algorithm --input {dir_name + os.sep + weight_type}_graph.in '
                          f'--output {dir_name + os.sep + weight_type}_genetic_partition.out '
                          f'--multi --algo genetic --seed 42 --gen louvain --mut separating extracting merging '
                          f'--engine girvan-newman --config genetic_config')
                print('partition was generated')

                with open(f'{dir_name}/real_partition.out', 'w') as real_partition:
                    real_partition.write('1\n')

                    for label in np.unique(hemisphere.zone):
                        real_partition.write(' '.join(list(map(str,
                                                               hemisphere[hemisphere.zone == label].id.values))) + '\n')
                print('real partition was saved')

                os.system(f'./metric -r {dir_name}/real_partition.out -f {dir_name + os.sep + weight_type}_genetic_partition.out '
                          f'-o {dir_name + os.sep + weight_type}_genetic_metric_value')


def test_minkowski():

    for name in names:
        brain = pd.read_csv(path_to_data + name)

        for hemi_ID in ['l', 'r']:
            dir_name = name.split('.')[0] + '-' + hemi_ID
            os.system(f'mkdir -p {dir_name}')
            hemisphere = brain[brain.h == hemi_ID]
            hemisphere['id'] = [i for i in range(len(hemisphere['id']))]
            ides = hemisphere.id.values + 1
            target = hemisphere.zone.values
            hemisphere_points = hemisphere.iloc[:, [1, 2, 3]].to_numpy()
            hemisphere.to_csv(f'{dir_name}/data.csv', index=False)
            param_dict = {'manhattan': 1,
                          'euclidian': 2,
                          'max': 'inf'}

            for param_name, param_value in param_dict.items():
                os.system(f'./graph_builder -i {dir_name}/data.csv -o {dir_name + os.sep + param_name}_graph.in '
                          f'-c graph_builder_config --distance minkowski --p {param_value}')
                print('graph was built')

                os.system(f'./algorithm --input {dir_name + os.sep + param_name}_graph.in '
                          f'--output {dir_name + os.sep + param_name}_genetic_partition.out '
                          f'--multi --algo genetic --seed 42 --gen louvain --mut separating extracting merging '
                          f'--engine girvan-newman --config genetic_config')
                print('partition was generated')

                with open(f'{dir_name}/real_partition.out', 'w') as real_partition:
                    real_partition.write('1\n')

                    for label in np.unique(hemisphere.zone):
                        real_partition.write(' '.join(list(map(str,
                                                               hemisphere[hemisphere.zone == label].id.values))) + '\n')
                print('real partition was saved')

                os.system(f'./metric -r {dir_name}/real_partition.out -f {dir_name + os.sep + param_name}_genetic_partition.out '
                          f'-o {dir_name + os.sep + param_name}_genetic_metric_value')


def test_louvain_leiden():

    for name in names:
        brain = pd.read_csv(path_to_data + name)

        for hemi_ID in ['l', 'r']:
            dir_name = name.split('.')[0] + '-' + hemi_ID
            os.system(f'mkdir -p {dir_name}')
            hemisphere = brain[brain.h == hemi_ID]
            hemisphere['id'] = [i for i in range(len(hemisphere['id']))]
            ides = hemisphere.id.values + 1
            target = hemisphere.zone.values
            hemisphere_points = hemisphere.iloc[:, [1, 2, 3]].to_numpy()
            hemisphere.to_csv(f'{dir_name}/data.csv', index=False)
            algorithms = ['louvain', 'leiden']

            for algorithm in algorithms:

                os.system(f'./algorithm --input {dir_name}/graph.in '
                          f'--output {dir_name + os.sep + algorithm}_partition.out '
                          f'--multi --algo {algorithm} --seed 42')
                os.system(f'./metric -r {dir_name}/real_partition.out -f {dir_name + os.sep + algorithm}_partition.out '
                          f'-o {dir_name + os.sep + algorithm}_metric_value')


def get_nonzero_indices(values):
    indices = []

    for i in range(len(values)):
        if values[i] != 0.0:
            indices.append(i)

    return indices


def plot_pies():
    data = pd.read_csv('../pipeline/zonetab-mel-animal-3-l/genetic_vs_real.csv')

    for indices, row in data.iterrows():
        gained_clusters_ides = row.keys().values[1:]
        tmp = np.array(row.array)
        real_cluster_id = tmp[0]
        data_structure = tmp[1:]
        ides = get_nonzero_indices(data_structure)
        ides = sorted(ides, key=lambda _id: data_structure[_id], reverse=True)
        plt.suptitle(str(int(real_cluster_id)), fontweight='bold', fontsize=20)
        plt.pie(data_structure[ides],
                labels=gained_clusters_ides[ides],
                autopct='%1.2f%%')
        _format = 'svg'
        plt.savefig('../pipeline/zonetab-mel-animal-3-l/pie_' + str(int(real_cluster_id)) + '.' + _format,
                    format=_format, transparent=True)
        plt.show()


def key_function(abscissa, x):
    key = 0

    if abscissa < x[0]:
        return 0
    elif abscissa > x[-1]:
        return -1
    else:
        while not (x[key] <= abscissa <= x[key + 1]):
            key += 1
        return key


def get_piecewise_function(x, y):
    functions = [lambda abscissa, i=i: (y[i] - y[i + 1])/(x[i] - x[i + 1]) * abscissa +
                                       (x[i] * y[i + 1] - x[i + 1] * y[i])/(x[i] - x[i + 1]) for i in range(len(x) - 1)]

    def piecewise_function(abscissa, funcs=functions, abscissas=x):
        return functions[key_function(abscissa, abscissas)](abscissa)

    return piecewise_function


def get_statistic():
    amount = defaultdict(int)
    seconds = defaultdict(list)

    for subdir, dirs, files in os.walk('./mel1'):
        for filename in files:
            filepath = subdir + os.sep + filename

            if filepath.endswith(".log"):
                number_of_threads = int(filename.split('-')[4])

                with open(filepath, 'r') as source:
                    data = source.read()

                    if data == '':
                        continue

                    data = data.split('\n')
                    seconds[number_of_threads].append(int(data[-2].split()[-1]))

                amount[number_of_threads] += 1

    for key in seconds.keys():
        seconds[key] = [min(seconds[key]), sum(seconds[key]) / len(seconds[key]), max(seconds[key])]

    x = sorted(seconds.keys())
    y = np.array([seconds[sorted_key] for sorted_key in x])
    _x = np.linspace(x[0], x[-1], num=100)
    plt.plot(x, y[:, 2], 'o-', color='tab:green')
    plt.fill_between(_x, [get_piecewise_function(x, y[:, 2])(abscissa) for abscissa in _x], 0,
                     alpha=0.5, color='tab:green')
    plt.plot(x, y[:, 1], 'o-', color='tab:blue')
    plt.fill_between(_x, [get_piecewise_function(x, y[:, 1])(abscissa) for abscissa in _x], 0,
                     alpha=0.5, color='tab:blue')
    plt.plot(x, y[:, 0], 'o-', color='tab:orange')
    plt.fill_between(_x, [get_piecewise_function(x, y[:, 0])(abscissa) for abscissa in _x], 0,
                     alpha=0.5, color='tab:orange')
    plt.legend(['max', 'mean', 'min'])
    plt.xlabel('количество потоков')
    plt.ylabel('время работы')
    plt.savefig('timedecrease')
    plt.show()

    numerator = seconds[1]

    for i in range(len(y)):
        y[i][1] = numerator[1] / y[i][1]

    plt.plot(x, y[:, 1], 'o-')
    plt.xlabel('количество потоков')
    plt.ylabel('относительное ускорение')
    plt.savefig('speedup')
    plt.show()


def new_plot():
    for name in names:
        brain = pd.read_csv(path_to_data + name)

        for hemi_ID in ['l', 'r']:
            dir_name = name.split('.')[0] + '-' + hemi_ID
            os.system(f'mkdir -p {dir_name}')
            hemisphere = brain[brain.h == hemi_ID]
            hemisphere['id'] = [i for i in range(len(hemisphere['id']))]
            ides = hemisphere.id.values + 1
            target = hemisphere.zone.values
            hemisphere_points = hemisphere.iloc[:, [1, 2, 3]].to_numpy()
            hemisphere.to_csv(f'{dir_name}/data.csv', index=False)
            genetic_labels = reorder_labels(get_labels_for_graph_partition(f'{dir_name}/k_means_partition.out'), ides)

            target_keys = np.unique(target)
            target_values = [i for i in range(len(target_keys))]
            target_map = {target_keys[i]: target_values[i] for i in range(len(target_keys))}

            labels_keys = np.unique(genetic_labels)
            labels_values = [i for i in range(len(labels_keys))]
            labels_map = {labels_keys[i]: labels_values[i] for i in range(len(labels_keys))}

            table = np.array([[0.] * len(target_keys)] * len(labels_keys))

            number_of_points = [0] * len(target_keys)

            for index in range(len(target)):
                target_index = target_map[target[index]]
                label_index = labels_map[genetic_labels[index]]
                number_of_points[target_index] += 1
                table[label_index][target_index] += 1

            table = np.transpose(table)
            data = pd.DataFrame(table, index=target_keys, columns=labels_keys)
            labels_to_target = defaultdict(int)

            for labels_key in labels_keys:
                current_max = target_keys[0]

                for target_key in target_keys[1:]:
                    if data.loc[current_max, labels_key] < data.loc[target_key, labels_key]:
                        current_max = target_key

                labels_to_target[labels_key] = current_max

            bool_vector = [target[i] == labels_to_target[genetic_labels[i]] for i in range(len(target))]
            plt.close('all')
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.view_init(elev=90, azim=90)
            ax.scatter(hemisphere_points[bool_vector][:, 0],
                       hemisphere_points[bool_vector][:, 1],
                       hemisphere_points[bool_vector][:, 2],
                       c=target[bool_vector])
            ax.scatter(hemisphere_points[[not elem for elem in bool_vector]][:, 0],
                       hemisphere_points[[not elem for elem in bool_vector]][:, 1],
                       hemisphere_points[[not elem for elem in bool_vector]][:, 2],
                       c='red')
            ax.set_title('match_not_match')
            plt.savefig(dir_name + os.sep + 'match_not_match_k_means')
            plt.show()

            #
            #draw(hemisphere_points[bool_vector], target[bool_vector], 'match', dir_name + os.sep + 'match_genetic')
            #draw(hemisphere_points[[not elem for elem in bool_vector]], target[[not elem for elem in bool_vector]],
            #     'not_match', dir_name + os.sep + 'not_match_genetic')


def generate_gif():
    for name in names:
        brain = pd.read_csv(path_to_data + name)

        for hemi_ID in ['l']:
            dir_name = name.split('.')[0] + '-' + hemi_ID
            os.system(f'mkdir -p {dir_name}')
            hemisphere = brain[brain.h == hemi_ID]
            hemisphere['id'] = [i for i in range(len(hemisphere['id']))]
            ides = hemisphere.id.values + 1
            target = hemisphere.zone.values
            hemisphere_points = hemisphere.iloc[:, [1, 2, 3]].to_numpy()
            hemisphere.to_csv(f'{dir_name}/data.csv', index=False)
            genetic_labels = reorder_labels(get_labels_for_graph_partition(f'{dir_name}/genetic_partition.out'), ides)

            filenames = []
            plt.grid(False)
            plt.axis('off')

            for azim in range(90, 450, 1):
                filename = f'gif/{azim}'
                filenames.append(filename)
                draw(hemisphere_points, target, 'MEL 1 l', filename, azim)

            with imageio.get_writer('gif/mygif.gif', mode='I') as writer:
                for filename in filenames:
                    image = imageio.imread(filename + '.png')
                    writer.append_data(image)


def main():
    os.system('cp ../graph_clusterisation/graph_builder/application/cmake-build-release/application graph_builder')
    os.system('cp ../graph_clusterisation/clusterizer/cmake-build-release/app/run_app algorithm')
    os.system('cp ../metric/cmake-build-release/metric metric')

    #test_genetic()
    #test_graph_builder()
    #test_minkowski()
    #test_louvain_leiden()
    #csv_to_tex('zonetab-mel-animal-3-l/genetic_vs_real.csv', 'zonetab-mel-animal-3-l/genetic_vs_real.tex', 'test', 'test')
    #plot_pies()
    #get_statistic()
    new_plot()
    #generate_gif()


main()
