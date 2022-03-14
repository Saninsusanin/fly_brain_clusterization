import pandas as pd


def get_new_id_column(data):
    tmp = data.h.value_counts()
    left_id = 0
    right_id = int(tmp['l'])
    middle_id = right_id + int(tmp['r'])
    new_id_column = []

    for index, row in data.iterrows():
        if row['h'] == 'l':
            new_id_column.append(left_id)
            left_id += 1
        elif row['h'] == 'r':
            new_id_column.append(right_id)
            right_id += 1
        else:
            new_id_column.append(middle_id)
            middle_id += 1

    return new_id_column


def data_preparation(path_to_file):
    data = pd.read_csv(path_to_file)
    data['h'] = data['h'].fillna('m')
    data_ides = data.n.value_counts()

    '''split data in different files'''
    for data_id in data_ides.keys():
        new_data = data[data.n == data_id]
        new_data['id'] = get_new_id_column(new_data)

        '''get new file names'''
        new_path_to_file = path_to_file.split('/')
        tmp = new_path_to_file[-1].split('.')
        tmp[0] += f'-{data_id}'
        tmp = '.'.join(tmp)
        new_path_to_file[-1] = tmp
        new_path_to_file = '/'.join(new_path_to_file)

        '''write data in new file'''
        new_data.to_csv(new_path_to_file, index=False)


