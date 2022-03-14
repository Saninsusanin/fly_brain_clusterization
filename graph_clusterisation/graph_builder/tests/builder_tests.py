import pytest
import pandas as pd
from pathlib import Path


class Case:
    def __init__(self, csv_filename: str, graph_filename: str, name: str):
        self.path_to_csv = Path('../../data') / csv_filename
        self.path_to_graph = Path('../classifier/data') / graph_filename
        self.name = name

    def __str__(self) -> str:
        return 'preparation_test_{}'.format(self.name)


CASES = [
    Case('zonetab-mel-animal-1.csv', 'brain_mel_1.in', 'test_1'),
    Case('zonetab-mel-animal-2.csv', 'brain_mel_2.in', 'test_2'),
    Case('zonetab-mel-animal-3.csv', 'brain_mel_3.in', 'test_3'),
    Case('zonetab-sim-animal-2.csv', 'brain_sim_2.in', 'test_4'),
    Case('zonetab-sim-animal-3.csv', 'brain_sim_3.in', 'test_5'),
    Case('zonetab-sim-animal-4.csv', 'brain_sim_4.in', 'test_6')
]


def check_file(path_to_csv, path_to_graph) -> bool:
    data = pd.read_csv(path_to_csv)

    ids = [-1 for _ in range(len(data))]
    for row in data.itertuples(index=False, name='Pandas'):
        if row.h == 'l':
            ids[row.id] = 1
        elif row.h == 'r':
            ids[row.id] = 2
        else:
            ids[row.id] = 3

    l_start, l_end = 0, 0
    r_start, r_end = 0, 0
    m_start, m_end = 0, len(data) - 1
    for ind in range(len(ids) - 1):
        if ids[ind] == 1 and ids[ind + 1] == 2:
            l_end, r_start = ind, ind + 1
        if ids[ind] == 2 and ids[ind + 1] == 3:
            r_end, m_start = ind, ind + 1

    def get_hemisphere(index):
        if l_start <= index <= l_end:
            return 'l'
        if r_start <= index <= r_end:
            return 'r'
        if m_start <= index <= m_end:
            return 'm'
        return None

    with open(path_to_graph) as graph_file:
        for line in graph_file:
            first, second, _ = list(map(int, line.split()))
            first_h = get_hemisphere(first)
            second_h = get_hemisphere(second)

            if not first_h or not second_h or not first_h == second_h:
                return False
    return True


@pytest.mark.parametrize('case', CASES, ids=str)
def test_builder(case) -> None:
    assert check_file(case.path_to_csv, case.path_to_graph)
