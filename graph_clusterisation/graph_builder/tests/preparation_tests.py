import pytest
import pandas as pd
from pathlib import Path


class Case:
    def __init__(self, filename: str, name: str):
        self.path_to_file = Path('../../data') / filename
        self.name = name

    def __str__(self) -> str:
        return 'preparation_test_{}'.format(self.name)


CASES = [
    Case('zonetab-mel-animal-1.csv', 'test_1'),
    Case('zonetab-mel-animal-2.csv', 'test_2'),
    Case('zonetab-mel-animal-3.csv', 'test_3'),
    Case('zonetab-sim-animal-2.csv', 'test_4'),
    Case('zonetab-sim-animal-3.csv', 'test_5'),
    Case('zonetab-sim-animal-4.csv', 'test_6')
]


def check_file(path_to_file) -> bool:
    data = pd.read_csv(path_to_file)

    ids = [-1 for _ in range(len(data))]
    for row in data.itertuples(index=False, name='Pandas'):
        if row.id >= len(data):
            return False
        if row.h == 'l':
            ids[row.id] = 1
        elif row.h == 'r':
            ids[row.id] = 2
        else:
            ids[row.id] = 3

    for ind in range(len(ids) - 1):
        if ids[ind] == -1:
            return False
        if ids[ind + 1] - ids[ind] < 0:
            return False
    return True


@pytest.mark.parametrize('case', CASES, ids=str)
def test_preparation(case) -> None:
    assert check_file(case.path_to_file)
