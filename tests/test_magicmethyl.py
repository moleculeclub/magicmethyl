import pytest
from magicmethyl import methylator

class TestMethylator:
    def test_benzene_to_toluene(self):
        assert methylator.generator('c1ccccc1') == ['Cc1ccccc1']
    
    @pytest.mark.parametrize('test_input,expected', [
        ('c1ccccc1', 1),
        ('Cc1ccccc1', 3),
    ])
    def test_len(self, test_input, expected):
        assert len(methylator.generator(test_input)) == expected
