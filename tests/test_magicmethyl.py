import pytest
from magicmethyl import magicmethyl

def test_convert(capsys):
    """Correct my_name argument prints"""
    magicmethyl.convert('Jill')
    captured = capsys.readouterr()
    assert 'Jill' in captured.out