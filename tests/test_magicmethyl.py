# import pytest
from magicmethyl import magicmethyl


def test_methylate():
    """Generates the correct methylated isomers"""
    analogs = magicmethyl.methylate('ONC(=O)C1CC=CC1')
    assert len(analogs) == 6
    assert 'C[C@@]1(C(=O)NO)CC=CC1' in analogs
    assert 'CONC(=O)C1CC=CC1' in analogs
    assert 'C[C@H]1C=CCC1C(=O)NO' in analogs
    assert 'CN(O)C(=O)C1CC=CC1' in analogs
    assert 'CC1=CCC(C(=O)NO)C1' in analogs
    assert 'C[C@@H]1C=CCC1C(=O)NO' in analogs

    analogs = magicmethyl.methylate('CN(CCN1CCCC1)C(=O)CC1=CC=C(Cl)C(Cl)=C1')
    assert len(analogs) == 14

    # this test fails due to inability to recognize symmetry in methyl groups
    analogs = magicmethyl.methylate('CCCC1=CC2=C(OC=CC2=O)C=C1OC(C)=O')
    assert len(analogs) == 9
