# import pytest
from magicmethyl import magicmethyl


def test_methylate():
    """Generates the correct methylated isomers"""
    analogs = magicmethyl.methylate('ONC(=O)C1CC=CC1')
    assert len(analogs) == 6

    analogs = magicmethyl.methylate('CN(CCN1CCCC1)C(=O)CC1=CC=C(Cl)C(Cl)=C1')
    assert len(analogs) == 14

    analogs = magicmethyl.methylate('CCCC1=CC2=C(OC=CC2=O)C=C1OC(C)=O')
    assert len(analogs) == 9

    analogs = magicmethyl.methylate('CC1=CC(C(=O)N2CCC[C@@H](COC3=CC(C)=C(F)C=C'
                                    '3)C2)=C(C=C1)N1N=CC=N1')
    assert len(analogs) == 21
