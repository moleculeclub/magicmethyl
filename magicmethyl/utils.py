from rdkit import Chem


def show_analogs(analogs):
    """Displays molecules from list of SMILES strings.

    Examples
    --------
    analogs = methylate('CC1CCCCC1')
    show_analogs(analogs)

    Parameters
    ----------
    analogs: list[str]
        A list of SMILES strings."""
    print(f'there are {len(analogs)} analogs\n')
    for analog in analogs:
        print(analog)
        display(Chem.MolFromSmiles(analog))  # noqa: F821
