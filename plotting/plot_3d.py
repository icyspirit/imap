
IDSs = [
    "charge_exchange",
    "coils_non_axisymmetric",
    "ec_launchers",
    "ece",
    "equilibrium",
    "interferometer",
    "magnetics",
    "nbi",
    "pf_active",
    "thomson_scattering",
    "wall",
]


if __name__ == "__main__":
    import os
    import sys
    import imas
    import importlib
    import matplotlib.pyplot as plt
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    try:
        path = sys.argv[1]
    except IndexError:
        print(f"Usage: {sys.argv[0]} $PATH")
        exit()

    plt.figure()
    plt.axes(projection="3d")

    with imas.DBEntry(f"imas:hdf5?path={path}", 'r') as db:
        for ids in IDSs:
            try:
                getattr(importlib.import_module('.'.join(["plotting", ids])), "plot_3d")(db.get(ids))
            except (AttributeError, imas.exception.DataEntryException):
                pass

    try:
        plt.axis("scaled")
    except:
        pass
    plt.show()
