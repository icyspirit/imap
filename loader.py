
from ruamel.yaml import YAML
import numpy as np
import functools
try:
    import MDSplus
    exceptions = MDSplus.mdsExceptions
except (ImportError, ModuleNotFoundError):
    import mdsthin as MDSplus
    exceptions = MDSplus.exceptions
from numpy_functions import numpy_functions
import logging
logger = logging.getLogger(__name__)


def debug(f):
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        logger.debug(f"{f.__name__} called with {args=}, {kwargs=}")
        ret = f(*args, **kwargs)
        logger.debug(f"{f.__name__} returns {ret}")
        return ret
    return wrapper


def add_unique_constructor(yaml):
    @debug
    def unique(data):
        assert(np.all(data[0] == data))
        return data[0]

    yaml.constructor.add_constructor(
        "!unique",
        lambda loader, node: unique(*loader.construct_sequence(node, deep=True)),
    )


def add_numpy_constructor(yaml):
    def _f(f, loader, node):
        return f(*loader.construct_sequence(node, deep=True))

    for name, f in numpy_functions.items():
        yaml.constructor.add_constructor(
            f"!np.{name}",
            functools.partial(_f, debug(f)),
        )


def add_mdsplus_constructor(yaml, connection):
    @debug
    def mdscheck(*exp):
        try:
            for exp_ in exp:
                connection.get(exp_)
            return True
        except (exceptions.TreeNNF, exceptions.TreeNODATA) as e:
            return False

    yaml.constructor.add_constructor(
        "!mdscheck",
        lambda loader, node: mdscheck(*loader.construct_sequence(node, deep=True)),
    )

    if MDSplus.__name__ == "mdsthin":
        # When the data is [], mdsthin returns None
        # This is a temporary fix for mdsthin in this case
        @debug
        def mdsget(exp, *slices):
            data = connection.get(exp).data()
            if data is None:
                return []
            if slices:
                return data[tuple(slices)]
            return data
    else:
        @debug
        def mdsget(exp, *slices):
            data = connection.get(exp).data()
            if slices:
                return data[tuple(slices)]
            return data

    yaml.constructor.add_constructor(
        "!mdsget",
        lambda loader, node: mdsget(*loader.construct_sequence(node, deep=True)),
    )


def add_custom_constructor(yaml, connection):
    add_unique_constructor(yaml)
    add_numpy_constructor(yaml)
    add_mdsplus_constructor(yaml, connection)

