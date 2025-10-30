import os
import types
from imas import ids_defs
import numpy as np
import datetime
import zoneinfo
import jinja2
from numpy_functions import numpy_functions


def get_global_config_template(env):
    return env.get_template("global_configs.yaml")


def get_local_config_template(env, ids_name):
    return env.get_template('/'.join([ids_name, "local_configs.yaml"]))


def get_mapping_template(env, ids_name):
    return env.get_template('/'.join([ids_name, "mappings.yaml"]))


def get_ids_defs():
    return {
        key: value for key, value in vars(ids_defs).items()
        if not key.startswith('_') and
           key.isupper() and
           not callable(value) and
           not isinstance(value, types.ModuleType)
    }


def get_env(device, dd_version="3.39.0", tree=""):
    env = jinja2.Environment(loader=jinja2.FileSystemLoader('/'.join([os.path.dirname(__file__), "mappings", device])))
    env.globals["zip"] = zip
    env.globals["np"] = numpy_functions

    env.globals["DD_VERSION"] = dd_version
    env.globals["TREE"] = tree.upper()
    env.globals["KSTNOW"] = datetime.datetime.now(zoneinfo.ZoneInfo("Asia/Seoul"))
    env.globals["PROVIDER"] = "Korea Institute of Fusion Energy"
    env.globals["FLOAT_TYPE"] = "float64"
    env.globals["PI"] = np.pi

    return env
