
from ruamel.yaml import YAML
import datetime
import zoneinfo
import jinja2
import imas
from loader import add_custom_constructor
from cached_connection import CachedConnection
from config import get_global_config_template, get_local_config_template, get_ids_defs
import postprocessing
import numpy as np
import functools
from numpy_functions import numpy_functions
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
from urllib.parse import urlparse, parse_qsl


def fortran_double(value):
    if isinstance(value, float):
        return f"{value:e}".replace('e', 'D')
    return '*'


def parse_url(ids_url):
    parsed = urlparse(ids_url)
    ids_name = parsed.path.strip("/")

    query_dict = dict(parse_qsl(parsed.query))
    postprocess = query_dict.get("postprocess", "").lower() == "true"
    tree = query_dict.get("tree", "")

    return ids_name, postprocess, tree


def exclude_private_keys(d):
    if d is None:
        return {}
    return dict(filter(lambda item: not item[0].startswith('_'), d.items()))


def update_object(obj, data):
    for key, value in data.items():
        if key.startswith('_'):
            continue

        if not hasattr(obj, key):
            raise KeyError(f"{obj} does not have attribute {key}")

        if isinstance(value, dict):
            update_object(getattr(obj, key), value)
        elif isinstance(value, list) and all(map(lambda v: isinstance(v, (dict, list)), value)):
            getattr(obj, key).resize(len(value))
            for o, d in zip(getattr(obj, key), value):
                update_object(o, d)
        else:
            try:
                setattr(obj, key, value)
            except TypeError as e:
                logger.error(f"Setting {key=} with {value=} failed in {obj}")
                raise e


def get_mapping_template(env, ids_name):
    return env.get_template(f"ids/{ids_name}/mappings.yaml")


def _generate_mapping(ids_name, shot, tbegin=None, tend=None, dd_version="3.39.0", tree=""):
    env = jinja2.Environment(loader=jinja2.FileSystemLoader('.'))
    env.globals["zip"] = zip
    env.globals["np"] = numpy_functions

    env.globals["DD_VERSION"] = dd_version
    env.globals["TREE"] = tree.upper()
    env.globals["KSTNOW"] = datetime.datetime.now(zoneinfo.ZoneInfo("Asia/Seoul"))
    env.globals["PROVIDER"] = "Korea Institute of Fusion Energy"
    env.globals["FLOAT_TYPE"] = "float64"
    env.globals["PI"] = np.pi

    connection = CachedConnection("tcp://172.17.250.21:8005")
    connection.openTree("KSTAR", shot)
    connection.get(f"SetTimeContext({fortran_double(tbegin)}, {fortran_double(tend)})")

    yaml = YAML(typ="safe")
    add_custom_constructor(yaml, connection)

    ids_defs = get_ids_defs()

    kwargs = {
        "TBEGIN": -abs(ids_defs["EMPTY_FLOAT"]) if tbegin is None else tbegin,
        "TEND": abs(ids_defs["EMPTY_FLOAT"]) if tend is None else tend,
    }

    global_config_template = get_global_config_template(env)
    global_configs = yaml.load(global_config_template.render(**ids_defs, **kwargs))
    global_configs = exclude_private_keys(global_configs)

    local_config_template = get_local_config_template(env, ids_name)
    local_configs = yaml.load(local_config_template.render(**ids_defs, **global_configs, **kwargs))
    local_configs = exclude_private_keys(local_configs)

    mapping_template = get_mapping_template(env, ids_name)
    mappings = yaml.load(mapping_template.render(**ids_defs, **global_configs, **local_configs, **kwargs))

    factory = imas.IDSFactory(dd_version)

    ids = getattr(factory, ids_name)()
    update_object(ids, mappings)

    return ids


def generate_mapping(ids_url, shot, tbegin, tend, dd_version):
    ids_name, postprocess, tree = parse_url(ids_url)
    f = getattr(postprocessing, ids_name).postprocess if postprocess else lambda x: None

    ids = _generate_mapping(ids_name, shot, tbegin, tend, dd_version, tree)
    f(ids)

    return ids


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--shot", type=int, required=True, help="Shot number")
    parser.add_argument("-b", "--tbegin", type=float, default=None, help="Start time")
    parser.add_argument("-e", "--tend", type=float, default=None, help="End time")
    parser.add_argument("-v", "--dd_version", type=str, default="3.39.0", help="Data dictionary version")
    parser.add_argument("-p", "--path", type=str, default=argparse.SUPPRESS, help="IDS data path")
    parser.add_argument("ids", nargs="+", help="List of IDSs")

    args = parser.parse_args()
    if not hasattr(args, "path"):
        args.path = str(args.shot)

    with imas.DBEntry(f"imas:hdf5?path={args.path}", 'w', dd_version=args.dd_version) as db:
        for ids_url in args.ids:
            db.put(generate_mapping(ids_url, args.shot, args.tbegin, args.tend, args.dd_version))
