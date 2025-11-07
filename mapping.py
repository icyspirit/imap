from ruamel.yaml import YAML
import imas
from loader import add_custom_constructor
from cached_connection import CachedConnection
from config import get_global_config_template, get_local_config_template, get_mapping_template, get_ids_defs, get_env
import postprocessing
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
from urllib.parse import urlparse, parse_qsl


connection = CachedConnection("tcp://172.17.250.21:8005")

yaml = YAML(typ="safe")
add_custom_constructor(yaml, connection)

ids_defs = get_ids_defs()


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


def get_mapping(device, ids_name, shot, tbegin=None, tend=None, dd_version="3.39.0", tree=""):
    env = get_env(device, dd_version, tree)

    connection.openTree("KSTAR", shot)
    connection.get(f"SetTimeContext({fortran_double(tbegin)}, {fortran_double(tend)})")

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

    return yaml.load(mapping_template.render(**ids_defs, **global_configs, **local_configs, **kwargs))


def generate_mapping(device, ids_name, shot, tbegin=None, tend=None, dd_version="3.39.0", tree=""):
    mapping = get_mapping(device, ids_name, shot, tbegin, tend, dd_version, tree)

    factory = imas.IDSFactory(dd_version)
    ids = getattr(factory, ids_name)()
    update_object(ids, mapping)

    return ids
