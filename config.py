
import types
from imas import ids_defs


def get_global_config_template(env):
    return env.get_template("ids/global_configs.yaml")


def get_local_config_template(env, ids_name):
    return env.get_template(f"ids/{ids_name}/local_configs.yaml")


def get_ids_defs():
    return {
        key: value for key, value in vars(ids_defs).items()
        if not key.startswith('_') and
           key.isupper() and
           not callable(value) and
           not isinstance(value, types.ModuleType)
    }

