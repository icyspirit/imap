import inspect
import imas
import postprocessings


def generate_required_ids_names(ids_name, dd_version="3.39.0"):
    return filter(
        imas.IDSFactory(dd_version).exists,
        inspect.signature(getattr(postprocessings, ids_name).postprocess).parameters.keys(),
    )


def postprocess(ids_name, ids, *ids_, **kwargs):
    getattr(postprocessings, ids_name).postprocess(ids, *ids_, **kwargs)
    return ids
