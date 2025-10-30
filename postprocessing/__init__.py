import importlib

__all__ = [
    "core_profiles",
    "equilibrium",
]


def __getattr__(name):
    if name in __all__:
        return importlib.import_module(f"{__name__}.{name}")
    raise AttributeError(f"Module '{__name__}' has no attribute {name}")
