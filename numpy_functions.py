import numpy as np

numpy_functions = {
    "abs": np.abs,
    "add": np.add,
    "arange": np.arange,
    "arccos": np.arccos,
    "arcsin": np.arcsin,
    "arctan2": np.arctan2,
    "array": np.array,
    "atleast_1d": np.atleast_1d,
    "cos": np.cos,
    "cross": np.cross,
    "divide": np.divide,
    "full": np.full,
    "full_like": np.full_like,
    "interp": np.interp,
    "linalg": np.linalg,
    "multiply": np.multiply,
    "ones": np.ones,
    "ones_like": np.ones_like,
    "searchsorted": np.searchsorted,
    "sign": np.sign,
    "sin": np.sin,
    "size": np.size,
    "subtract": np.subtract,
    "take": np.take,
    "transpose": np.transpose,
    "vstack": np.vstack,
    "zeros": np.zeros,
    "zeros_like": np.zeros_like,
}

assert(tuple(numpy_functions) == tuple(sorted(numpy_functions)))
