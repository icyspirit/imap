
try:
    import MDSplus

    def f32_to_f64(signal):
        if signal.dtype == "float32":
            return signal.astype("float64")
        return signal

except (ImportError, ModuleNotFoundError):
    import mdsthin as MDSplus

    def f32_to_f64(signal):
        if isinstance(signal, MDSplus.Float32):
            return MDSplus.Float64(signal)
        if isinstance(signal, MDSplus.Float32Array):
            return MDSplus.Float64Array(signal)
        return signal

import functools
import logging
logger = logging.getLogger(__name__)


class CachedConnection(MDSplus.Connection):
    CONVERT_TO_DOUBLE = True

    def openTree(self, tree, shot):
        super().openTree(tree, shot)
        self.get.cache_clear()

    @functools.cache
    def get(self, exp, *args, **kwargs):
        try:
            signal = super().get(exp, *args, **kwargs)
        except Exception as e:
            logger.error(f"Could not get signal {exp}")
            raise e

        if self.CONVERT_TO_DOUBLE:
            return f32_to_f64(signal)

        return signal
