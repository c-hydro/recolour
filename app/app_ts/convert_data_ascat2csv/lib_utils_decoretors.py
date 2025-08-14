"""
Library Features:

Name:          lib_utils_decoretors
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20250813'
Version:       '1.0.0'
"""
# ----------------------------------------------------------------------------------------------------------------------
# libraries
from __future__ import annotations
import logging
from typing import Any, Callable, Tuple, Type
import warnings

from functools import wraps
import pandas as pd

from lib_utils_info import logger_name

# set logger obj
logger_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to simplify a list by removing duplicates and returning a single element if only one exists
def simplify_list(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        result = func(*args, **kwargs)

        # Ensure we have a list-like result
        if not isinstance(result, list):
            return result

        # Remove duplicates (preserve order)
        seen = set()
        unique = []
        for item in result:
            if item not in seen:
                seen.add(item)
                unique.append(item)

        # Return single element if only one
        if len(unique) == 1:
            return unique[0]
        return unique
    return wrapper
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to iterate over items in collections
def iterate_items(
    *,
    iter_types: Tuple[Type[Any], ...] = (list, tuple, dict),
    strict_zip: bool = False,
    validate: Callable[[Any], Any] | None = None,
    warn_on_truncate: bool = False,
    dict_mode: str = "values",  # 'values', 'keys', or 'items'
    dict_key_source: str = None,  # 'first', 'second', or None
):
    """
    Decorator factory allowing a function to accept:
      - One positional arg: obj OR iterable[obj]
      - Two positional args: obj/iterable[obj] paired with obj/iterable[obj]
    'obj' can be ANY Python object.

    Args:
      iter_types: Types considered collections to iterate over.
      strict_zip: If True and both args are collections, lengths must match.
      validate: Optional callable(elem) -> elem for per-element validation.
      warn_on_truncate: If True and both args are collections with unequal lengths
                        while strict_zip=False, issue warning or raise.
      dict_mode: How to iterate single dicts ('values', 'keys', or 'items').
      dict_key_source: If both args are dicts, choose key source:
                       'first' → use keys from first dict,
                       'second' → use keys from second dict,
                       None → fall back to normal zip behavior.
    """
    def is_collection(x: Any) -> bool:
        return isinstance(x, iter_types) and not isinstance(x, (str, bytes))

    def to_iterable(x: Any):
        if isinstance(x, dict):
            if dict_mode == "values": return x.values()
            elif dict_mode == "keys": return x.keys()
            elif dict_mode == "items": return x.items()
            else:
                raise ValueError(f"Invalid dict_mode: {dict_mode}")
        return x

    def maybe_validate(x: Any) -> Any:
        return validate(x) if validate else x

    def decorator(func):
        @wraps(func)
        def wrapper(arg1, arg2=None, *args, **kwargs):
            # --- One-argument mode ---
            if arg2 is None:
                if is_collection(arg1):
                    if isinstance(arg1, dict):
                        result = {
                            k: func(maybe_validate(v), *args, **kwargs)
                            for k, v in arg1.items()
                        }
                        return result
                    return [func(maybe_validate(a1), *args, **kwargs)
                            for a1 in to_iterable(arg1)]
                return func(maybe_validate(arg1), *args, **kwargs)

            # --- Two-argument mode ---
            a1_is_col = is_collection(arg1)
            a2_is_col = is_collection(arg2)

            if not a1_is_col and not a2_is_col:
                return func(maybe_validate(arg1), maybe_validate(arg2), *args, **kwargs)

            if a1_is_col and not a2_is_col:
                if isinstance(arg1, dict):
                    return {
                        k: func(maybe_validate(v), maybe_validate(arg2), *args, **kwargs)
                        for k, v in arg1.items()
                    }
                return [func(maybe_validate(a1), maybe_validate(arg2), *args, **kwargs)
                        for a1 in to_iterable(arg1)]

            if not a1_is_col and a2_is_col:
                if isinstance(arg2, dict):
                    return {
                        k: func(maybe_validate(arg1), maybe_validate(v), *args, **kwargs)
                        for k, v in arg2.items()
                    }
                return [func(maybe_validate(arg1), maybe_validate(a2), *args, **kwargs)
                        for a2 in to_iterable(arg2)]

            # --- Both are collections ---
            if isinstance(arg1, dict) and isinstance(arg2, dict) and dict_key_source:
                if dict_key_source == "first":
                    keys = arg1.keys()
                    return {
                        k: func(
                            maybe_validate(arg1[k]),
                            maybe_validate(arg2.get(k)),
                            *args, **kwargs
                        )
                        for k in keys
                    }
                elif dict_key_source == "second":
                    keys = arg2.keys()
                    return {
                        k: func(
                            maybe_validate(arg1.get(k)),
                            maybe_validate(arg2[k]),
                            *args, **kwargs
                        )
                        for k in keys
                    }
                else:
                    raise ValueError(f"Invalid dict_key_source: {dict_key_source}")

            # --- Normal zip behavior ---
            col1, col2 = list(to_iterable(arg1)), list(to_iterable(arg2))
            if strict_zip and len(col1) != len(col2):
                raise ValueError(f"Length mismatch: l1={len(col1)} != l2={len(col2)}")
            elif warn_on_truncate and len(col1) != len(col2):
                warnings.warn(
                    f"Zip would truncate: l1={len(col1)}, l2={len(col2)}. "
                    f"Use strict_zip=True to require equality.",
                    UserWarning
                )

            return [
                func(maybe_validate(a1), maybe_validate(a2), *args, **kwargs)
                for a1, a2 in zip(col1, col2)
            ]

        return wrapper
    return decorator

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to iterate over time steps
def iterate_time_steps(func):
    """
    Decorator to iterate over time steps based on object's
    time_start, time_end, and frequency attributes.
    """

    @wraps(func)
    def wrapper(self, *args, **kwargs):
        if not all(hasattr(self, attr) for attr in ["time_start", "time_end", "frequency"]):
            raise AttributeError("Object must have 'time_start', 'time_end', and 'frequency' attributes.")

        time_range = pd.date_range(start=self.time_start, end=self.time_end, freq=self.frequency)

        results = []
        for time_step in time_range:
            result = func(self, time_step=time_step, *args, **kwargs)
            results.append(result)

        return results

    return wrapper
# ----------------------------------------------------------------------------------------------------------------------
