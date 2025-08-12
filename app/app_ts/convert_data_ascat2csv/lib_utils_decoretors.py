# ----------------------------------------------------------------------------------------------------------------------
# libraries
from functools import wraps
import pandas as pd
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
# decorator to iterate over a dictionary of DataFrames
def iterate_dict(func):
    """
    Decorator that applies a function to each DataFrame in a dict,
    preserving the same keys in the returned dict.
    """
    @wraps(func)
    def wrapper(obj, *args, **kwargs):
        if not isinstance(obj, dict):
            raise TypeError("Input must be a dictionary of DataFrames.")
        return {k: func(v, *args, **kwargs) for k, v in obj.items()}
    return wrapper
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# method to iterate over a file list or string
from functools import wraps
from collections.abc import Iterable

def iterate_file_list(func):
    """Decorator to allow a function to handle one or two path arguments,
    each being either a string or an iterable of strings.

    - If both args are strings → call once.
    - If one arg is iterable → iterate over it, pairing with the other arg.
    - If both args are iterable → iterate in parallel (zip).
    """

    @wraps(func)
    def wrapper(file_path1, file_path2=None, *args, **kwargs):

        def is_str(x):
            return isinstance(x, str)

        def is_str_iter(x):
            return isinstance(x, Iterable) and not isinstance(x, str) and all(isinstance(p, str) for p in x)

        # One argument case
        if file_path2 is None:
            if is_str(file_path1):
                return func(file_path1, *args, **kwargs)
            elif is_str_iter(file_path1):
                return [func(p, *args, **kwargs) for p in file_path1]
            else:
                raise TypeError("file_path1 must be a string or iterable of strings.")

        # Two arguments case
        else:
            if is_str(file_path1) and is_str(file_path2):
                return func(file_path1, file_path2, *args, **kwargs)
            elif is_str_iter(file_path1) and is_str(file_path2):
                return [func(p1, file_path2, *args, **kwargs) for p1 in file_path1]
            elif is_str(file_path1) and is_str_iter(file_path2):
                return [func(file_path1, p2, *args, **kwargs) for p2 in file_path2]
            elif is_str_iter(file_path1) and is_str_iter(file_path2):
                return [func(p1, p2, *args, **kwargs) for p1, p2 in zip(file_path1, file_path2)]
            else:
                raise TypeError("file_path1 and file_path2 must be strings or iterables of strings.")

    return wrapper
# ----------------------------------------------------------------------------------------------------------------------

from functools import wraps
from collections.abc import Iterable, Sequence

def iterate_items(*,
                  iter_types=(list, tuple),   # which container types should be iterated
                  strict_zip=False,           # raise if two iterables have different lengths
                  validate=None):             # optional per-element validator/coercer
    """
    Decorator factory that lets a function accept:
      - one arg: obj OR iterable[obj]
      - two args: obj/iterable[obj] paired with obj/iterable[obj]
    'obj' can be ANY Python object (not just strings/paths).

    Args:
      iter_types: containers that should be treated as "iterable-of-objs".
      strict_zip: if True and both args are iterables, lengths must match.
      validate: optional callable(elem) -> elem (e.g., type check or normalization).
    """
    def is_collection(x):
        # Only treat specific container types as collections to avoid
        # accidentally iterating over objects that are Iterable but "atomic".
        return isinstance(x, iter_types)

    def ensure_list(x):
        return list(x) if isinstance(x, Iterable) and not isinstance(x, (str, bytes)) else x

    def maybe_validate(x):
        return validate(x) if validate else x

    def decorator(func):
        @wraps(func)
        def wrapper(arg1, arg2=None, *args, **kwargs):
            # One-argument mode
            if arg2 is None:
                if is_collection(arg1):
                    return [func(maybe_validate(a1), *args, **kwargs) for a1 in arg1]
                else:
                    return func(maybe_validate(arg1), *args, **kwargs)

            # Two-argument mode
            a1_is_col = is_collection(arg1)
            a2_is_col = is_collection(arg2)

            if not a1_is_col and not a2_is_col:
                return func(maybe_validate(arg1), maybe_validate(arg2), *args, **kwargs)

            if a1_is_col and not a2_is_col:
                return [func(maybe_validate(a1), maybe_validate(arg2), *args, **kwargs) for a1 in arg1]

            if not a1_is_col and a2_is_col:
                return [func(maybe_validate(arg1), maybe_validate(a2), *args, **kwargs) for a2 in arg2]

            # both collections
            if strict_zip:
                l1, l2 = len(arg1), len(arg2)
                if l1 != l2:
                    raise ValueError(f"Length mismatch: {l1=} != {l2=}")
            return [func(maybe_validate(a1), maybe_validate(a2), *args, **kwargs) for a1, a2 in zip(arg1, arg2)]

        return wrapper
    return decorator


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
