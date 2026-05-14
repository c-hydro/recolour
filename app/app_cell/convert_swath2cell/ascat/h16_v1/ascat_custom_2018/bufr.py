import numpy as np

from eccodes import (
    codes_bufr_new_from_file,
    codes_set,
    codes_get,
    codes_get_array,
    codes_release,
)


class BUFRReader(object):
    """
    BUFR reader using ECMWF ecCodes.

    Compatible with code expecting:
        for message in bufr.messages():
            message[:, 12]
            message[:, 13]
            message[:, mid - 1]

    Each yielded message is a 2D NumPy array:
        shape = (number_of_subsets, number_of_numeric_values_per_subset)
    """

    def __init__(self, filename, kelem_guess=500, max_tries=10):
        self.filename = filename
        self.kelem_guess = kelem_guess
        self.max_tries = max_tries
        self.bufr = open(filename, "rb")

    def messages(self):
        while True:
            gid = codes_bufr_new_from_file(self.bufr)

            if gid is None:
                break

            try:
                codes_set(gid, "unpack", 1)

                n_subsets = int(codes_get(gid, "numberOfSubsets"))

                values = codes_get_array(gid, "numericValues")
                values = np.asarray(values, dtype=np.float64)

                if n_subsets <= 0:
                    n_subsets = 1

                if values.size == 0:
                    yield np.empty((0, 0), dtype=np.float64)
                    continue

                if values.size % n_subsets != 0:
                    raise IOError(
                        "Cannot reshape BUFR message. "
                        "numericValues size={} numberOfSubsets={}".format(
                            values.size,
                            n_subsets
                        )
                    )

                n_fields = values.size // n_subsets
                values = values.reshape((n_subsets, n_fields))

                yield values

            finally:
                codes_release(gid)

    def close(self):
        if self.bufr is not None:
            self.bufr.close()
            self.bufr = None

    def __enter__(self):
        return self

    def __exit__(self, exc, val, trace):
        self.close()