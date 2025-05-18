# import os
# import pytest
# import numpy as np
from dustysn.utils import (calc_distance)


def test_calc_distance():
    # Check that the get data returns None
    output = calc_distance(0.1).value
    # Check that the output is a float
    assert isinstance(output, float)
