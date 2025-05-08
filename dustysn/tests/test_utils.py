# import os
# import pytest
# import numpy as np
from dustyn.utils import (get_data)


def test_get_data():
    # Check that the get data returns None
    output = get_data('potato')
    assert output is None
