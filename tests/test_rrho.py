import pytest
import numpy as np
from src.rrho import calc_zpe


def test_calc_zpe():
    # Test case 1: frequency array is empty
    freq = np.array([])
    expected_zpe = 0
    assert np.isclose(calc_zpe(freq), expected_zpe)

    # Test case 2: frequency array has one value
    freq = np.array([1000])
    expected_zpe = 0.0022781676264559806
    assert np.isclose(calc_zpe(freq), expected_zpe)

    # Test case 3: frequency array has multiple values
    freq = np.array([1000, 2000, 3000])
    expected_zpe = 0.013669005758735883
    assert np.isclose(calc_zpe(freq), expected_zpe)
