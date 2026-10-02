r"""The script behind Magnus's points on the Earth speed-accuracy planes still runs.

``notebooks/gen_prem_plane_magnus.py`` passed the neutron-to-proton ratio to the three-flavor
Earth wrapper, which refuses it since 1.2.0, so the three-flavor panel could not be rebuilt.
Each panel is checked at one point against the error stored for the paper's figure.
"""
import pathlib
import sys
import warnings

import pytest

NOTEBOOKS = pathlib.Path(__file__).resolve().parents[1]/'notebooks'
sys.path.insert(0, str(NOTEBOOKS))


@pytest.mark.checkout_only     # reads notebooks/ and resources/benchmarks/
@pytest.mark.parametrize('panel, knob', [('3nu', 32), ('3+1', 128)])
def test_point_reproduces_stored_error(panel, knob):
    import gen_prem_plane_magnus as gpp
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        got = gpp.point(panel, knob, timed=False)['max_abs_error']
    stored = [q['max_abs_error'] for pts in gpp.stored_series(panel).values()
              for q in pts if q['knob'] == knob]
    assert stored, 'no stored Magnus point at knob %d' % knob
    assert got == pytest.approx(stored[0], rel=1.0e-3)
