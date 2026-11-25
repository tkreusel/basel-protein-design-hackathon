import sys
from pathlib import Path

import pandas as pd
import pytest

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent / "src"))
FIXTURES = HERE / "fixtures"
RUNS = ["hgnt__steer_tmol__b01", "hgwt__unsteered__b01", "hgnt__unsteered_sheet__b51"]


def read(*parts):
    """Read a fixture CSV (plain or gzipped)."""
    p = FIXTURES.joinpath(*parts)
    return pd.read_csv(p if p.exists() else str(p) + ".gz", low_memory=False)


@pytest.fixture(params=RUNS)
def run_name(request):
    return request.param
