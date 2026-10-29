# Shared setup for this repo's 01-07 pipeline + fit_refined_params.py. Defines
# an alternate, smaller, more symmetric PM(7)/AT(7) area grouping that mirrors
# the fit/simulate/plot/test logic of the manuscript's primary PM(7)/AT(13)
# analysis (a separate pipeline, not included in this repo) -- see
# ANALYSIS.md for the full area lists, data contract, and rationale.
#
# PM is the same 7 areas as common.py's PM (just a different order here --
# no scientific difference). AT is a different, smaller (7-area) curated
# set: GU, VISC, SS, ACA, MO are dropped outright (no substitute, genuinely
# excluded from this analysis); PL and ILA are NOT dropped -- their signal
# is folded into the single PFC area instead. PFC_length/PFC_endpoint are
# native columns already in the raw CSV (verified: PFC_length ==
# PL_length + ILA_length to float precision for all 160 subset rows), so
# PFC needs no special-case derivation code -- it loads through the exact
# same generic per-area comprehension as every other area below.
#
# Import before matplotlib.pyplot (this module does not import matplotlib).

from pathlib import Path
import numpy as np
import pandas as pd

_HERE     = Path(__file__).parent
_ROOT     = _HERE.parent
DATA_FILE = _ROOT / "data" / "master_detailed_comment.csv"
OUT_DIR   = _ROOT / "outputs"
OUT_DIR.mkdir(exist_ok=True, parents=True)

# Area lists (load-bearing -- see module docstring above)
PM    = ["VISpor", "RSP", "PAR", "PRE", "PTLp", "VIS", "POST"]
AT    = ["PERI", "ORB", "ECT", "TEa", "AI", "PFC", "OLF"]
AREAS = PM + AT
A     = len(AREAS)
ispm  = np.array([a in PM for a in AREAS])
nPM, nAT = len(PM), len(AT)

def sig(x):   return 1.0 / (1.0 + np.exp(-x))
def logit(p): return np.log(p / (1 - p))

# Load PMAT subset (160 neurons with null comment)
df  = pd.read_csv(DATA_FILE)
sub = df[df["comment"].isna()].reset_index(drop=True)
assert len(sub) == 160, f"Expected 160 PMAT neurons, got {len(sub)}"
N   = len(sub)

# Smoke check: PFC is a native column equal to PL+ILA (documented above),
# re-verified here on the loaded subset (not just the raw CSV) so a future
# data-file swap can't silently break this assumption.
assert np.allclose(sub["PFC_length"].values, sub["PL_length"].values + sub["ILA_length"].values, atol=1e-3), \
    "PFC_length no longer matches PL_length + ILA_length -- re-check the data contract"

# Projection binary matrix and axon-length matrix
Y = np.column_stack([(sub[f"{a}_length"].values > 0).astype(float) for a in AREAS])
L = np.column_stack([sub[f"{a}_length"].values for a in AREAS])

# ML coordinate: z-scored then scaled to [-1, 1]. Medial → PM, lateral → AT.
z   = sub["z"].values
z_z = (z - z.mean()) / z.std()
t   = z_z / np.max(np.abs(z_z))

# Per-area logit intercepts (marginal projection prevalence, clipped for stability)
a = logit(np.clip(Y.mean(axis=0), 0.02, 0.98))

# RSP/ORB group indices (used by broadcast-contrast statistics)
rsp_idx   = AREAS.index("RSP")
orb_idx   = AREAS.index("ORB")
other_idx = np.array([j for j in range(A) if j not in (rsp_idx, orb_idx)])

# Axon-length simulation constants (shared by scripts 01 and 02)
LMAX   = 10.0
KSHAPE = 3.0

# Structural sanity checks for this 14-area alternate definition.
assert A == 14 and nPM == 7 and nAT == 7 and len(other_idx) == 12
assert set(PM) == {"VISpor", "RSP", "PAR", "PRE", "PTLp", "VIS", "POST"}
assert set(AT) == {"PERI", "ORB", "ECT", "TEa", "AI", "PFC", "OLF"}
assert not ({"GU", "VISC", "SS", "ACA", "MO", "PL", "ILA"} & set(AREAS))
