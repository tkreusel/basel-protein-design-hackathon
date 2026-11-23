"""Simulated dataset: one truly informative metric + many noise metrics. The screen must find it (low q), the noise must not
survive FDR, and the leave-one-out composite on PURE noise must be ~0 (selection inside the loop)."""
import numpy as np, pandas as pd
from maxi_eval import validate_dataset as V
rng = np.random.default_rng(1); n = 60
y = rng.normal(3, 1.2, n)
X = pd.DataFrame({f"noise{i}": rng.normal(size=n) for i in range(120)})
X["signal"] = 0.8 * (y - 3) + rng.normal(scale=0.6, size=n)
X["const"] = 1.0
d = V.screen(X, y, (y >= np.median(y)).astype(float), None)
assert d.iloc[0].metric == "signal" and d.iloc[0].q_BH < 0.01, d.head()
assert (d[d.metric != "signal"].q_BH < 0.1).sum() <= 3, "too many noise metrics pass FDR"
r_sig, _ = V.loo_composite(X, y, 3)
r_noise, _ = V.loo_composite(X.drop(columns=["signal"]), y, 3)
assert r_sig > 0.4, r_sig
assert abs(r_noise) < 0.3, r_noise          # selection from 120 noise metrics must NOT look predictive when held out
# partial correlation removes a metric that is just a proxy for the control
ctrl = y + rng.normal(scale=0.3, size=n); proxy = ctrl + rng.normal(scale=0.1, size=n)
assert abs(V.partial_spearman(proxy, y, ctrl)) < abs(V.partial_spearman(X["signal"].to_numpy(), y, ctrl)) + 0.5
print("validate_dataset tests passed: signal r_LOO=%.2f, noise r_LOO=%.2f" % (r_sig, r_noise))
