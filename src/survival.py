"""Tree-plantation survival tracking: Kaplan-Meier estimates and log-rank test.

Field dataset generated with a fixed seed. Sites and species mix are modelled to be plausible for
northern Pakistan. Trees are checked at fixed visits; trees still alive at the last visit
(or not found at a visit) are right-censored.
"""
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
SEED = 3
VISITS = np.array([1, 3, 6, 12, 18, 24])          # months after planting
SPECIES_MEDIAN = {"Pinus roxburghii (chir pine)": 70, "Dalbergia sissoo (shisham)": 40,
                  "Morus alba (mulberry)": 55, "Eucalyptus camaldulensis": 90}   # median lifetime, months


def generate(seed: int = SEED, n: int = 2400) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    species = rng.choice(list(SPECIES_MEDIAN), n)
    site = rng.choice([f"Site {i}" for i in range(1, 9)], n)
    slope = rng.choice(["gentle", "steep"], n, p=[.6, .4])
    scale = np.vectorize(SPECIES_MEDIAN.get)(species) / np.log(2) * np.where(slope == "steep", 0.7, 1.0)
    true_death = rng.exponential(scale)
    lost = rng.random(n) < 0.05                                  # tag lost: censored at a random visit
    lost_at = rng.choice(VISITS, n)
    rows = []
    for i in range(n):
        last_alive, event_time, event = 0, None, 0
        for v in VISITS:
            if lost[i] and v > lost_at[i]:
                break
            if true_death[i] <= v:
                event_time, event = v, 1               # death observed at this visit
                break
            last_alive = v
        rows.append((i + 1, site[i], species[i], slope[i], event_time if event else last_alive, event))
    return pd.DataFrame(rows, columns=["tree_id", "site", "species", "slope", "time_months", "died"])


def kaplan_meier(t: np.ndarray, e: np.ndarray) -> pd.DataFrame:
    times = np.unique(t[e == 1]); s = 1.0; var = 0.0; out = [(0, 1.0, 1.0, 1.0)]
    for u in times:
        at_risk = (t >= u).sum(); d = ((t == u) & (e == 1)).sum()
        s *= 1 - d / at_risk
        var += d / (at_risk * (at_risk - d)) if at_risk > d else 0
        se = s * np.sqrt(var)                                   # Greenwood
        out.append((u, s, max(0, s - 1.96 * se), min(1, s + 1.96 * se)))
    return pd.DataFrame(out, columns=["time", "survival", "ci_low", "ci_high"])


def logrank(df: pd.DataFrame, group: str) -> tuple[float, float]:
    t, e, g = df["time_months"].to_numpy(), df["died"].to_numpy(), df[group].to_numpy()
    levels = np.unique(g); O = np.zeros(len(levels)); E = np.zeros(len(levels)); V = np.zeros((len(levels),) * 2)
    for u in np.unique(t[e == 1]):
        at = t >= u; n = at.sum(); d = ((t == u) & (e == 1)).sum()
        nk = np.array([(at & (g == l)).sum() for l in levels]); dk = np.array([((t == u) & (e == 1) & (g == l)).sum() for l in levels])
        O += dk; E += d * nk / n
        if n > 1:
            V += d * (n - d) / (n - 1) * (np.diag(nk / n) - np.outer(nk / n, nk / n))
    k = len(levels) - 1
    diff = (O - E)[:k]; chi2 = float(diff @ np.linalg.pinv(V[:k, :k]) @ diff)
    return chi2, float(stats.chi2.sf(chi2, k))


def main() -> dict:
    df = generate(); df.to_csv(ROOT / "data" / "plantation_checks.csv", index=False)
    res = {"trees": len(df), "deaths_observed": int(df["died"].sum()), "censored": int((df["died"] == 0).sum())}
    fig, ax = plt.subplots(1, 2, figsize=(14, 5))
    surv24 = {}
    for sp, g in df.groupby("species"):
        km = kaplan_meier(g["time_months"].to_numpy(), g["died"].to_numpy())
        ax[0].step(km["time"], km["survival"], where="post", label=sp)
        ax[0].fill_between(km["time"], km["ci_low"], km["ci_high"], step="post", alpha=0.15)
        surv24[sp] = round(100 * km["survival"].iloc[-1], 1)
    ax[0].set_ylim(0, 1.02); ax[0].set_xlabel("months after planting"); ax[0].set_ylabel("survival probability")
    ax[0].set_title("Kaplan-Meier survival by species (95% CI)"); ax[0].legend(fontsize=8)
    for sl, g in df.groupby("slope"):
        km = kaplan_meier(g["time_months"].to_numpy(), g["died"].to_numpy())
        ax[1].step(km["time"], km["survival"], where="post", label=f"{sl} slope")
    ax[1].set_ylim(0, 1.02); ax[1].set_xlabel("months after planting"); ax[1].set_title("Survival by slope"); ax[1].legend()
    fig.suptitle("Plantation survival tracker", fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.92]); fig.savefig(ROOT / "outputs" / "survival.png", dpi=130); plt.close(fig)
    res["survival_at_24_months_pct"] = surv24
    c, p = logrank(df, "species"); res["logrank_species"] = {"chi2": round(c, 1), "p_value": p}
    c, p = logrank(df, "slope"); res["logrank_slope"] = {"chi2": round(c, 1), "p_value": p}
    site = df.groupby("site").apply(lambda g: round(100 * kaplan_meier(g["time_months"].to_numpy(), g["died"].to_numpy())["survival"].iloc[-1], 1))
    res["survival_24m_by_site_pct"] = site.to_dict()
    (ROOT / "outputs" / "summary.json").write_text(json.dumps(res, indent=2))
    return res


if __name__ == "__main__":
    print(json.dumps(main(), indent=2))
