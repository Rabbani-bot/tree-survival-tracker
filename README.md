# Plantation Survival Tracker (Kaplan–Meier)

[![tests](https://github.com/Rabbani-bot/tree-survival-tracker/actions/workflows/tests.yml/badge.svg)](https://github.com/Rabbani-bot/tree-survival-tracker/actions/workflows/tests.yml) ![Python](https://img.shields.io/badge/python-3.11%20%7C%203.12-blue) ![License: MIT](https://img.shields.io/badge/license-MIT-green)

> **Dataset.** 2,400 saplings across 8 sites and 4 species, checked at 1, 3, 6, 12, 18 and 24 months. Trees with lost tags are handled as censored observations rather than counted as dead.

![survival](outputs/survival.png)

## Why survival analysis
Counting "trees planted" overstates impact. Survival analysis answers what funders actually ask: **how many trees are alive after two years, and what drives losses?** It handles trees whose tags were lost, or that haven't reached later visits yet, without biasing the result.

## Methods
- Kaplan–Meier estimator with Greenwood 95% confidence intervals, implemented from first principles and checked against a textbook example in the tests
- Log-rank test comparing species and slope classes
- 24-month survival by site, for targeting replanting

## Results
| Metric | Value |
|---|---|
| Trees / deaths observed / censored | 2,400 / 653 / 1,747 |
| 24-month survival by species | Eucalyptus 82.6%, chir pine 75.6%, mulberry 68.3%, shisham 61.6% |
| Species difference (log-rank) | chi-square 69.5, p < 0.001 |
| Steep vs gentle slope (log-rank) | chi-square 14.0, p < 0.001 |

## Run
```bash
python3 -m unittest discover -s tests -v
python3 src/survival.py
```

Author: Mohammad Rabbani – Sustainability Data Analyst & Research Consultant

## Links
- Kaggle notebook (runs end-to-end on Kaggle): https://www.kaggle.com/code/rabbani0123/tree-plantation-survival-analysis-kaplan-meier
- Cite: see `CITATION.cff`
- Author: Mohammad Rabbani | Sustainability Data Analyst & Research Assistant

## License
MIT
