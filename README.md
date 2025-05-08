# dustysn

<p align="center">
  <img src="docs/images/dustysn.png" width="200">
</p>

`dustysn` is a Python package designed to fit the SEDs of supernovae with dust models.

* Documentation: [https://dustysn.readthedocs.io/](https://dustysn.readthedocs.io/)
* Code: [https://github.com/gmzsebastian/dustysn](https://github.com/gmzsebastian/dustysn)
* License: MIT

![Tests](https://github.com/gmzsebastian/dustysn/actions/workflows/ci_tests.yml/badge.svg)
![License](http://img.shields.io/badge/license-MIT-blue.svg)
[![Coverage Status](https://coveralls.io/repos/github/gmzsebastian/dustysn/badge.svg?branch=main)](https://coveralls.io/github/gmzsebastian/dustysn?branch=main)

## Quick Start

The simplest way to use `dustysn` is:

```python
import dustysn 
```

<p align="center">
  <img src="docs/images/example_dustysn.png" width="800">
  <br>
  <em>Example output for a Supernova.</em>
</p>

## Installation

Install using pip:
```bash
pip install dustysn
```

Or install from source:
```bash
git clone https://github.com/gmzsebastian/dustysn.git
cd dustysn
pip install -e .
```

## Requirements

* Python 3.7 or later
* Having the ``emcee`` package installed

## Citation

If you use `dustysn` in your research, please cite:

```bibtex
@software{dustysn,
  author       = {Sebastian Gomez},
  title        = {dustysn: Fitting the Infrared SEDs of Dusty Supernovae},
  year         = {2025},
  publisher    = {GitHub},
  url          = {https://github.com/gmzsebastian/dustysn}
}
```

## License

Copyright 2025 Sebastian Gomez and contributors.

`dustysn` is free software made available under the MIT License. For details see the LICENSE file.
