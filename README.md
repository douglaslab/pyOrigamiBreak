# pyOrigamiBreak

## How to use

For an easy and interactive way to use Autobreak with your Cadnano designs, access the following [Colab notebook](https://colab.research.google.com/drive/1oavYcHN08N5lSMS1i5dtf-f3-j8dEe0-?usp=sharing). It provides a step-by-step guide to integrate Autobreak into your workflow.


## How to cite

We invite you to support the project by citing these references:

```
Design principles for accurate folding of DNA origami
Aksel et al. [TBD] (2024)
```
Upon publication, the complete citation details will be provided here. Please check back for updates or contact the authors for the most current information.

```
Rapid prototyping of 3D DNA-origami shapes with caDNAno
Douglas et al. Nucleic Acids Res: 37(15):5001–6 (2009)
https://doi.org/10.1093/nar/gkp436
```

## How to get help

- Before contacting us, please read the manuscript, this README, and all documentation in the Colab notebook.
- If you encounter a technical problem,  submit a [New Issue](https://github.com/douglaslab/pyOrigamiBreak/issues). Please describe the exact steps to reproduce it. Attach or email us the Cadnano file input, as well as any file outputs that resulted from an attempted run.


## Development

We use a pyproject.toml-based [build process](https://pip.pypa.io/en/stable/reference/build-system/pyproject-toml/) in pip. This workflow was tested with Python 3.12 on macOS and Linux distributions in late 2023.

**Setup a dev environment (Mac or Linux)**

* Create a virtualenv: `python3 -m venv ~/virtualenvs/abdev` 
* Activate virtualenv: `source ~/virtualenvs/abdev/bin/activate`
* Clone repo: `git clone git@github.com:douglaslab/pyOrigamiBreak.git`
* Change directory: `cd pyOrigamiBreak`
* Make desired code edits
* Build and install in [editable mode](https://pip.pypa.io/en/stable/cli/pip_install/#cmdoption-e): `pip install -e .` 
* Test: `autobreak -i cadnanofile.json`
* Repeat previous 3 steps as needed

**Build new dist and upload to PyPi**

* `pip install build twine` <- install [build](https://pypi.org/project/build/) and [twine](https://pypi.org/project/twine/)
* `cd /path/to/pyOrigamiBreak/` 
* `python3 -m build`  creates dist/pyOrigamiBreak-x.y.z.tar.gz and pyOrigamiBreak-x.y.z-py3-none-any.whl
* `python3 -m twine upload dist/pyOrigamiBreak-x.y.z`

## License

This version of pyOrigamiBreak is available under the MIT License.