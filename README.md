# pyOrigamiBreak

## How to use

To use this tool, we recommend using this [Colab notebook](https://colab.research.google.com/drive/1NbS7lE3wZkecC5-zXZ0KpfDWn_a3YFDE)

## How to cite

We invite you to support the project by citing these references:

```
Design principles for accurate folding of DNA origami
Aksel et al. [TBD] (2024)
```

```
Rapid prototyping of 3D DNA-origami shapes with caDNAno
Douglas et al. Nucleic Acids Res: 37(15):5001–6 (2009)
https://doi.org/10.1093/nar/gkp436
```

## Development

We use a pyproject.toml-based [build process](https://pip.pypa.io/en/stable/reference/build-system/pyproject-toml/) in pip. This workflow was tested with python 3.12 on macOS in late 2023.

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
* `python3 -m twine upload dist/pyOrigamiBreak-x.y.z

## License

This version of pyOrigamiBreak is available under the MIT License.