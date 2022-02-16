# Installation

OpenExplorer is distributed in its stable and testing version through the 'uibcdf' conda channel.
If there is no reason to install the library from the source code, we highly recommend working with
conda.

## Last stable version

There is no stable version yet

## Last developing version

If you want to work with the last testing version:

```bash
conda install -c uibcdf/label/dev openexplorer
```

To uninstall this library:

```bash
conda remove openexplorer
```

## The source code

The raw code fully alive can be cloned from the [github repository](https://github.com/uibcdf/OpenExplorer) as follows:

```bash
git clone https://github.com/uibcdf/OpenExplorer.git
```

Or with GitHub CLI:

```bash
gh repo clone uibcdf/OpenExplorer
```

Now, once you have cloned the repository. You need to install the required pakages to use it,
develope it, test it, or document it. Find all required libraries, depending on each usage case, in
the `OpenExplorer/devtools/conda-envs` directory. And if you want to create a conda environment to play
with OpenExplorer, feel free to make use of the Python scripts 'create\_conda\_env.py' and
'update\_conda\_env.py' in the same directory:

```
cd OpenExplorer/devtools/conda-envs
python create_conda_env.py -n OpenExplorer -p 3.7 production_env.yaml
```

You can now install the developing version of OpenExplorer from the source code:

```
conda activate OpenExplorer
cd OpenExplorer
python setup.py develop
```

