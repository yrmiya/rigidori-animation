# Rigid Origami Animation Tool

[![GitHub release (latest SemVer including pre-releases)](https://img.shields.io/github/v/release/yrmiya/rigidori-animation?color=0033&include_prereleases&sort=semver)](https://github.com/yrmiya/rigidori-animation/releases/latest)

The repository contains a package to create vtk files for rigid origami animation.

The current version supports:

- Single crease fold (a)
- Miura-ori (b)
- Single-orbit Hexagonal-triangular Resch (c)

![Supported origami types][origami]

[origami]: img/origami.png "Supported origami types"

## Requirements

- Python 3.10.x
  - NumPy
  - Matplotlib
  - SciPy
  - tqdm
- ParaView 5.8 (Download from [here](https://www.paraview.org/download/). The latest version is 5.11, but runs slow on laptop and requires powerful machine)

## Python installation

### Installation with installer

Download installer from [Python official website](https://www.python.org/downloads/release/python-3108/), and run installer.

After installation, open terminal (PowerShell for Windows, Terminal for Mac and Linux) and run

```sh
python --version
```

Make sure that it returns Python 3.10.\*.

### Installing packages

To install packages, run

```sh
python -m pip numpy matplotlib scipy tqdm PyQt5
```

## List of packages

- rigidori_animation: Exports animation of the single-orbit hexagonal-triangular Resch origami

## Usage

In root directory of the project, run with origami type specified (ori_type), e.g., for Miura-ori,

```sh
python -m rigidori_animation --ori_type miura
```

For single crease fold and Resch, give "crease" and "resch", respectively, instead of "miura".

Without other arguments, the code will generate vtk files with default parameters.
To define geometry, number of data points, and other options, see

```sh
python -m resch_animation -h
```

for optional arguments.

## Author

- Yasuhiro Miyazawa (LEMS, Dept. of Aero. & Astro., UW)
- Ted Chang (LEMS, Dept. of Aero. & Astro., UW)
