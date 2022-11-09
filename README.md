# SNU Guest Lecture

Some resources for the guest lecture held on 2022.11.15-17.

## Requirements

- Python 3.10.x
- ParaView (Can be downloaded from [here](https://www.paraview.org/download/))

## List of packages

- rigidori_animation: Exports animation of the single-orbit hexagonal-triangular Resch origami

## Usage

In root directory of the project, run with origami type specified (ori_type)

```sh
python -m rigidori_animation --ori_type [ori_type]
```

Without other arguments, the code will generate vtk files with default parameters.
To define geometry, number of data points, and other options, see

```sh
python -m resch_animation -h
```
