# SNU Guest Lecture

Some resources for the guest lecture held on 2022.11.15-17.

## Requirements

- Python 3.10.x
- ParaView (Can be downloaded from [here](https://www.paraview.org/download/))

## List of packages

- resch_animation: Exports animation of the single-orbit hexagonal-triangular Resch origami
- sc_animation: Exports single crease folding animation

## Usage

### Single-crease

In root directory of the project, run

```sh
python -m sc_animation
```

For optional arguments, see

```sh
python -m sc_animation -h
```

### Miura-ori animation

In root directory of the project, run

```sh
python -m miura_animation
```

Without arguments, the code generates single unit Miura-ori animation.
To define geometry, use the following options.

List of arguments can also be found by

```sh
python -m miura_animation -h
```
