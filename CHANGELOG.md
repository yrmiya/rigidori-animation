# Change Log

## v0.7.1 (2023-11-21)

- Updated `README.md`
- Added `LICENSE.md`

## v0.7.0 (2023-02-28)

### Added

### Changed

- Change Miura-ori model
  - Now support Miura-ori tessellation
  - Miura-ori unit definition is now single panel (i.e., 2x2 is conventional unit cell with 4 facets)
    - Default numbers of unit in x and y direction are set to 2 for both.
  - Orientation is changed

### Fixed

## v0.6.0 (2023-02-20)

### Added

- Add node number to projection
- Kresling geometry initializer

### Changed

- Change package name from `rigidori_animation` to `ori-anime`
- Update readme accordingly

### Fixed

- Remove empty entry in Resch vtk
- Add missing panels in Resch vtk
- Remove `./` from relative path import

## v0.5.2 (2022-11-27)

### Added

- Add CHANGELOG.md

### Changed

- Move plot_projection function to independent module
- Move calc_facet_area function to independent module

### Fixed

- Fix vtk file object name

## v0.5.1 (2022-11-24)

### Fixed

- Fix number of arguments passed onto packages

## v0.5.0 (2022-11-23)

### Added

- Add sub-commands "run" and "clean"
- Add triangular startuck geometry and example file

## v0.4.0 (2022-11-16)

### Changed

- Integrate all packages into rigidori_animation

## v0.3.1 (2022-11-08)

### Changed

- Modify plot_projection of resch_animation to include all edges and polygons
- Change boolean argument of argparser to store_true or store_false
- Update readme

## v0.3.0 (2022-11-08)

### Added

- Add animation package for Single-orbit Hexagon-triangle Resch origami
- Add example of vtk

## v0.2.0 (2022-11-05)

### Added

- Add Miura-ori animation package

## v0.1.1 (2022-11-03)

### Added

- Add optional arguments:
  - --theta0 for specifying initial fold angle
  - --thetaf for specifying final fold angle
  - --savezip for saving zip file containing vtk files
  - --figout for displaying plots

## v0.1.0 (2022-11-02)

### Added

- Add single crease fold animation package.
