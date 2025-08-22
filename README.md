# c3d: Real-time 3D Renderer

## Usage

To build the project:
```
make build
```

To run:
```
./c3d <path to model> <path to texture?
```

Some example models are included in the `models/` directory, and some example textures
are included in the 'textures/' directory.

As of now, the code has to be edited to change parameters such as FOV and screen resolution.

## Features

 - [x] Load models from file
 - [x] Custom normals
 - [x] Smooth shading
 - [x] Support for image textures (PPM only)
 - [ ] `.mtl` file handling
 - [ ] Config file format/command line arguments

## License

This project is licensed under the MIT License.
