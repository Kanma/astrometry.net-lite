# astrometry.net lite

This is a stripped-down version of [astrometry.net](https://github.com/dstndstn/astrometry.net)'s
source code (version 0.95.0), that only allows to:

  * detect stars in an image
  * perform plate solving to retrieve the celestial coordinates of the image

but has the following characteristics:

  * a single library that works on in-memory images (no need to write intermediate files anymore)
  * no reliance on external Python scripts
  * no GPL dependencies: ```astrometry.net``` is licensed under a ```BSD 3-Clause``` license,
  but relies on several GPL libraries. They were replaced by alternatives in
  ```astrometry.net lite```:
    * ```GSL``` by ```Eigen```
    * ```qfits``` by ```cfitsio```
  * Works on Linux, MacOS and Windows (tested on Visual Studio 2022)


An example program is provided, but you might want to look at how I use this library in
[astrophoto-toolbox](https://github.com/Kanma/astrophoto-toolbox) for a more complete usage.


## Compilation

```
$ mkdir build
$ cd build
$ cmake ..
```


## Running the example

To perform plate solving, you need to download index files. The recommended ones (unless you
know what you are doing) are available from the following URLs (you likely want all files at
both URLs, put them all in the same folder):

  * [http://data.astrometry.net/4100/](http://data.astrometry.net/4100/)
  * [https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/](https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/)

Then run the example program, giving it an image and the path to the folder containing the
index files, for example:

```
$ build/bin/example example/starfield.jpg /path/to/indexes
Loading the image 'example/starfield.jpg'...
Detecting stars...
   4646 stars found
Loading the index files from '/path/to/indexes'...
    150 index files used
Plate solving...
    282.654°, -12.9419°
    Orientation: 109.161°
    Pixel size: 1.17377 arcsec
```

Note that the variables ```minWidth``` and ```maxWidth``` are hard-coded to values suitable
for the example image, you might need to modify them to correspond to the parameters of your
image. They are respectively the minimum and maximum sizes (in degrees) of the image.


## Dependencies

CMake is required to compile ```astrometry.net lite```.

Other dependencies (*cfitsio*, *Eigen*, *zlib*) are automatically downloaded and compiled
alongside the library.

The example program use ```stb_image.h``` (a single header library) to load images.


## License

```astrometry.net lite``` is licensed under a BSD 3-Clause license.
