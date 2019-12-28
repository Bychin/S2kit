# S2kit - FFT and Convolution on Sphere

S2kit, version 1.1, is a library of C functions which compute forward and inverse discrete Fourier transforms, convolutions of functions defined on the sphere <img src="https://latex.codecogs.com/svg.latex?\Large&space;S^2" title="S^2"/>.

This is an updated version of the library S2kit 1.0, that was initially published there www.cs.dartmouth.edu/~geelong/sphere.

S2kit is free software and is distributed under the terms of the GNU General Public License.

## Documentation

Paper for S2kit 1.0 can be found [here](https://github.com/Bychin/S2kit/blob/master/dist/S2kitHowTo.pdf). I suggest you to read it first. More theoretical exposition can be found in the article "Computing Fourier Transforms and Convolutions on the 2-Sphere" by Driscoll J.R. and Healy D.M.

Full documentation of the library can be found [here](https://bychin.github.io/S2kit).

The list of core functions:

Name | Description
--- | ---
[`DLTNaive()`](https://bychin.github.io/S2kit/html/naive_8c.html) | Legendre transform using naive algorithm
[`InvDLTNaive()`](https://bychin.github.io/S2kit/html/naive_8c.html) | Inverse Legendre transform using naive algorithm
[`DLTSemi()`](https://bychin.github.io/S2kit/html/seminaive_8c.html) | Legendre transform using seminaive algorithm
[`InvDLTSemi()`](https://bychin.github.io/S2kit/html/seminaive_8c.html) | Inverse Legendre transform using seminaive algorithm
[`FSTSemiMemo()`](https://bychin.github.io/S2kit/html/_f_s_t__semi__memo_8c.html) | Spherical harmonic transform
[`InvFSTSemiMemo()`](https://bychin.github.io/S2kit/html/_f_s_t__semi__memo_8c.html) | Inverse spherical harmonic transform
[`FZTSemiMemo()`](https://bychin.github.io/S2kit/html/_f_s_t__semi__memo_8c.html) | Zonal harmonic transform
[`ConvOn2SphereSemiMemo()`](https://bychin.github.io/S2kit/html/_f_s_t__semi__memo_8c.html) | Convolution of two functions defined on the 2-sphere
[`FSTSemiFly()`](https://bychin.github.io/S2kit/html/_f_s_t__semi__fly_8c.html) | Spherical harmonic transform
[`InvFSTSemiFly()`](https://bychin.github.io/S2kit/html/_f_s_t__semi__fly_8c.html) | Inverse spherical harmonic transform
[`FZTSemiFly()`](https://bychin.github.io/S2kit/html/_f_s_t__semi__fly_8c.html) | Zonal harmonic transform
[`ConvOn2SphereSemiFly()`](https://bychin.github.io/S2kit/html/_f_s_t__semi__fly_8c.html) | Convolution of two functions defined on the 2-sphere

In most cases you should prefer to use `-Semi()` and `-Memo()` variants of the functions, since they are faster and accept reusing passed arguments (see more info in documentation).

## Examples

Complete examples can be found in documentation in `test` directory. It contains examples of call and use functions mentioned above.

### Changelog (from version 1.0)

* All files were formatted via clang-formatter
* Library's file structure was updated
  * `test_*` to `test/`
  * headers to `include/`
  * sources to `src/`
  * data for tests to `data/`
* All math consts and functions are imported from `<math.h>`
* All source and test files were refactored in prior to C11 standart and code style
* Due to huge amount of refactoring there were many optimization fixes

### Backward compatibility

// TODO create md file with history of renaming and moving files/functions
