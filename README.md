# Kuwahara Filters for MATLAB

MATLAB implementations of classic, planar and polynomial Kuwahara filters for grayscale and RGB images. The regression variants estimate image derivatives together with the filtered intensity.

This standalone package was extracted from a larger collection of experiments. It contains the three filters, a synthetic example and numerical checks. Historical benchmark scripts, generated figures, datasets and third-party implementations are excluded.

## Requirements

MATLAB is required. The filters and included tests use base MATLAB and do not require Image Processing Toolbox. The original workspace uses MATLAB R2025b; compatibility with earlier releases has not been established.

## Quick start

Download the repository, open MATLAB and change to its directory. Run:

```matlab
demoKuwaharaPlanar
```

The example creates a noisy piecewise affine image and displays the input, filtered image, gradient magnitude and residual variance. It requires no external dataset.

To process your own image:

```matlab
I = imread('your-image.png');
Jclassic = kuwaharaClassic(I, 3);
[Jplanar, gx, gy, info] = kuwaharaPlanar(I, 3);
[Jquadratic, qx, qy, polyInfo] = kuwaharaPolynomial(I, 3, 2);
imwrite(min(1, max(0, Jquadratic)), 'filtered.png');
```

## Methods

All methods examine four overlapping quadrants around each pixel. Each quadrant contains `(radius+1)^2` samples and includes the central pixel. Exact ties select the first quadrant in this order: top-left, top-right, bottom-left, bottom-right.

| Function | Local model | Selection score | Output at the center |
| --- | --- | --- | --- |
| `kuwaharaClassic` | Constant mean | Mean squared deviation from the mean | Quadrant mean |
| `kuwaharaPlanar` | `a + b*x + c*y` | `SSE/(N-3)` | `a`, with gradient `(b,c)` |
| `kuwaharaPolynomial` | All bivariate monomials of total degree at most `degree` | `SSE/(N-m)` | Polynomial value and derivatives |

Here `N=(radius+1)^2`, `m=(degree+1)*(degree+2)/2`, and SSE is the sum of squared residuals. The polynomial implementation requires `N>m` and full column rank in every quadrant. The sample-count condition alone is insufficient: degree 4 with radius 3 is rank deficient. Increase the radius or reduce the degree when a rank error occurs.

The polynomial degree defaults to 2. Degree 1 is equivalent to the planar model up to floating-point error. Degree 0 uses the residual score `SSE/(N-1)` and returns zero gradients; its variance values differ from the classic filter's normalization by `N`.

## Input and boundary conventions

- Inputs must be nonempty, real, full HxW grayscale or HxWx3 RGB arrays.
- Supported types: `logical`, `uint8`, `uint16`, `single`, `double`.
- Unsigned integers are divided by their type's maximum value. Floating-point values must be finite and in `[0,1]`.
- RGB is converted to grayscale using weights `[0.2989360213, 0.5870430745, 0.1140209042]`. Outputs are grayscale double arrays.
- `radius` must be a positive integer. The full window has width `2*radius+1`.
- Boundaries use symmetric reflection with repeated edge pixels. Affine and polynomial reproduction properties apply to interior pixels, away from reflected boundaries.
- Regression outputs are not clipped to `[0,1]`. Clip only when needed for display or image export.
- Coordinates point right (`x`) and down (`y`); derivatives are measured per pixel. Gradient orientation is `atan2(gy,gx)` in radians.

## Diagnostic outputs

The classic filter returns `info.Sector` and `info.MinimumVariance`. The regression filters return `info.Magnitude`, `info.Orientation`, `info.Sector` and `info.MinimumResidualVariance`.

The polynomial filter additionally returns `info.Degree`, `info.Exponents` and `info.Derivatives`. Each exponent pair `[alpha,beta]` identifies the corresponding mixed derivative page. The first page is the filtered intensity. For degree at least 1, the next two pages are `gx` and `gy`. For degree 0, only the intensity page exists.

## Validation

Run the supplied checks from the repository directory:

```matlab
testKuwaharaPlanar
testKuwaharaPolynomial
```

They cover grayscale/RGB equivalence, affine and quadratic reproduction, first and higher-order derivatives, degree-1 equivalence, constant images and invalid inputs/rank failures. The classic filter does not have a dedicated test function. See [VALIDATION.md](VALIDATION.md) for checks performed while preparing this package.

The implementations prioritize readability and use explicit pixel loops. Large images and high polynomial degrees may take substantial time. This package makes no claim of general superiority over other filters; benchmark results from the source workspace are not included.

## Names in the source workspace

| Original name | English package name |
| --- | --- |
| `kuwaharaClassico` | `kuwaharaClassic` |
| `kuwaharaPlanare` | `kuwaharaPlanar` |
| `kuwaharaPolinomiale` | `kuwaharaPolynomial` |
| `demoKuwaharaPlanare` | `demoKuwaharaPlanar` |
| `testKuwaharaPlanare` | `testKuwaharaPlanar` |
| `testKuwaharaPolinomiale` | `testKuwaharaPolynomial` |

Diagnostic field names were also translated, for example `Modulo` → `Magnitude`, `Settore` → `Sector`, and `Derivate` → `Derivatives`. Existing scripts that use the original names must be updated when using this package.

## License

No license has been selected for this package. Publication on GitHub does not by itself grant an open-source license. The copyright holder should choose and add a license before distributing the code under specific reuse terms.
