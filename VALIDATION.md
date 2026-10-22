# Validation record

Package preparation date: 2026-10-22.

- MATLAB R2025b Code Analyzer (`mlint.exe`) was run separately on all six MATLAB files and reported no diagnostics.
- English function names match their filenames, and the included example and tests use the renamed functions and diagnostic fields.
- SHA-256 checks confirmed that all 178 files in the original `matlab_kuwahara` directory were unchanged.
- Numerical MATLAB tests have not been confirmed. The first batch launch failed because MATLAB could not write to its default preferences directory. A second launch with a local preferences directory did not return test output during preparation.

To complete runtime validation locally, change to this directory in MATLAB and run:

```matlab
testKuwaharaPlanar
testKuwaharaPolynomial
assert(max(abs(kuwaharaClassic(0.4*ones(5),3)-0.4), [], 'all') < 1e-12)
demoKuwaharaPlanar
```

The last command opens a figure for visual inspection. Static analysis does not establish numerical correctness or successful execution.
