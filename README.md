# B-spline Basis Function Benchmark

## Overview

This benchmark measures and compares the computational performance of multiple implementations of B-spline basis function evaluations. Performance is assessed across a range of polynomial degrees and varying numbers of control points.

## B-spline Basis Function Definition

B-spline basis functions of degree $d$ are recursively defined by the Cox–de Boor formula. Starting with degree zero ($d = 0$), the definition is piecewise constant:

$$
B_{i,0}(t) = 
\begin{cases}
1 & \text{if } k_i \leq t < k_{i+1}, \\[5pt]
0 & \text{otherwise.}
\end{cases}
$$

For degrees $d \geq 1$, the functions are defined recursively as follows:

$$
B_{i,d}(t) = \frac{t - k_i}{k_{i+d} - k_i} B_{i,d-1}(t) 
\;+\; 
\frac{k_{i+d+1} - t}{k_{i+d+1} - k_{i+1}} B_{i+1,d-1}(t),
$$

where $t$ is the evaluation parameter, $k$ denotes the knot vector, and $i$ identifies the control point associated with the basis function.

## Benchmark Parameters

* Methods Evaluated: b1, b2, b3
* Polynomial Degrees: 1 to 5
* Number of Control Points: 0 to 1e7, incremented by 10000

Te above parameters can be adjusted in the `app/bench.f90` file.

## Compilation and Execution

Compile and run the benchmark using `fpm`:

**GNU Fortran (`gfortran`)**:

```bash
fpm @gf
```

**Intel Fortran (`ifx`)**:

```bash
fpm @if
```

**NVIDIA Fortran (`nvfortran`)**:

```bash
fpm @nv
```

**LLVM Fortran (`flang`)**:

```bash
fpm @fl
```

## Customizing Compiler Flags

Modify compiler flags by editing the file `fpm.rsp` located in the project's root directory.


## Adding New Methods

To add a new B-spline basis function implementation to the benchmark, follow these steps:

1. Implement the new method in the file `src/bspline_basis.f90`.
2. Update the array of method names in `app/bench.f90`.
3. Add a corresponding case to the `select case` statement in `app/bench.f90`.
