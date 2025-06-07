# B-spline Basis Function Benchmark

## Overview

This benchmark evaluates and compares the computational efficiency of multiple implementations of B-spline basis functions. Performance is measured as elapsed computation time across various polynomial degrees and numbers of control points, providing insights into scaling and efficiency.

## Definition of B-spline Basis Functions

B-spline basis functions of degree $d$ are defined using the recursive Cox–de Boor formula. Starting with degree zero $(d = 0)$, the basis functions are piecewise constants:

$$
B_{i,0}(t) = \begin{cases}
1, & k_i \leq t < k_{i+1}, \\
0, & \text{otherwise.}
\end{cases}
$$

For higher degrees $d \geq 1$, the basis functions are defined recursively as:

$$
B_{i,d}(t) = \frac{t - k_i}{k_{i+d} - k_i} B_{i,d-1}(t)
\; + \;
\frac{k_{i+d+1} - t}{k_{i+d+1} - k_{i+1}} B_{i+1,d-1}(t)
$$

Here, $t$ is the evaluation parameter, $k$ represents the knot vector, and $i$ indexes the control points associated with each basis function.

## Benchmark Parameters

* **Methods Evaluated:** `b1`, `b2`, `b3`
* **Polynomial Degrees:** 1 to 5
* **Number of Control Points:** 0 to 1,000,000 (incremented by 10,000)

Parameters can be adjusted in the file `app/bench.f90`.

## System Specifications

Benchmarks were performed on the following system:

| **Specification**    | **Details**                               |
| -------------------- | ----------------------------------------- |
| **Processor**        | Intel(R) Core(TM) i9-9980HK CPU @ 2.40GHz |
| **Memory**           | 64 GB                                     |
| **Operating System** | Ubuntu 24.04.2 LTS                        |
| **Fortran Compiler** | GNU Fortran 15.1.0                        |

## Elapsed Time Plots

The following plots illustrate the elapsed computation time for polynomial degrees 1 through 5:

<div align="center">

<img src="benchmark_degree_1.png" width="600">
<img src="benchmark_degree_2.png" width="600">
<img src="benchmark_degree_3.png" width="600">
<img src="benchmark_degree_4.png" width="600">
<img src="benchmark_degree_5.png" width="600">

</div>

## Compilation and Execution

Compile and run the benchmark using the Fortran Package Manager (`fpm`):

* **GNU Fortran (`gfortran`):**

```bash
fpm @gf
```

* **Intel Fortran (`ifx`):**

```bash
fpm @if
```

* **NVIDIA Fortran (`nvfortran`):**

```bash
fpm @nv
```

* **LLVM Fortran (`flang`):**

```bash
fpm @fl
```

## Customizing Compiler Flags

Modify compiler flags by editing `fpm.rsp` located in the project's root directory.

## Adding New Methods

To add a new B-spline basis function implementation:

1. Implement the method in `src/bspline_basis.f90`.
2. Update the method names array in `app/bench.f90`.
3. Extend the `select case` statement in `app/bench.f90` accordingly.

## Plotting Results

Generate the plots using the provided Python script located in `scripts/plot.py`. From the root directory, execute:

```bash
python scripts/plot.py
```

The resulting plots will be saved in the repository's root directory.
