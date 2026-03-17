# PCE and aPCE Toolbox (MATLAB)

This project provides MATLAB implementations of **Polynomial Chaos Expansion (PCE)** and **arbitrary Polynomial Chaos Expansion (aPCE)**.  
Given input-output samples, you can automatically or manually fit PCE/aPCE models for prediction, error analysis, and more.

---

## Table of Contents
- [Introduction](#introduction)
- [File Description](#file-description)
- [Installation](#installation)
- [Quick Start](#quick-start)
  - [1. Data Normalization](#1-data-normalization)
  - [2. Fitting an aPCE Model](#2-fitting-an-apce-model)
  - [3. Prediction](#3-prediction)
  - [4. Automatic Order Selection](#4-automatic-order-selection)
- [Function Documentation](#function-documentation)
  - [pce.m Functions](#pcem-functions)
  - [apce.m Functions](#apcem-functions)
- [Important Notes](#important-notes)
- [Contributing](#contributing)
- [License](#license)
- [Author](#author)

---

## Introduction

**Polynomial Chaos Expansion (PCE)** represents a model's output as a series expansion in terms of orthogonal polynomials of the input random variables. It is widely used in uncertainty quantification, sensitivity analysis, and related fields.  
**arbitrary Polynomial Chaos Expansion (aPCE)** generalizes PCE by constructing orthogonal polynomial bases that adapt to the data distribution using statistical moments of the input data, eliminating the need to assume a specific input distribution (e.g., uniform, Gaussian).

This toolbox provides:
- Standard PCE based on Legendre polynomials (inputs must be normalized to `[-1, 1]`)
- aPCE that constructs orthogonal polynomials from data moments (handles arbitrarily distributed inputs)
- Automatic selection of the expansion order
- Normalization / denormalization helper functions
- Error calculation (relative L2 norm)

---

## File Description

| File Name | Description                                           |
| --------- | ----------------------------------------------------- |
| `pce.m`   | Standard PCE implementation (basis: Legendre polynomials) |
| `apce.m`  | Arbitrary PCE implementation (basis constructed from input data moments) |

Both files share a consistent interface for easy switching.

---

## Installation

1. Download `pce.m` and `apce.m` to any local folder.
2. In MATLAB, add the folder to the search path:
   ```matlab
   addpath('folder_path');
   ```
   Or right-click the folder → **Add to Path** → **Selected Folders**.

---

## Quick Start

The following example demonstrates fitting an aPCE model to the simple one-dimensional function `y = x^2` and making predictions at new points.

### 1. Data Normalization
Inputs should typically be normalized to `[-1, 1]` to avoid numerical issues. `apce()` provides `norm` and `abnorm` functions.

```matlab
% Generate training data
x = rand(100, 1) * 4 - 2;        % x ∈ [-2, 2]
y = x.^2;

% Normalize input
x_min = min(x);
x_max = max(x);
x_norm = apce().norm(x, x_min, x_max);
```

### 2. Fitting an aPCE Model
Use `fit_auto` to automatically select an appropriate order (maximum allowed order = 10, tolerance = 1e-3):

```matlab
[degree, coeffs, c_matrix] = apce().fit_auto(10, x_norm, y, 1e-3);
```
- `degree`: Final expansion order used
- `coeffs`: PCE coefficient matrix (each column corresponds to one output)
- `c_matrix`: aPCE-specific basis coefficient matrix (required for prediction)

### 3. Prediction
Predict at new points (normalize new points first):

```matlab
x_new = [-1.5; 0; 1.8];
x_new_norm = apce().norm(x_new, x_min, x_max);
y_pred_norm = apce().pred(degree, coeffs, c_matrix, x_new_norm);
y_pred = apce().abnorm(y_pred_norm, min(y), max(y));   % Denormalize output
```

### 4. Automatic Order Selection
`fit_auto` increments the order from 0 upwards until the fitting error falls below `tol_allow` or `degree_max` is reached. It returns a positive order on success; on failure, it returns a negative value (its absolute value is the maximum order attempted).

---

## Function Documentation

### pce.m Functions

All functions are accessed via the struct returned by `pce()`, e.g., `pce().fit_auto(...)`.

| Function | Description |
|----------|-------------|
| `fit_auto(degree_max, x, y, tol_allow)` | Automatically selects the order and fits the model; returns `[fitted_deg, coeffs]` |
| `pred(degree, coeffs, x_new)` | Predicts using a trained PCE model |
| `norm(x, x_min, x_max)` | Linearly maps data to `[-1, 1]`; if `x_min, x_max` are omitted, uses min/max of the data |
| `abnorm(x, x_min, x_max)` | Maps `[-1, 1]` data back to the original interval |
| `legendre(n, x)` | Evaluates the n-th order Legendre polynomial at x (normalization handled internally) |
| `fit(degree, x, y)` | Fits a PCE model of specified order; returns `[coeffs, phi]` |
| `tol(degree, coeffs, phi, x, y)` | Computes the relative L2 error |
| `eval(degree, coeffs, phi, x)` | Predicts using an existing phi matrix (optimized for efficiency) |
| `test()` | Test function; displays an informational message only |

### apce.m Functions

Accessed via the `apce()` struct.

| Function | Description |
|----------|-------------|
| `fit_auto(degree_max, x, y, tol_allow)` | Automatically selects the order and fits the model; returns `[fitted_deg, coeffs, c_matrix]` |
| `pred(degree, coeffs, c_matrix, x_new)` | Predicts using a trained aPCE model |
| `norm(x, x_min, x_max)` | Same as in pce |
| `abnorm(x, x_min, x_max)` | Same as in pce |
| `fit(degree, x, y)` | Fits an aPCE model of specified order; returns `[coeffs, phi, c_matrix]` |
| `tol(degree, coeffs, c_matrix, phi, x, y)` | Computes the relative L2 error |
| `eval(degree, coeffs, c_matrix, phi, x)` | Predicts using an existing phi matrix (optimized for efficiency) |
| `test()` | Test function |

---

## Important Notes

1. **It is strongly recommended to normalize input data** (map to `[-1, 1]`). Failure to do so may cause numerical overflow or significant loss of precision, especially for higher-order terms. The toolbox provides `norm` / `abnorm` helper functions.
2. The aPCE basis functions are constructed from the first and second moments of the input data. Therefore, the training data should be representative of the true input distribution.
3. When the input dimension or expansion order is high, the number of basis functions grows rapidly (combinatorial explosion). Set the maximum order judiciously.
4. The `tol_allow` parameter in `fit_auto` is a threshold for the relative L2 error; adjust it according to your problem requirements.
5. The code uses `lsqminnorm` to solve the least-squares problem. If warnings appear, consider using `lscov` instead (a commented alternative is provided).

---

## Contributing

Contributions via Issues or Pull Requests for improvements, bug fixes, or new features are welcome. Please ensure code style is consistent with the existing files and add necessary comments.

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Author

Jensentsts Wang  
For questions, please contact: jensentsts.wang@example.com (replace with actual email)

---

*Last updated: 2026-03-18*
