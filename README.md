# linearDiagnostic

The goal of linearDiagnostic is to provide a set of tools to diagnose regression models,
focused on `lm` and `glm` models.

## Installation

You can install the development version of linearDiagnostic from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Alvaro-Kothe/linearDiagnostic")
```

## Example

This is a basic example to evaluate if the model assumption of a `glm` model is valid.

``` r
library(linearDiagnostic)

# Fit the glm model
plot_envelope(glm_model)
```

This will create a simulated envelope plot to verify if the model assumption is valid.

## Contribution and Bug Reports

Contributions to enhance the functionalities and usability of **LinearDiagnostic** are welcome.
If you encounter any bugs or issues while using the package, please file an [issue](https://github.com/Alvaro-Kothe/linearDiagnostic/issues).

## License

This package is open-source and distributed under the [MIT License](https://opensource.org/licenses/MIT). You are free to use, modify, and distribute the package, following the terms of the license.
