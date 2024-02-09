[![Coverage Status](https://app.travis-ci.com/bergio13/ConfSets.jl.svg?branch=main)](https://app.travis-ci.com/bergio13/ConfSets.jl.svg?branch=main)
[![codecov](https://codecov.io/gh/bergio13/ConfSets.jl/graph/badge.svg?token=2WDDX6XTIH)](https://codecov.io/gh/bergio13/ConfSets.jl)

# ConfSets.jl

ConfSets.jl is a Julia library for computing confidence intervals using various statistical methods including the Central Limit Theorem (CLT), Hoeffding's Inequality, and Chebyshev's Inequality and confidence sequences. In particular, if you are in a sequential setting, it is advisable to use confidence sequences as over time the confidence intervals will miscover the true mean in multiple instances. On the other hand, the confidence sequences are designed to simultaneously capture the true mean uniformly over time and asymptotically.

## Installation

You can install ConfSets.jl using the Julia package manager. From the Julia REPL, enter the package manager mode by pressing `]`, and then run:

```
add ConfSets
```

## Usage

Once installed, you can use the library to compute confidence intervals and confidence sequences using the following methods. In partiular, if you want to compute sequential sets, you can set the `sequential` parameter to `true`, otherwise sei it to `false`.

```julia
using ConfSets

# Compute confidence interval using CLT
data = [10.2, 12.3, 9.8, 11.5, 10.9]
clt_interval = clt_confidence_interval(data, alpha=0.05, sequential=true)

# Compute confidence interval using Hoeffding's Inequality
hoeffding_interval = hoeff_confidence_interval(data, alpha=0.05, range=3, sequential=false) # range is the range of the data (e.g. 10-13, range=3)

# Compute confidence interval using Chebyshev's Inequality
chebyshev_interval = cheby_confidence_interval(data, alpha=0.05, sequential=true)

# Compute Asymptotic Confidence Sequence
asymptotic_sequence = Asymp_conf_seq(data, alpha=0.05, sequential=true)
```

## Contributing

Contributions are welcome! If you encounter any issues, have suggestions for improvements, or would like to add new features, feel free to open an issue or submit a pull request on GitHub.

## License

The code is available as open source under the terms of the [MIT License](https://opensource.org/licenses/MIT).
