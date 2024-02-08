using ConfSets, Test

# Test the confidence intervals for the mean of a vector x, using the Central Limit Theorem.
x = randn(100)

lower = clt_confidence_interval(x, 0.05)["l"]
upper = clt_confidence_interval(x, 0.05)["u"]

@test lower < cumul_mean(x) < upper

lower = clt_confidence_interval(x, 0.01)["l"]
upper = clt_confidence_interval(x, 0.01)["u"]

@test lower < cumul_mean(x) < upper

@test clt_confidence_interval(x, 0.05) == clt_confidence_interval(x, 0.05)

# Test the confidence intervals for the mean of a vector x, using the Hoeffding inequality.
a, b = 1, 10
x = rand(a:b, 100)

range = b - a
float_range = float(range)

@test hoeff_confidence_interval(x, 0.1, range) == hoeff_confidence_interval(x, 0.1, float_range)
@test hoeff_confidence_interval(x, 0.05, range) == hoeff_confidence_interval(x, 0.05, float_range)
@test hoeff_confidence_interval(x, 0.01, range) == hoeff_confidence_interval(x, 0.01, float_range)

lower = hoeff_confidence_interval(x, 0.1, range)["l"]
upper = hoeff_confidence_interval(x, 0.1, range)["u"]

@test lower < cumul_mean(x) < upper

lower = hoeff_confidence_interval(x, 0.05, float_range)["l"]
upper = hoeff_confidence_interval(x, 0.05, float_range)["u"]

@test lower < cumul_mean(x) < upper

lower = hoeff_confidence_interval(x, 0.01, float_range)["l"]
upper = hoeff_confidence_interval(x, 0.01, float_range)["u"]

@test lower < cumul_mean(x) < upper

# Test the confidence intervals for the mean of a vector x, using the Chebyshev inequality.
x = randn(100)

lower = cheb_confidence_interval(x, 0.05)["l"]
upper = cheb_confidence_interval(x, 0.05)["u"]

@test lower < cumul_mean(x) < upper

lower = cheb_confidence_interval(x, 0.01)["l"]
upper = cheb_confidence_interval(x, 0.01)["u"]

@test lower < cumul_mean(x) < upper

@test cheb_confidence_interval(x, 0.05) == cheb_confidence_interval(x, 0.05)
@test cheb_confidence_interval(x, 0.01) == cheb_confidence_interval(x, 0.01)

# Test the confidence sets for the mean of a vector x, using the Asymptotic Confidence Sequence.
x = randn(100)

lower = Asymp_conf_seq(x, 0.05)["l"]
upper = Asymp_conf_seq(x, 0.05)["u"]

@test lower < cumul_mean(x) < upper

lower = Asymp_conf_seq(x, 0.01)["l"]
upper = Asymp_conf_seq(x, 0.01)["u"]

@test lower < cumul_mean(x) < upper

@test Asymp_conf_seq(x, 0.05) == Asymp_conf_seq(x, 0.05)
@test Asymp_conf_seq(x, 0.01) == Asymp_conf_seq(x, 0.01)