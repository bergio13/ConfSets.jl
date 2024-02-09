using ConfSets, Test, Statistics

function lower_upper(x, conf_set_fn, alpha)
    lower = conf_set_fn(x, alpha)["l"]
    upper = conf_set_fn(x, alpha)["u"]

    return lower, upper
end

function is_decreasing_on_average(vector)
    # Calculate the differences between consecutive elements
    diffs = diff(vector)

    # Calculate the average difference
    avg_diff = mean(diffs)

    # Check if the average difference is negative
    is_decreasing = avg_diff < 0

    return is_decreasing
end
# Test the confidence intervals for the mean of a vector x, using the Central Limit Theorem.
x = randn(1000)
true_mean = 0

l, u = lower_upper(x, clt_confidence_interval, 0.05)
@test l < cumul_mean(x) < u

diff1 = abs.(l .- true_mean)
diff2 = abs.(u .- true_mean)
@test is_decreasing_on_average(diff1) == true
@test is_decreasing_on_average(diff2) == true

l, u = lower_upper(x, clt_confidence_interval, 0.01)
@test l < cumul_mean(x) < u

diff1 = abs.(l .- true_mean)
diff2 = abs.(u .- true_mean)
@test is_decreasing_on_average(diff1) == true
@test is_decreasing_on_average(diff2) == true

@test clt_confidence_interval(x, 0.05) == clt_confidence_interval(x, 0.05)


# Test the confidence intervals for the mean of a vector x, using the Hoeffding inequality.
a, b = 1, 10
x = rand(a:b, 1000)

true_mean = (a + b) / 2

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


diff1 = abs.(hoeff_confidence_interval(x, 0.05, range)["l"] .- true_mean)
diff2 = abs.(hoeff_confidence_interval(x, 0.05, range)["u"] .- true_mean)
@test is_decreasing_on_average(diff1) == true
@test is_decreasing_on_average(diff2) == true


# Test the confidence intervals for the mean of a vector x, using the Chebyshev inequality.
x = randn(1000)
true_mean = 0

l, u = lower_upper(x, cheb_confidence_interval, 0.1)
@test l < cumul_mean(x) < u
diff1 = abs.(l .- true_mean)
diff2 = abs.(u .- true_mean)
@test is_decreasing_on_average(diff1) == true
@test is_decreasing_on_average(diff2) == true

l, u = lower_upper(x, cheb_confidence_interval, 0.05)
@test l < cumul_mean(x) < u

diff1 = abs.(l .- true_mean)
diff2 = abs.(u .- true_mean)
@test is_decreasing_on_average(diff1) == true
@test is_decreasing_on_average(diff2) == true

@test cheb_confidence_interval(x, 0.05) == cheb_confidence_interval(x, 0.05)
@test cheb_confidence_interval(x, 0.01) == cheb_confidence_interval(x, 0.01)


# Test the confidence sets for the mean of a vector x, using the Asymptotic Confidence Sequence.
x = randn(1000)
true_mean = 0

l, u = lower_upper(x, Asymp_conf_seq, 0.05)
@test l <= cumul_mean(x) <= u
diff1 = abs.(l .- true_mean)
diff2 = abs.(u .- true_mean)
@test is_decreasing_on_average(diff1) == true
@test is_decreasing_on_average(diff2) == true

l, u = lower_upper(x, Asymp_conf_seq, 0.01)
@test l <= cumul_mean(x) <= u
diff1 = abs.(l .- true_mean)
diff2 = abs.(u .- true_mean)
@test is_decreasing_on_average(diff1) == true
@test is_decreasing_on_average(diff2) == true

@test Asymp_conf_seq(x, 0.05) == Asymp_conf_seq(x, 0.05)
@test Asymp_conf_seq(x, 0.01) == Asymp_conf_seq(x, 0.01)