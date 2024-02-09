using ConfSets, Test, Statistics, Random


function lower_upper(x, conf_set_fn, alpha, sequential)
    lower = conf_set_fn(x, alpha, sequential)["l"]
    upper = conf_set_fn(x, alpha, sequential)["u"]

    return lower, upper
end

# Function to check if a vector is decreasing on average
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
Random.seed!(1234)
x = randn(1000)
true_mean = 0

@test isapprox(clt_confidence_interval(x, 0.1, false)["u"], clt_confidence_interval(x, 0.1, true)["u"][end])
@test isapprox(clt_confidence_interval(x, 0.1, false)["l"], clt_confidence_interval(x, 0.1, true)["l"][end])

l, u = lower_upper(x, clt_confidence_interval, 0.05, false)
@test l <= mean(x) <= u

l, u = lower_upper(x, clt_confidence_interval, 0.05, true)
@test l <= cumul_mean(x) <= u


diff1 = abs.(l .- true_mean)
diff2 = abs.(u .- true_mean)
@test is_decreasing_on_average(diff1) == true
@test is_decreasing_on_average(diff2) == true

l, u = lower_upper(x, clt_confidence_interval, 0.01, true)
@test l < cumul_mean(x) < u

diff1 = abs.(l .- true_mean)
diff2 = abs.(u .- true_mean)
@test is_decreasing_on_average(diff1) == true
@test is_decreasing_on_average(diff2) == true

@test clt_confidence_interval(x, 0.05, true) == clt_confidence_interval(x, 0.05, true)


# Test the confidence intervals for the mean of a vector x, using the Hoeffding inequality.
Random.seed!(1234)
a, b = 1, 10
x = rand(a:b, 1000)

true_mean = (a + b) / 2

range = b - a
float_range = float(range)

@test hoeff_confidence_interval(x, 0.1, range, true) == hoeff_confidence_interval(x, 0.1, float_range, true)
@test hoeff_confidence_interval(x, 0.05, range, true) == hoeff_confidence_interval(x, 0.05, float_range, true)
@test hoeff_confidence_interval(x, 0.01, range, true) == hoeff_confidence_interval(x, 0.01, float_range, true)

lower = hoeff_confidence_interval(x, 0.1, range, true)["l"]
upper = hoeff_confidence_interval(x, 0.1, range, true)["u"]
@test lower < cumul_mean(x) < upper

lower = hoeff_confidence_interval(x, 0.1, range, false)["l"]
upper = hoeff_confidence_interval(x, 0.1, range, false)["u"]
@test lower < mean(x) < upper

lower = hoeff_confidence_interval(x, 0.05, float_range, true)["l"]
upper = hoeff_confidence_interval(x, 0.05, float_range, true)["u"]

@test lower < cumul_mean(x) < upper

lower = hoeff_confidence_interval(x, 0.01, float_range, true)["l"]
upper = hoeff_confidence_interval(x, 0.01, float_range, true)["u"]

@test lower < cumul_mean(x) < upper


diff1 = abs.(hoeff_confidence_interval(x, 0.05, range, true)["l"] .- true_mean)
diff2 = abs.(hoeff_confidence_interval(x, 0.05, range, true)["u"] .- true_mean)
@test is_decreasing_on_average(diff1) == true
@test is_decreasing_on_average(diff2) == true

@test isapprox(hoeff_confidence_interval(x, 0.1, range, false)["u"], hoeff_confidence_interval(x, 0.1, range, true)["u"][end])
@test isapprox(hoeff_confidence_interval(x, 0.1, range, false)["l"], hoeff_confidence_interval(x, 0.1, range, true)["l"][end])


# Test the confidence intervals for the mean of a vector x, using the Chebyshev inequality.
Random.seed!(1234)
x = randn(1000)
true_mean = 0

l, u = lower_upper(x, cheb_confidence_interval, 0.1, true)
@test l < cumul_mean(x) < u
diff1 = abs.(l .- true_mean)
diff2 = abs.(u .- true_mean)
@test is_decreasing_on_average(diff1) == true
@test is_decreasing_on_average(diff2) == true

l, u = lower_upper(x, cheb_confidence_interval, 0.05, true)
@test l < cumul_mean(x) < u

diff1 = abs.(l .- true_mean)
diff2 = abs.(u .- true_mean)
@test is_decreasing_on_average(diff1) == true
@test is_decreasing_on_average(diff2) == true

@test cheb_confidence_interval(x, 0.05, true) == cheb_confidence_interval(x, 0.05, true)
@test cheb_confidence_interval(x, 0.01, true) == cheb_confidence_interval(x, 0.01, true)

@test isapprox(cheb_confidence_interval(x, 0.1, false)["u"], cheb_confidence_interval(x, 0.1, true)["u"][end])
@test isapprox(cheb_confidence_interval(x, 0.1, false)["l"], cheb_confidence_interval(x, 0.1, true)["l"][end])

l, u = lower_upper(x, cheb_confidence_interval, 0.05, false)
@test l <= mean(x) <= u


# Test the confidence sets for the mean of a vector x, using the Asymptotic Confidence Sequence.
Random.seed!(1234)
x = randn(1000)
true_mean = 0

l, u = lower_upper(x, Asymp_conf_seq, 0.05, true)
@test l <= cumul_mean(x) <= u
diff1 = abs.(l .- true_mean)
diff2 = abs.(u .- true_mean)
@test is_decreasing_on_average(diff1) == true
@test is_decreasing_on_average(diff2) == true

l, u = lower_upper(x, Asymp_conf_seq, 0.01, true)
@test l <= cumul_mean(x) <= u
diff1 = abs.(l .- true_mean)
diff2 = abs.(u .- true_mean)
@test is_decreasing_on_average(diff1) == true
@test is_decreasing_on_average(diff2) == true

@test Asymp_conf_seq(x, 0.05, true) == Asymp_conf_seq(x, 0.05, true)
@test Asymp_conf_seq(x, 0.01, true) == Asymp_conf_seq(x, 0.01, true)

@test isapprox(Asymp_conf_seq(x, 0.1, false)["u"], Asymp_conf_seq(x, 0.1, true)["u"][end])
@test isapprox(Asymp_conf_seq(x, 0.1, false)["l"], Asymp_conf_seq(x, 0.1, true)["l"][end])

l, u = lower_upper(x, Asymp_conf_seq, 0.05, false)
@test l <= mean(x) <= u