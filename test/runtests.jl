using ConfSets, Test, Statistics, Random

function lower_upper(x, method, alpha, sequential)
    if method == "clt"
        lower = confint(x, alpha, "clt", sequential=sequential).l
        upper = confint(x, alpha, "clt", sequential=sequential).u
    elseif method == "hoeffding"
        lower = confint(x, alpha, "hoeffding", 10, sequential).l
        upper = confint(x, alpha, "hoeffding", 10, sequential).u
    elseif method == "chebyshev"
        lower = confint(x, alpha, "chebyshev", sequential=sequential).l
        upper = confint(x, alpha, "chebyshev", sequential=sequential).u
    elseif method == "asymptotic_cs"
        lower = confseq(x, alpha, "asymptotic", sequential=sequential).l
        upper = confseq(x, alpha, "asymptotic", sequential=sequential).u
    end
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

#-----------------------------------------------------------------------------------------------------------------
# Test the confidence intervals for the mean of a vector x, using the Central Limit Theorem.
Random.seed!(1234)
x = randn(1000)
true_mean = 0

@test isapprox(confint(x, 0.1, "clt", sequential=false).u, confint(x, 0.1, "clt", sequential=true).u[end])
@test isapprox(confint(x, 0.1, "clt", sequential=false).l, confint(x, 0.1, "clt", sequential=true).l[end])

l, u = lower_upper(x, "clt", 0.05, false)
@test l <= mean(x) <= u

l, u = lower_upper(x, "clt", 0.05, true)
@test l <= cumul_mean(x) <= u


diff1 = abs.(l .- true_mean)
diff2 = abs.(u .- true_mean)
@test is_decreasing_on_average(diff1) == true
@test is_decreasing_on_average(diff2) == true

l, u = lower_upper(x, "clt", 0.01, true)
@test l < cumul_mean(x) < u

diff1 = abs.(l .- true_mean)
diff2 = abs.(u .- true_mean)
@test is_decreasing_on_average(diff1) == true
@test is_decreasing_on_average(diff2) == true

@test confint(x, 0.05, "clt", sequential=true) == confint(x, 0.05, "clt", sequential=true)


# Test the confidence intervals for the mean of a vector x, using the Hoeffding inequality.
Random.seed!(1234)
a, b = 1, 10
x = rand(a:b, 1000)

true_mean = (a + b) / 2

range = b - a
float_range = float(range)

@test confint(x, 0.1, "hoeffding", rang=range, sequential=true) == confint(x, 0.1, "hoeffding", rang=float_range, sequential=true)
@test confint(x, 0.05, "hoeffding", rang=range, sequential=true) == confint(x, 0.05, "hoeffding", rang=float_range, sequential=true)
@test confint(x, 0.01, "hoeffding", rang=range, sequential=true) == confint(x, 0.01, "hoeffding", rang=float_range, sequential=true)

lower = confint(x, 0.1, "hoeffding", rang=range, sequential=true).l
upper = confint(x, 0.1, "hoeffding", rang=range, sequential=true).u
@test lower < cumul_mean(x) < upper

lower = confint(x, 0.1, "hoeffding", rang=range, sequential=false).l
upper = confint(x, 0.1, "hoeffding", rang=range, sequential=false).u
@test lower < mean(x) < upper

lower = confint(x, 0.05, "hoeffding", rang=float_range, sequential=true).l
upper = confint(x, 0.05, "hoeffding", rang=float_range, sequential=true).u
@test lower < cumul_mean(x) < upper

lower = confint(x, 0.01, "hoeffding", rang=float_range, sequential=true).l
upper = confint(x, 0.01, "hoeffding", rang=float_range, sequential=true).u
@test lower < cumul_mean(x) < upper


diff1 = abs.(confint(x, 0.05, "hoeffding", rang=range, sequential=true).l .- true_mean)
diff2 = abs.(confint(x, 0.05, "hoeffding", rang=range, sequential=true).u .- true_mean)
@test is_decreasing_on_average(diff1) == true
@test is_decreasing_on_average(diff2) == true

@test isapprox(confint(x, 0.1, "hoeffding", rang=range, sequential=false).u, confint(x, 0.1, "hoeffding", rang=range, sequential=true).u[end])
@test isapprox(confint(x, 0.1, "hoeffding", rang=range, sequential=false).l, confint(x, 0.1, "hoeffding", rang=range, sequential=true).l[end])


# Test the confidence intervals for the mean of a vector x, using the Chebyshev inequality.
Random.seed!(1234)
x = randn(1000)
true_mean = 0

l, u = lower_upper(x, "chebyshev", 0.1, true)
@test l < cumul_mean(x) < u
diff1 = abs.(l .- true_mean)
diff2 = abs.(u .- true_mean)
@test is_decreasing_on_average(diff1) == true
@test is_decreasing_on_average(diff2) == true

l, u = lower_upper(x, "chebyshev", 0.05, true)
@test l < cumul_mean(x) < u

diff1 = abs.(l .- true_mean)
diff2 = abs.(u .- true_mean)
@test is_decreasing_on_average(diff1) == true
@test is_decreasing_on_average(diff2) == true

@test confint(x, 0.05, "chebyshev", sequential=true) == confint(x, 0.05, "chebyshev", sequential=true)
@test confint(x, 0.01, "chebyshev", sequential=true) == confint(x, 0.01, "chebyshev", sequential=true)

@test isapprox(confint(x, 0.1, "chebyshev", sequential=false).u, confint(x, 0.1, "chebyshev", sequential=true).u[end])
@test isapprox(confint(x, 0.1, "chebyshev", sequential=false).l, confint(x, 0.1, "chebyshev", sequential=true).l[end])

l, u = lower_upper(x, "chebyshev", 0.05, false)
@test l <= mean(x) <= u


# Test the confidence sets for the mean of a vector x, using the Asymptotic Confidence Sequence.
Random.seed!(1234)
x = randn(1000)
true_mean = 0

l, u = lower_upper(x, "asymptotic_cs", 0.05, true)
@test l <= cumul_mean(x) <= u
diff1 = abs.(l .- true_mean)
diff2 = abs.(u .- true_mean)
@test is_decreasing_on_average(diff1) == true
@test is_decreasing_on_average(diff2) == true

l, u = lower_upper(x, "asymptotic_cs", 0.01, true)
@test l <= cumul_mean(x) <= u
diff1 = abs.(l .- true_mean)
diff2 = abs.(u .- true_mean)
@test is_decreasing_on_average(diff1) == true
@test is_decreasing_on_average(diff2) == true

@test confseq(x, 0.05, "asymptotic", sequential=true) == confseq(x, 0.05, "asymptotic", sequential=true)
@test confseq(x, 0.01, "asymptotic", sequential=true) == confseq(x, 0.01, "asymptotic", sequential=true)

@test isapprox(confseq(x, 0.1, "asymptotic", sequential=false).u, confseq(x, 0.1, "asymptotic", sequential=true).u[end])
@test isapprox(confseq(x, 0.1, "asymptotic", sequential=false).l, confseq(x, 0.1, "asymptotic", sequential=true).l[end])

l, u = lower_upper(x, "asymptotic_cs", 0.05, false)
@test l <= mean(x) <= u

#---------------------------------------------------------------------------------------------------------

@test_throws ArgumentError confint(x, 0.05, "cheb", sequential=true)
@test_throws ArgumentError confint(x, 0.05, "ard", sequential=true)