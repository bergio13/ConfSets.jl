module ConfSets

using Statistics
using Distributions

export cumul_mean, cumul_var, clt_std_margin, hoeff_margin, cheb_margin, ACS_margin, clt_confidence_interval, hoeff_confidence_interval, cheb_confidence_interval, Asymp_conf_seq

# This function computes the cumulative mean of a vector x, with regularizers for the observations and the mean.
function cumul_mean(x::Vector, regularizer_obs=0, regularizer_mean=1 / 2)
    t = 1:length(x)
    return ((cumsum(x) .+ regularizer_obs * regularizer_mean) ./ (t .+ regularizer_obs))
end

# This function computes the cumulative variance of a vector x
function cumul_var(x::Vector)
    t = 1:length(x)
    return ((cumul_mean(x .^ 2) .- cumul_mean(x) .^ 2) .* t ./ (t .- 1))
end


# Gaussian margin
function clt_std_margin(t::Vector, alpha::Float64)
    return quantile(Normal(), 1 - alpha / 2) ./ (t .^ 0.5)
end

# Hoeffding margin
function hoeff_margin(t::Vector, alpha::Float64, rang::Union{Float64,Int})
    return (((rang ./ (2 .* t)) .* log(2 / alpha)) .^ 0.5)
end

# Chebyshev margin
function cheb_margin(t::Vector, alpha::Float64)
    return (1 ./ ((t .* alpha) .^ 0.5))
end

# Asympt Conf Sequence margin
function ACS_margin(t::Vector, alpha::Float64)
    nom = ((log.(log.(2 .* t)) .+ 0.72 * log(10.4 / alpha))) .^ 0.5
    denom = t .^ 0.5
    return 1.7 .* nom ./ denom
end

# This function computes the confidence intervals for the mean of a vector x, using the Central Limit Theorem.
function clt_confidence_interval(x::Vector, alpha::Float64)
    t = collect(1:length(x))
    mu_hat = cumul_mean(x)
    std = cumul_var(x) .^ 0.5
    margin = clt_std_margin(t, alpha)
    std_margin = std .* margin

    std_margin[isnan.(std_margin)] .= Inf

    return Dict("l" => mu_hat .- std_margin, "u" => mu_hat .+ std_margin)

end

# This function computes the confidence intervals for the mean of a vector x, using the Hoeffding inequality.
function hoeff_confidence_interval(x::Vector, alpha::Float64, rang::Union{Float64,Int,Nothing})
    t = collect(1:length(x))
    mu_hat = cumul_mean(x)
    if isnothing(rang)
        rang = maximum(x) - minimum(x)
    end
    margin = hoeff_margin(t, alpha, rang)

    margin[isnan.(margin)] .= Inf

    return Dict("l" => mu_hat .- margin, "u" => mu_hat .+ margin)

end

# This function computes the confidence intervals for the mean of a vector x, using the Chebyshev inequality.
function cheb_confidence_interval(x::Vector, alpha::Float64)
    t = collect(1:length(x))
    mu_hat = cumul_mean(x)
    std = cumul_var(x) .^ 0.5
    margin = cheb_margin(t, alpha) .* std

    margin[isnan.(margin)] .= Inf

    return Dict("l" => mu_hat .- margin, "u" => mu_hat .+ margin)

end

# This function computes the confidence sets for the mean of a vector x, using the Asymptotic Confidence Sequence.
function Asymp_conf_seq(x::Vector, alpha::Float64)
    t = collect(1:length(x))
    mu_hat = cumul_mean(x)
    std = cumul_var(x) .^ 0.5

    margin = ACS_margin(t, alpha) .* std

    margin[isnan.(margin)] .= Inf

    return Dict("l" => mu_hat .- margin, "u" => mu_hat .+ margin)

end



end # module ConfSets
