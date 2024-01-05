#based on:
#Handbook of Parametric and Nonparametric Statistical Procedures, 5th edition, Chapter 2
using Random;
Random.seed!(783376894512);
using Pkg;
Pkg.activate(".")
Pkg.instantiate()
using Distributions
using Intervals
using Plots

n_sim = 1_000_000

n = 10
σₜ = 2.0
μₜ = 3.0
Mₜ = Product([Normal(μₜ, σₜ) for _ in 1:n]) # assumed known to be Normal and i.i.d.
generate_data(Mₜ) = rand(Mₜ)

function test_statistic(y, μ₀)
    n = length(y)
    μ̂ = sum(y) / n
    σ̂² = sum((y .- μ̂) .^ 2) / (n - 1)
    return (μ̂ - μ₀) / √(σ̂² / n)
end

# Null Hypothesis is true

μ₀ = 3.0
μₜ == μ₀

simulated_statistics = [test_statistic(generate_data(Mₜ), μ₀) for _ in 1:n_sim]
distribution_under_H₀ = t -> pdf(TDist(n - 1), t)

histogram(simulated_statistics; normalize=true, label="Simulated")
tmin, tmax = minimum(simulated_statistics), maximum(simulated_statistics)
t_plot = tmin:((tmax - tmin) / 100):tmax
plot!(distribution_under_H₀, t_plot; lw=5, label="Theoretical")
plot!(;
    title="Distribution of single sample T test statistic, \n when Null Hypothesis is correct",
    xlimits=(tmin, tmax),
    xlabel="t",
    ylabel="pdf",
    grid=false,
    legendfontsize=12,
    tickfontsize=12,
    guidefontsize=16,
)

# Null Hypothesis is false

μ₀ = 5.0
μₜ != μ₀

simulated_statistics = [test_statistic(generate_data(Mₜ), μ₀) for _ in 1:n_sim]
distribution_under_H₁ = tnc -> pdf(NoncentralT(n - 1, (μₜ - μ₀) / √(σₜ^2 / n)), tnc)

histogram(simulated_statistics; normalize=true, label="Simulated")
tmin, tmax = minimum(simulated_statistics), maximum(simulated_statistics)
t_plot = tmin:((tmax - tmin) / 100):tmax
plot!(distribution_under_H₁, t_plot; lw=5, label="Theoretical")
plot!(;
    title="Distribution of single sample T test statistic, \n when Null Hypothesis is incorrect",
    xlimits=(tmin, tmax),
    xlabel="t",
    ylabel="pdf",
    grid=false,
    legendfontsize=12,
    tickfontsize=12,
    guidefontsize=16,
)
