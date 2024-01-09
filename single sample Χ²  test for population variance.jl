#based on:
#Handbook of Parametric and Nonparametric Statistical Procedures, 5th edition, Chapter 3
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

function test_statistic(y, σ₀)
    n = length(y)
    μ̂ = sum(y) / n
    σ̂² = sum((y .- μ̂) .^ 2) / (n - 1)
    return σ̂² / (σ₀^2 / (n - 1))
end

# Null Hypothesis is true

σ₀ = 2.0
σₜ == σ₀

simulated_statistics = [test_statistic(generate_data(Mₜ), σ₀) for _ in 1:n_sim]
distribution_under_H₀ = χ² -> pdf(Chisq(n - 1), χ²)

histogram(simulated_statistics; normalize=true, label="Simulated")
χmin, χmax = minimum(simulated_statistics), maximum(simulated_statistics)
χ_plot = χmin:((χmax - χmin) / 100):χmax
plot!(distribution_under_H₀, χ_plot; lw=5, label="Theoretical")
plot!(;
    title="Distribution of Χ² test statistic for a pop. var., \n when Null Hypothesis is correct",
    xlimits=(χmin, χmax),
    xlabel="χ",
    ylabel="pdf",
    grid=false,
    legendfontsize=12,
    tickfontsize=12,
    guidefontsize=16,
)

# Null Hypothesis is false

σ₀ = 3.0
σₜ != σ₀

simulated_statistics = [test_statistic(generate_data(Mₜ), σ₀) for _ in 1:n_sim]
distribution_under_H₁ = χnc -> pdf(Chisq(n - 1) * σₜ^2 / σ₀^2, χnc)

histogram(simulated_statistics; normalize=true, label="Simulated")
χmin, χmax = minimum(simulated_statistics), maximum(simulated_statistics)
χ_plot = χmin:((χmax - χmin) / 100):χmax
plot!(distribution_under_H₁, χ_plot; lw=5, label="Theoretical")
plot!(;
    title="Distribution of Χ² test statistic for a pop. var., \n when Null Hypothesis is incorrect",
    xlimits=(χmin, χmax),
    xlabel="χ",
    ylabel="pdf",
    grid=false,
    legendfontsize=12,
    tickfontsize=12,
    guidefontsize=16,
)

# Confidence Interval

function confidence_interval(y, α)
    n = length(y)
    μ̂ = sum(y) / n
    σ̂² = sum((y .- μ̂) .^ 2) / (n - 1)
    L = (n - 1) * σ̂² / quantile(Chisq(n - 1), 1 - α / 2)
    U = (n - 1) * σ̂² / quantile(Chisq(n - 1), α / 2)
    return Interval(L, U)
end

α = 0.05

simulated_CIs = [confidence_interval(generate_data(Mₜ), α) for _ in 1:n_sim]
coverage = count([σₜ^2 in CI for CI in simulated_CIs]) / n_sim

bar(["inside CI", "outside CI"], [coverage, 1 - coverage])
plot!(;
    title="Simulated Coverage of 95% confidence interval \n corresponding to Χ² test statistic for a pop. var.",
    ylimits=(0.0, 1.0),
    yticks=[0.05, 0.95],
    ylabel="proportion",
    grid=false,
    legend=false,
    legendfontsize=12,
    tickfontsize=12,
    guidefontsize=16,
)
