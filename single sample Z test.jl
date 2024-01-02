#based on:
#Handbook of Parametric and Nonparametric Statistical Procedures, 5th edition, Chapter 1
using Random;
Random.seed!(783376894512);
using Pkg;
Pkg.activate("."); Pkg.instantiate()
using Distributions
using Plots

n_sim = 1_000_000

n = 5
σₜ = 2.0 # assumed known
μₜ = 3.0
Mₜ = Product([Normal(μₜ, σₜ) for _ in 1:n]) # assumed known to be Normal and i.i.d.
generate_data(Mₜ) = rand(Mₜ)

function test_statistic(y, μ₀, σₜ)
    n = length(y)
    μ̂ = sum(y) / n
    return (μ̂ - μ₀) / (σₜ / √n)
end

# Null Hypothesis is true

μ₀ = 3.0
μₜ == μ₀

simulated_statistics = [test_statistic(generate_data(Mₜ), μ₀, σₜ) for _ in 1:n_sim]
distribution_under_H₀ = z -> pdf(Normal(), z)

histogram(simulated_statistics; normalize=true, label="Simulated")
zmin, zmax = minimum(simulated_statistics), maximum(simulated_statistics)
z_plot = zmin:((zmax - zmin) / 100):zmax
plot!(distribution_under_H₀, z_plot; lw=5, label="Theoretical")
plot!(;
    title="Distribution of single sample Z test statistic, \n when Null Hypothesis is correct",
    xlimits=(zmin, zmax),
    xlabel="z",
    ylabel="pdf",
    grid=false,
    legendfontsize=12,
    tickfontsize=12,
    guidefontsize=16,
)

# Null Hypothesis is false

μ₀ = 5.0
μₜ != μ₀

simulated_statistics = [test_statistic(generate_data(Mₜ), μ₀, σₜ) for _ in 1:n_sim]
distribution_under_H₁ = znc -> pdf(Normal((μₜ - μ₀) / (σₜ / √n)), znc)

histogram(simulated_statistics; normalize=true, label="Simulated")
zmin, zmax = minimum(simulated_statistics), maximum(simulated_statistics)
z_plot = zmin:((zmax - zmin) / 100):zmax
plot!(distribution_under_H₁, z_plot; lw=5, label="Theoretical")
plot!(;
    title="Distribution of single sample Z test statistic, \n when Null Hypothesis is incorrect",
    xlimits=(zmin, zmax),
    xlabel="z",
    ylabel="pdf",
    grid=false,
    legendfontsize=12,
    tickfontsize=12,
    guidefontsize=16,
)
