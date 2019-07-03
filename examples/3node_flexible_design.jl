# include("C:/Users/puls446/.julia/dev/InfOpt/examples/3node_flexible_design.jl")

using Revise, InfOpt, JuMP, Clp, Distributions

# Set the covariance matrix for the uncertain parameters
θ_nom = [0.; 60.; 10.]
covar = [80. 0 0; 0 80. 0; 0 0 120.]

# Set the dimensions
n_z = 3
n_θ = 3
n_d = 3

# Set the problem parameters
c = ones(n_d) / sqrt(n_d)
c_max = 5
U = 10000
num_samples = 10000

# Initialize the model
m = InfiniteModel(with_optimizer(Clp.Optimizer))

# Set the uncertainty parameters
dist = MvNormal(θ_nom, covar)
θs = rand(dist, num_samples)
@infinite_parameter(m, θ[i = 1:n_θ] in dist, supports = θs[i, :])

# Initialize the variables
@infinite_variable(m, 0 <= y(θ) <= 1)
# @infinite_variable(m, y(θ), Bin)
@infinite_variable(m, z[1:n_z](θ))
@global_variable(m, d[1:n_d] >= 0)

# Set objective function
expect_data = DiscreteMeasureData(θ, ones(num_samples) / num_samples, supports(θ), name = "expect")
@objective(m, Max, measure(1 - y, expect_data))

# Set the line capacity constraints
@constraint(m, f1, -z[1] - 35 - d[1] <= y * U)
@constraint(m, f2, z[1] - 35 - d[1] <= y * U)
@constraint(m, f3, -z[2] - 50 - d[2] <= y * U)
@constraint(m, f4, z[1] - 50 - d[2] <= y * U)

# Set the generator capacity constraints
@constraint(m, f5, -z[3] <= y * U)
@constraint(m, f6, z[3] - 100 - d[3] <= y * U)

# Set the node balance constraints
@constraint(m, h1, z[1] - θ[1] == 0)
@constraint(m, h2, -z[1] -z[2] + z[3] - θ[2] == 0)
@constraint(m, h3, z[2] - θ[3] == 0)

# Enforce the minimum SF
@constraint(m, max_cost, sum(c[i] * d[i] for i = 1:n_d) <= c_max)

# Solve and and obtain results
optimize!(m)
if has_values(m)
    opt_y = value(y)
    opt_d = value.(d)
    opt_obj = objective_value(m)
end

# Estimate the value of SF
SF = 1 - sum(opt_y .>= 1e-8) / num_samples

# Print the results
print("------------------RESULTS------------------\n")
print("Optimal Objective:     ", opt_obj, "\n")
print("Optimal Cost:          ", sum(c[i] * opt_d[i] for i = 1:n_d), "\n")
print("Maximum Cost:          ", c_max, "\n")
print("Predicted SF:          ", 100 * SF, "%\n")
print("Optimal Design Values: ", opt_d, "\n\n")
