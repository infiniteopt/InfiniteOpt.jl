# # Quadcopter Trajectory and Set Point Tracking

# The quadcopter is an unmanned aerial vehicle with 4 propellors and in this case study,
# we seek to determine an optimal control policy for the 
# 4 control inputs that must be adjusted over a specified time period to enable the 
# quadcopter to closely follow given set points/trajectory while minimizing the difference between 
# the set point and actual position and propellor input

# ## Background

# Modelling the quadcopter trajectory is an infinite dimensional optimization problem.
# The variables are dependent on time. 
# Using a sinusoidal setpoint trajectory and 5 minute time horizon, the problem formulation is as follows:

# ## Formulation

# ```math 
#\begin{aligned}
#	&&\underset{x(t), u(t)}{\text{min}} &&& \int_{0}^{T} \left( (x_1 - \bar{x}_1)^2 + (x_3 - \bar{x}_3)^2 + (x_5 - \bar{x}_5)^2 + x_7^2 + x_8^2 + x_9^2 + 0.1 \sum_{i=1}^{4} u_i^2 \right) dt  \\
#	&&\text{s.t.} &&& x_j(0) = 0, \quad j \in \{1,2,3,\dots,9\} \\
#	&&&&& \frac{\partial x_j}{\partial t} = x_{j+1}, \quad j \in \{1,3,5\} \\
#	&&&&& \frac{\partial x_2}{\partial t} = u_1 \cos(x_7) \sin(x_8) \cos(x_9) + u_1 \sin(x_7) \sin(x_9) \\
#	&&&&& \frac{\partial x_4}{\partial t} = u_1 \cos(x_7) \sin(x_8) \cos(x_9) - u_1 \sin(x_7) \sin(x_9) \\
#	&&&&& \frac{\partial x_6}{\partial t} = u_1 \cos(x_7) \cos(x_8) - 9.8 \\
#	&&&&& \frac{\partial x_7}{\partial t} = u_2 \frac{\cos(x_7)}{\cos(x_8)} + u_3 \frac{\sin(x_7)}{\cos(x_8)} \\
#	&&&&& \frac{\partial x_8}{\partial t} = - u_2 \sin(x_7) + u_3 \cos(x_7) \\
#	&&&&& \frac{\partial x_9}{\partial t} = u_2 \cos(x_7) \tan(x_8) + u_3 \sin(x_7) \tan(x_8) + u_4 \\
#	&&&&& \bar{x}_1 = \sin\left(\frac{2\pi t}{T}\right), \quad \bar{x}_3 = \sin\left(\frac{2\pi t}{T}\right), \quad \text{and} \quad \bar{x}_5 = \frac{2t}{T}
#\end{aligned}
# ```

# ## Model Definition

# Import the required packages

using InfiniteOpt, Ipopt, Plots

# Specify the time horizon and control indices

time_horizon = 60 * 5 # 5 Minute Time Horizon

I = 1:4 # control indices

# Initialize the infinite model

model = InfiniteModel(optimizer_with_attributes(Ipopt.Optimizer)); 

# Specify infinite parameter

@infinite_parameter(model, t ∈ [0, 60 * 5], num_supports = (60 * 4))

# Specify decision variables

@variables(model, begin
    x[1:9], Infinite(t)
    u[1:4], Infinite(t), (start = 0)
end)


# Specify parameter functions
## t needs to be explicitly defined as an infinite parameter because time horizon introduces ambiguity
@parameter_function(model, x1bar == t -> sin(2 * pi * t/time_horizon))
@parameter_function(model, x3bar == t -> 2*sin(4 * pi * t/time_horizon))
@parameter_function(model, x5bar == t -> 2 * t/time_horizon)

# Specify the objective

@objective(model, Min, ∫((x[1]-x1bar)^2 + (x[3]-x3bar)^2 + (x[5]-x5bar)^2 + x[7]^2 + x[8]^2 + x[9]^2+ 0.1 * sum(u[i]^2 for i in I) , t))

# Set the initial condition for the states

@constraint(model, [j = 1:9], x[j](0) == 0)

# Define the point physics ODEs
## Note that states 1, 3, 5 correspond to the x,y, and z spatial positions
@constraint(model, [i=[1,3,5]],  ∂(x[i], t) == x[i+1]) 
@constraint(model, ∂(x[2], t) == u[1]*cos(x[7])*sin(x[8])*cos(x[9])+u[1]*sin(x[7])*sin(x[9]))
@constraint(model, ∂(x[4], t) == u[1]*cos(x[7])*sin(x[8])*sin(x[9])-u[1]*sin(x[7])*cos(x[9]))
@constraint(model, ∂(x[6], t) == u[1]*cos(x[7])*cos(x[8])-9.8)
@constraint(model, ∂(x[7], t) == u[2]*(cos(x[7])/cos(x[8]))+u[3]*(sin(x[7])/cos(x[8])))
@constraint(model, ∂(x[8], t) == -u[2]*sin(x[7])+u[3]*cos(x[7]))
@constraint(model, ∂(x[9], t) == u[2]*cos(x[7])*tan(x[8])+u[3]*sin(x[7])*tan(x[8])+u[4])

# ## Problem Solution

# Solve the optimization problem and check the status of the solution
optimize!(model)

if is_solved_and_feasible(model)
    println("Optimal solution found!")
    
    ## Extract solution

    t_values = value.(t)  
    x_values = value.(x)            
    u_values = value.(u)          

    ## Plot the Quadcopter Trajectory with respect to time (3D plot)

    x = x_values[1]

    xset = sin.(2 * pi.* t_values./time_horizon)

    y = x_values[3]

    yset = 2 * sin.(4 * pi.* t_values./time_horizon)

    z = x_values[5]

    zset = 2 * t_values./time_horizon

    traj_plot = scatter(x,y,z, label = "Quadcopter Trajectory")

    plot!(xset, yset, zset, label = "Setpoint", lw = 3, ls=:dot)

    savefig("Quadcopter Trajectory.png")


    ## Plot the 4 control states with respect to time

    l = @layout [a ; b ; c ; d]
    p1 = plot(t_values, u_values[1], label = "Propeller 1", ylabel = "u1(t)")
    p2 = plot(t_values, u_values[2], label = "Propeller 2", ylabel = "u2(t)")
    p3 = plot(t_values, u_values[3], label = "Propeller 3", ylabel = "u3(t)")
    p4 = plot(t_values, u_values[4], label = "Propeller 4", ylabel = "u4(t)")
    xlabel!("Time (t)")

    plot(p1, p2, p3, p4, layout = l)
    plot!(size=(800,600))

    savefig("Control Input (U).png")

else
    println("Optimization did not converge: ", termination_status(model))
end





