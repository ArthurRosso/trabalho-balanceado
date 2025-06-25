import Pkg; Pkg.add("JuMP"); Pkg.add("GLPK"); Pkg.add("DelimitedFiles"); Pkg.add("Random");

using JuMP          # Modeling language for optimization
using GLPK          # Open-source solver for linear programming
using Printf        # For formatted printing of time
using Random
using Plots

# Start timing
start_time = time()

# Dados do problema
n = 27  # número de operações
m = 11  # número de trabalhadores

production_times = [
    1.313962 0.597045 1.009159 0.612462 0.449617 0.837132 0.772335 0.516805 0.578815 0.100000 0.685093 0.100000 0.436276 0.738441 0.118670 0.277524 0.493666 0.885146 0.203521 0.100000 0.135818 0.629099 0.680649 0.481319 0.180385 0.944709 0.747208
    1.125384 0.549770 0.746493 0.337872 0.001165 1.077978 0.538118 0.945970 0.635612 0.592962 0.907819 0.127168 1.033690 0.837369 0.100000 0.731133 0.610901 0.507580 0.524510 0.341730 0.373138 0.538946 0.662680 0.278364 0.622943 0.873428 0.344741
    0.974094 0.518205 0.773904 0.279570 0.000150 0.440662 0.339608 0.537872 0.962631 0.100000 0.358857 0.332600 0.100000 0.820566 0.529051 0.100000 0.571007 0.234868 0.843442 1.192587 0.640188 1.160046 0.984348 0.356233 0.836932 0.956617 0.983841
    1.572660 0.665359 0.856818 0.842129 0.150582 1.009583 0.484141 1.232283 0.640873 0.449888 1.265831 0.582249 0.460993 0.341157 0.339155 0.058748 0.680368 1.142364 0.100000 0.348297 0.333913 1.137994 0.987456 0.339435 0.284941 0.929942 0.707277
    0.898069 0.125949 0.467786 0.546594 0.142375 0.207885 0.210388 0.425543 0.237767 0.179150 0.094757 0.100000 0.068130 0.100000 0.100000 0.602273 0.514534 0.391336 0.182386 0.100000 0.546251 0.445065 0.748578 0.731170 0.564640 0.415168 0.100000
    1.379605 0.130188 0.486603 0.299697 0.100000 0.100000 0.100000 1.085317 0.762615 0.100000 0.068504 0.933243 0.100000 0.100000 0.336508 0.100000 1.275376 0.213094 0.674059 1.230953 0.100000 0.863346 1.075218 1.411747 1.376286 1.162985 1.114326
    1.109833 0.100000 0.524150 0.156218 0.238192 0.100000 1.252440 0.943218 1.579160 0.100000 0.734170 0.465704 0.204304 0.344996 0.100000 0.869373 0.946586 0.695260 0.831999 1.157619 0.539516 0.082603 0.206743 0.350447 0.309646 0.848950 0.091500
    0.845763 0.961070 0.511301 0.100000 0.100000 0.290111 1.075047 1.092602 0.824097 0.100000 0.521137 0.765702 0.010663 0.292768 0.187448 0.096697 0.692770 0.901463 0.660722 0.290358 0.965068 0.497296 1.450147 0.377958 0.314187 0.886801 0.984891
    1.512588 1.298165 0.763043 0.755011 0.100000 1.029159 0.345728 0.734917 1.009231 0.146111 0.077315 0.523648 0.100000 0.652456 0.515476 0.777121 0.100000 0.935739 0.620474 0.750382 0.538626 0.570449 0.540436 0.535995 1.090752 0.869258 0.121960
    1.091883 0.959550 0.772693 0.299005 0.279401 0.999538 0.704504 0.610980 0.452057 0.146066 0.564167 0.100000 0.102006 0.741028 0.105038 0.483031 0.282629 1.010945 0.749318 0.816155 0.227229 1.419271 0.820718 0.288646 0.608495 0.404634 0.493236
    1.370025 0.351032 0.739575 0.100000 0.821125 0.303201 0.581145 0.595318 0.333740 0.100000 0.460903 0.148374 0.424300 0.351381 0.014628 0.540198 0.364985 0.671686 0.869736 0.100000 0.100000 0.877114 1.269992 0.549942 0.227411 0.321557 0.226957
]

model_setup_time = @elapsed begin
    model = Model(GLPK.Optimizer)

    # Variables
    @variable(model, x[1:n, 1:m], Bin)  # Assignment of tasks to workers
    @variable(model, b[1:n, 1:m], Bin)  # Start of a block for a worker
    @variable(model, e[1:n, 1:m], Bin)  # End of a block for a worker
    @variable(model, T >= 0)  # Maximum workload

    @objective(model, Min, T)

    for i in 1:n
        @constraint(model, sum(x[i, j] for j in 1:m) == 1)
    end

    # Contiguity constraints: Each worker has exactly one contiguous block
    for j in 1:m
        # Exactly one start and one end per worker
        @constraint(model, sum(b[i, j] for i in 1:n) == 1)
        @constraint(model, sum(e[i, j] for i in 1:n) == 1)

        # Start position definition
        for i in 1:n
            # Handle boundary cases (i=1 and i=num_operations)
            prev_x = (i == 1) ? 0 : x[i-1, j]
            next_x = (i == n) ? 0 : x[i+1, j]

            # d[i,j] = 1 if x[i,j] = 1 and x[i-1,j] = 0
            @constraint(model, b[i, j] >= x[i, j] - prev_x)
            @constraint(model, b[i, j] <= x[i, j])
            @constraint(model, b[i, j] <= 1 - prev_x)

            # e[i,j] = 1 if x[i,j] = 1 and x[i+1,j] = 0
            @constraint(model, e[i, j] >= x[i, j] - next_x)
            @constraint(model, e[i, j] <= x[i, j])
            @constraint(model, e[i, j] <= 1 - next_x)
        end
    end

    # Workload constraints: T >= sum of production times for each worker
    for j in 1:m
        @constraint(model, T >= sum(production_times[j, i] * x[i, j] for i in 1:n))
    end
end

# Solve the model and measure optimization time
optimization_time = @elapsed optimize!(model)

total_time = time() - start_time

println("\nTiming Results:")
println("Model setup time: ", @sprintf("%.4f", model_setup_time), " seconds")
println("Optimization time: ", @sprintf("%.4f", optimization_time), " seconds")
println("Total execution time: ", @sprintf("%.4f", total_time), " seconds\n")

println("Maximum workload (T): ", objective_value(model))

# Print assignments
for j in 1:m
    assigned_tasks = [i for i in 1:n if value(x[i, j]) > 0.5]
    if !isempty(assigned_tasks)
        println("Worker $j: Tasks ", assigned_tasks, " | Workload: ", sum(production_times[j, i] for i in assigned_tasks))
    end
end

# Start timing for the new section
new_section_start_time = time()

# Funções auxiliares
function calculate_max_workload(assignment, production_times, n, m)
    workloads = zeros(m)
    for j in 1:m
        for i in 1:n
            if assignment[i] == j
                workloads[j] += production_times[j, i]
            end
        end
    end
    return maximum(workloads)
end

function is_valid_assignment(assignment, m)
    for j in 1:m
        assigned_indices = findall(x -> x == j, assignment)
        if !isempty(assigned_indices)
            if maximum(assigned_indices) - minimum(assigned_indices) + 1 != length(assigned_indices)
                return false
            end
        end
    end
    return true
end
# Geração de solução inicial
function generate_initial_solution(n, m)
    # Divide as tarefas em m blocos contíguos
    split_points = sort([rand(1:n-1) for _ in 1:(m-1)])
    split_points = unique(split_points)
    while length(split_points) < m-1
        new_point = rand(1:n-1)
        if !(new_point in split_points)
            push!(split_points, new_point)
        end
    end
    split_points = sort(unique(split_points))
    push!(split_points, n)

    assignment = zeros(Int, n)
    start = 1
    for (j, stop) in enumerate(split_points)
        assignment[start:stop] .= j
        start = stop + 1
    end

    # Alguns workers podem ficar sem tarefas, então redistribuímos
    used_workers = unique(assignment)
    if length(used_workers) < m
        unused_workers = setdiff(1:m, used_workers)
        for worker in unused_workers
            if length(used_workers) >= 1
                donor = rand(used_workers)
                donor_tasks = findall(x -> x == donor, assignment)
                if length(donor_tasks) > 1
                    split_point = rand(1:(length(donor_tasks)-1))
                    assignment[donor_tasks[1:split_point]] .= worker
                    push!(used_workers, worker)
                end
            end
        end
    end

    return assignment
end
# Geração de vizinhos
function generate_neighbor(assignment, m)
    n = length(assignment)
    new_assignment = copy(assignment)

    i = rand(1:n)
    current_worker = assignment[i]

    start = i
    while start > 1 && assignment[start-1] == current_worker
        start -= 1
    end
    stop = i
    while stop < n && assignment[stop+1] == current_worker
        stop += 1
    end

    if rand() < 0.5 && (stop - start + 1) > 1
        if i - start > 0 && stop - i > 0
            split_point = i
            new_worker = rand(setdiff(1:m, [current_worker]))
            new_assignment[split_point+1:stop] .= new_worker
        else
            new_worker = rand(setdiff(1:m, [current_worker]))
            new_assignment[start:stop] .= new_worker
        end
    else
        new_worker = rand(setdiff(1:m, [current_worker]))
        new_assignment[start:stop] .= new_worker
    end

    if !is_valid_assignment(new_assignment, m)
        return assignment
    end

    return new_assignment
end
# Algoritmo Simulated Annealing
function simulated_annealing(production_times, n, m; max_iter=10000, initial_temp=100.0, cooling_rate=0.995)
    current_solution = generate_initial_solution(n, m)
    current_cost = calculate_max_workload(current_solution, production_times, n, m)

    best_solution = copy(current_solution)
    best_cost = current_cost

    temp = initial_temp
    costs = Float64[]

    for iter in 1:max_iter
        new_solution = generate_neighbor(current_solution, m)
        new_cost = calculate_max_workload(new_solution, production_times, n, m)

        delta_cost = new_cost - current_cost

        if delta_cost < 0 || rand() < exp(-delta_cost / temp)
            current_solution = new_solution
            current_cost = new_cost

            if current_cost < best_cost
                best_solution = copy(current_solution)
                best_cost = current_cost
            end
        end

        push!(costs, best_cost)
        temp *= cooling_rate

        if iter % 1000 == 0
            temp = initial_temp * 0.5
        end
    end

    return best_solution, best_cost, costs
end

# Execute and measure the simulated annealing part
sa_time = @elapsed begin
    Random.seed!(123)  # For reproducibility
    best_solution, best_cost, cost_history = simulated_annealing(
        production_times, n, m,
        max_iter=1200000,
        initial_temp=250.0,
        cooling_rate=0.8
    )
end

# Execute and measure the warm start model setup
model_warm_setup_time = @elapsed begin
    # Create complete LP model
    model_warm = Model(GLPK.Optimizer)

    @variable(model_warm, x[1:n, 1:m], Bin)
    @variable(model_warm, b[1:n, 1:m], Bin)  # Block start
    @variable(model_warm, e[1:n, 1:m], Bin)  # Block end
    @variable(model_warm, T >= 0)

    @objective(model_warm, Min, T)

    # Original constraints (contiguity)
    for i in 1:n
        @constraint(model_warm, sum(x[i, j] for j in 1:m) == 1)
    end

    for j in 1:m
        @constraint(model_warm, sum(b[i, j] for i in 1:n) == 1)
        @constraint(model_warm, sum(e[i, j] for i in 1:n) == 1)

        for i in 1:n
            prev_x = (i == 1) ? 0 : x[i-1, j]
            next_x = (i == n) ? 0 : x[i+1, j]

            @constraint(model_warm, b[i, j] >= x[i, j] - prev_x)
            @constraint(model_warm, b[i, j] <= x[i, j])
            @constraint(model_warm, b[i, j] <= 1 - prev_x)

            @constraint(model_warm, e[i, j] >= x[i, j] - next_x)
            @constraint(model_warm, e[i, j] <= x[i, j])
            @constraint(model_warm, e[i, j] <= 1 - next_x)
        end
    end

    for j in 1:m
        @constraint(model_warm, T >= sum(production_times[j, i] * x[i, j] for i in 1:n))
    end

    # Warm Start using SA solution
    for i in 1:n, j in 1:m
        set_start_value(x[i,j], best_solution[i] == j ? 1.0 : 0.0)

        # Set initial values for b and e (block start/end)
        if best_solution[i] == j
            prev_val = (i > 1) ? best_solution[i-1] : 0
            next_val = (i < n) ? best_solution[i+1] : 0
            set_start_value(b[i,j], (best_solution[i] == j) && (prev_val != j))
            set_start_value(e[i,j], (best_solution[i] == j) && (next_val != j))
        else
            set_start_value(b[i,j], 0.0)
            set_start_value(e[i,j], 0.0)
        end
    end
    set_start_value(T, best_cost)

    # Set solver parameters
    set_optimizer_attribute(model_warm, "presolve", 1)
    set_optimizer_attribute(model_warm, "mip_gap", 0.01)  # Accept 1% gap
    #set_optimizer_attribute(model_warm, "tm_lim", 300000)  # Time limit (5 minutes)
end

# Measure optimization time with warm start
warm_optimization_time = @elapsed optimize!(model_warm)

total_new_section_time = time() - new_section_start_time

println("\nAdditional Timing Results:")
println("Simulated Annealing execution time: ", @sprintf("%.4f", sa_time), " seconds")
println("Warm start model setup time: ", @sprintf("%.4f", model_warm_setup_time), " seconds")
println("Warm start optimization time: ", @sprintf("%.4f", warm_optimization_time), " seconds")
println("Total time for new section: ", @sprintf("%.4f", total_new_section_time), " seconds\n")

println("Solution with Warm Start:")
println("Maximum workload: ", objective_value(model_warm))
println("Optimality gap: ", relative_gap(model_warm))

# Comparative visualization
sa_assignment = best_solution
hybrid_assignment = [findfirst(j -> value(x[i,j]) > 0.5, 1:m) for i in 1:n]

p1 = bar(sa_assignment, title="Pure SA", ylim=(0,m+1))
p2 = bar(hybrid_assignment, title="SA + LP with Warm Start", ylim=(0,m+1))
plot(p1, p2, layout=(2,1), size=(800,600))
