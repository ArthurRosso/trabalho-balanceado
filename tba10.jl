import Pkg; Pkg.add("JuMP"); Pkg.add("GLPK"); Pkg.add("DelimitedFiles"); Pkg.add("Random");

using JuMP          # Modeling language for optimization
using GLPK          # Open-source solver for linear programming
using Printf        # For formatted printing of time
using Random
using Plots

# Start timing
start_time = time()

# Dados do problema
n = 42  # número de operações
m = 10  # número de trabalhadores
production_times = [
    0.100000 0.157916 1.262799 1.290917 0.109419 0.956469 1.048012 0.267745 0.616112 0.381447 0.142017 0.100000 0.100000 0.275619 1.067472 0.116311 0.043384 0.100000 0.753886 0.477781 0.100000 1.404547 0.127177 1.373362 0.714355 0.217598 0.309765 1.207631 1.205433 0.130243 1.219805 0.335577 0.589780 1.709745 0.139897 0.626678 1.334373 0.803798 0.700504 0.466742 0.722569 0.557902
    0.988762 0.026033 0.344655 0.603306 0.587490 0.652222 0.421817 0.088780 0.852775 0.218461 0.801375 0.222353 0.549195 0.600392 0.637823 0.544367 0.610101 0.234584 0.665633 0.691910 0.571449 0.589709 0.100000 0.286885 0.100000 0.436417 0.730759 1.207174 0.985311 0.151022 0.792664 0.817342 0.778123 0.994229 0.494307 0.405089 0.169122 0.752759 0.663545 0.496304 0.122936 0.876286
    0.572378 0.100000 0.011968 0.950888 0.072999 0.812977 0.084497 0.427929 0.772320 0.616919 0.185853 0.573961 0.339514 0.353842 0.722161 0.100000 0.100000 0.069952 0.628917 0.005640 0.100000 0.367631 0.100000 1.223211 0.645359 0.100000 0.281969 0.281988 0.699690 0.100000 0.850426 0.498823 1.022732 0.998399 0.284923 1.037996 0.756077 0.412065 0.914034 0.520254 0.387817 0.610792
    0.068620 0.014490 0.055620 0.100000 0.583422 0.199094 0.657882 0.348501 0.614044 0.754357 0.147757 0.100000 0.015072 0.435123 1.035354 0.100000 0.194400 0.172988 0.100000 0.343841 0.100000 0.022228 0.100000 1.509331 0.815410 0.322853 0.571434 0.676112 0.945204 0.386271 0.394986 0.983189 0.283517 0.628565 0.145346 0.191645 0.412015 0.882287 0.877377 1.424677 0.385175 0.505852
    0.797531 0.399259 0.203260 0.814548 0.465033 0.966060 0.072941 0.749886 0.531740 1.046584 0.100000 0.100000 0.465523 0.909956 0.815424 0.431637 0.387499 0.037826 0.891858 0.025489 0.100000 0.996827 0.177246 0.673288 0.100000 0.352558 0.450197 1.442929 0.001378 0.150140 0.674537 0.100000 0.124727 1.084218 0.100000 0.757676 0.752393 0.100000 1.530032 0.682603 1.550781 0.572838
    0.844141 0.100000 0.261215 0.306577 0.983179 0.275689 0.715033 0.299600 0.063611 0.678545 0.057644 0.100000 0.399239 1.384731 0.445457 0.314536 0.100000 0.428565 0.175239 0.591247 0.152156 0.653230 0.202500 0.674880 0.652734 0.625663 0.986547 0.543733 0.474079 0.412397 1.076239 0.669097 0.650080 0.927165 0.333914 1.090906 0.202114 0.303397 1.474339 0.683823 0.221622 0.575815
    0.535806 0.221485 0.100000 0.402019 0.235406 0.618712 0.345675 0.010342 1.003911 0.374969 0.420879 0.100000 0.267348 1.104489 0.100000 0.143550 0.100000 0.100000 0.106673 0.243619 0.100000 0.511684 0.269377 0.495636 0.710855 0.142100 0.597025 0.304712 0.280579 0.502231 0.774837 0.945471 0.985205 0.226248 0.593024 0.404204 0.100000 0.877897 0.315674 0.570285 0.964077 0.897042
    0.183328 0.736867 0.100000 0.100000 1.325092 0.100000 0.105104 1.338662 0.669275 1.209392 0.527939 0.018967 0.108079 0.544422 1.794255 0.697402 0.472250 0.228696 1.347705 0.306829 0.100000 0.885774 0.774910 0.430882 0.100000 0.100000 0.100000 0.100000 0.438878 0.100000 0.781502 1.475757 0.100000 0.484092 0.297532 0.374247 0.100000 1.056941 0.972823 0.352911 0.100000 0.575560
    0.466190 0.509117 0.325959 0.372250 0.100000 0.249220 0.100000 0.382986 0.505143 0.408174 0.205062 0.162168 0.257930 0.743567 0.672855 0.346580 0.820722 0.164462 0.454791 0.552017 0.252859 1.004156 0.280069 1.028491 0.646793 0.353609 0.290439 0.263874 0.656616 0.100000 0.885402 0.725558 0.824257 0.522019 0.675642 0.684066 0.865765 0.408710 1.021669 0.632090 0.575717 0.691220
    0.873186 0.882416 0.188929 0.295729 0.359003 0.597654 1.232744 0.001473 0.694448 1.028069 0.464547 0.577962 0.100000 0.770422 0.311185 0.638233 0.100000 0.678897 0.550846 0.358884 0.513845 0.400113 0.278348 0.528436 1.246539 0.100000 0.855153 0.377221 0.178270 0.138585 1.657865 0.966718 0.219303 0.100000 0.100000 0.961428 0.895742 0.453777 1.075138 0.175247 0.387570 0.191486
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
