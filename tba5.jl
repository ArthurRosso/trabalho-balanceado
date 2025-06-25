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
    0.824972 0.748855 0.985319 0.368152 0.494423 1.352047 0.365010 1.477899 0.117156 0.164884 1.168871 0.286475 0.299580 0.506642 1.019126 0.156562 1.550568 0.663989 0.653186 0.529661 0.100000 0.611770 0.495251 0.869660 0.041371 0.558254 0.107564 0.783252 0.953264 0.567330 0.342755 0.197970 1.259813 0.681452 0.861444 0.917369 0.484125 0.085160 0.088678 0.730138 0.100000 0.754547
    0.100000 0.100000 0.850774 0.872260 0.581743 0.298839 1.144109 0.650960 1.123835 0.222335 0.880949 0.898076 0.215019 0.453726 0.865229 0.689777 0.637151 0.773566 1.314217 1.093563 0.100000 0.689749 0.324512 0.100000 0.821920 0.124092 0.509316 0.476336 1.650988 0.100000 0.419885 0.100000 1.372962 0.766863 0.100000 0.243731 0.297435 1.045586 1.572023 1.144180 0.270404 0.100000
    0.894116 0.143745 0.042561 0.609144 0.100000 0.100000 0.798137 1.095015 0.102034 0.682849 0.854231 0.753674 0.131011 0.292548 0.471175 0.614766 0.655073 0.561503 0.100000 0.434639 0.100000 0.638959 0.299925 0.937692 0.384593 0.100000 0.100000 0.417338 0.575980 0.100000 0.299497 0.100000 1.002387 0.850871 0.185761 0.134252 0.100000 0.899535 0.698389 0.244300 0.559528 0.352028
    0.535348 0.100000 0.088750 0.550408 1.022337 0.552360 0.326791 0.911376 0.069824 1.133790 0.037736 0.542586 1.032420 0.254768 0.852856 0.514424 0.604272 0.848203 0.100000 0.613280 0.100000 0.629088 0.050542 0.536960 0.004763 0.860313 0.395858 0.283070 0.100000 0.299935 0.554562 0.065988 1.263242 0.100000 1.060334 0.540585 0.100000 0.100000 0.399530 0.659336 0.383587 0.100000
    0.375510 0.100000 0.100000 0.973676 0.100000 0.076413 0.230734 0.992201 0.100000 0.100000 1.470864 0.952046 1.223138 0.562517 0.325127 1.200134 0.604143 0.174244 0.339153 1.235328 0.100000 0.100000 0.747576 0.100000 0.100000 0.324085 0.100000 0.832647 0.534109 0.556055 0.397129 0.313009 1.010865 0.122663 0.296786 0.711442 0.569859 0.480454 0.424974 0.735372 0.212469 0.807730
    0.253857 0.445623 0.505591 0.634472 0.100000 0.294177 0.387532 0.100000 0.100000 0.607924 0.835065 0.265243 0.100000 0.409799 0.637082 0.100000 0.848454 0.100000 0.578866 0.100000 0.232893 0.498980 0.308330 0.812449 0.079669 0.884524 0.100000 0.337672 0.447057 0.156672 0.857158 0.100000 0.568384 0.100000 0.662403 0.566868 2.061965 0.869289 0.100000 0.820106 0.499168 0.100000
    1.241708 0.411287 0.238223 0.395653 0.100000 0.100000 0.936979 0.220080 0.205308 0.572478 1.352376 0.552238 0.220584 1.024534 0.574103 0.218287 0.100000 0.100000 0.100000 1.200176 0.100000 0.933059 0.536454 2.045580 0.041700 0.789185 1.295710 0.620364 0.100000 0.100000 0.458513 0.100000 0.874109 0.478668 0.594506 0.229283 0.692941 0.476497 0.953956 0.100000 0.869887 0.977500
    0.426756 0.100000 0.165860 0.456450 0.100000 0.100000 0.666598 0.100000 0.475684 0.527811 0.171359 0.261725 0.100000 1.014117 0.986371 0.100000 0.714449 0.100000 0.731452 0.263199 0.211229 0.100000 0.836487 0.134902 0.691623 0.100000 0.567231 0.507526 1.216192 0.184252 0.482263 0.100000 1.721344 0.100000 0.042550 1.371902 0.100000 0.116663 1.287267 0.001099 0.100000 0.374444
    0.100000 0.566423 0.484050 0.679989 0.100000 0.100000 0.833016 0.100000 0.173963 0.100000 0.723740 0.388600 0.163106 0.503850 1.425759 0.210710 0.838725 0.100000 0.983289 0.500409 0.100000 0.886672 0.841904 0.110103 0.688876 0.824044 0.598258 0.445340 0.592687 0.403196 0.572549 0.732552 0.749755 0.862013 0.100000 0.100000 0.115333 1.015375 0.898168 0.722237 0.387101 0.069993
    0.100000 0.073973 0.500032 0.084874 0.482922 0.056343 0.268271 0.677015 0.342105 0.104008 0.100000 0.589934 0.867905 0.100000 0.780094 0.859836 0.823127 0.483041 0.100000 0.489803 0.100000 1.092856 0.442385 0.433423 0.448439 0.100000 0.100000 0.369570 0.948891 0.178382 0.100000 0.100000 1.153027 1.263044 0.100000 0.148087 0.769704 0.392272 0.185973 0.146923 0.100000 0.299229
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
