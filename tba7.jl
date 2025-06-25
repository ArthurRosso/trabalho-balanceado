import Pkg; Pkg.add("JuMP"); Pkg.add("GLPK"); Pkg.add("DelimitedFiles"); Pkg.add("Random");

using JuMP          # Modeling language for optimization
using GLPK          # Open-source solver for linear programming
using Printf        # For formatted printing of time
using Random
using Plots

# Start timing
start_time = time()

# Dados do problema
n = 32  # número de operações
m = 14  # número de trabalhadores
production_times = [
    0.394455 0.887469 0.059056 0.630017 0.192982 0.309118 0.347316 0.356712 1.480903 1.244382 0.593374 0.804714 0.185590 0.472639 0.526727 0.669363 0.175841 0.631347 0.425798 0.901230 0.550198 0.287530 0.748237 0.073036 0.557304 0.218358 0.375037 0.215329 0.184584 0.290624 0.844265 0.075421
    0.175642 0.138501 0.989093 0.100000 0.046554 0.100000 0.100000 0.520904 1.317292 1.336659 0.899442 0.307156 0.478133 0.100000 0.771518 0.788452 0.814638 0.535061 0.306833 0.100000 0.372234 0.384788 0.100000 0.786044 0.634378 0.433438 1.071173 0.100000 0.117785 0.469208 0.444131 0.060578
    0.491818 1.041722 0.366432 0.463663 0.503861 0.100000 0.334810 0.466481 0.335455 0.210806 0.100000 0.808876 0.140211 0.311972 0.679312 0.509171 0.100000 0.339150 0.890086 0.419357 0.576498 0.135279 0.169215 0.236398 0.289106 0.051084 0.093069 0.454653 0.040544 0.044223 0.495571 0.275243
    0.152060 0.564769 0.609526 0.641453 0.471718 0.361538 0.596131 0.275872 0.123922 1.257859 0.230196 0.339581 0.479209 0.793445 0.161461 0.661456 0.864926 0.538613 0.616115 0.637031 0.179539 1.062452 0.533542 0.194476 0.118978 0.178660 0.954042 0.461457 0.912611 0.287839 0.167807 0.849978
    0.702860 0.100000 0.582303 0.311126 0.246378 0.015724 0.103491 0.100000 1.185220 0.787430 0.100000 0.100000 0.100000 0.891934 1.224721 0.496911 0.631752 0.947047 0.100000 0.531682 0.705449 0.828490 0.265224 0.100000 0.079252 0.619900 1.466792 0.419110 0.100000 0.035892 0.100000 0.083834
    0.100000 0.772797 0.625066 0.467431 0.254559 0.100000 0.357131 0.647448 0.345371 0.575543 0.100000 0.369325 0.100000 0.100000 0.301792 0.480311 0.196490 0.055209 0.397203 0.627999 0.481369 0.016959 0.297014 0.117415 0.773319 0.100000 0.618367 0.469299 0.231496 0.349081 0.379292 0.100000
    1.323624 0.283758 0.762827 0.100000 0.487584 0.479639 0.280230 0.635129 1.074527 0.912477 0.899896 1.001921 0.548699 0.100000 0.100000 0.789036 0.100000 1.061424 0.289565 0.590118 0.732408 0.100000 0.354167 0.638481 1.293346 0.100000 0.219919 0.188593 0.688738 0.250313 1.829285 0.560429
    0.968449 0.005319 1.502469 0.017305 0.581375 0.100000 0.772587 0.428938 0.337353 0.742352 0.217709 0.510060 0.132199 0.615168 0.846812 0.184614 0.617952 0.637964 0.071788 0.797156 0.303474 0.116420 0.100000 0.051611 0.730077 0.449382 1.039248 0.015841 0.100000 0.428449 0.351001 0.100000
    0.501568 0.831582 1.052152 0.891741 0.194699 0.019053 0.397246 0.401026 1.173963 0.497298 0.342927 0.117926 0.919601 0.689348 1.079065 1.046354 0.100000 0.312415 0.489932 0.693381 0.577590 0.155443 0.100000 0.696266 0.887702 0.506390 0.106709 0.363311 0.581033 0.699732 1.084987 0.391246
    0.416857 1.105155 0.502641 0.501636 0.192490 0.065014 1.140042 0.905681 0.709029 0.438350 0.100000 0.279268 0.781009 0.293687 0.601738 0.713376 0.221829 0.884528 0.852070 0.639786 0.100000 0.182462 0.899737 0.558000 0.338634 0.338584 0.216920 0.490104 0.121888 0.459100 0.605081 0.100000
    0.692023 0.703051 0.100000 0.493539 0.504104 0.100000 0.603327 0.100000 0.063652 0.871094 0.735658 0.429101 0.470496 0.014040 0.340105 0.364202 0.100000 0.386374 0.657068 0.602427 0.451308 0.100000 0.078635 0.649230 0.494453 0.100000 0.915078 0.100000 1.221627 0.054939 0.241175 0.479627
    0.425164 1.257863 0.100000 0.740668 0.100000 1.157179 0.888621 0.521177 1.169033 0.100000 0.100000 0.493194 0.990943 0.627544 0.607544 1.114142 0.100000 0.795785 0.100000 0.941198 1.105178 0.373278 0.930544 0.536590 0.100000 0.704848 1.309884 0.831666 0.559446 0.747816 0.138316 0.138855
    0.240873 0.591146 0.711213 0.339020 0.936425 0.394629 0.338113 0.100000 0.285101 0.617951 0.100000 0.554995 0.304907 0.018462 0.442478 0.100000 0.082952 0.188861 0.517139 0.960084 0.187193 0.100000 0.012574 0.040193 0.033629 0.621127 0.745115 0.044945 0.643159 0.194942 0.161233 0.658101
    0.390067 0.455697 1.353690 0.623134 0.803194 0.100000 0.433356 0.262178 0.133002 0.308981 0.423883 0.382316 0.353663 0.283399 0.745939 0.887463 0.366587 0.285726 0.632973 0.384544 1.271791 0.599024 0.141999 0.100000 0.327378 0.766905 1.076373 0.546487 1.032184 0.397710 0.880639 0.823773
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
