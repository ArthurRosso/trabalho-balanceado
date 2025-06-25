import Pkg; Pkg.add("JuMP"); Pkg.add("GLPK"); Pkg.add("DelimitedFiles"); Pkg.add("Random");

using JuMP          # Modeling language for optimization
using GLPK          # Open-source solver for linear programming
using Printf        # For formatted printing of time
using Random
using Plots

# Start timing
start_time = time()

# Dados do problema
n = 37  # número de operações
m = 13  # número de trabalhadores
production_times = [
    0.423654 0.283088 0.401749 0.537266 0.659383 0.187169 0.726688 0.100000 1.501611 0.100000 0.299022 0.265373 0.408561 0.012306 0.100000 1.357203 0.893916 1.045484 0.100000 0.447577 0.100000 0.100000 0.642908 0.161556 0.222289 1.049944 1.222469 0.100000 0.580845 0.985271 1.864559 0.124676 0.366306 0.498465 1.570959 0.692720 0.607364
    0.869255 0.726535 0.139632 0.100000 0.406066 0.424875 1.391138 0.819067 0.792095 0.458505 0.569375 0.480167 0.808194 0.100000 1.074828 0.923875 0.100000 0.221318 0.647856 1.036310 0.354233 0.523441 0.803523 0.527607 0.364845 0.415553 0.100000 0.057583 0.727878 0.737818 0.652523 0.784421 0.841741 1.164389 0.267074 0.100000 1.062979
    0.572249 0.949738 0.586244 0.100000 0.804295 0.518671 0.078673 0.774149 0.385562 0.859574 0.138032 0.106480 1.527984 0.415206 0.100000 0.133900 0.289637 0.970430 0.761422 1.243011 0.219962 0.122468 0.673184 0.100000 0.089253 0.828612 0.690065 1.023505 0.223692 0.559279 1.125790 0.363079 0.357954 0.377114 0.087226 0.002081 0.077594
    0.773920 0.228788 0.115829 0.363470 2.006115 0.767712 0.332513 1.066157 0.587722 1.212666 0.100000 0.100000 0.182856 0.100000 0.947802 0.100000 0.405775 1.018406 0.676851 0.552494 0.119828 0.605882 0.100000 0.436218 0.100000 0.962646 0.100000 0.976239 0.100000 1.223247 0.585397 0.257549 0.289149 0.833516 0.426481 0.441588 0.100000
    0.100000 0.250595 0.589650 0.436220 0.690579 1.396425 0.100000 0.830534 0.644896 0.100000 0.100000 0.226418 0.100000 0.767794 0.956787 0.414027 0.100000 0.564348 1.210821 0.883664 0.558584 0.552752 0.622719 0.100000 0.115932 0.069594 1.551684 0.279594 0.100000 1.097523 0.100000 1.329288 1.280115 0.098221 0.442129 0.473180 0.100000
    0.766733 0.081720 0.794455 0.100000 1.075815 1.521521 0.132180 0.395028 0.100000 0.100000 0.100000 0.100000 0.951540 0.100000 0.980411 0.369661 0.373938 0.662027 1.083990 0.472567 0.720664 0.563668 0.638337 0.470541 0.013354 0.570563 0.100000 0.298356 0.510770 1.140620 0.778996 0.196933 1.292942 0.100000 0.402502 0.569648 0.980543
    0.357801 0.236615 0.473655 0.262793 0.100000 0.100000 0.100000 0.428089 1.694532 0.100000 0.583089 0.327736 0.100000 0.100000 0.100000 0.070174 0.997952 1.211436 1.214322 0.100000 0.452432 0.565002 0.351045 0.629129 0.100000 0.356196 0.114589 1.360129 0.100000 1.028605 0.886417 0.100000 0.997328 0.627691 0.585411 0.404629 0.412620
    0.100000 0.547307 0.615576 0.100000 0.566854 0.672745 0.099588 0.100000 0.100000 1.031771 0.100000 0.143497 0.152893 0.375519 0.100000 0.088041 0.394794 0.100000 0.526791 0.100000 0.100000 0.449829 1.199344 0.797251 1.510302 0.592457 0.667199 0.553175 0.225428 1.185482 0.516557 0.589015 0.100000 0.490431 0.835605 0.100000 0.057046
    0.983934 0.261846 1.231931 0.884645 0.423723 0.455395 0.365500 0.100000 1.151994 0.605645 0.392618 0.747435 0.664120 0.790678 0.539219 0.413247 0.100000 0.797657 0.328296 0.321891 0.606485 0.683920 0.340310 0.619162 0.405280 0.022599 1.070977 0.113888 0.100000 0.100000 1.238311 0.100000 0.799720 0.777482 0.099338 0.098271 0.230929
    0.789876 0.166174 0.100000 0.223350 1.259682 0.415396 0.100000 0.268183 0.180231 0.328989 0.100000 0.384930 0.831347 0.607636 0.100000 0.706883 0.100000 1.114479 0.228185 0.495110 0.195735 0.271490 0.126910 0.311468 0.783895 0.390367 0.731954 0.372603 0.754047 0.853063 0.100000 0.100000 0.185096 0.565500 0.536702 0.686774 0.100000
    0.749664 0.291782 0.404903 0.444787 0.100000 0.183903 0.100000 0.562395 1.088673 1.194003 0.100000 0.100000 0.388049 0.100000 0.238763 0.281265 0.308395 0.858590 0.606299 0.016527 0.575453 0.436060 0.741274 0.994668 0.670449 0.868317 0.099896 0.863107 0.168342 0.100000 0.615704 0.100000 0.042430 0.100000 0.267699 0.130255 0.424713
    0.100000 0.100000 0.314032 0.507438 0.362484 0.345224 0.797920 0.100000 0.770617 0.688759 0.165577 0.935507 1.060342 0.269227 0.514151 0.695363 0.839036 0.100000 0.100000 0.011191 1.292022 0.100000 0.100000 1.259202 0.100000 0.407771 0.994697 0.658166 0.349541 0.894319 0.848602 0.100000 1.862418 0.534835 0.584303 0.847889 0.945105
    1.306493 0.065853 0.091090 0.310499 0.794966 0.129894 0.630383 0.100000 0.443096 0.245461 0.052134 0.750293 0.156868 0.925495 0.100000 0.100000 0.100000 0.223969 0.100000 0.100000 0.404933 1.124190 0.100000 0.626568 0.100000 0.100000 0.100000 0.100000 0.415234 1.021791 0.359019 0.100000 0.748070 0.100000 0.890559 0.233397 0.100000
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
