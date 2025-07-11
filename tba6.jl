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
    0.183139 0.605582 0.305518 0.090858 0.370164 0.724743 0.590451 0.379136 0.601041 0.817994 0.767640 0.452775 0.637407 0.645494 0.616733 0.592551 0.391614 0.658357 1.004405 0.098053 0.220737 0.192181 0.023187 0.590805 0.748719 0.100000 0.744637 0.210445 0.477795 0.290140 0.100000 0.080564 0.410204 0.857599 0.100000 0.575311 0.431696
    0.722704 0.379281 0.217166 0.823136 0.265645 0.578194 1.118839 0.405270 0.289360 0.565096 0.906647 0.290319 0.442410 0.536681 1.049764 0.687675 0.938496 0.600024 0.600196 0.501004 0.675099 0.204400 0.290622 1.042604 0.731260 0.415747 1.095737 0.390081 0.768204 0.589006 0.486085 0.490498 0.429869 0.516831 0.348893 0.352515 0.423103
    0.549786 0.968937 0.073298 0.436735 0.005966 0.384408 0.183142 0.158692 0.849123 0.300749 1.151384 0.140483 0.901022 0.397100 0.572663 0.334458 0.549263 0.112161 0.818907 0.101201 0.100000 0.100000 0.100000 0.105655 0.054419 0.100000 0.192317 0.299575 0.100000 0.364829 0.119001 0.348591 0.526368 0.264619 0.333689 0.447103 0.709307
    1.065906 0.100000 0.100000 0.100000 0.100000 0.475433 0.965144 0.965669 0.831711 0.100000 0.214345 1.247420 0.868751 0.790319 0.579308 0.351324 0.100000 0.100000 0.659217 0.630332 0.100000 0.100000 0.447301 1.438199 1.333385 0.100000 0.798105 0.100000 0.943303 0.027273 0.083881 0.604483 0.157030 0.197905 0.958297 1.372184 1.122930
    0.402125 0.342582 0.100000 0.523577 0.189882 0.283262 0.461351 0.100000 0.830919 0.100000 1.251473 0.252590 0.653480 0.183139 0.465052 0.537518 0.100000 0.773098 0.943565 0.143015 0.282882 0.100000 0.184242 0.335046 0.382780 0.076320 0.528994 0.036794 0.268950 0.010768 0.190597 0.175687 0.071637 0.100000 0.100000 0.240933 0.530501
    0.283069 0.685761 0.085950 0.101148 0.451844 0.587780 0.528450 0.100000 0.450506 0.100000 0.369422 0.100000 0.431462 0.100000 0.100000 1.270666 1.271223 1.013762 0.451800 0.100000 0.100000 0.100000 0.100000 0.514748 0.899151 0.198009 0.733316 0.100000 0.056750 0.483165 0.100000 0.314797 0.858601 0.188334 0.109919 0.100000 1.347325
    0.100000 0.100000 0.100000 0.100000 0.463117 1.300689 1.268158 0.328727 0.768743 0.881979 0.685140 0.415193 0.100000 1.145187 0.792196 0.100000 0.467790 0.818631 0.977748 0.811144 0.156274 0.095643 0.926747 1.654110 1.314019 0.100000 0.614464 0.325087 0.141961 0.094953 0.100000 0.100000 0.000613 1.049692 0.026766 0.100000 0.959405
    0.100000 0.492369 0.100000 0.248441 0.169909 0.560805 0.244093 0.156967 0.236739 0.609297 1.290406 0.606156 0.411470 0.769899 0.100000 0.443961 0.289807 0.341442 0.399150 0.100000 0.227023 0.100000 0.041254 1.020021 0.327186 0.100000 0.748342 0.423610 0.457473 0.100000 0.691022 0.100000 0.613031 0.619714 0.100000 0.281635 0.260174
    0.436333 0.642633 0.051075 0.749964 0.100000 0.299087 0.849088 0.938026 0.423404 0.182445 1.059520 0.100000 0.666932 0.752961 0.218387 0.362307 0.100000 0.461497 0.649591 0.500226 0.174386 0.152401 0.118339 0.571553 0.440500 0.228265 1.177162 0.100000 0.100000 0.875976 0.196459 0.100000 0.334072 0.359768 0.100000 0.226761 0.808430
    1.071027 0.100000 0.128019 0.100000 0.665422 0.823401 0.762405 0.855959 0.100000 0.224711 0.871897 0.181989 0.262065 1.178418 0.436375 0.037698 1.309908 0.393916 0.783464 0.312663 0.513756 0.943532 0.058844 0.758983 1.272479 0.049405 1.142029 0.247801 0.677035 0.023073 0.100000 1.043186 0.785654 1.014401 0.100000 1.201406 0.100000
    0.100000 0.277295 0.100000 0.354871 0.165165 0.353344 0.773199 0.505180 0.170796 0.705505 1.012005 0.100000 0.893118 0.238793 0.100000 0.012517 0.033618 0.360886 1.241041 0.021858 0.100000 0.389864 0.032311 0.480881 0.435830 0.426154 0.991971 0.593925 0.130469 0.210255 0.525162 0.187549 0.295723 0.100000 0.100000 0.100000 0.100000
    0.415755 0.422200 0.356851 0.220622 0.100000 0.100000 0.560276 0.083672 0.136295 0.027899 0.174859 0.177206 0.100000 0.124122 0.237650 0.009147 0.159341 0.661184 0.668758 0.100000 0.534783 0.282027 0.266042 0.690345 0.831951 0.100000 0.663226 0.100000 0.522760 0.359334 0.140553 0.182167 0.342812 0.804803 0.100000 0.768090 0.364154
    0.986747 0.097920 0.100000 0.399361 0.149589 0.471358 0.579629 0.246008 0.950071 0.518191 0.780711 0.383945 0.383532 0.636100 0.446225 0.134626 0.175217 0.466751 0.267180 0.099315 0.200185 0.394893 0.399592 0.752007 1.025858 0.345202 1.065698 0.100000 0.116381 0.117044 0.068059 0.526282 0.307963 0.421333 0.741523 0.181912 0.821711
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
