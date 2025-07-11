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
m = 15  # número de trabalhadores
production_times = [
    0.890661 0.100000 0.177132 0.100000 0.529553 0.390666 0.100000 0.100000 0.821462 0.785439 1.617629 0.100000 0.852795 0.560338 0.060270 0.261902 1.130127 0.909629 0.422074 0.055287 1.057074 0.313438 0.100000 0.165839 0.232135 0.425175 0.627948 0.100000 0.100000 0.385133 0.349847 0.100000
    0.100000 0.201730 1.328912 0.100000 0.100000 0.978650 0.780160 0.562959 0.583223 1.397648 0.100000 0.412369 0.486314 0.397564 0.100000 0.052441 0.957283 0.386480 0.474806 0.335646 0.036059 0.020055 0.947573 0.081915 0.659697 0.100000 1.039202 0.278111 0.377961 0.100000 0.100000 0.374428
    0.382417 0.100000 1.726658 0.358937 0.100000 0.712248 0.361727 0.100000 0.561401 0.986316 1.727977 0.100000 0.579869 1.200957 0.196316 0.725943 0.100000 1.241777 0.098095 1.112497 0.100000 0.371163 1.290880 0.582376 0.343506 0.855024 0.275685 1.560870 0.100000 0.563130 0.708329 0.100000
    0.258192 0.655838 0.872570 0.289637 0.493474 0.701791 1.384474 0.280318 0.192337 0.706908 0.873142 1.147647 0.455163 0.408700 1.571699 0.492463 0.100000 1.816261 0.100000 0.030963 0.136302 0.861470 0.687186 0.336406 0.590922 0.111183 0.122248 0.373282 0.332134 0.100000 0.615357 0.732453
    0.495037 0.100000 0.758145 0.280339 0.100000 0.633850 0.100000 0.248400 0.303857 0.552071 0.410473 1.162012 0.100000 0.821449 0.198540 1.012698 0.100000 0.606464 0.100000 0.027242 0.100000 0.276056 0.310901 0.449378 1.293435 0.850291 0.429150 0.100000 0.100000 0.100000 0.477724 0.932241
    0.321022 0.100000 1.711453 0.100000 0.100000 1.158414 0.186435 0.781424 0.692581 0.423235 0.771024 0.384682 0.100000 0.968224 0.100000 1.888774 0.659884 0.882312 0.884424 0.369410 0.583516 0.092775 0.140019 1.014709 0.100000 0.100000 0.423446 0.910685 0.013962 0.126919 0.110869 0.100000
    0.007691 1.044083 0.807365 0.043492 0.100000 1.616439 0.501842 0.100000 0.305458 0.352975 0.299552 0.100000 0.321675 0.674008 1.185942 0.100000 0.186941 0.830165 0.100000 1.000456 0.100000 0.027622 0.692752 0.245957 0.258509 0.501406 0.100000 1.107569 0.086911 1.057563 0.100000 0.100000
    0.695813 0.100000 0.190341 0.823744 0.095337 1.093666 0.320435 0.977830 0.440740 0.708548 0.402224 0.259841 0.966075 0.100000 0.045052 0.316506 0.376642 0.744349 1.140295 0.681118 1.691094 0.100000 0.554755 0.647983 0.817335 0.100000 0.880458 0.115274 0.100000 0.900057 1.054724 1.260020
    0.253746 0.116859 0.488320 0.848778 0.100000 0.322532 0.100000 0.100000 0.100000 0.843201 0.388036 0.100000 0.611051 1.040706 0.750740 0.494748 0.100000 0.794476 0.100000 0.795443 0.431126 0.100000 1.435126 0.516085 0.747728 1.066573 0.772020 0.795012 0.306435 0.100000 0.558718 1.119541
    0.656893 0.416342 1.650894 0.984406 0.742756 0.169307 0.727764 1.083716 0.100000 0.687584 0.737152 0.100000 0.631042 0.549619 0.472649 0.781801 0.897168 1.518833 1.032926 0.100000 0.323902 0.459780 0.660436 0.695000 0.012182 0.297469 0.070472 0.331113 0.100000 0.100000 0.100000 0.728237
    0.499228 0.904172 1.185095 0.534964 0.100000 0.406320 0.585389 0.100000 0.100000 0.984580 1.247980 0.583254 0.104195 1.159520 0.748554 0.531281 0.100000 0.906721 0.695308 0.773314 1.123692 0.775907 0.167090 1.083280 0.472975 0.100000 0.100000 1.263128 0.587620 0.113164 0.837569 0.208604
    0.142952 0.686266 0.762997 0.588483 0.100000 0.404907 0.100000 1.377906 0.525836 0.560117 0.133794 0.269903 0.731784 1.547233 1.441720 0.824900 0.750518 0.935396 0.344110 0.249128 0.100000 0.170307 0.908666 1.051160 0.857423 0.100000 0.514229 0.741925 0.100000 1.599592 1.014286 1.443868
    1.437351 0.100000 1.650064 0.043201 0.216370 0.994453 1.308924 0.196658 0.925381 0.649797 0.625702 0.100000 0.794193 0.630061 0.507755 1.234356 0.705815 1.411497 0.990627 1.066452 0.908234 0.100000 0.985587 0.008307 0.385334 0.793140 0.906863 0.277619 0.025984 0.231853 0.682592 0.768844
    0.649845 0.100000 0.100000 0.608595 0.445009 0.137827 1.080639 0.689533 0.168804 0.841617 0.596912 0.641956 0.567373 0.881226 0.326193 0.708751 0.711494 0.441111 0.245572 0.340515 0.854952 0.329267 0.100000 0.961670 1.254379 1.216361 1.230060 0.311179 0.184767 0.250508 1.363192 0.224355
    0.963625 0.137698 0.567678 0.660414 0.143778 0.441318 0.100000 0.697878 0.275290 0.066416 1.165624 0.344407 0.100000 0.797215 0.603002 1.956298 0.100000 0.629862 0.217097 0.201358 0.765382 0.100000 0.970146 0.870645 0.803379 0.633148 0.875797 0.753701 0.232835 0.100000 0.761570 0.171585
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
