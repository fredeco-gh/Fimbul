"""
    btes(; <keyword arguments>)

Setup function for borehole thermal energy storage (BTES) system.

# Keyword arguments
- `field = missing`: Explicit wells for the system. If given, this overrides
  `num_wells`, `num_sectors` and `well_spacing`. The wells are given using
  the following structure:
  - A full field is represented as `[sector_1, sector_2, ..., sector_n]`.
  - `sector_k` contains all wells in sector `k`, represented as
    `[well_1, well_2, ..., well_nk]`.
  - Each `well_l` is a `3 x m` matrix containing the coordinates of the well
    trajectory.
  Wells within a sector are coupled in series, in the order they appear in
  `sector_k`: `well_1` is always charged first and `well_nk` charged last.
  Discharge order depends on `reversed_discharge`. A single well (a `3 x m`
  matrix) can also be passed directly as `field`, representing the special
  case of a field with a single sector containing a single well.
- `pattern = :sunflower`: Well placement pattern used when `field` is not
  given. One of `:sunflower`, `:rectangular`, `:circular` or `:polygonal`.
- `num_sides = 6`: Number of sides of the polygon when `pattern = :polygonal`.
- `num_wells = 48`: Number of wells in the BTES system.
- `num_sections = 6`: Number of sections in the BTES system. The system is
- `num_sectors = 6`: Number of sections in the BTES system. The system is
  divided into equal circle sectors, and all wells in each sector are coupled in series.
- `well_spacing = 5.0`: Horizontal spacing between wells [m].
- `well_depth = 50.0`: Well depth [m]. Can be either a scalar (all wells), one
    value per sector, or one value per well.
- `layer_depths = [0.0, 0.5, 50, 65]`: Depths delineating geological layers [m].
- `density = [30, 2580, 2580]: Rock density in the layers [kg/m³].
- `thermal_conductivity = [0.034, 3.7, 3.7]: Thermal conductivity in the layers [W/(m⋅K)].
- `heat_capacity = [1500, 900, 900]`: Heat capacity in the layers [J/(kg⋅K)].
- `geothermal_gradient = 0.03 K/m`: Geothermal gradient [K/m].
- `temperature_charge = 90 °C/363.15 K`: Injection temperature during charging [K].
- `temperature_discharge = 10 °C/283.15 K`: Injection temperature during discharging [K].
- `rate_charge = 0.5 l/s`: Injection rate during charging [m³/s].
- `rate_discharge = rate_charge`: Injection rate during discharging [m³/s].
- `reversed_discharge = false`: All sectors are operated in parallel. During
  charging, flow runs from the first to the last well in each sector. If
  `reversed_discharge = false`, discharge uses the same direction as charge.
  If `true`, discharge flow is reversed, so it runs from the last well to the
  first well in each sector.
- `temperature_surface = 10 °C/283.15 K`: Temperature at the surface [K].
- `num_years = 5`: Number of years to run the simulation.
- `charge_period = ["June", "September"]`: Period during which the system is charged.
- `discharge_period = ["December", "March"]`: Period during which the system is discharged.
- `report_interval = 14 day`: Reporting interval for the simulation.
- `utes_schedule_args = NamedTuple()`: Additional arguments for the UTES schedule.
- `n_z = [3, 8, 3]`: Number of layers in the vertical direction for each layer.
- `n_xy = 3`: Number of layers in the horizontal direction for each layer.
- `mesh_args = NamedTuple()`: Additional arguments for the mesh generation.
"""
function btes(
    pattern::Symbol = :sunflower;
    num_sides = 6,
    num_wells = 48,
    num_sectors = 6,
    well_spacing = 5.0,
    well_depth = 50.0,
    kwargs...
    )

    if pattern == :sunflower
        field = sunflower_pattern(num_wells, well_spacing; num_sectors = num_sectors, well_depth = well_depth)
    elseif pattern == :rectangular
        field = rectangular_pattern(num_wells, well_spacing; num_sectors = num_sectors, well_depth = well_depth)
    elseif pattern == :circular
        field = circular_pattern(num_wells, well_spacing; num_sectors = num_sectors, well_depth = well_depth)
    elseif pattern == :polygonal
        field = polygonal_pattern(num_wells, well_spacing, num_sides; num_sectors = num_sectors, well_depth = well_depth)
    else
        error("Unknown pattern: $pattern. Supported patterns are :sunflower, :rectangular, :circular and :polygonal.")
    end

    return btes(field; kwargs...)

end

function btes(
    field::Vector{Vector{Matrix{Float64}}};
    depths = [0.0, 0.5, 50, 65],
    density = [30, 2580, 2580]*si_unit(:kilogram)/si_unit(:meter)^3,
    thermal_conductivity = [0.034, 3.7, 3.7]*si_unit(:watt)/si_unit(:meter)/si_unit(:Kelvin),
    heat_capacity = [1500, 900, 900]*si_unit(:joule)/si_unit(:kilogram)/si_unit(:Kelvin),
    geothermal_gradient = 0.03*si_unit(:Kelvin)/si_unit(:meter),
    temperature_charge = convert_to_si(90.0, :Celsius),
    temperature_discharge = convert_to_si(10.0, :Celsius),
    rate_charge = 0.5*si_unit(:litre)/si_unit(:second),
    rate_discharge = rate_charge,
    reversed_discharge = false,
    temperature_surface = convert_to_si(10.0, :Celsius),
    num_years = 4,
    charge_period = ["June", "September"],
    discharge_period = ["December", "March"],
    report_interval = 14*si_unit(:day),
    utes_schedule_args = NamedTuple(),
    n_z = [3, 8, 3],
    n_xy = 3,
    mesh_args = NamedTuple(),
    )

    if field isa AbstractMatrix
        # Special case: a single well, given as a single 3 x m matrix
        field = [[field]]
    end

    well_coords_3d = vcat(field...)
    num_wells = length(well_coords_3d)
    well_depths = [maximum(wc[3, :]) for wc in well_coords_3d]

    depths, updated_layer_data = insert_well_depths_in_layers(
        depths,
        well_depths,
        (
            density = density,
            thermal_conductivity = thermal_conductivity,
            heat_capacity = heat_capacity,
            n_z = n_z,
        ),
    )
    density = updated_layer_data.density
    thermal_conductivity = updated_layer_data.thermal_conductivity
    heat_capacity = updated_layer_data.heat_capacity
    n_z = Int.(updated_layer_data.n_z)

    # ## Create mesh
    # Use the (x,y) projection of each well trajectory as a mesh constraint
    cell_constraints = [unique(wc[1:2, :], dims=2) for wc in well_coords_3d]
    well_spacing = min_distance(hcat(cell_constraints...))
    hz = diff(depths)./n_z
    hxy = well_spacing/n_xy

    # ## Set up model
    # Set up reservoir domain with rock properties similar to that of granite,
    # with a styrofoam layer on top                
    domain, layers, metrics = layered_reservoir_domain(cell_constraints, depths,
        (
            rock_density = density,
            rock_thermal_conductivity = thermal_conductivity,
            rock_heat_capacity = heat_capacity
        );
        mesh_args = (; hxy_min = hxy, hz = hz, mesh_args...),
        permeability = 1e-6*si_unit(:darcy),
        porosity = 0.01,
        component_heat_capacity = 4.278e3*si_unit(:joule)/si_unit(:kilogram)/si_unit(:Kelvin),
    )
    mesh = physical_representation(domain)
    # Set up BTES wells
    hxy_min = metrics.hxy_min
    well_models = []
    well_names = Symbol[]
    nl = length(layers)
    geo = tpfv_geometry(mesh)

    for (wno, wc) in enumerate(well_coords_3d)
        name = Symbol("B$wno")
        println("Adding well $name ($wno/$num_wells)")
        cells = Jutul.find_enclosing_cells(mesh, permutedims(wc), n = 100)
        w_sup, w_ret = setup_btes_well(domain, cells, name=name, closed_loop_type=:u1)
        push!(well_models, w_sup, w_ret)
        push!(well_names, name)
    end

    # Make the model
    model = setup_reservoir_model(
        domain, :geothermal,
        wells = well_models,
    );

    # ## Set up initial state and boundary conditions
    geo = tpfv_geometry(mesh)
    z_bc = geo.boundary_centroids[3, :]
    bottom = map(v -> isapprox(v, maximum(z_bc)), z_bc)
    # Define pressure and temperature profiles
    rho = reservoir_model(model).system.rho_ref[1]
    dpdz = rho*gravity_constant
    dTdz = geothermal_gradient
    T = z -> temperature_surface .+ dTdz*z
    p = z -> 5atm .+ dpdz.*z
    # Set initial conditions
    z_cells = geo.cell_centroids[3, :]
    z_hat = z_cells .- minimum(z_cells)
    state0 = setup_reservoir_state(model,
        Pressure = p(z_hat),
        Temperature = T(z_hat)
    );
    # Set boundary conditions
    z_bc = z_bc[.!bottom]
    z_hat = z_bc .- minimum(z_bc)
    bc_cells = geo.boundary_neighbors[.!bottom]
    bc = flow_boundary_condition(bc_cells, domain, p(z_hat), T(z_hat));

    # Group supply well names by sector, preserving the order given in `field`
    wells_per_sector = Vector{Vector{Symbol}}()
    wtot = 0
    for sector in field
        sw = [Symbol(well_names[wtot + l], "_supply") for l in 1:length(sector)]
        wtot += length(sector)
        push!(wells_per_sector, sw)
    end
    control_charge, control_discharge, sectors = setup_controls(model, wells_per_sector,
        rate_charge, rate_discharge, temperature_charge, temperature_discharge;
        reversed_discharge = reversed_discharge);
    
    forces_charge = setup_reservoir_forces(model, control=control_charge, bc=bc)
    forces_discharge = setup_reservoir_forces(model, control=control_discharge, bc=bc);
    forces_rest = setup_reservoir_forces(model, bc=bc)
    # Make schedule
    dt, forces, timestamps = make_utes_schedule(
        forces_charge, forces_discharge, forces_rest;
        charge_period = charge_period,
        discharge_period = discharge_period,
        num_years = num_years,
        report_interval = report_interval,
        utes_schedule_args...,
    )

    # ## Useful case info
    info = Dict()
    info[:description] = "Borehole thermal energy storage (BTES) case set up using Fimbul.btes()"
    info[:sectors] = sectors
    info[:timestamps] = timestamps

    # ## Assemble and return model
    case = JutulCase(model, dt, forces, state0 = state0, input_data = info)
    return case

end

function btes(field::AbstractMatrix; kwargs...)
    return btes([[Matrix{Float64}(field)]]; kwargs...)
end

function expand_layer_values(value, num_layers::Int, name::Symbol)
    if value isa Number
        return fill(value, num_layers)
    elseif value isa AbstractVector
        if length(value) == 1
            return fill(value[1], num_layers)
        elseif length(value) == num_layers
            return collect(value)
        else
            error("Length of $name ($(length(value))) does not match number of layers ($num_layers)")
        end
    else
        error("$name must be a scalar or an AbstractVector")
    end
end

function insert_well_depths_in_layers(layer_depths, well_depths, layer_data::NamedTuple; tol = 1e-6)
    layer_depths = sort(collect(layer_depths))
    @assert length(layer_depths) >= 2 "layer_depths must contain at least two depth values"

    num_layers = length(layer_depths) - 1
    values = Dict{Symbol, Any}()
    for (name, value) in pairs(layer_data)
        values[name] = expand_layer_values(value, num_layers, name)
    end

    too_deep = sort(unique([d for d in well_depths if d > layer_depths[end] + tol]))
    if !isempty(too_deep)
        @warn "One or more wells are deeper than the deepest layer depth ($(layer_depths[end]) m)." max_well_depth = maximum(too_deep)
    end

    insertion_depths = sort(unique([d for d in well_depths if d > layer_depths[1] + tol && d < layer_depths[end] - tol]))
    for d in insertion_depths
        if any(isapprox.(layer_depths, d; atol = tol))
            continue
        end
        upper_idx = findfirst(z -> z > d, layer_depths)
        isnothing(upper_idx) && continue

        layer_idx = upper_idx - 1
        insert!(layer_depths, upper_idx, d)
        for name in keys(layer_data)
            insert!(values[name], layer_idx, values[name][layer_idx])
        end
    end

    expanded_values = (; (name => values[name] for name in keys(layer_data))...)
    return layer_depths, expanded_values
end

# ## Patterns
#
# Each pattern function takes the number of wells and the approximate
# spacing between neighboring wells, and returns a full field

function sunflower_pattern(num_wells, spacing; num_sectors, well_depth = 50.0)
    xy = fibonacci_pattern_2d(num_wells; spacing = spacing)
    return field_from_points(:angular,xy, num_sectors, well_depth)
end

function rectangular_pattern(num_wells, spacing; num_sectors, well_depth = 50.0)
    nx = max(1, round(Int, sqrt(num_wells)))
    ny = ceil(Int, num_wells/nx)
    xs = ((0:nx-1) .- (nx-1)/2).*spacing
    ys = ((0:ny-1) .- (ny-1)/2).*spacing
    xy = Matrix{Float64}(undef, 2, nx*ny)
    k = 0
    for y in ys, x in xs
        k += 1
        xy[:, k] = [x, y]
    end
    xy = xy[:, 1:num_wells]
    return field_from_points(:cartesian, xy, num_sectors, well_depth)
end

function circular_pattern(num_wells, spacing; num_sectors, well_depth = 50.0)
    points = Vector{Vector{Float64}}()
    push!(points, [0.0, 0.0])
    ring = 1
    while length(points) < num_wells
        r = ring*spacing
        n_ring = min(max(1, round(Int, 2π*r/spacing)), num_wells - length(points))
        for k in 0:n_ring-1
            θ = 2π*k/n_ring
            push!(points, [r*cos(θ), r*sin(θ)])
        end
        ring += 1
    end
    xy = hcat(points...)
    return field_from_points(:angular, xy, num_sectors, well_depth)
end

function polygonal_pattern(num_wells, spacing, num_sides; num_sectors, well_depth = 50.0)
    # Regular polygon with an area roughly matching num_wells points at the
    # given spacing
    R = sqrt(2*num_wells*spacing^2/(num_sides*sin(2π/num_sides)))
    θ = range(0.0, 2π, length = num_sides + 1)[1:num_sides]
    polygon = vcat((R*cos.(θ))', (R*sin.(θ))')

    # Build a rectangular grid covering the polygon with some margin, and
    # keep only the points that fall inside it
    margin = 1.2
    nxy = 2*ceil(Int, margin*R/spacing) + 1
    xs = ((0:nxy-1) .- (nxy-1)/2).*spacing
    xy = Matrix{Float64}(undef, 2, nxy^2)
    k = 0
    for y in xs, x in xs
        k += 1
        xy[:, k] = [x, y]
    end
    xy = xy[:, points_in_polygon(xy, polygon)]

    # Keep the num_wells points closest to the center
    r = sqrt.(xy[1,:].^2 .+ xy[2,:].^2)
    order = sortperm(r)[1:min(num_wells, size(xy, 2))]
    xy = xy[:, order]

    return field_from_points(:angular,xy, num_sectors, well_depth)
end

function setup_controls(model, wells_per_sector::AbstractVector{<:AbstractVector{Symbol}},
    rate_charge, rate_discharge, temperature_charge, temperature_discharge;
    reversed_discharge::Bool = false)

    rho = reservoir_model(model).system.rho_ref[1]
    rate_target = TotalRateTarget(rate_charge)
    ctrl_charge = InjectorControl(rate_target, [1.0],
        density=rho, temperature=temperature_charge)
    rate_target = TotalRateTarget(rate_discharge)
    ctrl_discharge = InjectorControl(rate_target, [1.0],
        density=rho, temperature=temperature_discharge);
    # BHP control for return side
    bhp_target = BottomHolePressureTarget(1.0si_unit(:atm))
    ctrl_ret = ProducerControl(bhp_target);
    # Set up forces
    control_charge = Dict()
    control_discharge = Dict()
    assigned = []
    get_return = (well) -> Symbol(replace(String(well), "_supply" => "_return"))
    sectors = Dict()

    for (sno, sw) in enumerate(wells_per_sector)
        sec_wells = Symbol[]
        for (k, well_sup) in enumerate(sw)
            well_ret = get_return(well_sup)
            @assert well_sup ∉ assigned
            @assert well_ret ∉ assigned
            if length(sw) == 1
                # Single well in sector: charge and discharge it directly
                control_charge[well_sup] = ctrl_charge
                control_discharge[well_sup] = ctrl_discharge
            else
                # Charging always runs from the first to the last well in the sector
                if k == 1
                    control_charge[well_sup] = ctrl_charge
                else
                    well_prev = get_return(sw[k-1])
                    target = JutulDarcy.ReinjectionTarget([well_prev])
                    control_charge[well_sup] = InjectorControl(target, [1.0],
                        density=rho, temperature=NaN; check=false)
                end
                # Discharging runs from first to last (reversed_discharge = false)
                # or from last to first (reversed_discharge = true)
                discharge_first = reversed_discharge ? length(sw) : 1
                if k == discharge_first
                    control_discharge[well_sup] = ctrl_discharge
                else
                    well_prev = get_return(sw[reversed_discharge ? k+1 : k-1])
                    target = JutulDarcy.ReinjectionTarget([well_prev])
                    control_discharge[well_sup] = InjectorControl(target, [1.0],
                        density=rho, temperature=NaN; check=false)
                end
            end
            control_charge[well_ret] = ctrl_ret
            control_discharge[well_ret] = ctrl_ret
            push!(assigned, well_sup, well_ret)
            push!(sec_wells, well_sup, well_ret)
        end
        sectors[Symbol("S$sno")] = sec_wells
    end

    @assert sort(assigned) == sort(well_symbols(model))

    return control_charge, control_discharge, sectors

end
