include("../riemann_solvers/RiemannSolvers.jl")

using JSON
using StaticArrays
using .RiemannSolvers

const R_SHOCK = 1.0
const R_RIGHT = 5.0
const R_INFTY = 100.0

function geometric_progression_space(
  x_left::Float64, 
  x_right::Float64, 
  initial_step::Float64, 
  factor::Float64
)::Vector{Float64}
  x = x_left
  step = initial_step
  res::Vector{Float64} = []
  while true
    res = append!(res, x)
    if x >= x_right
      break
    end

    step *= factor
    x = min(x + step, x_right)
  end
  return res
end 

function main()
  M_E = parse(Float64, ARGS[1])
  K = parse(Float64, ARGS[2])
  A = parse(Float64, ARGS[3])
  St = parse(Float64, ARGS[4])
  t_end = parse(Float64, ARGS[5])
  frames_count::Int = 10000

  r_e = sqrt(2.0 * (GasFlow.GAMMA + 1.0) / (GasFlow.GAMMA + 3.0) * K)
  after_termination_shock = LinRange(r_e, R_SHOCK + 0.5, 7500 + 1)
  before_temination_shock = LinRange(R_SHOCK + 0.5, R_RIGHT, 5000 + 1)
  infty = geometric_progression_space(R_RIGHT, R_INFTY, step(before_temination_shock), 1.5)
  r = unique(cat(
    collect(after_termination_shock),
    collect(before_temination_shock),
    infty;
    dims=1
  ))

  flow_e = GasFlow.Params(
    1.0 / (GasFlow.GAMMA * M_E^2),
    1.0,
    1.0
  )
  function supersonic_source(r::Float64)::GasFlow.Params
    return GasFlow.Params(
      flow_e.pressure * (r_e / r)^(2 * GasFlow.GAMMA),
      flow_e.density * (r_e / r)^2,
      flow_e.velocity
    )
  end
  function subsonic_source(r::Float64)::GasFlow.Params
    density = (GasFlow.GAMMA + 1.0) / (GasFlow.GAMMA - 1.0) * r_e^2
    velocity = (GasFlow.GAMMA - 1.0) / (GasFlow.GAMMA + 1.0) * (1.0 / r)^2
    pressure = K - 0.5 * density * velocity^2
    return GasFlow.Params(
      pressure,
      density,
      velocity
    )
  end
  function initial_flow(r::Float64)::GasFlow.Params
    if r < R_SHOCK
      return supersonic_source(r)
    else
      return subsonic_source(r)
    end
  end
  state = Godunov.State{HLLC.State}(
    r,
    # [
    #   GasFlow.Params(
    #     payload["pressure"], 
    #     payload["density"], 
    #     payload["velocity"]
    #   ) 
    #   for payload in data["solution"]
    # ]
    [initial_flow(0.5 * (r[i+1] + r[i])) for i in 1:length(r) - 1]
  )

  function left_boundary_condition(state::Godunov.State)::SVector{3,Float64}
    flow = GasFlow.Params(
      flow_e.pressure,
      flow_e.density * (1.0 + A * sin(state.t / St)),
      flow_e.velocity
    )
    return GasFlow.to_flux(flow)
  end

  function right_boundary_condition(state::Godunov.State)::SVector{3,Float64}
    bound_cell = GasFlow.from_conservative(last(state.u))
    bound_cell_r = 0.5 * (last(state.x) + state.x[lastindex(state.x) - 1])
    bound_cell_dr = last(state.x) - state.x[lastindex(state.x) - 1]
    fake_cell_r = bound_cell_r + bound_cell_dr
    fake_cell_density = bound_cell.density
    fake_cell_velocity = bound_cell.velocity * (bound_cell_r / fake_cell_r)^2
    fake_cell_pressure = K - 0.5 * fake_cell_density * fake_cell_velocity^2
    fake_cell = GasFlow.Params(
      fake_cell_pressure,
      fake_cell_density,
      fake_cell_velocity
    )
    solution = HLLC.State(bound_cell, fake_cell)
    return FluxSolver.get_flux(solution)
  end

  for t in LinRange(0.0, t_end, frames_count)
    print(stderr, t, " / ", t_end, "\n")
    Godunov.run!(
      state,
      t,
      Godunov.spherical_difference_schema,
      left_boundary_condition,
      right_boundary_condition
    )
    print(json(state))
    println()
  end
end

main()
