include("../base/gas_flow.jl")
include("../base/flux_solver.jl")
include("../base/hllc.jl")
include("../base/godunov.jl")

using JSON
using StaticArrays

import Main.Godunov
import Main.GasFlow

const R_SHOCK = 1.0
const R_INFTY = 5.0

function main()
  M_E = parse(Float64, ARGS[1])
  K = parse(Float64, ARGS[2])
  A = parse(Float64, ARGS[3])
  St = parse(Float64, ARGS[4])
  t_end = parse(Float64, ARGS[5])
  cells_count = parse(Int, ARGS[6])

  r_e = sqrt(2.0 * (GasFlow.GAMMA + 1.0) / (GasFlow.GAMMA + 3.0) * K)
  r = collect(LinRange(r_e, R_INFTY, cells_count + 1))

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

  state = Godunov.State(
    r,
    [initial_flow(0.5 * (r[i+1] + r[i])) for _ in range(1, cells_count)]
  )

  function left_boundary_condition(state::Godunov.State)::SVector{3,Float64}
    return GasFlow.to_flux(flow_e)
  end

  function right_boundary_condition(state::Godunov.State)::SVector{3,Float64}
    bound_cell_flow = GasFlow.from_conservative(last(state.u))
    bound_cell_r = 0.5 * (last(state.x) + state.x[lastindex(state.x) - 1])
    density = bound_cell_flow.density
    velocity = bound_cell_flow.velocity * (bound_cell_r / last(x))^2
    pressure = K - 0.5 * density * velocity^2
    return GasFlow.to_flux(GasFlow.Params(pressure, density, velocity))
  end

  Godunov.run!(
    state,
    t_end,
    Godunov.spherical_difference_schema,
    left_boundary_condition,
    right_boundary_condition
  )

  print(json(state))
  println()
end

main()
