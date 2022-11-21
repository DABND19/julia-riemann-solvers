include("gas_flow.jl")
include("flux_solver.jl")
include("hll.jl")

import JSON
using Printf

import Main.FluxSolver
import Main.GasFlow
import Main.HLL

const CFL::Float64 = 0.2

mutable struct State
  x::Vector{Float64}
  u::Vector{Vector{Float64}}
  F::Vector{Vector{Float64}}
  t::Float64
  dt::Float64

  function State(x::Vector{Float64}, flow::Vector{GasFlow.Params})
    @assert length(x) == length(flow) + 1
    return new(
      x,
      map(GasFlow.to_conservative, flow),
      [zeros(Float64, 3) for _ in range(1, length(x))],
      0.0,
      0.0
    )
  end
end

function get_left_params(state::State, i_knot::Int)::GasFlow.Params
  @assert i_knot != firstindex(state.F)
  return GasFlow.from_conservative(state.u[i_knot-1])
end

function get_right_params(state::State, i_knot::Int)::GasFlow.Params
  @assert i_knot != lastindex(state.F)
  return GasFlow.from_conservative(state.u[i_knot])
end

function get_left_boundary_condition(state::State)::Vector{Float64}
  flow = GasFlow.from_conservative(first(state.u))
  return GasFlow.to_flux(flow)
end

function get_right_boundary_condition(state::State)::Vector{Float64}
  flow = GasFlow.from_conservative(last(state.u))
  return GasFlow.to_flux(flow)
end

function calculate_fluxes!(state::State)
  state.F[firstindex(state.F)] = get_left_boundary_condition(state)
  state.F[lastindex(state.F)] = get_right_boundary_condition(state)

  state.dt = Inf64
  for i_knot in range(firstindex(state.F) + 1, lastindex(state.F) - 1)
    left = get_left_params(state, i_knot)
    right = get_right_params(state, i_knot)

    solution = HLL.State(left, right)
    state.F[i_knot] = FluxSolver.get_flux(solution)

    S_l, S_r = FluxSolver.get_wave_velocities(solution)
    if S_l != 0.0
      dx = state.x[i_knot] - state.x[i_knot-1]
      state.dt = min(state.dt, CFL * dx / abs(S_l))
    end
    if S_r != 0.0
      dx = state.x[i_knot+1] - state.x[i_knot]
      state.dt = min(state.dt, CFL * dx / abs(S_r))
    end
  end
end

function difference_schema(state::State, i_cell::Int)::Vector{Float64}
  u = state.u[i_cell]
  F_l = state.F[i_cell]
  F_r = state.F[i_cell+1]
  dx = state.x[i_cell+1] - state.x[i_cell]
  dt = state.dt
  return u - dt / dx * (F_r - F_l)
end

function run!(state::State, t_end::Float64)
  while state.t < t_end
    calculate_fluxes!(state)
    state.dt = min(state.dt, t_end - state.t)

    @printf(stderr, "Current time: %f. Current time step: %e.\n", state.t, state.dt)

    for (i_cell, _) in enumerate(state.u)
      state.u[i_cell] = difference_schema(state, i_cell)
    end

    state.t = state.t + state.dt
  end
end

const CELLS_COUNT = 5000
const X_LEFT::Float64 = -1.0
const X_RIGHT::Float64 = 1.0
const X_DIAPH::Float64 = 0.0

function JSON.lower(state::State)
  x = [0.5 * (state.x[i+1] + state.x[i]) for (i, _) in enumerate(state.u)]
  return Dict(
    "time" => state.t,
    "x" => x,
    "solution" => map(GasFlow.from_conservative, state.u)
  )
end

function main()
  if length(ARGS) != 7
    print(stderr, "You must provide 7 positional arguments.")
  end
  left = GasFlow.Params(
    parse(Float64, ARGS[1]),
    parse(Float64, ARGS[2]),
    parse(Float64, ARGS[3])
  )
  right = GasFlow.Params(
    parse(Float64, ARGS[4]),
    parse(Float64, ARGS[5]),
    parse(Float64, ARGS[6])
  )
  t_end = parse(Float64, ARGS[7])

  x = collect(LinRange(X_LEFT, X_RIGHT, CELLS_COUNT + 1))
  initial_flow = [0.5 * (x[i] + x[i+1]) < X_DIAPH ? left : right for i in range(1, CELLS_COUNT)]
  state = State(x, initial_flow)
  run!(state, t_end)

  print(JSON.json(state))
  println()
end

main()
