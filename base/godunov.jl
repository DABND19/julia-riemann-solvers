module Godunov

import JSON
using Printf
using StaticArrays

import Main.FluxSolver
import Main.GasFlow
import Main.HLLC

const CFL::Float64 = 0.2

mutable struct State
  x::Vector{Float64}
  u::Vector{SVector{3,Float64}}
  F::Vector{SVector{3,Float64}}
  t::Float64
  dt::Float64

  function State(x::Vector{Float64}, flow::Vector{GasFlow.Params})
    @assert length(x) == length(flow) + 1
    return new(
      x,
      map(GasFlow.to_conservative, flow),
      [@SVector zeros(Float64, 3) for _ in x],
      0.0,
      0.0
    )
  end
end

function minmod(x::Float64, y::Float64)::Float64
  if sign(x) != sign(y)
    return 0.0
  end
  return sign(x) * min(abs(x), abs(y))
end

function get_left_params(state::State, i_knot::Int)::GasFlow.Params
  @assert i_knot != firstindex(state.F)
  current = GasFlow.from_conservative(state.u[i_knot-1])

  if i_knot == firstindex(state.F) + 1 || i_knot == lastindex(state.F) - 1
    return current
  end
  prev = GasFlow.from_conservative(state.u[i_knot-2])
  next = GasFlow.from_conservative(state.u[i_knot])

  return GasFlow.Params(
    current.pressure + 0.5 * minmod(next.pressure - current.pressure, current.pressure - prev.pressure),
    current.density + 0.5 * minmod(next.density - current.density, current.density - prev.density),
    current.velocity + 0.5 * minmod(next.velocity - current.velocity, current.velocity - prev.velocity)
  )
end

function get_right_params(state::State, i_knot::Int)::GasFlow.Params
  @assert i_knot != lastindex(state.F)
  current = GasFlow.from_conservative(state.u[i_knot])

  if i_knot == firstindex(state.F) + 1 || i_knot == lastindex(state.F) - 1
    return current
  end
  prev = GasFlow.from_conservative(state.u[i_knot-1])
  next = GasFlow.from_conservative(state.u[i_knot+1])

  return GasFlow.Params(
    current.pressure - 0.5 * minmod(next.pressure - current.pressure, current.pressure - prev.pressure),
    current.density - 0.5 * minmod(next.density - current.density, current.density - prev.density),
    current.velocity - 0.5 * minmod(next.velocity - current.velocity, current.velocity - prev.velocity)
  )
end

function calculate_fluxes!(state::State)
  state.dt = Inf64
  for i_knot in range(firstindex(state.F) + 1, lastindex(state.F) - 1)
    left = get_left_params(state, i_knot)
    right = get_right_params(state, i_knot)

    solution = HLLC.State(left, right)
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

function plain_difference_schema(state::State, i_cell::Int)::SVector{3,Float64}
  u = state.u[i_cell]
  F_l = state.F[i_cell]
  F_r = state.F[i_cell+1]
  dx = state.x[i_cell+1] - state.x[i_cell]
  dt = state.dt
  return u - dt / dx * (F_r - F_l)
end

function left_soft_boundary_condition(state::Godunov.State)::SVector{3,Float64}
  flow = GasFlow.from_conservative(first(state.u))
  return GasFlow.to_flux(flow)
end

function right_soft_boundary_condition(state::Godunov.State)::SVector{3,Float64}
  flow = GasFlow.from_conservative(last(state.u))
  return GasFlow.to_flux(flow)
end

function update_left_bound_flux!(state::State, boundary_condition::Function)
  state.F[firstindex(state.F)] = boundary_condition(state)
end

function update_right_bound_flux!(state::State, boundary_condition::Function)
  state.F[lastindex(state.F)] = boundary_condition(state)
end

function run!(state::State, t_end::Float64, difference_schema::Function, left_boundary_condition::Function, right_boundary_condition::Function)
  while state.t < t_end
    calculate_fluxes!(state)
    state.dt = min(state.dt, t_end - state.t)

    update_left_bound_flux!(state, left_boundary_condition)
    update_right_bound_flux!(state, right_boundary_condition)

    for (i_cell, _) in enumerate(state.u)
      state.u[i_cell] = difference_schema(state, i_cell)
    end

    @printf(stderr, "Current time: %f. Current time step: %e.\n", state.t, state.dt)
    state.t = state.t + state.dt
  end
end

function JSON.lower(state::State)
  x = [0.5 * (state.x[i+1] + state.x[i]) for (i, _) in enumerate(state.u)]
  return Dict(
    "time" => state.t,
    "x" => x,
    "solution" => map(GasFlow.from_conservative, state.u)
  )
end

end
