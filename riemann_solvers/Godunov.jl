module Godunov

using Base.Threads

import JSON
using StaticArrays
import Polyester

import ..RiemannSolvers.FluxSolver
import ..RiemannSolvers.GasFlow
import ..RiemannSolvers.HLLC

const CFL::Float64 = 0.2

mutable struct State{SolverStateT<:FluxSolver.State}
  x::Vector{Float64}
  u::Vector{SVector{3,Float64}}
  F::Vector{SVector{3,Float64}}
  t::Float64
  dt::Float64
  time_steps::Vector{Float64}

  function State{SolverStateT}(
    x::Vector{Float64},
    flow::Vector{GasFlow.Params}
  ) where {SolverStateT<:FluxSolver.State}
    @assert length(x) == length(flow) + 1
    return new(
      x,
      map(GasFlow.to_conservative, flow),
      [@SVector zeros(Float64, 3) for _ in x],
      0.0,
      0.0,
      zeros(Float64, size(x)),
    )
  end
end

@inline function minmod(x::Float64, y::Float64)::Float64
  if sign(x) != sign(y)
    return 0.0
  end
  return sign(x) * min(abs(x), abs(y))
end

@inline function left_minmod(prev::Float64, current::Float64, next::Float64)::Float64
  return current + 0.5 * minmod(next - current, current - prev)
end

@inline function right_minmod(prev::Float64, current::Float64, next::Float64)::Float64
  return current - 0.5 * minmod(next - current, current - prev)
end

@inline function get_left_params(state::State, i_knot::Int)::GasFlow.Params
  @assert i_knot != firstindex(state.F)
  current = GasFlow.from_conservative(state.u[i_knot-1])
  @inbounds current_dx = 0.5 * (state.x[i_knot] + state.x[i_knot-1])

  if i_knot == firstindex(state.F) + 1 || i_knot == lastindex(state.F) - 1
    return current
  end

  prev = GasFlow.from_conservative(state.u[i_knot-2])
  @inbounds prev_dx = 0.5 * (state.x[i_knot-1] + state.x[i_knot-2])

  next = GasFlow.from_conservative(state.u[i_knot])
  @inbounds next_dx = 0.5 * (state.x[i_knot+1] + state.x[i_knot])

  return GasFlow.Params(
    current.pressure + current_dx * minmod(
      (current.pressure - prev.pressure) / (current_dx + prev_dx),
      (next.pressure - current.pressure) / (current_dx + next_dx)
    ),
    current.density + current_dx * minmod(
      (current.density - prev.density) / (current_dx + prev_dx),
      (next.density - current.density) / (current_dx + next_dx)
    ),
    current.velocity + current_dx * minmod(
      (current.velocity - prev.velocity) / (current_dx + prev_dx),
      (next.velocity - current.velocity) / (current_dx + next_dx)
    ),
  )
end

@inline function get_right_params(state::State, i_knot::Int)::GasFlow.Params
  @assert i_knot != lastindex(state.F)
  current = GasFlow.from_conservative(state.u[i_knot])
  @inbounds current_dx = 0.5 * (state.x[i_knot] + state.x[i_knot+1])

  if i_knot == firstindex(state.F) + 1 || i_knot == lastindex(state.F) - 1
    return current
  end

  prev = GasFlow.from_conservative(state.u[i_knot-1])
  @inbounds prev_dx = 0.5 * (state.x[i_knot] + state.x[i_knot-1])

  next = GasFlow.from_conservative(state.u[i_knot+1])
  @inbounds next_dx = 0.5 * (state.x[i_knot+1] + state.x[i_knot+2])

  return GasFlow.Params(
    current.pressure - current_dx * minmod(
      (current.pressure - prev.pressure) / (current_dx + prev_dx),
      (next.pressure - current.pressure) / (current_dx + next_dx)
    ),
    current.density - current_dx * minmod(
      (current.density - prev.density) / (current_dx + prev_dx),
      (next.density - current.density) / (current_dx + next_dx)
    ),
    current.velocity - current_dx * minmod(
      (current.velocity - prev.velocity) / (current_dx + prev_dx),
      (next.velocity - current.velocity) / (current_dx + next_dx)
    ),
  )
end

function calculate_fluxes!(state::State{SolverStateT}) where {SolverStateT<:FluxSolver.State}
  @inbounds Polyester.@batch per=thread threadlocal=Inf64::Float64 for i_knot in range(firstindex(state.F) + 1, lastindex(state.F) - 1)
    left = get_left_params(state, i_knot)
    right = get_right_params(state, i_knot)

    solution = SolverStateT(left, right)
    state.F[i_knot] = FluxSolver.get_flux(solution)

    S_l, S_r = FluxSolver.get_wave_velocities(solution)
    if S_l != 0.0
      dx = state.x[i_knot] - state.x[i_knot-1]
      threadlocal = min(threadlocal, CFL * dx / abs(S_l))
    end
    if S_r != 0.0
      dx = state.x[i_knot+1] - state.x[i_knot]
      threadlocal = min(threadlocal, CFL * dx / abs(S_r))
    end
  end
  state.dt = minimum(threadlocal)
end

function plain_difference_schema(state::State, i_cell::Int)::SVector{3,Float64}
  @inbounds begin
    u = state.u[i_cell]
    F_l = state.F[i_cell]
    F_r = state.F[i_cell+1]
    dx = state.x[i_cell+1] - state.x[i_cell]
    dt = state.dt
    return u - dt / dx * (F_r - F_l)
  end
end

function spherical_difference_schema(state::State, i_cell::Int)::SVector{3,Float64}
  @inbounds begin
    u = state.u[i_cell]
    F_l = state.F[i_cell]
    F_r = state.F[i_cell+1]
    dx = state.x[i_cell+1] - state.x[i_cell]
    x = 0.5 * (state.x[i_cell+1] + state.x[i_cell])
    dt = state.dt
    flow = GasFlow.from_conservative(u)
    Q = GasFlow.to_flux(flow)
    q = @SVector [0.0, flow.pressure, 0.0]
    return u - dt / dx * (F_r - F_l) + 2.0 / x * dt * (q - Q)
  end
end

function left_soft_boundary_condition(state::Godunov.State)::SVector{3,Float64}
  @inbounds flow = GasFlow.from_conservative(first(state.u))
  return GasFlow.to_flux(flow)
end

function right_soft_boundary_condition(state::Godunov.State)::SVector{3,Float64}
  @inbounds flow = GasFlow.from_conservative(last(state.u))
  return GasFlow.to_flux(flow)
end

function update_left_bound_flux!(state::State, boundary_condition::Function)
  @inbounds state.F[firstindex(state.F)] = boundary_condition(state)
end

function update_right_bound_flux!(state::State, boundary_condition::Function)
  @inbounds state.F[lastindex(state.F)] = boundary_condition(state)
end

function run!(state::State, t_end::Float64, difference_schema::Function, left_boundary_condition::Function, right_boundary_condition::Function)
  while state.t < t_end
    calculate_fluxes!(state)
    state.dt = min(state.dt, t_end - state.t)

    update_left_bound_flux!(state, left_boundary_condition)
    update_right_bound_flux!(state, right_boundary_condition)

    @inbounds Polyester.@batch per=thread for i_cell in eachindex(state.u)
      state.u[i_cell] = difference_schema(state, i_cell)
    end

    state.t = state.t + state.dt
  end
end

function JSON.lower(state::State)
  x = [0.5 * (state.x[i+1] + state.x[i]) for i in eachindex(state.u)]
  return Dict(
    "time" => state.t,
    "x" => x,
    "solution" => map(GasFlow.from_conservative, state.u)
  )
end

end
