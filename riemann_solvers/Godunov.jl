module Godunov

using Base.Threads

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
  @inbounds begin
    current = GasFlow.from_conservative(state.u[i_knot-1])

    if i_knot == firstindex(state.F) + 1 || i_knot == lastindex(state.F) - 1
      return current
    end
    prev = GasFlow.from_conservative(state.u[i_knot-2])
    next = GasFlow.from_conservative(state.u[i_knot])

    return GasFlow.Params(
      left_minmod(prev.pressure, current.pressure, next.pressure),
      left_minmod(prev.density, current.density, next.density),
      left_minmod(prev.velocity, current.velocity, next.velocity)
    )
  end
end

@inline function get_right_params(state::State, i_knot::Int)::GasFlow.Params
  @assert i_knot != lastindex(state.F)
  @inbounds begin
    current = GasFlow.from_conservative(state.u[i_knot])

    if i_knot == firstindex(state.F) + 1 || i_knot == lastindex(state.F) - 1
      return current
    end
    prev = GasFlow.from_conservative(state.u[i_knot-1])
    next = GasFlow.from_conservative(state.u[i_knot+1])

    return GasFlow.Params(
      right_minmod(prev.pressure, current.pressure, next.pressure),
      right_minmod(prev.density, current.density, next.density),
      right_minmod(prev.velocity, current.velocity, next.velocity)
    )
  end
end

function calculate_fluxes!(state::State)
  time_steps = zeros(Float64, size(state.F))
  fill!(time_steps, Inf64)
  @inbounds @threads for i_knot in range(firstindex(state.F) + 1, lastindex(state.F) - 1)
    left = get_left_params(state, i_knot)
    right = get_right_params(state, i_knot)

    solution = HLLC.State(left, right)
    state.F[i_knot] = FluxSolver.get_flux(solution)

    dt = Inf64
    S_l, S_r = FluxSolver.get_wave_velocities(solution)
    if S_l != 0.0
      dx = state.x[i_knot] - state.x[i_knot-1]
      dt = min(dt, CFL * dx / abs(S_l))
    end
    if S_r != 0.0
      dx = state.x[i_knot+1] - state.x[i_knot]
      dt = min(dt, CFL * dx / abs(S_r))
    end
    time_steps[i_knot] = dt
  end
  state.dt = minimum(time_steps)
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

    @threads for i_cell in eachindex(state.u)
      @inbounds state.u[i_cell] = difference_schema(state, i_cell)
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
