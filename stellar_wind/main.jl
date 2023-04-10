include("../riemann_solvers/RiemannSolvers.jl")

using ArgParse
using JSON
using StaticArrays

import .RiemannSolvers.FluxSolver
import .RiemannSolvers.Godunov
import .RiemannSolvers.GasFlow
import .RiemannSolvers.HLLC

const R_SHOCK::Float64 = 1.0
const R_RIGHT::Float64 = 5.0
const R_INFTY::Float64 = 100.0

function parse_shell_args()
  s = ArgParseSettings()

  @add_arg_table s begin
    "--A"
      help = "Solar wind oscillation amplitude"
      arg_type = Float64
      default = 0.0
    "--K"
      help = "Pressure at infinity"
      arg_type = Float64
      default = 6.76e-5
    "--St"
      help = "Strukhal number"
      arg_type = Float64
      default = 1.365
    "--Me"
      help = "Mach number of stellar wind on Earth orbit"
      arg_type = Float64
      default = 10.0
    "--nu-ph"
      help = "Photoionization frequency"
      arg_type = Float64
      default = 34.05
    "--Kn"
      help = "Knudsen number"
      arg_type = Float64
      default = 1.02e-3
    "--Nh"
      help = "Hydrogen atom concentration ratio"
      arg_type = Float64
      default = 1.5e-2
    "--t-end"
      help = "End time"
      arg_type = Float64
      required = true
  end

  return parse_args(s)
end

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
  args = parse_shell_args()

  M_E = args["Me"]
  K = args["K"]
  A = args["A"]
  St = args["St"]
  t_end = args["t-end"]
  nu_ph_e = args["nu-ph"]
  Kn = args["Kn"]
  Nh = args["Nh"]
  r_0 = 1.332

  function n_h(r::Float64)::Float64
    return Nh * exp(-r_0 / r)
  end

  function nu_ph(r::Float64)::Float64
    return K * nu_ph_e / (r^2)
  end

  function nu_ce(flow::GasFlow.Params)::Float64
    return flow.density * flow.velocity / Kn
  end

  function h_difference_schema(state::Godunov.State, i_cell::Int)::SVector{3,Float64}
    @inbounds begin
      u = state.u[i_cell]
      F_l = state.F[i_cell]
      F_r = state.F[i_cell+1]
      dx = state.x[i_cell+1] - state.x[i_cell]
      x = 0.5 * (state.x[i_cell+1] + state.x[i_cell])
      dt = state.dt
      flow = GasFlow.from_conservative(u)
      Q = GasFlow.to_flux(flow)

      nu_ph_x = nu_ph(x)
      n_h_x = n_h(x)
      freqs_sum = (nu_ph_x + nu_ce(flow))
      Q_h = @SVector [
        n_h_x * nu_ph_x,
        -n_h_x * flow.velocity * freqs_sum,
        -0.5 * n_h_x * flow.velocity^2 * freqs_sum,
      ]

      q = @SVector [0.0, flow.pressure, 0.0]
      return u - dt / dx * (F_r - F_l) + 2.0 / x * dt * (q - Q) + Q_h * dt
    end
  end

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
    [
      GasFlow.Params(K, (GasFlow.GAMMA + 1.0) / (GasFlow.GAMMA - 1.0) * r_e^2, 0.0) 
      for i in range(1, lastindex(r)-1)
    ]
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

  Godunov.run!(
    state,
    t_end,
    h_difference_schema,
    left_boundary_condition,
    Godunov.right_soft_boundary_condition
  )

  print(json(state))
  println()
end

main()
