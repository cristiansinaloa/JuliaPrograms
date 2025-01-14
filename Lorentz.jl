#!/usr/bin/env julia

using OrdinaryDiffEq
using Plots

function lorentz!(du,u,p,t)
    x, y, z = u
    sigma, rho, beta = p

    du[1] = sigma*(y - x)
    du[2] = x*(rho - z) - y
    du[3] = x*y - beta*z

    return du
end

u0 = [1.0,0.0,0.0]
p = [10,28,8/3]
tspan = (0.0,100.0)

prob = ODEProblem(lorentz!,u0,tspan,p)

sol = solve(prob,Tsit5())
sol = convert(Array,sol)

plt = plot(sol[1,:],sol[2,:],sol[3,:],
           title="Lorentz Equations",
           legend=false,
           denseplot=false)
display(plt)
readline()
