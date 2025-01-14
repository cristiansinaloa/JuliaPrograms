#!/usr/bin/env julia

using DifferentialEquations
using Plots

function lorentz(u,p,t)
    x, y, z = u
    sigma, rho, beta = p

    dx = sigma*(y - x)
    dy = x*(rho - z) - y
    dz = x*y - beta*z

    return [dx,dy,dz]
end

u0 = [1.0,0.0,0.0]
p = [10,28,8/3]
t = (0.0,100.0)

prob = ODEProblem(lorentz,u0,t,p)

sol = solve(prob,Tsit5())
sol = convert(Array,sol)

plt = plot(sol[1,:],sol[2,:],sol[3,:],
           title="Lorentz Equations",
           legend=false,
           denseplot=false)
display(plt)
readline()
