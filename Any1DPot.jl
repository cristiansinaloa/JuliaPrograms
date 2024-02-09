#!/usr/bin/env julia
using Pkg
Pkg.activate(".")
Pkg.add(Pkg.PackageSpec(name="Plots",version="1.39.0"))

using LinearAlgebra
using Plots

L=1 #Length of the box from 0 to L
dx = 0.001 #with of the spacing 
N = round(Int,(L-0)/dx) # num intervals = b-a/dx 1000 intvls

x = range(0,L,length=N+1) # Since we are including 0 we need N+1 points

function V(x)
    return 0*x #Infinite square well potential
    # return  10*(x .- 0.5) .^ 2 # Harmonic Potential Centered at x=0.5
    # return 0.5*exp.(-0.5*((x .- 0.7)/(.01)) .^ 2) / sqrt((0.1)^2) # Gaussian Potential at x=0.75
end

#create diags
d = 2*ones(N-1) + V(x)[2:end-1]
e = -ones(N-2)

#Create the Tridiagonal Sparse Hamiltonian
H = (1/(2*dx^2))*SymTridiagonal(d,e)

w,Psi = eigen(H)

#Normalize the wavefunctions 
for i in range(1,N-1)
   A = 1 / sqrt(sum(Psi[:,i].^2*dx))
   Psi[:,i] = A*Psi[:,i]
end

#Wavefunction Plot
plot(size=(1000,800))
plot!(title="Wavefunctions",titlefontsize=24)
plot!(xlabel="x",xguidefontsize=24,xtickfontsize=24)
plot!(ylabel="ψₙ(x)",yguidefontsize=24,ytickfontsize=24)
#plot!(xlims=(0,1))
#plot!(ylims=(-5,5))
hline!(0,color=:black)
vline!(0,color=:black)
plot!(range(dx,1-dx,step=dx),Psi[:,1],label="Ψ₁(x)",linecolor=:blue,linewidth=4)
#plot!(range(dx,1-dx,step=dx),(sqrt(2)*sin.(pi.*x[2:end-1])),label="Analytic ψ₁(x)",linecolor=:black,linewidth=4)
plot!(range(dx,1-dx,step=dx),Psi[:,2],label="Ψ₂(x)",linewidth=4,linecolor=:orange)
plot!(range(dx,1-dx,step=dx),Psi[:,3],label="Ψ₃(x)",linewidth=4,linecolor=:green)
plot!(range(dx,1-dx,step=dx),V(x)[2:end-1],linecolor=:red,linestyle=:dash,linewidth=4,label="Potential") #Plot of the potential
plot!(legend=true,legendfontsize=14,gridlinewidth=2)
display(plot!())
readline()


# savefig("plot.png")


#Eigenvalue Plot
bar(size=(1000,800))
bar!(title="H.O. Eigen Values")
bar!(xlabel="Eigen State")
bar!(ylabel="Eₙ(x)")
bar!(range(1,10),w[1:10],label="Eigen Energies")
bar!(grid=true)
display(bar!())
readline()

# savefig("plot2.png")

# Print Evals
# for i in range(1,10)
#     println(w[i])
# end
