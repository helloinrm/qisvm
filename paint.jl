using Plots
using DelimitedFiles
using Statistics
gr()

function paint(times::Int,x::Array{Float64},y::Array{Float64},y2::Array{Float64})
    plt1=plot(x,y)
    k_bar=[ones(times,1) x] \ y
    b_bar=k_bar[1]
    k_bar=k_bar[2]
    fx=k_bar*x.+b_bar
    plt2=plot!(x,fx)
    ave=mean(y)
    ave2=mean(y2)
    println("qisvm average acc",ave)
    println("svm average acc",ave2)
    println("acc sum divided",sum(ave)/sum(ave2))
    display(plt1)
    display(plt2)
end


times=100
s1=readdlm("experiment_eta_s1.csv",',')
s2=readdlm("experiment_eta_s2.csv",',')
eta=[i/times for i=1:times]
paint(times,eta,s1,s2)
