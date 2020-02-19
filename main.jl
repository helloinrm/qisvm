include("qisvm.jl")
using LinearAlgebra
using Statistics
using Plots
using Dates
using Notifier
using DelimitedFiles
#gr()

function experiment_eta()
    n,m,k=100,100,1
    e=5.0
    max_round=5000
    subsize=1
    times=10
    times2=10
    s1=zeros(times)
    s2=zeros(times)
    for i=1:times
        eta=i/times
        for j in 1:times2
            s1[i],s2[i]+=work(n,m,k,subsize,e,eta,max_round)
        end
    end
    writedlm("experiment_eta_s1.csv",s1./times2,',')
    writedlm("experiment_eta_s2.csv",s2./times2,',')
end

function experiment_e()
    n,m,k=100,100,1
    eta=0.1
    max_round=5000
    subsize=1
    times=19
    s1=zeros(times)
    s2=zeros(times)
    for i=2:20
        e=i/2
        s1[i],s2[i]=work(n,m,k,subsize,e,eta,max_round)
    end
    writedlm("experiment_e_s1.csv",s1,',')
    writedlm("experiment_e_s2.csv",s2,',')
end

function experiment_mr()
    n,m,k=100,100,1
    e,eta=3.0,0.1
    subsize=1
    times=100
    s1=zeros(times)
    s2=zeros(times)
    for i=1:times
        max_round=i*100
        s1[i],s2[i]=work(n,m,k,subsize,e,eta,max_round)
    end
    writedlm("experiment_mr_s1.csv",s1,',')
    writedlm("experiment_mr_s2.csv",s2,',')
end

function experiment_ss()
    n,m,k=100,100,1
    e,eta=5.0,0.1
    max_round=5000
    times=100
    s1=zeros(times)
    s2=zeros(times)
    for subsize=1:times
        s1[subsize],s2[subsize]=work(n,m,k,subsize,e,eta,max_round)
    end
    writedlm("experiment_ss_s1.csv",s1,',')
    writedlm("experiment_ss_s2.csv",s2,',')
end

function experiment_rank()
    n,m=2,2
    e,eta=0.5,0.1
    max_round=5000
    subsize=1
    times=m
    s1=zeros(times)
    s2=zeros(times)
    for k=1:times
        s1[i],s2[i]=work(n,m,k,subsize,e,eta,max_round)
    end
    writedlm("experiment_rank_s1.csv",s1,',')
    writedlm("experiment_rank_s2.csv",s2,',')
end

function experiment_time()
    k=1
    e,eta=1.0,0.5
    max_round=5000
    times=100
    t=zeros(times)
    for i=1:times
        n=m=5*times
        r=Int(ceil(4*log2(n/eta)/e^2))
        c=Int(ceil(4*log2(r/eta)/e^2))
        println("n,m,k,r,c,e,eta,max_round ",n," ",m," ",k," ",r," ",c," ",e," ",eta," ",max_round)
        X,y=generate(n,m,k)
        tat=mat2tatrix(X)
        ranj=rand(1:m)
        t[i]=@elapsed svm(tat,y,X[:,ranj],r,c,e,eta,max_round)
        println("running time:  ",t[i])
    end
    writedlm("experiment_time_t.csv",t,',')
end




experiment_eta()
#experiment_e()
#experiment_mr()
#experiment_ss()
notify("Task completed")
