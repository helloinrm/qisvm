using Plots
using CSV
using DataFrames
using Statistics
using Printf
using RecipesBase
using StatsPlots
using StatsBase
gr()



function paint(file::String,x::Array,work_times::Int,xname::String,savename::String,ylowest::Float64=0.4)
    data=CSV.read(file)
    times=length(x)
    mean=zeros(times)
    sqvar=zeros(times)
    for i in 1:times
        grep=data[1+(i-1)*times : work_times+(i-1)*times,8]
        mean[i]=Statistics.mean(grep)
        sqvar[i]=sqrt(var(grep))
    end
    plt=plot(x,mean-sqvar,guidefontsize=20,tickfontsize=14,line=(0,:transparent),ylims=(ylowest,1),xlabel=xname,ylabel="average success rate",legend=:no,frame=:box,fillrange=mean+sqvar,fillalpha=0.3)
    plot!(x,mean,lw=2,linecolor=:black,marker=:*,yerror=sqvar)
    display(plt)
    savefig(plt,savename)
end

# eta=[i/10 for i=1:10]
# paint("experiment_eta.csv",eta,50,"\\eta","experiment_eta.pdf",0.6)
# e=[i for i=1:10]
# paint("experiment_e.csv",e,50,"\\epsilon","experiment_e.pdf",0.6)
# mr=[i*50 for i=1:10]
# paint("experiment_mr.csv",mr,50,"max round","experiment_mr.pdf",0.6)
# ss=[i for i=1:10]
# paint("experiment_ss.csv",ss,50,"subsampling size","experiment_ss.pdf",0.6)
# rank=[i for i=1:10]
# paint("experiment_rank_e_2.csv",rank,10,"rank","experiment_rank_e_2.pdf")
# paint("experiment_rank_base.csv",rank,10,"rank","experiment_rank_base.pdf")
# paint("experiment_rank_eta_0.01.csv",rank,10,"rank","experiment_rank_eta_001.pdf")
# paint("experiment_rank_mr_10000.csv",rank,10,"rank","experiment_rank_mr_10000.pdf")
# paint("experiment_rank_ss_2.csv",rank,10,"rank","experiment_rank_ss_2.pdf")

function paint2(file::String,sorc::String,title::String)
    data=CSV.read(file)
    s1,s2,s3,c1,c2,c3=[data[:,i] for i=1:6]
    # for i=1:6
    #     println(i)
    #     @printf "%0.2f\n" mean(data[:,i])*100
    # end
    if sorc=="s"
        plt1=histogram(Any[s1,s2,s3],legend=:topleft,xlabel="success rate",ylabel="frequency",label=["qiSVM" "SVM" "LIBSVM"],bins=0:0.1:1.1, line=(0,:transparent), fillcolor=[:red :green :blue], fillalpha=0.2)
    elseif sorc=="c"
        plt1=histogram(Any[c1,c2,c3],legend=:topleft,xlabel="success rate",ylabel="frequency",label=["qiSVM" "SVM" "LIBSVM"],bins=0:0.1:1.1, line=(0,:transparent), fillcolor=[:red :green :blue], fillalpha=0.2)
    end
    #plt2=histogram(Any[c1,c2,c3],legend=:topleft,xlabel="success rate",ylabel="frequency",label=["qiSVM" "SVM" "LIBSVM"],bins=16, line=(1,:black), fillcolor=[:red :green :blue], fillalpha=0.2)
    display(plt1)
    savefig(plt1,title*".pdf")
end

function paint3(file::String,sorc::String,title::String)
    function newcounts(s::Array{Float64,2})
        n,m=size(s)
        val=zeros(5,m)
        for j=1:m
            for x=5:9
                val[x-4,j]=count(i->(x/10<i<=(x+1)/10),s[:,j])
            end
        end
        val[:]
    end
    data=CSV.read(file)
    s1,s2,s3,c1,c2,c3=[data[:,i] for i=1:6]
    if sorc=="s"
        c=newcounts([s1 s2 s3])/1000
        x=repeat(["("*string(i/10)*","*string((i+1)/10)*"]" for i=5:9],outer=3)
        g=repeat(["qiSVM", "SVM","LIBSVM"], inner = 5)
        plt=groupedbar(x,c,guidefontsize=20,tickfontsize=14,legendfontsize=14,group=g, bar_position = :dodge, bar_width=0.7,legend=:topleft,xlabel="success rate", line=(0,:transparent),ylabel="frequency",fillcolor=[:red :green :blue], fillalpha=0.2)
    elseif sorc=="c"
        c=newcounts([c1 c2 c3])/1000
        x=repeat(["("*string(i/10)*","*string((i+1)/10)*"]" for i=5:9],outer=3)
        g=repeat(["qiSVM", "SVM","LIBSVM"], inner = 5)
        plt=groupedbar(x,c,guidefontsize=20,tickfontsize=14,legendfontsize=14,group=g, bar_position = :dodge, bar_width=0.7,legend=:topleft,xlabel="success rate", line=(0,:transparent),ylabel="frequency",fillcolor=[:red :green :blue], fillalpha=0.2)
    end
    display(plt)
    savefig(plt,title*".pdf")
end



function paint4(file::String,sorc::String,title::String)
    function newcounts(s::Array{Float64,2})
        n,m=size(s)
        val=zeros(5,m)
        for j=1:m
            for x=5:9
                val[x-4,j]=count(i->(x/10<i<=(x+1)/10),s[:,j])
            end
        end
        val[:]
    end
    data=CSV.read(file)
    s1,s2,s3,c1,c2,c3=[data[:,i] for i=1:6]
    plt=plot()
    if sorc=="s"
        c=newcounts([s1 s3])/1000
        x=repeat(["("*string(i/10)*","*string((i+1)/10)*"]" for i=5:9],outer=2)
        g=repeat(["qiSVM","LIBSVM"], inner = 5)
        groupedbar!(x,c,guidefontsize=20,tickfontsize=14,legendfontsize=14,group=g, bar_position = :dodge, bar_width=0.7,legend=:topleft,xlabel="success rate", line=(0,:transparent),ylabel="frequency",fillcolor=[:red :blue], fillalpha=0.2,frame=:box,grid=(:y, :olivedrab, :dot, 1, 0.9))
    elseif sorc=="c"
        c=newcounts([c1 c3])/1000
        x=repeat(["("*string(i/10)*","*string((i+1)/10)*"]" for i=5:9],outer=2)
        g=repeat(["qiSVM","LIBSVM"], inner = 5)
        groupedbar!(x,c,guidefontsize=20,tickfontsize=14,legendfontsize=14,group=g, bar_position = :dodge, bar_width=0.7,legend=:topleft,xlabel="success rate", line=(0,:transparent),ylabel="frequency",fillcolor=[:red :blue], fillalpha=0.2,frame=:box,grid=(:y, :olivedrab, :dot, 1, 0.9))
    end
    display(plt)
    savefig(plt,title*".pdf")
end

# paint4("competition_0.csv","s","ding12")
# paint4("competition_0.csv","c","ding15")
# paint4("competition_2.csv","s","ding13")
# paint4("competition_2.csv","c","ding16")
# paint4("competition_1.csv","s","ding14")
# paint4("competition_1.csv","c","ding17")

function simpledata(file::String,x::Array,work_times::Int,savename::String)
    data=CSV.read(file)
    times=10
    me=zeros(10)
    sqvar=zeros(10)
    for i in 1:times
        grep=data[1+(i-1)*times : work_times+(i-1)*times,8]
        me[i]=Statistics.mean(grep)
        sqvar[i]=sqrt(var(grep))
    end
    data=zeros(10,3)
    data[:,1]=x
    data[:,2]=me
    data[:,3]=sqvar
    data=DataFrame(data)
    CSV.write(savename,data)
end


# eta=[i/10 for i=1:10]
# simpledata("experiment_eta.csv",eta,50,"simple_eta.csv")
# e=[i for i=1:10]
# simpledata("experiment_e.csv",e,50,"simple_e.csv")
# mr=[i*50 for i=1:10]
# simpledata("experiment_mr.csv",mr,50,"simple_mr.csv")
# ss=[i for i=1:10]
# simpledata("experiment_ss.csv",ss,50,"simple_ss.csv")
# rank=[i for i=1:10]
# simpledata("experiment_rank_base.csv",rank,10,"simple_rank_base.csv")
# simpledata("experiment_rank_e_2.csv",rank,10,"simple_rank_e.csv")
# simpledata("experiment_rank_eta_0.01.csv",rank,10,"simple_rank_eta.csv")
# simpledata("experiment_rank_mr_10000.csv",rank,10,"simple_rank_mr.csv")
# simpledata("experiment_rank_ss_2.csv",rank,10,"simple_rank_ss.csv")
function paint5()
    zu=5
    data=CSV.read("0702.csv")
    s1=zeros(8)
    s3=zeros(8)
    c1=zeros(8)
    c3=zeros(8)
    for i=1:8
        s1[i]=mean(data[1+5(i-1):5+5(i-1),1])
        s3[i]=mean(data[1+5(i-1):5+5(i-1),3])
        c1[i]=mean(data[1+5(i-1):5+5(i-1),4])
        c3[i]=mean(data[1+5(i-1):5+5(i-1),6])
        df=DataFrame([s1[i] s3[i] c1[i] c3[i]])
        if i%2==1
            CSV.write("0702_notur.csv",df,append=true)
        else
            CSV.write("0702_1tur.csv",df,append=true)
        end
    end
end
paint5()
