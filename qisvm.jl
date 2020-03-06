using LinearAlgebra
using Statistics
using Dates
using LIBSVM
using Plots
using DataFrames
using CSV
using DelimitedFiles
using Notifier
using Convex
using SCS
using ECOS


struct Tector
    array::Array{Float64}
    height::Int
end

struct Tatrix
    mat::Array{Float64,2}
    col_norm_tector::Tector
    col_tector_list::Array{Tector}
end

function vec2tector(vec::Array{Float64})
    n=length(vec)
    k::Int=ceil(log2(n))
    s=2^k
    array=zeros(3s-1)
    array[2s:2s+n-1]=sign.(vec)
    array[s:s+n-1]=abs.(vec).^2
    point=s-1
    for i in (k-1):-1:0
        for j in point-2^i+1:point
            array[j]=array[2j]+array[2j+1]
        end
        point=point-2^i
    end
    return Tector(array,k)
end

function samtector(tec::Tector)
    point=1
    for i in 1:(tec.height)
        p_list=[tec.array[2point],tec.array[2point+1]]
        p_list=p_list./sum(p_list)
        x=rand()
        point=x<p_list[1] ? 2point : 2point+1
    end
    return point-2^tec.height+1, tec.array[point+2^tec.height]*sqrt(tec.array[point])
end

function mat2tatrix(mat::Array{Float64,2})
    col_norm_vec=[norm(mat[:,i]) for i in 1:size(mat)[2]]
    col_norm_tector=vec2tector(col_norm_vec)
    col_tector_list=[vec2tector(mat[:,i]) for i in 1:size(mat)[2]]
    return Tatrix(mat,col_norm_tector,col_tector_list)
end

function samtatrix1(tat::Tatrix)
    j,col_norm=samtector(tat.col_norm_tector)
    i,val=samtector(tat.col_tector_list[j])
    return i,j,col_norm,val
end

function samtatrix2(X::Tatrix,r::Int,c::Int)
    Fro_div_r=sqrt(X.col_norm_tector.array[1]/r)
    j_list=zeros(Int,r)
    vj_list=zeros(r)
    i_list=zeros(Int,c)
    vi_list=zeros(c)
    X_prime_s=fill(NaN,(size(X.mat)[1],r))# this operation is in complexity O(n*r), however it takes only half seconds allocating 10^8 NaNs, we choose to ingore it
    X_dprime=zeros(c,r)
    for j in 1:r
        j_list[j],vj_list[j]=samtector(X.col_norm_tector)
    end
    for i in 1:c
        j=rand(j_list)
        i_list[i],vi_list[i]=samtector(X.col_tector_list[j])
    end
    X_dprime=X.mat[i_list[1:c],j_list[1:r]]
    for j in 1:r
        X_dprime[:,j]=X_dprime[:,j]/vj_list[j]
    end
    X_prime_s[i_list[1:c],:]=X_dprime*Fro_div_r
    for i in 1:c
        X_dprime[i,:]=X_dprime[i,:]/norm(X_dprime[i,:])
    end
    X_dprime=X_dprime*sqrt(X.col_norm_tector.array[1]/c)
    return Fro_div_r,j_list,vj_list,X_prime_s,X_dprime
end

function trace_ab(a::Tatrix,b,xi::Float64,eta::Float64,max_round::Int)
    n1=round_control(6*log2(2/eta),max_round)
    n2=round_control(9/xi^2,max_round)
    outcome=zeros(n1)
    for p in 1:n1
        for q in 1:n2
            i,j,col_norm,val=samtatrix1(a)
            outcome[p]=outcome[p]+a.col_norm_tector.array[1]*b(j,i)/val
        end
    end
    return median(outcome)/n2
end

function dot_ab(a::Tector,b,xi::Float64,eta::Float64,max_round::Int)
    n1=round_control(6*log2(2/eta),max_round)
    n2=round_control(9/xi^2,max_round)
    outcome=zeros(n1)
    for p in 1:n1
        for q in 1:n2
            index,val=samtector(a)
            outcome[p]=outcome[p]+a.array[1]*b(index)/val
        end
    end
    return median(outcome)/n2
end

function dot_ab(a::Tector,b::Array{Float64},xi::Float64,eta::Float64,max_round::Int)
    n1=round_control(6*log2(2/eta),max_round)
    n2=round_control(9/xi^2,max_round)
    outcome=zeros(n1)
    for p in 1:n1
        for q in 1:n2
            index,val=samtector(a)
            outcome[p]=outcome[p]+a.array[1]*b[index]/val
        end
    end
    return median(outcome)/n2
end

function round_control(x::Float64,mr::Int)
    if x<0
        return 2
    elseif mr==-1||x<mr
        if x>2^63-2
            return 2^63-1
        else
            return Int(ceil(x))
        end
    else
        return mr
    end
end

function getlam(X::Tatrix,y::Array{Float64},k::Int,r::Int,e::Float64,eta::Float64,max_round::Int,Fro_div_r::Float64,j_list::Array{Int},vj_list::Array{Float64},X_prime_s::Array{Float64,2},X_dprime::Array{Float64,2})
    sig2,V=eigen(X_dprime'*X_dprime)
    #tol=sig2[end]*size(X.mat)[1]*size(X.mat)[2]*eps(Float64)
    tol=1e-10
    comp=sig2.>tol
    if k==-1
        k=sum(comp)
    end
    if r-k+1<=0
        k=r
    end
    sig2=sig2[end-k+1:end]
    V=V[:,end-k+1:end]
    lam=zeros(k)
    function Xprime(i,j)
        if isnan(X_prime_s[i,j])
            X_prime_s[i,j]=X.mat[i,j_list[j]]/vj_list[j]*Fro_div_r
        end
        return X_prime_s[i,j]
    end
    for l in 1:k
        Bforlam(i,j)=y[i]*Xprime.(j,1:r)'*V[:,l]
        lam[l]=trace_ab(X,Bforlam,3e*sig2[l]/8/sqrt(k),eta/4/k,max_round)/sig2[l]
    end
    ls=lam./sig2.^2
    return V[:,end-k+1:end]*ls
end

function getclas(X::Tatrix,y::Array{Float64},x::Array{Float64},r::Int,e::Float64,eta::Float64,max_round::Int,Fro_div_r::Float64,j_list::Array{Int},vj_list::Array{Float64},u::Array{Float64})
    times=ceil(36/e^2)*ceil(6*log2(16/eta))
    relax_parameter=10
    #in practice we don't think e1 and eta1 being that small is both effective and efficiency, so we set a relax_parameter
    e1=e/2/sqrt(r)/times*relax_parameter
    eta1=eta/8/r/times*relax_parameter
    R_s=fill(NaN,(r,size(X.mat)[2]))
    alpha_s=fill(NaN,size(X.mat)[2])
    function R(i,j)
        if isnan(R_s[i,j])
            R_s[i,j]=dot_ab(X.col_tector_list[j],X.mat[:,j_list[i]],e1,eta1,max_round)*Fro_div_r/vj_list[i]
        end
        return R_s[i,j]
    end
    function alpha(p)
        if isnan(alpha_s[p])
            alpha_s[p]=u'*R.(1:r,p)
        end
        return alpha_s[p]
    end
    j_flag=-1
    if length(size(x))==1
        num_h=1
    else
        num_h=size(x)[2]
    end
    clas=zeros(num_h)
    for h in 1:num_h
        function Bforclas(p,q)
            alpha_p=alpha(p)
            if j_flag==-1
                if abs(alpha_p)<=size(X.mat)[1]*size(X.mat)[2]*eps(Float64)
                    return 0
                else
                    j_flag=p
                end
            end
            return alpha_p*(x[q,h]-X.mat[q,j_flag])
        end
        cla=trace_ab(X,Bforclas,e/2,eta/8,max_round)
        try
            clas[h]=y[j_flag]+cla
        catch err
            if h==1 && isa(err,BoundsError)
                #println(alpha_s)
                j_flag=rand(1:size(X.mat)[2])
            end
        end
    end
    return clas
end

function svm(X::Tatrix,y::Array{Float64},x::Array{Float64},r::Int,c::Int,e::Float64,eta::Float64,max_round::Int,k::Int=-1)
    Fro_div_r,j_list,vj_list,X_prime_s,X_dprime=samtatrix2(X,r,c)
    u=getlam(X,y,k,r,e,eta,max_round,Fro_div_r,j_list,vj_list,X_prime_s,X_dprime)
    clas=getclas(X,y,x,r,e,eta,max_round,Fro_div_r,j_list,vj_list,u)
    return sign.(clas)
end

function generate(n::Int,m::Int,k::Int,turbulence::Int=0)
    #a,b,c=svd(rand(n,m))
    #b[end-k+1,end]=0
    #X=a*Diagonal(b)*c
    X=(rand(n,k).-0.5)*(rand(k,m).-0.5)
    if turbulence==1
        X=X+0.1*mean(abs.(X))*2*(rand(n,m).-0.5)
    elseif turbulence==2
        X[1:10,:]=X[1:10,:]+0.1*mean(abs.(X))*2*(rand(10,m).-0.5)
    end
    #the randomness might bring a X with rank less than k, but that is very rare case
    X=X./opnorm(X)
    flag=0
    y=ones(m)
    while flag==0
        a=rand(n).-0.5
        for i in 1:m
            if a'*X[:,i]<0
                flag=1
                y[i]=-1
            end
        end
    end
    return X,y
end

function generate(n::Int,m::Int,k::Int,test_rate::Float64,turbulence::Int=0)
    if test_rate>=1 || test_rate<=0
        error("test rate not good!")
    end
    X,y=generate(n,m,k,turbulence)
    test_m=Int(ceil(test_rate*m))
    train_m=m-Int(ceil(test_rate*m))
    flag=0
    train_col_index=rand(1:m,train_m)
    while flag==0
        for i in train_col_index
            if y[i]*y[train_col_index[1]]<0
                flag=1
            end
        end
        train_col_index=rand(1:m,train_m)
    end
    test_col_index=setdiff(1:m,train_col_index)
    train_set=X[:,train_col_index]
    test_set=X[:,test_col_index]
    train_label=y[train_col_index]
    test_label=y[test_col_index]
    return train_set,test_set,train_label,test_label
end

function svm2(X::Array{Float64,2},y::Array{Float64},x::Array{Float64})
    a=pinv(X'*X)*y
    a=a[:]
    j=rand((1:size(X)[2])[abs.(a).>size(X)[1]*size(X)[2]*eps(Float64)])
    if length(size(x))==1
        num_h=1
    else
        num_h=size(x)[2]
    end
    clas=zeros(num_h)
    for h in 1:num_h
        clas[h]=y[j]+(x[:,h]-X[:,j])'*X*a
    end
    return sign.(clas)
end

function svm3(X::Array{Float64,2},y::Array{Float64},x::Array{Float64})
    model = svmtrain(X,y,kernel=Kernel.Linear)
    #println(model)
    (predicted_labels, decision_values) = svmpredict(model, x)
    return predicted_labels
end

function svm4(X::Array{Float64,2},y::Array{Float64},x::Array{Float64})
    alpha=Variable(size(X)[2])
    problem=maximize(sum(alpha)-1/2*sumsquares(X*Diagonal(y)*alpha),[alpha>=0,alpha'*y==0])
    solve!(problem,SCSSolver())
    a=alpha.value.*y
    a=a[:]
    anothera=pinv(X'*X)*y
    println(norm(a-anothera))
    j=rand((1:size(X)[2])[abs.(a).>size(X)[1]*size(X)[2]*eps(Float64)])
    if length(size(x))==1
        num_h=1
    else
        num_h=size(x)[2]
    end
    clas=zeros(num_h)
    for h in 1:num_h
        clas[h]=y[j]+(x[:,h]-X[:,j])'*X*a
    end
    return sign.(clas)
end

function svm5(X::Array{Float64,2},y::Array{Float64},x::Array{Float64})
    alpha=Variable(size(X)[2])
    #problem=maximize(alpha'*y-1/2*sumsquares(X*alpha),[alpha'*y>=0,alpha==0])
    problem=maximize(alpha'*y-1/2*sumsquares(X*alpha))
    solve!(problem,ECOSSolver(verbose=false))
    #solve!(problem,GurobiSolver(verbose=false))
    a=alpha.value[:]
    anothera=pinv(X'*X)*y
    #println(X*a-y)
    println(norm(a-anothera))
    j=rand((1:size(X)[2])[abs.(a).>size(X)[1]*size(X)[2]*eps(Float64)])
    if length(size(x))==1
        num_h=1
    else
        num_h=size(x)[2]
    end
    clas=zeros(num_h)
    for h in 1:num_h
        clas[h]=y[j]+(x[:,h]-X[:,j])'*X*a
    end
    return sign.(clas)
end

acc(y1::Array{Float64},y2::Array{Float64})=1-sum(abs.(y1-y2))/2length(y1)

function find_good_matrix(n::Int,m::Int,k::Int)
    score=0
    X_s=zeros(n,m)
    y_s=zeros(m)
    for i in 1:10
        X,y=generate(n,m,k)
        label2=svm3(X,y,X)
        ac=acc(label2,y)
        if ac>score
            score=ac
            X_s=X
            y_s=y
        end
    end
    #writedlm("X.csv",X_s)
    #writedlm("y.csv",y_s)
    return X_s,y_s
end

function maxroundcal(r::Int,e::Float64,eta::Float64,relax_parameter::Int)
    times=ceil(36/e^2)*ceil(6*log2(16/eta))
    e1=e/2/sqrt(r)/times*relax_parameter
    eta1=eta/8/r/times*relax_parameter
    n1=6*log2(2/eta1)
    n2=9/e1^2
    return Int(ceil(max(n1,n2)))
end
#224,22303

function work(X::Array{Float64,2},y::Array{Float64},tat::Tatrix,n::Int,m::Int,k::Int,subsize::Int,e::Float64,eta::Float64,max_round::Int,work_times::Int)
    r=Int(ceil(4*log2(n/eta)/e^2))*subsize
    c=Int(ceil(4*log2(r/eta)/e^2))*subsize
    s1=zeros(work_times)
    for i in 1:work_times
        println(Dates.now())
        println("n,m,k,r,c,e,eta,max_round ",n," ",m," ",k," ",r," ",c," ",e," ",eta," ",max_round)
        label1=svm(tat,y,X,r,c,e,eta,max_round)
        s1[i]=acc(label1,y)
        println("acc of qisvm:  ",s1[i])
        df=DataFrame([n m k r c e eta max_round subsize s1[i]])
        CSV.write("data.csv",df,append=true)
    end
    return s1
end

function experiment_eta()
    n,m,k=100,100,1
    X=readdlm("X.csv")
    y=readdlm("y.csv")
    tat=mat2tatrix(X)
    e=5.0
    max_round=5000
    subsize=1
    times=10
    work_times=50
    s1=zeros(work_times)
    df=DataFrame()
    for i=1:times
        eta=i/times
        s1=work(X,y,tat,n,m,k,subsize,e,eta,max_round,work_times)
        df=vcat(df,array2df(n,m,k,subsize,e,eta,max_round,work_times,s1))
    end
    CSV.write("experiment_eta.csv",df)
end

function experiment_e()
    n,m,k=100,100,1
    X=readdlm("X.csv")
    y=readdlm("y.csv")
    tat=mat2tatrix(X)
    eta=0.1
    max_round=5000
    subsize=1
    times=10
    work_times=50
    s1=zeros(work_times)
    df=DataFrame()
    for i=1:times
        e=11-i/times*10
        s1=work(X,y,tat,n,m,k,subsize,e,eta,max_round,work_times)
        df=vcat(df,array2df(n,m,k,subsize,e,eta,max_round,work_times,s1))
    end
    CSV.write("experiment_e.csv",df)
end

function experiment_mr()
    n,m,k=100,100,1
    X=readdlm("X.csv")
    y=readdlm("y.csv")
    tat=mat2tatrix(X)
    e,eta=5.0,0.1
    subsize=1
    times=10
    work_times=50
    df=DataFrame()
    for i=1:times
        max_round=i*50
        s1=work(X,y,tat,n,m,k,subsize,e,eta,max_round,work_times)
        df=vcat(df,array2df(n,m,k,subsize,e,eta,max_round,work_times,s1))
    end
    CSV.write("experiment_mr.csv",df)
end

function experiment_ss()
    n,m,k=100,100,1
    X=readdlm("X.csv")
    y=readdlm("y.csv")
    tat=mat2tatrix(X)
    e,eta=5.0,0.1
    max_round=5000
    times=10
    work_times=50
    df=DataFrame()
    for i=1:times
        subsize=i
        s1=work(X,y,tat,n,m,k,subsize,e,eta,max_round,work_times)
        df=vcat(df,array2df(n,m,k,subsize,e,eta,max_round,work_times,s1))
    end
    CSV.write("experiment_ss.csv",df)
end

function experiment_rank(e::Float64,eta::Float64,max_round::Int,subsize::Int)
    n,m=100,100
    times=5
    work_times=10
    for k=1:times
        X,y=find_good_matrix(n,m,k)
        tat=mat2tatrix(X)
        s=work(X,y,tat,n,m,k,subsize,e,eta,max_round,work_times)
        # df=array2df(n,m,k,subsize,e,eta,max_round,work_times,s)
        # CSV.write("experiment_rank.csv",df,append=true)
    end
end

function competition(turbulence::Int=0)
    n,m=100,100
    k=1
    e,eta=5.0,0.1
    max_round=5000
    r=Int(ceil(4*log2(n/eta)/e^2))
    c=Int(ceil(4*log2(r/eta)/e^2))
    for i=1:5
        println(Dates.now())
        println("competition: ",i)
        train_set,test_set,train_label,test_label=generate(n,m,k,0.2,turbulence)
        tat=mat2tatrix(train_set)
        label1=svm(tat,train_label,test_set,r,c,e,eta,max_round)
        label2=svm2(train_set,train_label,test_set)
        label3=svm3(train_set,train_label,test_set)
        s1=acc(label1,test_label)
        s2=acc(label2,test_label)
        s3=acc(label3,test_label)
        label1=svm(tat,train_label,train_set,r,c,e,eta,max_round)
        label2=svm2(train_set,train_label,train_set)
        label3=svm3(train_set,train_label,train_set)
        c1=acc(label1,train_label)
        c2=acc(label2,train_label)
        c3=acc(label3,train_label)
        rank=LinearAlgebra.rank(train_set)
        df=DataFrame([s1 s2 s3 c1 c2 c3 rank])
        CSV.write("competition.csv",df,append=true)
        #println(s1," ",s2," ",s3," ",c1," ",c2," ",c3," ",rank)
        println(c2)
    end
end


function asamlltest()
    sco=zeros(100)
    for n=2:100
        X,y=generate(n,100,1,1)
        label1=svm2(X,y,X)
        println(acc(label1,y))
        sco[n]=acc(label1,y)
    end
    plot(sco)

    #label1=svm5(X,y,X)
    #s1=acc(label1,y)
    #println("acc of svm:  ",s1)
    #println(LinearAlgebra.rank(X))
end

function showdataset(n::Int,m::Int,k::Int)
    test_rate=0.2
    X=(rand(n,k).-0.5)*(rand(k,m).-0.5)
    X=X+0.1*mean(abs.(X))*2*(rand(n,m).-0.5)
    X=X./opnorm(X)
    flag=0
    y=ones(m)
    a=zeros(n)
    while flag==0
        a=rand(n).-0.5
        for i in 1:m
            if a'*X[:,i]<0
                flag=1
                y[i]=-1
            end
        end
    end
    test_m=Int(ceil(test_rate*m))
    train_m=m-Int(ceil(test_rate*m))
    flag=0
    train_col_index=rand(1:m,train_m)
    while flag==0
        for i in train_col_index
            if y[i]*y[train_col_index[1]]<0
                flag=1
            end
        end
        train_col_index=rand(1:m,train_m)
    end
    test_col_index=setdiff(1:m,train_col_index)
    train_set=X[:,train_col_index]
    test_set=X[:,test_col_index]
    train_label=y[train_col_index]
    test_label=y[test_col_index]



    yscale1=X'*a
    al=pinv(train_set'*train_set)*train_label
    al=al[:]
    j=rand((1:size(train_set)[2])[abs.(al).>size(train_set)[1]*size(train_set)[2]*eps(Float64)])
    w=train_set*al
    b=train_label[j]-train_set[:,j]'*w
    yscale2=X'*w.+b
    plt=plot()
    for i in test_col_index
        if yscale1[i]>0
            plot!([0,1],[yscale1[i],yscale2[i]],color=:blue,legend=:none)
        else
            plot!([0,1],[yscale1[i],yscale2[i]],color=:red)
        end
    end

    #
    # num_h=size(x)[2]
    # clas=zeros(num_h)
    # for h in 1:num_h
    #     clas[h]=y[j]+(x[:,h]-X[:,j])'*X*a
    # end
    # label=svm2(X,y,X)
    # println(acc(label,y))
    #
    # mi,ma=extrema(X[1,:])
    # x=range(mi,ma,length=100)
    # fx=-(w[1].*x.+b)./w[2]
    # plt=plot(x,fx,color=:green)
    # for i=1:100
    #     if y[i]==1
    #         scatter!((X[1,i],X[2,i]),color=:red,legend=:none)
    #     else
    #         scatter!((X[1,i],X[2,i]),color=:blue)
    #     end
    # end
    display(plt)
end

#showdataset(100,100,1)
#asamlltest()
#competition(1)
#competition(2)
#experiment_rank(5.0,0.01,5000,1)
#experiment_rank(5.0,0.1,10000,1)
#experiment_rank(5.0,0.1,5000,2)
#alarm()
#say("Finish calculation!")
#asamlltest()
