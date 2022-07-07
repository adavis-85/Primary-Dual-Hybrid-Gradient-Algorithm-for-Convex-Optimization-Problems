using LinearAlgebra

function proj(y)

    N=length(y)
    ##Sorting Y into U
    u=sort(y,rev=true)
    x_n=zeros(N)
    K=0
         
    for i in 1:N
        check_sum=((sum(u[j] for j in 1:i)-1)/i)
        if check_sum<u[i]
            K=i
        else
            break
        end
    end
        
    τ=(sum(u[i] for i in 1:K)-1)/K
    ##Loop finds and returns max{yn-τ,0} for entire vector.  This is the probability distrubution.     
    for i in 1:N
        x_n[i]=maximum([y[i]-τ,0.0])
    end
        
    return x_n

end

function PDHG(A)

        θ=1
        
        x_L=norm(A*A')

        N=length(A[:,1])
        ##Step sizes for each incremented step.  Follow τσL^2<1 .
        σ=1
        τ=1/(x_L^2)
        ##Tolerance for possible stopping condition. 
        tol=1e-8
        check_tol=Inf
        ##Initializing the vectors
        x_n=ones(N) ./N
        y_n=ones(N) ./N
        x_bar_zero=x_n
        prob_x=0
        prob_y=0
    
        count=1
        
        ##Running algorithm 100,000 times if convergence isn't achieved
        while count<100000 
                ##Projection of Y onto the simplex to have probabilies evenly distributed.
                prob_y=proj(y_n+(A*(x_bar_zero*σ)))
                ##Projection of X onto the simplex to have probabilies evenly distributed.
                prob_x=proj(x_n-(A'*(prob_y*τ)))
                ##Incrementing for projection of Y.
                x_bar_zero=prob_x+θ*(prob_x-x_n)*τ
                ##Checking if tolerance is met during increments
                check_tol=((norm(prob_y-y_n)^2)/(2σ))+((norm(prob_x-x_n)^2)/2τ)
        
                x_n=prob_x
                y_n=prob_y
               ##View performance each 1,000 iterations
               if count%1000==0
                   println(prob_y)
                   println(prob_x)
                end
        
                count+=1
                ##Stopping condition
                if check_tol<tol
                    break
                end
                
        end
            
        return prob_x,prob_y

end

##equilibrium matrix
function eq_prob_matrix(x,y)
    
    N=length(x)
    
    mate=zeros(N,N)
    ##Multiplying probability vectors to achieve a Nash equilibrium
    ##Each player will choose the most advantageous strategy.
    for i in 1:N
        mate[:,i]=x[i]*y
    end
    
    return  mate
end
