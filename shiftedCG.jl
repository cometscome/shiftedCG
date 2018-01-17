module shiftedcg
    function solve(mat_A,n,b,vec_σ,N,eps=1e-12)
        println("--------------------------------------------------------")
        println("Solving linear equations with the shifted CG method...")
        tb = typeof(b[1]) 
        x = zeros(tb,n)
        r = zeros(tb,n)
        p = zeros(tb,n)
        Ap = zeros(tb,n)
        r[:] = b[:]
        p[:] = b[:]
        αm = 1.0
        βm = 0.0
        tσ = typeof(vec_σ[1])
        vec_x = zeros(tσ,n,N)
        vec_r = zeros(tσ,n,N)
        vec_p = zeros(tσ,n,N)
        ρm = ones(tσ,N)
        ρ0 = ones(tσ,N)
        ρp = ones(tσ,N)
        for j in 1:N
            vec_r[:,j] = b[:] 
            vec_p[:,j] = b[:] 
        end

        k = 0
        hi = 1.0
#        println(eps)
#        println(hi)
        while abs(hi) > eps
            A_mul_B!(Ap,mat_A,p)
            #Ap = mat_A*p
            pAp = p'*Ap
            rr = r'*r
            αk = rr/pAp

            for i=1:n
                x[i] = αk*p[i]
                r[i] += -αk*Ap[i]
            end

            #x += αk*p
            #r += -αk*Ap
            βk = r'*r/rr

            for i=1:n
                p[i] = r[i] + βk*p[i]
            end
            #p = r + βk*p
            for j in 1:N
                ρkj = ρ0[j]
                ρkmj =ρm[j]
                ρp[j] = ifelse(abs(ρkj) > eps,ρkj*ρkmj*αm/(ρkmj*αm*(1.0+αk*vec_σ[j])+αk*βm*(ρkmj-ρkj)),ρkj)
                αkj = ifelse(abs(ρkj) > eps,(ρp[j]/ρkj)*αk,0.0)

                for i=1:n
                    vec_x[i,j] += αkj*vec_p[i,j]
                end
                
                βkj = (ρp[j]/ρkj)^2*βk
                if abs(ρkj) > eps
                    for i=1:n
                        vec_p[i,j] = ρp[j]*r[i]+βkj*vec_p[i,j]
                    end   
                end
                

            end

            for i=1:N
                ρm[i] = ρ0[i] 
                ρ0[i] = ρp[i]
            end  

            
            
            αm = αk
            βm = βk
    
            ρMAX = maximum(abs.(ρp))^2
            hi = rr*ρMAX
            #println(hi," ",k)
            k +=1
            if k > 10*n
                println(hi," ",k)
                println("not converged! there might be zero eigenvalue")
                conv = false
                return vec_x,x,conv
            end
        end
        println("done.")
        println("error = ",hi)
        println("Num. of iterations = ", k)
        conv = true

        println("--------------------------------------------------------")
        

        return vec_x,x,conv
    end



end