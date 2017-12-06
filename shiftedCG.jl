module shiftedcg
    function solve(mat_A,n,b,vec_σ,N,eps=1e-12)
        x = zeros(Float64,n)
        r = zeros(Float64,n)
        p = zeros(Float64,n)
        Ap = zeros(Float64,n)
        r[:] = b[:]
        p[:] = b[:]
        αm = 1.0
        βm = 0.0
        vec_x = zeros(Float64,n,N)
        vec_r = zeros(Float64,n,N)
        vec_p = zeros(Float64,n,N)
        ρm = ones(Float64,N)
        ρ0 = ones(Float64,N)
        ρp = ones(Float64,N)
        for j in 1:N
            vec_r[:,j] = b[:] 
            vec_p[:,j] = b[:] 
        end

        k = 0
        hi = 1.0
        println(eps)
        println(hi)
        while abs(hi) > eps
            Ap = mat_A*p
            pAp = p'*Ap
            rr = r'*r
            αk = rr/pAp
            x += αk*p
            r += -αk*Ap
            βk = r'*r/rr
            p = r + βk*p
            for j in 1:N
                ρkj = ρ0[j]
                ρkmj =ρm[j]
                if abs(ρkj) > eps
                    ρp[j] = ρkj*ρkmj*αm/(ρkmj*αm*(1.0+αk*vec_σ[j])+αk*βm*(ρkmj-ρkj))
                    αkj = (ρp[j]/ρkj)*αk
                else
                    ρp[j] = ρkj
                    αkj = 0.0
                end

                vec_x[:,j] = vec_x[:,j]+αkj*vec_p[:,j]
                βkj = (ρp[j]/ρkj)^2*βk
                if abs(ρkj) > eps 
                    vec_p[:,j] = ρp[j]*r+βkj*vec_p[:,j]
                end
                #println("j=",j," ",ρp[j])
            end

            ρm[:] = ρ0[:]
            ρ0[:] = ρp[:]
            αm = αk
            βm = βk

            hi = rr
            #println(hi," ",k)
            k +=1
        end
        println("error = ",hi)
        println("Num. of iterations = ", k)

        return vec_x,x
    end



end