using Printf
using Dates

function gridlookup2(x0,xgrid)

    nx = size(xgrid,1)

    ix = 0
    @inbounds for jx in 1:nx
        if (x0<=xgrid[jx]); break; end;
        ix = ix+1
    end

    ix = min(max(1,ix),nx-1)

    return ix

end

function driver(r0,vmat0,mu0,Ge,Pe,knotsb,β,α,δ,znow,lnow,critin,critmu)

    ne = size(Ge,1)
    nb = size(vmat0,1)

    vmat1 = zeros(Float64,nb,ne)
    gmat0 = zeros(Float64,nb,ne)

    # Factor Price
    m0 = lnow*(r0/α/znow)^(1/(α-1))
    w0 = (1-α)*znow*m0^(α)*lnow^(-α)

    # instanteneous utility
    umat0 = zeros(nb,nb,ne)
    @inbounds for ie in 1:ne
        @inbounds for ib in 1:nb
            @inbounds for jb in 1:nb

                enow = Ge[ie]
                know = knotsb[ib]
                kp = knotsb[jb]
                cnow = w0*enow + (1.0+r0-δ)*know - kp

                if (cnow>0)
                    umat0[jb,ib,ie] = log(cnow)
                else
                    umat0[jb,ib,ie] = -Inf
                end

            end
        end
    end

    # value function iteration
    diff = 1e+4
    critin = 1e-4
    iterin = 0

    utemp = zeros(Float64,nb,nb)
    vtemp = zeros(Float64,nb,nb)
    # utemp = zeros(Float64,nb)
    # vtemp = zeros(Float64,nb)
    vcond = zeros(Float64,nb)

    while (diff>critin)

        @inbounds for ie in 1:ne

            utemp[:,:] = umat0[:,:,ie]
            vcond[:] = Pe[ie,1]*vmat0[:,1] + Pe[ie,2]*vmat0[:,2]
            vtemp[:,:] = utemp .+ β*vcond
            vmat1[:,ie] = maximum(vtemp,dims=1)

        # @inbounds for ie in 1:ne
        #
        #     @inbounds for ib in 1:nb
        #
        #         utemp[:] = umat0[:,ib,ie]
        #         vcond[:] = Pe[ie,1]*vmat0[:,1] + Pe[ie,2]*vmat0[:,2]
        #         vtemp[:] = utemp .+ β*vcond
        #         vmat1[ib,ie] = maximum(vtemp)
        #
        #     end

        end

        diff = maximum(abs.(vmat1-vmat0))
        iterin = iterin+1
        vmat0[:,:] = copy(vmat1)
#         println([iterin diff])

    end

    # policy function
    utemp = zeros(Float64,nb)
    vtemp = zeros(Float64,nb)

    @inbounds for ie in 1:ne

        @inbounds for ib in 1:nb

            utemp[:] = umat0[:,ib,ie]
            vcond[:] = Pe[ie,1]*vmat0[:,1] + Pe[ie,2]*vmat0[:,2]
            vtemp[:] = utemp .+ β*vcond
            jb = findmax(vtemp)[2] # argmax(vtemp)
            gmat0[ib,ie] = knotsb[jb]

        end

    end

    # transition matrix
    #AA = sparse(nb*ne,nb*ne);
    AA = zeros(Float64,nb*ne,nb*ne)
    wb = zeros(Float64,nb,ne)

    @inbounds for ie in 1:ne

        @inbounds for ib in 1:nb

            know = knotsb[ib]
            kp = gmat0[ib,ie]
            kb = gridlookup2(kp,knotsb)
            wb[ib,ie] = (knotsb[kb+1]-kp)/(knotsb[kb+1]-knotsb[kb])

            @inbounds for je in 1:ne

                ia = nb*(ie-1)+ib
                ja = nb*(je-1)+kb
                AA[ia,ja]   = wb[ib,ie]*Pe[ie,je]
                AA[ia,ja+1] = (1.0-wb[ib,ie])*Pe[ie,je]

            end

        end

    end

    # distribution
    diffmu = 1e+4
    dist0 = reshape(mu0,2*nb)
    dist1 = zeros(Float64,2*nb)

    while (diffmu>critmu)

        dist1[:] = AA'*dist0
        diffmu = maximum(abs.(dist1-dist0))
#         println(diffmu)
        dist0[:] = dist1/sum(dist1)

    end

    mu0[:,:] = reshape(dist0,(nb,2))

    # Calculate K
    m1 = 0.0;
    @inbounds for ie in 1:ne

        m1 = m1 + mu0[:,ie]'*knotsb
#         m1 = m1 + dot(mu0[:,ie],knotsb)

    end

    r1 = (α)*znow*m1^(α-1)*lnow^(1-α)

    return r1, vmat0, mu0

end

function main()

    β  = 0.96
    α  = 0.36
    δ  = 1-0.92
    σ  = 1.5
    μ  = 0.05
    pee = 0.925
    puu = 0.5

    # grid points
    ne = 2
    Ge = [1.0,μ]
    Pe = zeros(Float64,2,2)
    Pe[1,:] = [pee 1-pee]
    Pe[2,:] = [1-puu puu]

    nb = 1001
    kmin = 0
    kmax = 20.0
    knotsb = collect(LinRange(kmin,kmax,nb))

    mue = Pe^10000

    mue = mue[1,:]
    znow = 1.0
    # efficiency unit of labor
    lnow = Ge'*mue

    # initial distribution
    vmat0 = zeros(Float64,nb,ne)
    mu0   = ones(Float64,nb,ne)/(nb*ne)

    # initial value of r
    mnow = lnow*(α*β/(1.0-β*(1.0-δ)))^(1.0/(1.0-α))
    mnow = 5.2074 # from My_Aiyagari.m
    m0 = mnow
    r0 = (α)*znow*mnow^(α-1)*lnow^(1-α)


    start = now()
    critin = 1e-4
    critmu = 1e-8
    critout = 1e-3
    diffout = 1e+3
    damp = 0.01
    iter = 0
    maxiter = 1000

    # driver(r0,vmat0,mu0,Ge,Pe,knotsb,β,α,δ,znow,lnow,critin,critmu)

    while (diffout>critout && iter<maxiter)

        r1,vmat0,mu0 = driver(r0,vmat0,mu0,Ge,Pe,knotsb,β,α,δ,znow,lnow,critin,critmu)

        diffout = abs(log(r1)-log(r0))
        iter = iter+1
        println(@sprintf("iter = %4d, diff = %5.6f, oldr = %5.6f, newr = %5.6f",iter,diffout,r0,r1))

        # Update K
        r0 = damp*r1 + (1.0-damp)*r0

    end

    elapsed = now()-start
    println(@sprintf("Elapsed time is %s", canonicalize(Dates.CompoundPeriod(elapsed))))

end

main()
