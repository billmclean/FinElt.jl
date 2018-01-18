
# FIXME: This doesn't really work if x is anything but a Vector or a scalar
function fdjacobian(F, x::Number, t)
    ftx = F(t, x)

    # The 100 below is heuristic
    dx = (x .+ (x==0))./100
    dFdx = (F(t,x+dx)-ftx)./dx

    return dFdx
end

function fdjacobian(F, x, t)
    ftx = F(t, x)
    lx = max(length(x),1)
    dFdx = zeros(eltype(x), lx, lx)
    for j = 1:lx
        # The 100 below is heuristic
        dx = zeros(eltype(x), lx)
        dx[j] = (x[j] .+ (x[j]==0))./100
        dFdx[:,j] = (F(t,x+dx)-ftx)./dx[j]
    end
    return dFdx
end

# ODE23S  Solve stiff systems based on a modified Rosenbrock triple
# (also used by MATLABS ODE23s); see Sec. 4.1 in
#
# [SR97] L.F. Shampine and M.W. Reichelt: "The MATLAB ODE Suite," SIAM Journal on Scientific Computing, Vol. 18, 1997, pp. 1–22
#
# supports keywords: points = :all | :specified (using dense output)
#                    jacobian = G(t,y)::Function | nothing (FD)
#                    mass = M::Matrix | nothing (I)
function ode23s(F, y0, tspan; reltol = 1.0e-5, abstol = 1.0e-8,
                                                jacobian=nothing,
                                                mass=nothing,
                                                points=:all,
                                                norm=Base.norm)

    # select method for computing the Jacobian
    if typeof(jacobian) == Function
        jac = jacobian
    else
        # fallback finite-difference
        jac = (t, y)->fdjacobian(F, y, t)
    end
    if mass == nothing
        M = I
    else
        M = mass
    end

    # constants
    const d = 1/(2 + sqrt(2))
    const e32 = 6 + sqrt(2)


    # initialization
    t = tspan[1]
    tfinal = tspan[end]
    tdir = sign(tfinal - t)

    hmax = abs(tfinal - t)/10
    hmin = abs(tfinal - t)/1e9
    h = tdir*abs(tfinal - t)/100  # initial step size

    y = y0
    tout = Array{typeof(t)}(1)
    tout[1] = t         # first output time
    yout = Array{typeof(y0)}(1)
    yout[1] = copy(y)         # first output solution

    F0 = F(t,y)

    J = jac(t,y)    # get Jacobian of F wrt y

    while abs(t) < abs(tfinal) && hmin < abs(h)
        if abs(t-tfinal) < abs(h)
            h = tfinal - t
        end

        if size(J,1) == 1
            W = M - h*d*J
        else
            # (see Sec. 5 in [SR97])
            W = lufact( M - h*d*J )
        end

        # approximate time-derivative of F
        T = h*d*(F(t + h/100, y) - F0)/(h/100)

        # modified Rosenbrock formula
        k1 = W\(F0 + T)
        F1 = F(t + 0.5*h, y + 0.5*h*k1)
        k2 = W\(F1 - M*k1) + k1
        ynew = y + h*k2
        F2 = F(t + h, ynew)
        k3 = W\(F2 - e32*(M*k2 - F1) - 2*(M*k1 - F0) + T )

        err = (h/6)*norm(k1 - 2*k2 + k3) # error estimate
        delta = max(reltol*max(norm(y),norm(ynew)), abstol) # allowable error

        # check if new solution is acceptable
        if  err <= delta

            if points==:specified || points==:all
                # only points in tspan are requested
                # -> find relevant points in (t,t+h]
                for toi in tspan[(tspan.>t) .& (tspan.<=t+h)]
                    # rescale to (0,1]
                    s = (toi-t)/h

                    # use interpolation formula to get solutions at t=toi
                    push!(tout, toi)
                    push!(yout, y + h*( k1*s*(1-s)/(1-2*d) + k2*s*(s-2*d)/(1-2*d)))
                end
            end
            if (points==:all) && (tout[end]!=t+h)
                # add the intermediate points
                push!(tout, t + h)
                push!(yout, ynew)
            end

            # update solution
            t = t + h
            y = ynew

            F0 = F2         # use FSAL property
            J = jac(t,y)    # get Jacobian of F wrt y
                            # for new solution
        end

        # update of the step size
        h = tdir*min( hmax, abs(h)*0.8*(delta/err)^(1/3) )
    end

    return tout, yout
end


