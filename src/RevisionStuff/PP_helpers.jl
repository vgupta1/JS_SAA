##Some helpers for hte predictive/prescrive experiment

#Vector version of zbar with multiple p0s
#assumes lambdas = 1
function zbar_vec(xs, cs, p0s, alpha, mhats, ps)
    K = size(cs, 2)
    out = 0.
    for k = 1:K
        out += JS.z_k(xs[k], cs[:, k], p0s[:, k], alpha, mhats[:, k], ps[:, k], 1, 1)
    end
    out/K
end

#Vector version of oracle with multiple anchors
#assumes lambdas = 1
function oracle_alpha_vec(xs, cs, mhats, ps, p0s, alpha_grid)
    alphaOR = 0.
    jstar = -1
    best_val = Inf
    out = zeros(length(alpha_grid))
    for (j, alpha) in enumerate(alpha_grid)
        out[j] = zbar_vec(xs, cs, p0s, alpha, mhats, ps)
        if out[j] < best_val
            jstar = j
            alphaOR = alpha
            best_val = out[j]
        end
    end
    return alphaOR, jstar, out
end

#Oracle that optimizes over alpha and h
#assumes lambdas = 1
#fun_anchor(phats, Xs, h) -> anchor matrix
#
function oracle_alpha_h(xs, cs, mhats, Xs, ps, fun_anchor, alpha_grid, h_grid)
    alphaOR = 0.
    hOR = 0.
    jstar = -1
    kstar = -1
    best_val = Inf
    out = zeros(length(alpha_grid), length(h_grid))
    for (k, h) in enumerate(h_grid)
        p0s = fun_anchor(mhats ./ sum(mhats, dims=1), Xs, h)
        for (j, alpha) in enumerate(alpha_grid)
            out[j, k] = zbar_vec(xs, cs, p0s, alpha, mhats, ps)
            if out[j, k] < best_val
                jstar = j
                kstar = k
                alphaOR = alpha
                hOR = h
                best_val = out[j, k]
            end
        end
    end
    return alphaOR, jstar, kstar, out
end




###Functions to compute the predictive/prescriptivie anchors.  
k(x, y, h) = exp( -LinearAlgebra.norm(x-y)^2 / h^2) / h
function get_NW_weights(xs, Xs, h)
    out = map(ix -> k(xs, Xs[:, ix], h), 1:size(Xs, 2))
    out ./ sum(out)
end

function get_anchors(phats, Xs, h)
    K = size(Xs, 2)
    ws = zeros(K)
    p0s = zeros(size(phats))
    for k = 1:K
        ws[:] = get_NW_weights(Xs[:, k], Xs, h)
        p0s[:, k] = phats * ws
    end
    p0s
end

#an oracle tuning of h
function tune_h_mse(phats, ps, Xs, h_grid)
    p0s = zeros(size(phats))
    err(h) = maximum(abs.(get_anchors(phats, Xs, h) .- ps))
    out = map(err, h_grid)
    hstar = h_grid[argmin(out)]
    return hstar, out, get_anchors(phats, Xs, hstar)
end

####
# High dimensional variants
####
#kernel using high-dimensional featbar
#idx is the unit labeling of x and y (categorical)
function k(x, idx, y, idy, h)
    t = LinearAlgebra.norm(x - y)^2
    if idx != idy
        t += 2 
    end
    exp( -t / h^2 ) / h
end
function get_NW_weights(xs, idx, Xs, h)
    out = map(ix -> k(xs, idx, Xs[:, ix], ix, h), 1:size(Xs, 2))
    out ./ sum(out)
end

function get_anchors_bar(phats, Xs, h)
    K = size(Xs, 2)
    ws = zeros(K)
    p0s = zeros(size(phats))
    for k = 1:K
        ws[:] = get_NW_weights(Xs[:, k], k, Xs, h)
        p0s[:, k] = phats * ws
    end
    p0s
end

#an oracle tuning of h
function tune_h_mse_bar(phats, ps, Xs, h_grid)
    p0s = zeros(size(phats))
    err(h) = maximum(abs.(get_anchors_bar(phats, Xs, h) .- ps))
    out = map(err, h_grid)
    hstar = h_grid[argmin(out)]
    return hstar, out, get_anchors_bar(phats, Xs, hstar)
end

