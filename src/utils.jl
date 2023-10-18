function devectorizel2normsq(x, y::Vector{T})::T where T <: AbstractFloat
    
    out = zero(T)
    for i in eachindex(x)
        out += (x[i] - y[i])^2
    end
    return out
end

function devectorizel2norm(x, y::Vector{T})::T where T <: AbstractFloat
    return sqrt(devectorizel2norm(x,y))
end

# given a length r, find the radius of the D-dimensional ball with radius h, such that:
# such that sqrt(r^2 + r^2 + ... + r^2) = sqrt(h). The LHS sum has D times.
function getradius(r::T, D::Int)::T where T
    
    # based on:
    # h^2 = D*r^2 # strive to get the D-dim distance with the given 1D length.
    
    h = sqrt(D*r^2)
    return h
end