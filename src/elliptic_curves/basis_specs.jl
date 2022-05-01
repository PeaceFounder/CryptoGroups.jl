### X9.62 spec talks about the way optimal basis is selected as for polynomial and for gaussian normal basis. Currently, selection rule is implemented only for gaussian normal basis 


function _isprime(p)
    
    for i in 2:(p - 1)
        mod(p, i) == 0 && return false
    end

    return true
end


function gn_basis_exist(m, T)

    @assert mod(m, 8) != 0 # May as well return false?
    
    p = T*m + 1

    !_isprime(p) && return false

    k = compute_integer_order(2, p) 

    h = div(T*m, k)

    d = gcd(h, m)

    if d == 1
        return true
    else
        return false
    end

end


function gn_basis_representation_rule(m)

    @assert mod(m, 8) != 0    

    if gn_basis_exist(m, 2)
        return 2
    elseif gn_basis_exist(m, 1)
        return 1
    else
        T = 3

        while true
            if gn_basis_exist(m, T)
                return T
            end
            T += 1            
        end
    end
end




function specialize(::Type{F2GNB}, N::Int)

    @assert div(N, 8) != 0 "Out of X9.62 spec"
    T = gn_basis_representation_rule(m)

    return F2GNB{N, T}
end

function F2GNB(x::BitVector)

    N = length(x)
    F = specialize(F2GNB, N)

    return F(x)
end


