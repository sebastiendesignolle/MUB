using Combinatorics
using LinearAlgebra
using Nemo
using Primes

# Construction of the standard complete set of MUBs
# The dimension d can be any integer greater than two
# The output contains min_i p_i^r_i+1 bases where d = p_1^r_1*...*p_n^r_n
# Reference: arXiv:1004.3348
# Contact sebastien.designolle@gmail.com for questions
function mub(d; T=Float64)
    # Auxiliary function to compute the trace in finite fields
    function tr_ff(a)
        parse(Int, string(tr(a)))
    end
    f = collect(Primes.factor(d))
    p = f[1][1]
    r = f[1][2]
    if length(f) > 1
        B_aux1 = mub(p^r; T=T)
        B_aux2 = mub(d÷p^r; T=T)
        k = min(size(B_aux1, 3), size(B_aux2, 3))
        B = zeros(Complex{T}, d, d, k)
        for j in 1:k
            B[:, :, j] = kron(B_aux1[:, :, j], B_aux2[:, :, j])
        end
    else
        B = zeros(Complex{T}, d, d, d+1)
        B[:, :, 1] = Matrix(I, d, d)
        #  f, x = finite_field(p, r, "x") # syntax for newer versions of Nemo
        f, x = FiniteField(p, r, "x") # syntax for older versions of Nemo
        pow = [x^i for i in 0:r-1]
        el = [sum(digits(i; base=p, pad=r) .* pow) for i in 0:d-1]
        if p == 2
            for i in 1:d, k in 0:d-1, q in 0:d-1
                aux = one(Complex{T})
                q_bin = digits(q; base=2, pad=r)
                for m in 0:r-1, n in 0:r-1
                    aux *= conj(im^tr_ff(el[i]*el[q_bin[m+1]*2^m+1]*el[q_bin[n+1]*2^n+1]))
                end
                B[:, k+1, i+1] += (-1)^tr_ff(el[q+1]*el[k+1]) * aux * B[:, q+1, 1] / √T(d)
            end
        else
            γ = exp(2*im*T(π)/p)
            inv_two = inv(2*one(f))
            for i in 1:d, k in 0:d-1, q in 0:d-1
                B[:, k+1, i+1] += γ^tr_ff(-el[q+1]*el[k+1]) * γ^tr_ff(el[i]*el[q+1]*el[q+1]*inv_two) * B[:, q+1, 1] / √T(d)
            end
        end
    end
    return B
end

# Select a specific subset with k bases
function mub(d::Int, k::Int, s::Int)
    B = mub(d)
    subs = collect(combinations(1:size(B, 3), k))
    sub = subs[s]
    return B[:, :, sub]
end

# Select the first subset with k bases
function mub(d::Int, k::Int)
    return mub(d, k, 1)
end

# Check whether the input is indeed mutually unbiased
function is_mu(B::Array{Complex{T}, 3}; tol=Base.rtoldefault(T)) where {T <: AbstractFloat}
    d = size(B, 1)
    k = size(B, 3)
    for x in 1:k, y in x:k, a in 1:d, b in 1:d
        if x == y
            aux = T(a == b)
        else
            aux = 1 / √T(d)
        end
        if abs(dot(B[:, a, x], B[:, b, y])) - aux > tol
            #  println([x, y, a, b])
            #  println(abs(dot(B[:, a, x], B[:, b, y])) - aux)
            return false
        end
    end
    return true
end
