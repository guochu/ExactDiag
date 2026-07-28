

struct Lindbladian{M1<:TimeDependentOperator, M2<:AbstractMatrix, M3<:AbstractMatrix, M4<:AbstractMatrix}
	ht::Vector{M1}
	hc::M2
	dissipators::Vector{M3}
	workspace::M4
end

Base.eltype(::Type{Lindbladian{M1, M2, M3, M4}}) where {M1, M2, M3, M4} = eltype(M4)
isconstant(h::Lindbladian) = isempty(h.ht)

function Lindbladian(hc::AbstractMatrix, ht::AbstractVector{<:TimeDependentOperator}, dissipators::Vector) 
	for item in ht
		(size(item.op) == size(hc)) || throw(DimensionMismatch("ht hc matrix size mismatch"))
	end
	for jump in dissipators
		(size(jump) == size(hc)) || throw(DimensionMismatch("jump operator size mismatch with Hamiltonian size"))
	end
	# ht = Vector{TimeDependentOperator{typeof(hc)}}()

	m = -im*hc
	for jump in dissipators
		(size(jump) == size(hc)) || throw(DimensionMismatch("jump operator size mismatch with Hamiltonian size"))
	    m -= jump' * jump
	end
	workspace = zeros(eltype(m), size(m))
	return Lindbladian(complex.(ht), m, complex.(dissipators), workspace)
end

function Lindbladian(hc::AbstractMatrix, dissipators::Vector)
	ht = Vector{TimeDependentOperator{typeof(hc)}}()
	return Lindbladian(hc, ht, dissipators)
end
function Lindbladian(ht::AbstractVector{<:TimeDependentOperator}, dissipators::Vector)
	isempty(ht) && throw(ArgumentError("ht can not be empty"))
	return Lindbladian(zero(ht[1]), ht, dissipators)
end
Lindbladian(h::Hamiltonian, dissipators::Vector) = Lindbladian(h.hc, h.ht, dissipators)

Base.size(s::Lindbladian, i::Int) = begin
    l = size(s.hc, i)
    return l * l
end
Base.size(s::Lindbladian) = size(s, 1), size(s, 2)


# function Lindbladian(h::Hamiltonian, dissipators::Vector) 
# 	m = -im*h.hc
# 	for jump in dissipators
# 		(size(jump) == size(h)) || throw(DimensionMismatch("jump operator size mismatch with Hamiltonian size"))
# 	    m -= jump' * jump
# 	end
# 	workspace = zeros(eltype(m), size(m))
# 	return Lindbladian(h.ht, m, dissipators, workspace)
# end
Lindbladian(h::Hamiltonian) = Lindbladian(h, [])

# function (m::Lindbladian)(t::Real)
# 	h2 = Hamiltonian(m.hc, m.ht)
# 	h2c = h2(t)
# 	return Lindbladian(h2c.ht, h2c.hc, m.dissipators, m.workspace)
# end
function (h::Lindbladian)(t::Real)
	m = copy(h.hc)
	for ht in h.ht
		op, coef = ht.op, ht.coef(t)
		axpy!(-im*coef, op, m)
		# m += tmp * item
	end	
	ht = typeof(h.ht)()
	return Lindbladian(ht, m, h.dissipators, h.workspace)
end


function apply!(b::AbstractMatrix, m::Lindbladian, a::AbstractMatrix)
	isconstant(m) || error("only constant Lindbladian supported")
	hc = m.hc
	mul!(b, hc, a)
	# mul!(b, hc, a)
	mul!(b, a, hc', 1., 1.)
	# b =  * a + a * get_hc(m)'
	workspace = m.workspace
	for item in m.dissipators
	    # b += 2. * (item * a * item')
	    mul!(workspace, item, a)
	    mul!(b, workspace, item', 2., 1.)
	end
	return b
end

function apply!(b::AbstractMatrix, m::Lindbladian, t::Real, a::AbstractMatrix)
	hc = m.hc
	mul!(b, hc, a)
	mul!(b, a, hc', 1., 1.)
	for ht in m.ht
		item, coef = ht.op, ht.coef(t)
	    v = -im * coef
	    # b += v * (item * a - a * item)
	    mul!(b, item, a, v, 1.)
	    mul!(b, a, item, -v, 1.)
	end
	workspace = m.workspace
	for item in m.dissipators
	    # b += 2. * (item * a * item')
	    mul!(workspace, item, a)
	    mul!(b, workspace, item', 2., 1.)	    
	end
	return b
end

Base.:*(m::Lindbladian, a::AbstractMatrix) = apply!(similar(a), m, a)
(m::Lindbladian)(a::AbstractMatrix) = m * a

function (m::Lindbladian)(t::Real, a::AbstractMatrix)
	return apply!(similar(a), m, t, a)	
end

function matrix(h::Lindbladian)
    isconstant(h) || error("can not convert non-constant Lindbladian into a constant matrix")
    hc = h.hc
    iden = one(hc)
    b = kron(iden, hc) + kron(conj(hc), iden)
	for item in h.dissipators
	    b .+= 2. .* kron(conj(item), item) 
	end
	return b   
end

