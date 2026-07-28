


struct TimeDependentOperator{M<:AbstractMatrix}
	op::M
	coef::Function
end

Base.size(m::TimeDependentOperator, args...) = size(m.op, args...)
Base.eltype(::Type{TimeDependentOperator{M}}) where {M<:AbstractMatrix} = eltype(M)

(h::TimeDependentOperator)(t::Real) = h.op .* h.coef(t)

Base.:*(h::TimeDependentOperator, r::Number) = TimeDependentOperator(r .* h.op, h.coef)
Base.:*(r::Number, h::TimeDependentOperator) = h * r

struct Hamiltonian{M1<:AbstractMatrix, M2<:AbstractMatrix}
	ht::Vector{TimeDependentOperator{M2}}
	hc::M1
end
function Hamiltonian(hc::AbstractMatrix, ht::AbstractVector{<:TimeDependentOperator})
	for item in ht
		(size(item.op) == size(hc)) || throw(DimensionMismatch("ht hc matrix size mismatch"))
	end
	return Hamiltonian(ht, hc)
end

function Hamiltonian(hc::AbstractMatrix) 
	ht = Vector{TimeDependentOperator{typeof(hc)}}()
	return Hamiltonian(hc, ht)
end

Base.size(h::Hamiltonian, args...) = size(h.hc, args...)
Base.eltype(::Type{Hamiltonian{M1, M2}}) where {M1<:AbstractMatrix, M2<:AbstractMatrix} = promote_type(eltype(M1), eltype(M2))

function Base.:*(h::Hamiltonian, r::Number) 
	hc = r .* h.hc
	if isconstant(h)
		return Hamiltonian(hc)
	else
		return Hamiltonian([item * r for item in h.ht], hc)
	end
end
Base.:*(r::Number, h::Hamiltonian) = h * r

isconstant(h::Hamiltonian) = isempty(h.ht)





function (h::Hamiltonian)(t::Real)
	m = copy(h.hc)
	for ht in h.ht
		op, coef = ht.op, ht.coef(t)
		axpy!(coef, op, m)
		# m += tmp * item
	end	
	return Hamiltonian(m)
end

function apply!(b::AbstractVector, m::Hamiltonian, a::AbstractVector)
	isconstant(m) || error("only constant Hamiltonian support mv")
	return mul!(b, m.hc, a)
end

function apply!(y::AbstractVector, h::Hamiltonian, t::Real, x::AbstractVector)
	mul!(y, h.hc, x)
	for ht in h.ht
	    op, coef = ht.op, ht.coef(t)
	    # y += v * (item * x)
	    # axpy!(v, item * x, y)
	    mul!(y, op, x, coef, one(coef))
	end
	return y	
end

Base.:*(m::Hamiltonian, a::AbstractVector) = apply!(similar(a), m, a)
(h::Hamiltonian)(x::AbstractVector) = h * x


function (h::Hamiltonian)(t::Real, x::AbstractVector)
	return apply!(similar(x), h, t, x)
end


function matrix(h::Hamiltonian)
    isconstant(h) || error("can not convert non-constant Hamiltonian into a constant matrix")
    return h.hc
end


