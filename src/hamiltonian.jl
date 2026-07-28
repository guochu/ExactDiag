

struct TimeDependentOperator{M<:AbstractMatrix} 
	coef::Function
	op::M
end

function TimeDependentOperator(m::AbstractMatrix, f::Function)
	r = f(0)
	if (!isreal(r)) && isreal(eltype(m))
		m = complex(m)
	end
	return TimeDependentOperator(f, m)
end

Base.size(m::TimeDependentOperator, args...) = size(m.op, args...)
Base.eltype(::Type{TimeDependentOperator{M}}) where {M<:AbstractMatrix} = eltype(M)
Base.adjoint(h::TimeDependentOperator) = TimeDependentOperator(h.op', t->conj(h.coef(t)))
Base.complex(h::TimeDependentOperator) = TimeDependentOperator(h.coef, complex(h.op))
Base.isreal(h::TimeDependentOperator) = isreal(h.op)

(h::TimeDependentOperator)(t::Real) = h.op .* h.coef(t)

Base.:*(h::TimeDependentOperator, r::Number) = TimeDependentOperator(r .* h.op, h.coef)
Base.:*(r::Number, h::TimeDependentOperator) = h * r
Base.zero(h::TimeDependentOperator) = zero(h.op)

struct Hamiltonian{M1<:TimeDependentOperator, M2<:Union{AbstractMatrix, Nothing}}
	ht::Vector{M1}
	hc::M2
end
function Hamiltonian(hc::AbstractMatrix, ht::AbstractVector{<:TimeDependentOperator}) 
	for item in ht
		(size(item.op) == size(hc)) || throw(DimensionMismatch("ht hc matrix size mismatch"))
	end
	r = isreal(hc) && all(isreal, ht)
	if !r
		hc = complex(hc)
		ht = complex.(ht)
	end
	return Hamiltonian(ht, hc)
end

function Hamiltonian(hc::AbstractMatrix) 
	ht = Vector{TimeDependentOperator{typeof(hc)}}()
	return Hamiltonian(hc, ht)
end
function Hamiltonian(ht::AbstractVector{<:TimeDependentOperator})
	isempty(ht) && throw(ArgumentError("ht can not be empty"))
	return Hamiltonian(ht, nothing)
end 


function Base.size(h::Hamiltonian, args...)
	if isnothing(h.hc)
		return size(h.ht[1], args...)
	else
		return size(h.hc, args...)
	end
end 
Base.eltype(::Type{Hamiltonian{M1, M2}}) where {M1, M2} = eltype(M1)
function Base.complex(h::Hamiltonian)
	if isnothing(h.hc)
		return Hamiltonian(complex.(h.ht))
	else
		return Hamiltonian(complex(h.hc), complex.(h.ht))
	end
end
function Base.zero(h::TimeDependentOperator) 
	if isnothing(h.hc)
		return zero(h.ht[1])
	else
		return zero(h.hc)
	end
end

function Base.:*(h::Hamiltonian, r::Number) 
	if isconstant(h)
		return Hamiltonian(r .* h.hc)
	else
		if isnothing(h.hc)
			return Hamiltonian([item * r for item in h.ht])
		else
			return Hamiltonian([item * r for item in h.ht], r .* h.hc)
		end
	end
end
Base.:*(r::Number, h::Hamiltonian) = h * r

isconstant(h::Hamiltonian) = isempty(h.ht)





function (h::Hamiltonian)(t::Real)
	if isnothing(h.hc)
		m = zero(h)
	else
		m = copy(h.hc)
	end
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
	if !isnothing(h.hc)
		mul!(y, h.hc, x)
	end
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


