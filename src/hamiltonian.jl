

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

struct Hamiltonian{M1<:TimeDependentOperator, M2<:AbstractMatrix}
	ht::Vector{M1}
	hc::M2
end
function Hamiltonian(hc::AbstractMatrix, ht::AbstractVector)
	if isempty(ht)
		ht2 =  Vector{TimeDependentOperator{typeof(hc)}}()
		return Hamiltonian(ht2, hc)
	end
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
Hamiltonian(hc::AbstractMatrix) = Hamiltonian(hc, [])
function Hamiltonian(ht::AbstractVector{<:TimeDependentOperator})
	isempty(ht) && throw(ArgumentError("ht can not be empty"))
	return Hamiltonian(zero(ht[1]), ht)
end 


Base.size(h::Hamiltonian, args...) = size(h.hc, args...)

Base.eltype(::Type{Hamiltonian{M1, M2}}) where {M1, M2} = eltype(M2)
Base.complex(h::Hamiltonian) = Hamiltonian(complex(h.hc), complex.(h.ht))

Base.zero(h::Hamiltonian) = zero(h.hc)

Base.:*(h::Hamiltonian, r::Number) = Hamiltonian(r .* h.hc, [item * r for item in h.ht])
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

expectation(m::AbstractMatrix, v::AbstractVector) = dot(v, m, v)
expectation(m::AbstractMatrix, v::AbstractMatrix) = tr(m * v)
function expectation(h::Hamiltonian, v::AbstractVector)
	isconstant(h) || throw(ArgumentError("expectation only applies for constant Hamiltonian"))
	return expectation(h.hc, v)
end


function matrix(h::Hamiltonian)
    isconstant(h) || error("can not convert non-constant Hamiltonian into a constant matrix")
    return h.hc
end


