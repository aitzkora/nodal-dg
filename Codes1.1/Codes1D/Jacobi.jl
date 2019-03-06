"""
JacobiGQ(α::Float64,β::Float64,N::Int)

return the N'th Gauß quadrature points
"""
function JacobiGQ(α::Float64,β::Float64,N::Int)
    if (N == 0) 
        return QuadratureFormula(-(α-β)/(α+β+2.), 2.)
    end
    J = zeros(N+1)
    h₁ = 2*(0:N) .+ α .+ β
    J = diagm(0 => -1/2*(α^2-β^2)./(h₁.+2)./h₁)
    J += 
    diagm(1 => 2. ./(h₁[1:N].+2.) .* .√((1:N).*((1:N) .+ α .+ β) .* ((1:N) .+ α) .* ((1:N) .+ β) ./ (h₁[1:N] .+ 1.) ./ (h₁[1:N] .+ 3.)))
    if (α + β < 10 * eps(1.) ) 
        J[1,1] = 0.0
    end 
    J = J + J';
    F = eigen!(J)
    points = F.values
    weights = (F.vectors[1,:]').^2*2^(α+β+1)/(α+β+1)*gamma(α+1)*gamma(β+1)/gamma(β+α+1)
    return QuadratureFormula(points, weights)
end


"""
Compute the N'th order Gauß Lobatto quadrature formula points
"""
function JacobiGL(α::Float64, β::Float64, N::Int)
    x = zeros(N + 1, 1)
    if (N == 1) 
        return  [-1.,1.]
    end
    gq= JacobiGQ(α + 1, β + 1, N - 2)
    x = [-1, gq.points, 1]
end

let x = [ -8.385864502265894e-01, -5.857254804559920e-01, -2.613290131006463e-01, 9.639069017068973e-02, 4.452559470863178e-01,  7.449277954410395e-01], 
    w = [ 3.640915066793222e-02, 2.148950338971699e-01, 3.935534533149237e-01, 3.151164618818836e-01,1.085195495207685e-01, 1.100558293610547e-02]
    @test JacobiGQ(3.14, 2, 5).points ≈ x atol=1e-16
    @test JacobiGQ(3.14, 2, 5).weights ≈ w atol=1e-16
end

"""
JacobiP(x::Float64,α::Float64, β::Float64, N::Int)
Evaluate Jacobi Polynomial of type(α,β) > -1
"""
function JacobiP(x::Float64, α::Float64, β::Float64, N::Int)
    xp = copy(x)
    PL = zeros(N + 1, length(xp))
    γ₀ = 2^(α + β + 1)/(α + β + 1) * gamma(α + 1)*gamma(β + 1)/gamma(α + β + 1);
    PL[1,:] = 1.0/√gamma₀
    if (N==0) 
        return PL'
    end
    γ₁ = (α + 1)*(β + 1)/(α + β + 3) * γ₀
    PL[2,:] = ((α + β + 2) *xp/2 + (α - β)/2)/√γ₁;
    if ( N==1 ) 
      return PL(N+1, :)'
    end
# Repeat value in recurrence.
aold = 2/(2+alpha+beta)*sqrt((alpha+1)*(beta+1)/(alpha+beta+3));

% Forward recurrence using the symmetry of the recurrence.
for i=1:N-1
  h1 = 2*i+alpha+beta;
  anew = 2/(h1+2)*sqrt( (i+1)*(i+1+alpha+beta)*(i+1+alpha)*...
      (i+1+beta)/(h1+1)/(h1+3));
  bnew = - (alpha^2-beta^2)/h1/(h1+2);
  PL(i+2,:) = 1/anew*( -aold*PL(i,:) + (xp-bnew).*PL(i+1,:));
  aold =anew;
end;

P = PL(N+1,:)';
return
end
