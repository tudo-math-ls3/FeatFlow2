function testelem(x,y)

m_el=zeros(4,4);
for i=1:4,
    for j=1:4,
        m(i,j) = quadRTRZ(@(xi1,xi2) phi_e031(i,xi1,xi2) .* ...
                                   phi_e031(j,xi1,xi2) .* ...
                                   detDPhi(x,y,xi1,xi2), ...
                        -1,1,-1,1);
    end
end
m

end

function p = phi_e030(k,xi1,xi2)
% BASIS FUNCTION phi_k for parametric nonconforming integral mean value
% based rotated bilinear finite elements

switch(k)
    case{1}
        p = -3/8*(xi1.^2-xi2.^2)-1/2*xi2+1/4;
    case{2}
        p =  3/8*(xi1.^2-xi2.^2)+1/2*xi1+1/4;
    case{3}
        p = -3/8*(xi1.^2-xi2.^2)+1/2*xi2+1/4;
    case{4}
        p =  3/8*(xi1.^2-xi2.^2)-1/2*xi1+1/4;
    otherwise
        error('Invalid basis function!');
end
end

function p = phi_e031(k,xi1,xi2)
% BASIS FUNCTION phi_k for parametric nonconforming midpoint-value
% based rotated bilinear finite elements

switch(k)
    case{1}
        p = -1/4*(xi1.^2-xi2.^2)-1/2*xi2+1/4;
    case{2}
        p =  1/4*(xi1.^2-xi2.^2)+1/2*xi1+1/4;
    case{3}
        p = -1/4*(xi1.^2-xi2.^2)+1/2*xi2+1/4;
    case{4}
        p =  1/4*(xi1.^2-xi2.^2)-1/2*xi1+1/4;
    otherwise
        error('Invalid basis function!');
end
end

function d = detDPhi(x,y,xi1,xi2)
% DETERMINENT of the Jacobian of the bilinear mapping
a1 = 1/4 * ( x(1) + x(2) + x(3) + x(4));
b1 = 1/4 * ( y(1) + y(2) + y(3) + y(4));
a2 = 1/4 * (-x(1) + x(2) + x(3) - x(4));
b2 = 1/4 * (-y(1) + y(2) + y(3) - y(4));
a3 = 1/4 * (-x(1) - x(2) + x(3) + x(4));
b3 = 1/4 * (-y(1) - y(2) + y(3) + y(4));
a4 = 1/4 * ( x(1) - x(2) + x(3) - x(4));
b4 = 1/4 * ( y(1) - y(2) + y(3) - y(4));

d =  (a2+a4.*xi2).*(b3+b4.*xi1)...
    -(b2+b4.*xi2).*(a3+a4.*xi1);
end

function q = quadRTRZ(f,x1,x2,y1,y2)
% Rotated trapezoidal quadrature rule
q = 0.25 * ( f(0.5*(x1+x2),y1) + f(x2,0.5*(y1+y2)) +...
             f(0.5*(x1+x2),y2) + f(x1,0.5*(y1+y2)) );
end