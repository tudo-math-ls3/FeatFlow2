function e = error(A,b,p)
    
    A_gg = full(A(1:2,1:2));
    A_go = full(A(1:2,3:end));
    A_og = full(A(3:end,1:2));
    A_oo = full(A(3:end,3:end));
    
    b_g  = full(b(1:2))';
    b_o  = full(b(3:end))';
    
    u_D  = [1 0]';
    
    % Compute matrix B and its inverse
    B = A_gg-A_go*inv(A_oo)*A_og;
    invB = inv(B);
    
    % Compute norms    
    %    e = norm(invB,p) * ( norm(b_g-A_gg*u_D,p) + ...
    %                     norm(A_go*inv(A_oo),p) * norm(b_o-A_og*u_D,p)
    %                     );
    
    e = invB*(b_g-A_go*inv(A_oo)*b_o)-u_D;
end