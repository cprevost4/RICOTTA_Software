function [U,V,W,G] = ricotta(Y1,Y2,Y3,P1,P2,P3,R,l1,l2,l3,alpha)


[U, ~, ~] = svds([tens2mat(Y2,1,[]) tens2mat(Y3,1,[])],R(1));
[V, ~, ~] = svds([tens2mat(Y1,2,[]) tens2mat(Y3,2,[])],R(2));
[W, ~, ~] = svds([tens2mat(Y1,3,[]) tens2mat(Y2,3,[])], R(3));


U_tilde = P1*U; V_tilde = P2*V; W_tilde = P3*W;

h = size(Y1,1); l = size(Y2,2); p = size(Y3,3);

if (R(1)>h && R(2)>l && R(3)>p)
    fprintf("Out of the identifiability region !")
    SRI_hat = NaN; info = "Non-identifiable";
elseif R(1)>R(3) || R(2)||R(3)
    A = l1*(U_tilde'*U_tilde);
    A = A + alpha*eye(size(A));
    B = eye(R(2)*R(3));
    C = eye(R(1));
    D = kron(l2*eye(R(3)),V_tilde'*V_tilde)+ kron(W_tilde'*W_tilde,l3*eye(R(2)));
    D = D + alpha*eye(size(D));
    E = l1*tmprod(Y1,{U_tilde', V', W'},[1,2,3]) + l2*tmprod(Y2,{U', V_tilde', W'},[1,2,3]) ...
        + l3*tmprod(Y3,{U', V', W_tilde'},[1,2,3]);
    E = reshape(E,R(1),R(2)*R(3));
    G = reshape(bartelsStewart(A,[],[],D,E), R);
else
    A = l1*kron(eye(R(2),U_tilde'*U_tilde)) + l2 * kron(V_tilde'*V_tilde,eye(R(1)));
    A = A + alpha*eye(size(A));
    B = eye(R(3));
    C = eye(R(1)*R(2));
    D = l3*(W_tilde'*W_tilde);
    D = D + alpha*eye(size(D));
    E = l1*tmprod(Y1,{U_tilde', V', W'},[1,2,3]) + l2*tmprod(Y2,{U', V_tilde', W'},[1,2,3]) ...
        + l3*tmprod(Y3,{U', V', W_tilde'},[1,2,3]);
    E = reshape(E,R(3),R(1)*R(2));
    G = reshape(bartelsStewart(A,[],[],D,E), R);
    
end


end

