function x = admm(gamma, H, d, Nv, P)

[m,n] = size(H);
numElem = nnz(H);

%param
mu = 3;
alpha = 0.8;

% init
x = zeros(n,1);
z = zeros(numElem,1);
lambda = zeros(numElem,1);

limit = 200;
for counter=1:limit
    % x update
    t = P'*(z-lambda/mu) - gamma/mu;
    x = (t - alpha/mu)./(Nv-2*alpha/mu);
    x = min(max(x,0),1);
 
    % z update 
    v = P*x + lambda/mu;
    z = PolytopeProjectionPacket(v,d);
    
    %lambda update
    lambda = lambda + mu*(P*x - z);
    
    %check
    if (counter>10)
        xhat = x >= 0.5;
        synd = mod(H*xhat,2);
        if (isempty(find(synd, 1)))
            break;
        end
    end
end
end