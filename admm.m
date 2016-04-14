function x = admm(y)

    n = length(y);
    % questa H di prova
    m = n;
    H = eye(n);

    % init
    x = zeros(n);
    d = zeros(n);

    % parametri
    mu = 1;
    alpha = 1;
    epsilon = 1;

    for j=1:m
        for i=1:n
           if H(j,i) == 1 
               d(j) = d(j)+1;
           end
        end
    end

    P = zeros(1,m);

    % costruisco sparsa mettendo come indice j la Pj e 
    % k come indice di riga e i di colonna
    for j=1:m
        k=0;
        P(1,j) = sparse(d(j),n);
        current = P(1,j);
        for i=1:n
            if H(j,i) == 1
                current(k,i) = 1;
                k = k+1;
            end
        end
    end

    % manca /2*sigma_w^2
    gamma = -y;

    % seleziono il j-esimo attraverso la colonna
    lambda = zeros(n,m);
    z = 0.5*ones(n,m);

    counter = 0;
    limit = 100;

    while counter<=limit

        % x-update
        x = inv(mu*(P')*P - 2*alpha)*(P'*(mu*z-lambda)-gamma-alpha);
        for i=1:n
            if x(i)<0
                x(i) = 0;
            elseif x(i)>1
                x(i) = 1;
            end
        end
        z_old = z;
        % lambda&zeta-update
        v = zeros(n,m); % seleziono il j-esimo attraverso la colonna
        for j=1:m 
           v(:,j) = P(j)*x + lambda(:,j)/mu;
           z(:,j) = ProjPolytope(v(:,j));
           lambda(:,j) = lambda(:,j) + mu*(P(j)*x - z(:,j));
        end

        counter = counter + 1;
        
        % check
        sum1 = zeros(m);
        sum2 = zeros(m);
        for j=1:m
            sum1 = sum1 + (P(j)*x-z(:,j))^2;
            sum2 = sum2 + (z(:,j)-z_old(:,j))^2;
        end
        if sum1<epsilon^2*m*d
            break;
        elseif sum2<epsilon^2*m*d
            break;
        end
    end

    return x;

end