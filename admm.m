function x = admm(y)

n = length(y);

% H di prova
H = [1 0 0 1 0 0; 
     0 1 0 0 1 0; 
     0 0 1 0 0 1];
[m,n] = size(H);

% init
x = zeros(n,1);
d = zeros(m,1);   % numero di check
Nv = zeros(n,1);  % numero di vicini

% parametri
mu = 4;  
alpha = 0.6;     
epsilon = 1^-5;
sigma_square = 10^-16;

%costruzione d, Nv
for j=1:m
    for i=1:n
        if H(j,i) == 1
            d(j) = d(j)+1;
            Nv(i) = Nv(i) +1;
        end
    end
end

%costruzione P
P = zeros(n, n, m); %P(x,y,z) x riga, y colonna, z la z-esima matrice

for j=1:m
    k = 1; %tiene conto della riga della Pj
    for i=1:n
        if H(j,i) == 1
            P(k, i, j) = 1;
            k = k + 1;
        end
    end
end


% log-likelihood ratio
gamma = -y/(2*sigma_square);

% seleziono il j-esimo attraverso la colonna
lambda = zeros(n,m);
z = 0.5*ones(n,m);

counter = 0;
limit = 100;

while counter<=limit
    
    % x-update
    
    %paper version
    
    t = zeros(n,1); %vettore ausilario
    for i=1:n
        for j=1:m
            t(i) = t(i) + z(i,j) - lambda(i,j)/mu ;
        end
        t(i) = t(i) - gamma(i)/mu;
    end
    for i=1:n
        x(i) = (1/(Nv(i) - 2*alpha/mu))*(t(i)-alpha/mu);
    end
    
    
    %fast version
    %x = inv(mu*(P')*P - 2*alpha)*(P'*(mu*z-lambda)-gamma-alpha);
    
    
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
        v(:,j) = P(:,:,j)*x + lambda(:,j)/mu;
        z(:,j) = ProjPolytope(v(:,j));
        lambda(:,j) = lambda(:,j) + mu*(P(:,:,j)*x - z(:,j));
    end
    
    counter = counter + 1;
    
%     % check
%     sum1 = zeros(m);
%     sum2 = zeros(m);
%     for j=1:m
%         sum1 = sum1 + (P(1:d(j),:,j)*x-z(:,j))^2;
%         sum2 = sum2 + (z(:,j)-z_old(:,j))^2;
%     end
%     if sum1<epsilon^2*m*d
%         break;
%     elseif sum2<epsilon^2*m*d
%         break;
%     end
end
end