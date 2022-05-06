function [Zn] = LSGMC(X,lambda,mu,rho,p)
% lambda: λ
V = size(X,2);
n = size(X{1,1},2);
mu_max = 1e8;
iter_max = 30;
%% initialization
C = zeros(n,n);
L = cell(V,1);R = cell(V,1);
Z = cell(V,1);J = cell(V,1);
E = cell(V,1);
Q1 = cell(V,1);Q2 = cell(V,1);Q3 = cell(V,1);
tempX = cell(V,1);
for v = 1:V
    R{v}=zeros(n,n);
    Z{v}=zeros(n,n);J{v}=zeros(n,n);
    Q1{v} = zeros(size(X{v}));Q2{v} = zeros(n,n);Q3{v} = zeros(n,n);
    tempX{v} = X{v}'*X{v};
end
%% optimization
oldstop=zeros(1,V);
for iter = 1:iter_max
    tempC = zeros(n,n);
    for v = 1:V
        %update E
        tempE = Q1{v}/mu+X{v}-X{v}*Z{v};
        for i=1:n
            nw = norm(tempE(:,i));
            if nw>1/mu
                x = (nw-1/mu)*tempE(:,i)/nw;
            else
                x = zeros(length(tempE(:,i)),1);
            end
            tempE(:,i) = x;
        end
        E{v} = tempE;
        %update L、R
        tempL = (Q2{v}/mu+Z{v})*R{v}*C';
        [u,~,vv] = svd(tempL,'econ');
        L{v} = u*vv';
        tempR = (Q2{v}/mu+Z{v})'*L{v}*C;
        [u,~,vv] = svd(tempR,'econ');
        R{v} = u*vv';       
        % update Z
        tempZ = X{v}'*(Q1{v}/mu-E{v})+tempX{v}+L{v}*C*R{v}'-Q2{v}/mu+J{v}-Q3{v}/mu;
        Z{v} = (2*eye(n)+tempX{v})\tempZ;
        tempC = tempC + L{v}'*(Q2{v}/mu+Z{v})*R{v};
        % update J
        temp = Z{v} + Q3{v}/mu;
        J{v} = (temp + temp') / 2;
    end
    %update C
    [UUU,sigma,VVV] = svd(tempC/V,'econ');
    sigma = diag(sigma);
    xi = spw(sigma,lambda/(V*mu),p);
    C = UUU*diag(xi)*VVV';
    % update Q1,Q2,Q3 and mu
    for v = 1:V
        tempQ1 = X{v}-X{v}*Z{v}-E{v};Q1{v} = Q1{v}+mu*tempQ1;
        tempQ2 = Z{v}-L{v}*C*R{v}';Q2{v} = Q2{v}+mu*tempQ2;
        Q3{v} = Q3{v}+mu*(Z{v}-J{v});
    end
    mu = min(rho*mu,mu_max);
    % stop
    stop=0;
    for v=1:V
        tempstop=norm(Z{v}-C,inf);
        stop=stop+(tempstop-oldstop(1,v))/oldstop(1,v);
        oldstop(1,v)=tempstop;
    end
    disp(['iter' num2str(iter) 'stop' num2str(stop)])
    if abs(stop)<0.01||p==0.5
        break
    end
end
Zn = zeros(n, n);
for v = 1 : V
    Zn = Zn + Z{v};
end
end