function Q = Denoise3(P,Omega)
[M, N]= size(P);
A= Omega.*P;
maxiter=50;
q=min(30, N);
tol= 10^-5;
X= rand(M,q);
Y= rand(q, N);
Z=A;
U= zeros(M,q);
V= zeros(q, N);
GAM= zeros(M,q);
PI= zeros(q, N);
fac= 2.5*10^5/ norm(A,'fro');
A= A*fac;
alpha= 2.0*10^-4 * (norm(A,'fro'))* max(M,N)/q;
gamma= 1.618;
beta= N*alpha/M;
f0= norm(Omega.*(X*Y-A), 'fro')/norm(A, 'fro');

for i=1:maxiter
    [X,Y,Z,U,V,GAM,PI]=update(X,Y,Z,U,V,GAM, PI, alpha, beta, gamma, A, Omega);
    f1=norm(Omega.*(X*Y-A), 'fro')/norm(A, 'fro');
    if(abs(f0-f1)/max(1,f0) <=tol || f0<=tol)
        break;
    end
end
Q= X*Y/fac;
end

function [X1, Y1, Z1, U1,V1, GAM1, PI1]= update(X, Y, Z, U, V, GAM, PI, alpha, beta, gamma, A, Omega)
    Yup= Y*Y';
    
    X1=((Yup+ alpha*eye(size(Yup)))')\(Z*Y'+ alpha*U- GAM)';
    X1= X1';
    Xup= X1'*X1;
    Y1= (Xup+ beta*eye(size(Xup)))\(X1'*Z+ beta*V- PI);
    Z1= X1*Y1+ Omega.*(A- X1*Y1);
    U1= max(X1+ GAM/alpha, 0);
    V1= max(Y1+ PI/beta, 0);
    GAM1= GAM+ gamma*alpha*(X1- U1);
    PI1= PI+ gamma*beta*(Y1- V1);
    
end

