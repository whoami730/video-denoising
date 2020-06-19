function X = Denoise2(P,Omega)

[N,M] = size(P);

Pomega= P.*Omega;
SumVec= sum(Pomega, 2);
CountMat= Omega;
CountVec= sum(CountMat, 2);
CountVec= max(CountVec,1);
AvVec= SumVec./CountVec;
Pomega= Pomega-(AvVec.*Omega);
Pomega= Pomega.^2;
SumVec= sum(Pomega, 2);
AvVec= SumVec./CountVec;
sigma= mean(AvVec, 'all','omitnan');
rgl = 0.05;
if isnan(sigma)
    warning("sigma is NaN, not a good thing");
    sigma = rgl;
else
    sigma = sigma+rgl;
end
sigma = sqrt(sigma);
p= sum(Omega, 'all')/ (M*N);
t= (sqrt(M)+ sqrt(N))*sqrt(p)*sigma;


Y0= zeros(N,M);
del= 1.2*N*M/sum(Omega,'all')/2;
ep = 10^-5;
it=1;

X= softShrink(Y0, t);
while it~=100 && norm(P-X, 'fro')> ep
    Y1= Y0 + del*Omega.*(P-X);
    %Y1= iterate(Y1, Omega, P, t, del);
    X= softShrink(Y1, t);
    it= it+1;
end
%X= softShrink(Y1,t);    
end

function D= softShrink(R, t)
[U,S,V]= svd(R);
S= S-t;
S= max(S, 0);
D= U*S* V';

end

function Y1= iterate(Y0, Omega, P, t, del)
X = softShrink(Y0, t);
Y1= Y0 + del*Omega.*(P - X);
end