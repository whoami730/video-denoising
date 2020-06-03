function Q = Denoise(P, Omega)
N= size(P,1);
M= size(P,2);
Q0= zeros(N,M);
Pomega= P.*Omega;
SumVec= sum(Pomega, 2);
CountMat= Omega;
CountVec= sum(CountMat, 2);
AvVec= SumVec./CountVec;
Pomega= Pomega-AvVec;
Pomega= Pomega.*Omega;
Pomega= Pomega.*Pomega;
CountVec= CountVec-1;
SumVec= sum(Pomega, 2);
AvVec= SumVec./CountVec;
sigma= sum(AvVec, 'all');
sigma= sigma/N;
sigma= sqrt(sigma);
p= sum(Omega, 'all')/ (M*N);
u= (sqrt(M)+ sqrt(N))*sqrt(p)*sigma;
t=1;
it=1;
Q1= iterate(Q0, Omega, P, t, u);
ep= 10^(-5);
while it~=30 && sqrt(sum((Q1-Q0).*(Q1-Q0), 'all'))> ep
    Q0=Q1;
    Q1= iterate(Q1, Omega, P, t, u);
    it= it+1;
end

Q=Q1;

end


function D= softShrink(R, t)
[U,S,V]= svd(R);
S= S-t;
S= max(S, 0);
D= U*S* V';

end

function Q1= iterate(Q0, Omega, P, t, u)
R= Q0- t*(Q0-P).*Omega;
Q1= softShrink(R, t*u);
end


