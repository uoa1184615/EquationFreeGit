function f = gen_linear_system(DIM1,DIM2,fastband,slowband)

DIM=DIM1+DIM2; %total dimension
cD1 = floor(DIM1/2); %half the number of complex eigenvalues in first band
cD2 = floor(DIM2/2); %half the number of complex eigenvalues in second band

E = [];
V = [];

for n=1:cD1 %complex eigenvalues
    tmp=rand;
    e1 = [(1-tmp) tmp]*fastband+rand*1i;
    v1 = randn(DIM,1) + randn(DIM,1)*1i;
    E=diag([diag(E); e1; conj(e1)]);
    V = [V v1 conj(v1)];
end
for n = 2*cD1+1:DIM1 %real eigenvalues
    tmp=rand;
    e1 = [(1-tmp) tmp]*fastband;%+randn*1i;
    v1 = randn(DIM,1);% + rand(DIM,1)*1i;
    E=diag([diag(E); e1]);
    V = [V v1];
end


for n=1:cD2 %complex eigenvalues
    tmp=rand;
    e1 = [(1-tmp) tmp]*slowband+rand*1i;
    v1 = randn(DIM,1) + randn(DIM,1)*1i;
    E=diag([diag(E); e1; conj(e1)]);
    V = [V v1 conj(v1)];
end
for n = 2*cD2+1:DIM2 %real eigenvalues
    tmp=rand;
    e1 = [(1-tmp) tmp]*slowband;%+randn*1i;
    v1 = randn(DIM,1);% + rand(DIM,1)*1i;
    E=diag([diag(E); e1]);
    V = [V v1];
end

A=real(V*E/V);
B=randn(DIM); B=B/norm(B);
A=B*A/B;

b=randn(DIM,1);

f=@(t,x)A*x+b;