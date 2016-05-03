load out.txt
Nx=401;
Ny=326;
psi=reshape(out,[Ny,Nx]);
figure
mesh(psi)
%%
fileID=fopen('Self.txt','w');
mech_i=13;
for i=1:10000
    fprintf(fileID,'%f\t%f\t%f\n',Scatter_Table(1,i,mech_i),Scatter_Table(2,i,mech_i),Scatter_Table(3,i,mech_i));
end

%%
fileID=fopen('GammaM.txt','w');
fprintf(fileID,'%f\t%f\t%f\n',log(GammaM(1)),log(GammaM(2)),log(GammaM(3)));

%%
out=load('out.txt');
figure
quiver(out(:,1),out(:,2),out(:,3),out(:,4))

%%
figure
scatter(out(:,1),out(:,2),1)

%%
Vt=1.38e-23*300/(1.6e-19);
opsi=load('data2/opsi23000.txt');
Nx=401;
Ny=326;
opsi=reshape(opsi,[Ny,Nx]);
figure
mesh(opsi*Vt)
%%
Vt=1.38e-23*300/(1.6e-19);
opsi=load('opsi.txt');
Nx=401;
Ny=326;
opsi=reshape(opsi,[Ny,Nx]);
figure
mesh(opsi*Vt)
%%
omesh=load('omesh.txt');
Nx=401;
Ny=326;
omesh=reshape(omesh,[Ny,Nx]);
figure
mesh(omesh)

%%
Result=load('data2/Result.txt');

accu_s=Result(:,1);
accu_d=Result(:,2);

for i=2:length(accu_s)
    accu_s(i)=accu_s(i-1)+accu_s(i);
    accu_d(i)=accu_d(i-1)+accu_d(i);
end

figure
plot((1:length(accu_s))*0.5e-3,accu_s,'b')
hold on
plot((1:length(accu_s))*0.5e-3,-accu_d,'r')