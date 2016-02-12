x=[ 
     NaN 0 1 NaN % NaN= not a number 
       0 0 1 1 
       0 0 NaN 1 
     NaN 0 1 NaN 
     NaN 0 1 NaN 
     NaN NaN NaN NaN 
     ]; 
y=[ 
     NaN 0 0 NaN 
       0 0 0 0 
       1 1 1 1 
     NaN 1 1 NaN 
     NaN 0 0 NaN 
     NaN NaN NaN NaN 
     ]; 
z=[ 
     NaN 0 0 NaN 
       0 1 1 0 
       0 1 1 0 
     NaN 0 0 NaN 
     NaN 0 0 NaN 
     NaN NaN NaN NaN 
     ];

     cc=zeros(8,3); % zeros crea un arreglo de 8(Colores) por 3(x,y,z) 
     cc(1,:)=[0 0 0]; % black 
     cc(2,:)=[1 0 0]; % red 
     cc(3,:)=[0 1 0]; % green 
     cc(4,:)=[0 0 1]; % blue 
     cc(5,:)=[1 0 1]; % magenta 
     cc(6,:)=[0 1 1]; % cyan 
     cc(7,:)=[1 1 0]; % yellow 
     cc(8,:)=[1 1 1]; % white 
     cs=size(x); 
     c=repmat(zeros(cs),[1 1 3]);%repmat crea una nueva matriz donde cada reglon contiene el arreglo zeros 
for i=1:size(cc,1) 
     ix=find(x==cc(i,1) &... 
             y==cc(i,2) &... 
             z==cc(i,3)); 
     [ir,ic]=ind2sub(cs,ix); 
for k=1:3 
for m=1:length(ir) 
     c(ir(m),ic(m),k)=cc(i,k); 
end 
end 
end 
     s=surf(x,y,z,c); 
     shading interp; % Para que difumine el color 
     xlabel('x'); 
     ylabel('y'); 
     zlabel('z'); 
     title('CuboRGB'); 
     axis equal; 
     axis on; 
rotate3d on