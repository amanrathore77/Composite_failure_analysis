clear all;
clc;

t=input('Enter thickness of lamina (m)\n');%thickness
num=input('Enter number of even lamina (m)\n');%number of plies
arry = input('Enter the stacking sequence of the symmetric laminate (for eg.[0,45,-45,90])\n');
E1=input('Enter E1 value (N/m^2)\n');
E2=input('Enter E2 value (N/m^2)\n');
G12=input('Enter G12 value (N/m^2)\n');
V12=input('Enter V12 value\n');
V21=V12*E2/E1;
delt=input('Enter temperature change value \n');

sig1=input('Enter longitudinal tensile and compressive(negative) strength in N/m^2(for eg. [1062,-610]*10^6\n');
sig2=input('Enter transverse tensile and compressive(negative) strength in N/m^2(for eg. [1062,-610]*10^6\n');
tau=input('Enter shear strength in N/m^2(for eg. 1062*10^6\n');
a1 = input('Enter coefficient of thermal expansion value in longitudinal direction m per degree\n');
a2 = input('Enter coefficient of thermal expansion value in transverse direction m per degree\n');

a=[a1;a2;0];
NX=input('Enter Nx value\n');
NY=input('Enter Ny value\n');
NZ=input('Enter Nxy value\n');
MX=input('Enter Mx value\n');
MY=input('Enter My value\n');
MZ=input('Enter Mxy value\n');

L_mat=[NX NY NZ MX MY MZ]';



z=zeros(num+1,1);

x=(pi/180)*arry;

%z_k values for each lamina
%----------------------------------------------------------------
in=(num/2)+1;%index of neutral for even number of lamina 
z(in,1)=0;
for i=1:num/2
    z(i,1)=(in-i)*t;
    z(num+2-i,1)=-z(i,1);
end
%-------------------------------------------------------------------

Q11=E1/(1-V12*V21);
Q12=V12*E2/(1-V21*V12);
Q22=E2/(1-V12*V21);
Q66=G12;
Q=[Q11 Q12 0;Q12 Q22 0;0 0 Q66];
A=zeros(3,3);
B=zeros(3,3);
D=zeros(3,3);
global_ax_str=zeros(3,num); %global axis strain for each lamina
mat_ax_strain=zeros(3,num); %material axis strain for each lamina
mat_ax_stress=zeros(3,num); %material axis stress for each lamina
str=zeros(3,num);
Nt=zeros(3,1);
Mt=zeros(3,1);
res_stress=zeros(3,num);
%----------------------------------------------------------------
for i=1:num/2   %4 as number of laminas are 8 but it is symmetric
C=cos(x(1,i)); S=sin(x(1,i));

q11=Q11*C^4+Q22*S^4+2*(Q12+2*Q66)*S^2*C^2;
q12=(Q11+Q22-4*Q66)*C^2*S^2+Q12*(C^4+S^2);
q22=Q11*S^4+Q22*C^4+2*(Q12+2*Q66)*S^2*C^2;
q16=(Q11-Q12-2*Q66)*C^3*S-(Q22-Q12-2*Q66)*S^3*C;
q26=(Q11-Q12-2*Q66)*S^3*C-(Q22-Q12-2*Q66)*C^3*S;
q66=(Q11+Q22-2*Q12-2*Q66)*S^2*C^2+Q66*(S^4+C^4);

%reduced transformation stiffness matrix

q(:,:,i)=[q11 q12 q16;q12 q22 q26;q16 q26 q66];
A=A+q(:,:,i)*((z(i,1)-z(i+1,1))+(z(num+1-i,1)-z(num+2-i,1)));
B=B+0.5*((z(i,1)^2-z(i+1,1)^2)+(z(num+1-i,1)^2-z(num+2-i,1)^2))*q(:,:,i);
D=D+(1/3)*((z(i,1)^3-z(i+1,1)^3)+(z(num+1-i,1)^3-z(num+2-i,1)^3))*q(:,:,i);

end
%-------------------------------------------------------------
inv_A=inv(A);
E_x=1/(inv_A(1,1)*t*num); %effective modulus of laminate

ABBD=[A,B;B,D];
mid_plane_strain=inv(ABBD)*L_mat; %mid plane strain in global axis of laminate(epsilon) and curvatures(K) 

% now calculating global axis strains in each lamina
for i=1:num
    global_ax_str(:,i)=mid_plane_strain(1:3,1)+(z(i,1)+z(i+1,1))*mid_plane_strain(4:6,1)/2;
end
global_ax_str_reduced=global_ax_str;
global_ax_str_reduced(3,num)=global_ax_str(3,num)/2;

%material axis strain in each lamina

  for i=1:num/2
      mat_ax_strain(:,i)=transf(x(1,i))*global_ax_str_reduced(:,i);
      mat_ax_strain(:,num+1-i)= mat_ax_strain(:,i); %remember gamma is 1/2 should multiply  by 2
  end
   mat_ax_strain(3,:)=mat_ax_strain(3,:)*2;
   
 %material axis stress in each lamina
 
 for i=1:num/2
     mat_ax_stress(:,i)=Q*mat_ax_strain(:,i);
     mat_ax_stress(:,num+1-i)=mat_ax_stress(:,i);
 end
%strength ratio 
    for i=1:num
        
     if mat_ax_stress(1,i)<0
         str(1,i)= mat_ax_stress(1,i)/sig1(1,2);
     else 
         str(1,i)=mat_ax_stress(1,i)/sig1(1,1);
     end
     
     if mat_ax_stress(2,i)<0
         str(2,i)= mat_ax_stress(2,i)/sig2(1,2);
     else 
         str(2,i)=mat_ax_stress(2,i)/sig2(1,1);
     end
     str(3,i)=mat_ax_stress(3,i)/tau;
    end
    
 %first ply failure load without hygrothermal effect
    fpf=NX/max(str,[],'all');
    
    for i=1:num/2
alpg(:,i)=inv(transf(x(1,i)))*a;
alpg(:,num+1-i)=alpg(:,i);
end
alpg(3,:)=2*alpg(3,:);
for i=1:num/2
    Nt=Nt+ delt*q(:,:,i)*alpg(:,i)*((z(i,1)-z(i+1,1))+(z(num+1-i,1)-z(num+2-i,1)));
    Mt=Mt+delt*q(:,:,i)*alpg(:,i)*((z(i,1)^2-z(i+1,1)^2)+(z(num+1-i,1)^2-z(num+2-i,1)^2))/2;
end
nt=[Nt;Mt];
%mid surface strain due to hygrothermal change
mid_str_ht=inv(ABBD)*nt;
for i=1:num
    gb_str(:,i)=mid_str_ht(1:3,1)+(z(i,1)+z(i+1,1))*mid_str_ht(4:6,1)/2;
end
%free thermal strain in each ply
    ft_str=delt*alpg;
 %residual strain
 res_str=gb_str-ft_str;
 %residual stress in plies in global axis
 for i=1:num/2
     res_stress(:,i)=q(:,:,i)*res_str(:,i);
     res_stress(:,num+1-i)=res_stress(:,i);
 end
 %residual stress in each plies in material axis
 for i=1:num/2
     res_stress_ma(:,i)=transf(x(1,i))*res_stress(:,i);
     res_stress_ma(:,num+1-i)=res_stress_ma(:,i);
 end
    U=max(str);
    V=max(str');
    [k u]=max(U); %coloumn
    [K v]=max(V); %row
   % total_str=mat_ax_stress(v,u)+res_stress_ma(v,u);
   if v==1
       sig=sig1;
   else 
       sig=sig2;
   end
   if mat_ax_stress(v,u)>0
       s=sig(1,1);
   else 
       s=sig(1,2);
   end
    str1=mat_ax_stress(v,u)/(s-res_stress_ma(v,u));
    FPF=NX/str1;

 