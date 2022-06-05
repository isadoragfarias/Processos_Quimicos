clear
clc
function [dy]=reacao_str1(t,y)
    CA = y(1); // s^-1
    dCA = -(k1 + k2) * CA;
    dCQ = k1 * CA;
    dCS = k2 * CA;
    dy=[dCA; dCQ; dCS];
endfunction

function [dy]=reacao_str2(t,y)
    CA = y(1);
    CQ = y(2);
    dCA = - k1 * CA;
    dCQ = k1 * CA - k2*CQ;
    dCS = k2 * CQ;
    dy=[dCA; dCQ; dCS];
endfunction

function [dy]=reacao_str3(t,y)
    CA=y(1);
    CQ=y(2);
    CS=y(3);
    CT=y(4);
 
    dCA = -k0*CA;
    dCQ = k0*CA - k1*CQ*CS^2 - k2*CQ;
    dCS = k1*CQ*CS^2 + k2*CQ - k3*CS;
    dCT = k3*CS;
    dy=[dCA; dCQ; dCS; dCT];

endfunction

// Processo 1
t0 = 0;  
h = 0.02;
tn = 30;
t = [t0:h:tn];
t2 = [t0:h:tn + 20];    
y0=[0.1; 0.0; 0.0];  
inter = [0.1 1 2 5 10]
n = size(inter,2)

for i = 1:n
    k1 = sample(1,inter)
    k2 = sample(1,inter)
    k11(i) = k1
    k22(i) = k2
    y1 = ode(y0,t0,t,reacao_str1);
    y2 = ode(y0,t0,t2,reacao_str2);
    val(:,:,i) = y1 // valores processo 1
    val2(:,:,i) = y2 // valores processo 2
end

for i = 1:n
    k(i) = 'k1 = '+string(k11(i)) + '  k2 = '+string(k22(i))
end

scf(0)
plot(t,val(1,:,1),'-',t,val(1,:,2),'--',t,val(1,:,3),':',t,val(1,:,4),'-.',t,val(1,:,5))
xtitle(["Concentração de A no processo 1"],boxed = %t)
ylabel('Concentração molar (mol/L)')
xlabel('Tempo (s)')
legend(k,1)
xs2eps(0, 'R11')



scf(1)
plot(t,val(2,:,1),'-',t,val(2,:,2),'--',t,val(2,:,3),':',t,val(2,:,4),'-.',t,val(2,:,5))
xtitle(["Concentração de Q no processo 1"],boxed = %t)
ylabel('Concentração molar (mol/L)')
xlabel('Tempo (s)')
legend(k,1)
xs2eps(1, 'R12')


scf(2)
plot(t,val(3,:,1),'-',t,val(3,:,2),'--',t,val(3,:,3),':',t,val(3,:,4),'-.',t,val(3,:,5))
xtitle(["Concentração de S no processo 1"],boxed = %t)
ylabel('Concentração molar (mol/L)')
xlabel('Tempo (s)')
legend(k,1)
xs2eps(2, 'R13')



// Processo 2

t = t2
scf(3)
plot(t,val2(1,:,1),'-',t,val2(1,:,2),'--',t,val2(1,:,3),':',t,val2(1,:,4),'-.',t,val2(1,:,5))
xtitle(["Concentração de A no processo 2"],boxed = %t)
ylabel('Concentração molar (mol/L)')
xlabel('Tempo (s)')
legend(k,1)
xs2eps(3, 'R21')


[a,b] = max(val2(2,:,1))
[c,d] = max(val2(2,:,2))
[e,f] = max(val2(2,:,3))
[g,h] = max(val2(2,:,4))
[i,j] = max(val2(2,:,5))
tempos = [t(b) t(d) t(f) t(h) t(j)]
concmax =[a c e g i]

scf(4)
plot(t,val2(2,:,1),'-',t,val2(2,:,2),'--',t,val2(2,:,3),':',t,val2(2,:,4),'-.',t,val2(2,:,5),t(b),a,'.r',t(d),c,'.r',t(f),e,'.r',t(h),g,'.r',t(j),i,'.r')
xtitle(["Concentração de Q no processo 2"],boxed = %t)
ylabel('Concentração molar (mol/L)')
xlabel('Tempo (s)')
legend(k,1)
xs2eps(4, 'R22')

mprintf('Interrupção da batelada no processo 2:\n')
for i = 1:size(k,1)
    disp(k(i))
    mprintf('Tempo(s) = %f',tempos(i))
    mprintf('\nConcentração máxima(mol/L) = %f\n',concmax(i))
end

scf(5)
plot(t,val2(3,:,1),'-',t,val2(3,:,2),'--',t,val2(3,:,3),':',t,val2(3,:,4),'-.',t,val2(3,:,5))
xtitle(["Concentração de S no processo 2"],boxed = %t)
ylabel('Concentração molar (mol/L)')
legend(k,1)
xlabel('Tempo (s)')
xs2eps(5, 'R23')


y0=[0.1; 0.0; 0.0; 0.0];
t = [t0:1:10000];
k0 = 1e-3; 
k1 = 2.5e9; 
k2 = 1e-2;
k3 = 1.0;
y3 = ode(y0,t0,t,reacao_str3);
scf(6)     
subplot(2,2,1)
plot(t,y3(1,:),'k:')
xtitle("CA em função do tempo")
subplot(2,2,2)
plot(t,y3(2,:),'k:')
xtitle("CQ em função do tempo")
subplot(2,2,3)
plot(t,y3(3,:),'k:')
xtitle("CS em função do tempo")
subplot(2,2,4)
plot(t,y3(4,:),'k:')
xtitle("CT em função do tempo")
xlabel('Tempo (s)')
xs2eps(6, 'R3')

J3 = [-k0 0 0; k0 -k1*y3(3,$)^2-k2 -2*k1*y3(3,$)*y3(2,$);0 k1*y3(3,$)^2+k2 2*k1*y3(2,$)*y3(3,$)-k3 ]
evals=spec(J3)
