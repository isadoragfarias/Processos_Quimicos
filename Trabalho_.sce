clear
clc
function [dy]=reacao_str(t,y)

CA=y(1);
CQ=y(2);
CS=y(3);
CT=y(4);

k0 = 2;
k1 = 5;

//k = [k0; k1]; // s^-1
//    dCA = -( k(1) + k(2) ) * CA;
//    dCQ = k(1) * CA;
//    dCS = k(2) * CA;
    
//k = [k0 ; k1]; // s^-1
//   dCA = - k(1) * CA;
//    dCQ = k(1) * CA - k(2)*CQ;
//    dCS = k(2) * CQ;

k0 = 1e-3; 
k1 = 2.5e-2; 
k2 = 1e-2;
k3 = 1.0; 
  dCA = -k0*CA;
  dCQ = k0*CA - k1*CQ*CS^2 - k2*CQ;
  dCS = k1*CQ*CS^2 + k2*CQ - k3*CS;
  dCT = k3*CS;


//J = [-(k(1)+k(2)) 0 0; k(1) 0 0; k(2) 0 0];
//disp(J)
//evals=spec(J);

//dy=[dCA; dCQ; dCS]
dy=[dCA; dCQ; dCS; dCT];

endfunction


t0=0;  
h=1;
tn=10000;
t=[t0:h:tn];     // Vetor da variavel independente

y0=[0.1; 0.0; 0.0; 0];   // Vetor de condicoes iniciais
//y0=[0.1; 0.0; 0.0]


// Chamada do pacote ode para integracao do sistema de EDOs
y=ode(y0,t0,t,reacao_str); 
disp("O tempo no qual a concentracao de //Q no processo 2 eh maxima: ")
tmax = t(y(2,:) == max(y(2,:)) )
disp(tmax)

// Apresentando a solucao na forma grafica
clf();
n = 4;
subplot(n,1,1)
plot(t,y(1,:),'k:')
xtitle("CA em funcao do tempo")
subplot(n,1,2)
plot(t,y(2,:),'k:')
xtitle("CQ em funcao do tempo")
subplot(n,1,3)
plot(t,y(3,:),'k:')
xtitle("CS em funcao do tempo")
subplot(n,1,4)
plot(t,y(4,:),'k:')
xtitle("CT em funcao do tempo")

//soma=(y(1,:)+y(2,:)+y(3,:)+y(4,:))'
//soma=(y(1,:)+y(2,:)+y(3,:))'
//disp(soma)

