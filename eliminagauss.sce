// Script que implementa o algoritmo de substituição de Gauss

// Matriz A e vetor b  (exemplo)
A = [1 -1 -2 -1; 2 -2 3 -3; 1 1 1 0; 1 -1 4 3]
b = [-8; -20; 4; 4]

// Verificando se a matriz A eh quadrada (nl=nc)
[nl,nc]=size(A);
if (nl <> nc) then
    mprintf('Matriz A nao eh quadrada');
    abort
end
n=nl;

// Montando a matriz aumentada A=[A b]
A=[A b];

// Passo 1
for i=1:(n-1),

// Passo 2
    k=i;
    while (A(k,i) == 0),
      k=k+1;
      if (k > n) then
        mprintf ('matriz eh singular');
        abort
      end,
    end,

// Passo 3
    if (k <> i) then
      for j=1:(n+1),
        aprov=A(i,j);
        A(i,j)=A(k,j);
        A(k,j)=aprov;
      end,
    end,

// Passo 4

    for j=i+1:n,
      if (A(j,i) <> 0) then

// Passo 5
        coef=A(j,i)/A(i,i);

// Passo 6
         for k=i+1:(n+1),
           A(j,k)=A(j,k) - coef*A(i,k);
         end,
      end,
    end,
end,

// Passo 7
  if (A(n,n) == 0) then
    mprintf ('Matriz A singular');
    abort
  end,

// Passo 8
x(n)=A(n,n+1)/A(n,n);

// Passo 9
for i=(n-1):-1:1,
  soma=0.0;
    for j=i+1:n,
      soma=soma+A(i,j)*x(j);
    end,
  x(i)=(A(i,n+1)-soma)/A(i,i);
end

mprintf('Valor de x:\n x(1)=%f \n x(2)=%f \n x(3)=%f',x(1),x(2),x(3))

