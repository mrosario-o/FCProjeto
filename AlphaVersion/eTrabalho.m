clear;
clc;


% Campo el�trico

L=1;
h=0.05; %Passo

N=(2*L/h)+1;

V_old=zeros(N,N);

for xIndex=1:N
    for yIndex = 1:N
        if xIndex == N || xIndex == 1
            V_old(xIndex,yIndex) = (yIndex-((N+1)/2))*h;
        end
        
        if yIndex == N
            V_old(xIndex,yIndex) = 1;
        end
        if yIndex == 1
            V_old(xIndex,yIndex) = -1;
        end
    end
end

V_new=V_old;

iteracao_max=10000;
tolerancia=1e-7;
alpha=(2/(1+(pi/N)));

for ite=1:iteracao_max
    
    for xIndex=2:N-1 % Fronteira externa nao alterada
        for yIndex=2:N-1 % Fronteira externa nao alterada
            V_new(xIndex,yIndex)=(1-alpha)*V_old(xIndex,yIndex)+alpha*(V_new(xIndex,yIndex+1)+V_new(xIndex,yIndex-1)+V_new(xIndex+1,yIndex)+V_new(xIndex-1,yIndex))/4;
        end
    end
    
    
    if (sqrt(sum(sum((V_new-V_old).^2)))/sqrt(sum(sum(V_new.^2)))) < tolerancia % Sum faz soma ao longo das colunas
        
        n_iter=ite;
        xy= -L:h:L;
        mesh(xy,xy,V_new)
        break
        
    end
    
    V_old=V_new;
    
end

[E_x,E_y]=gradient(V_new,h,h); %E_x e E_y sao os parametros do campo eletrico

%   Fazemos isto por causa do menos da expressão, nao podemos meter na função
% gradient
E_x=-E_x;
E_y=-E_y;

figure(2)
quiver(xy,xy,E_x,E_y)
grid on

E=sqrt(E_x(1,1).^2+E_y(1,1).^2);
