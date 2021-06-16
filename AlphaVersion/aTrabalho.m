clear;
clc;


% Método de Jacobi

L=1;
h=0.05; %Passo

N=(2*L/h)+1;

V_old=zeros(N,N);

for xIndex=1:N
    for yIndex = 1:N
        if xIndex == N || xIndex == 1
            V_old(xIndex,yIndex) = ((yIndex-((N+1)/2))*h)/L;
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

for ite=1:iteracao_max
    
    for xIndex=2:N-1 % Fronteira externa nao alterada
        for yIndex=2:N-1 % Fronteira externa nao alterada
            V_new(xIndex,yIndex)=(V_old(xIndex,yIndex+1)+V_old(xIndex,yIndex-1)+V_old(xIndex+1,yIndex)+V_old(xIndex-1,yIndex))/4;
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

fprintf('Número de iterações realizadas foi %d.\n',n_iter)
