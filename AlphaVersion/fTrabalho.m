clear;
clc;


% Campo el�trico

L=1;
h=0.1; %Passo

N=(2*L/h)+1;

V_old=zeros(N,N);

for xIndex=1:N
    for yIndex = 1:N
        if xIndex == N || xIndex == 1
            V_old(xIndex,yIndex) = (yIndex-((N+1)/2))*h;
            V_old(xIndex,yIndex) = (h*(yIndex-1)/L-L)/L;
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
        if ~(((xIndex >= (N-1)/4 && xIndex <= 3*(N-1)/4) && (yIndex == (N-1)/4 || yIndex == 3*(N-1)/4)) || ((yIndex >= (N-1)/4 && yIndex <= 3*(N-1)/4) && (xIndex == (N-1)/4 || xIndex == 3*(N-1)/4))) 
            V_new(xIndex,yIndex)=(1-alpha)*V_old(xIndex,yIndex)+alpha*(V_new(xIndex,yIndex+1)+V_new(xIndex,yIndex-1)+V_new(xIndex+1,yIndex)+V_new(xIndex-1,yIndex))/4;
        end
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
