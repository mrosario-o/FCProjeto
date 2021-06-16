clear;
clc;


% Método de Jacobi

L=1;

h = [0.01 0.02 0.05 0.1 0.2 0.5];
for hIndex = 1:(length(h))

    N=(2*L/h(hIndex))+1;
    M(hIndex) = N;

    V_old=zeros(N,N);

    for xIndex=1:N
        for yIndex = 1:N
            if xIndex == N || xIndex == 1
                V_old(xIndex,yIndex) = (yIndex-((N+1)/2))*h(hIndex);
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

    iteracao_max=1000000;
    tolerancia=1e-7;
    
    tic
    for ite=1:iteracao_max

        for xIndex=2:N-1 % Fronteira externa nao alterada
            for yIndex=2:N-1 % Fronteira externa nao alterada
                V_new(xIndex,yIndex)=(V_old(xIndex,yIndex+1)+V_old(xIndex,yIndex-1)+V_old(xIndex+1,yIndex)+V_old(xIndex-1,yIndex))/4;
            end
        end


        if (sqrt(sum(sum((V_new-V_old).^2)))/sqrt(sum(sum(V_new.^2)))) < tolerancia % Sum faz soma ao longo das colunas

            n_iter(hIndex)=ite;
            break

        end

        V_old=V_new;

    end
    time(hIndex) = toc;
end

plot(log(n_iter),log(M));
figure;
plot(log(time),log(M));