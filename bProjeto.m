clear;
clc;

%Análise do método de Jacobi B1) e B2)

%Properties:
L=1; %Width
tol=1e-7; %Tolerance

h = [0.05 0.1 0.2 0.5]; %Diferent steps
%........................................

for hIndex = 1:length(h)

    N=(2*L/h(hIndex))+1; % Or N = length(-L:h(hIndex):L);
    M(hIndex) = N;
    V_old=zeros(N,N);


    %{
    Mathematicaly:
    Coords(index) = (2*L*(index-1))/(N-1)-L
    IndexN(coordinate) = (coordinate+L)*(N-1)/(2*L)+1

    Como N = (2*L/h(hIndex))+1
    Coords(index) = h(hIndex)*(index-1) - L
    IndexN(coordinate) = (coordinate+L)/h(hIndex) + 1
    %}


    for xIndex=1:N
        for yIndex = 1:N
            if xIndex == N || xIndex == 1
                V_old(xIndex,yIndex) = (h(hIndex)*(yIndex-1)-L)/L; %Or V_old(xIndex,yIndex) = h(hIndex)*(yIndex-((N+1)/2))/L;
            end

            if yIndex == N
                V_old(xIndex,yIndex) = 1;
            end
            if yIndex == 1
                V_old(xIndex,yIndex) = -1;
            end
        end
    end

    %{
    De um modo alternativo:
    (Tem a vantagem de usar coordenadas em vez de indices como o anterior)
    for x = -L:h(hIndex):L
        for y = -L:h(hIndex):L
            if x == L || x == -L
                V_old((x+L)/h(hIndex) + 1,(y+L)/h(hIndex) + 1) = y/L;
            end

            if y == L
                V_old((x+L)/h(hIndex) + 1,(y+L)/h(hIndex) + 1) = 1;
            end
            if y == -L
                V_old((x+L)/h(hIndex) + 1,(y+L)/h(hIndex) + 1) = -1;
            end
        end
    end
    %}

    V_new=V_old;

    nIte = 0;
    tic;
    while true
        nIte = nIte + 1; %Incremento de número de iterações

        for xIndex = 2:N-1 %Fronteiras não alteradas
            for yIndex = 2:N-1
                V_new(xIndex,yIndex)=(V_old(xIndex,yIndex+1)+V_old(xIndex,yIndex-1)+V_old(xIndex+1,yIndex)+V_old(xIndex-1,yIndex))/4;
            end
        end

        if (sqrt(sum(sum((V_new-V_old).^2)))/sqrt(sum(sum(V_new.^2)))) < tol %Condição de tolerância
            n_Iter(hIndex) = nIte;
            break
        end

        V_old=V_new;
    end
    
    time(hIndex) = toc;
end

plot(log(n_Iter),log(M));
figure;
plot(log(time),log(M));

%AJUSTE POLINOMIAL!!!
%PLOTS
%STEPS!!!