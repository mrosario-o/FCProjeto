clear;
clc;

%Prisma condutor D)

%Properties:
L=1; %Width
h=0.1; %Step
tol=1e-7; %Tolerance
%....................


N=(2*L/h)+1; % Or N = length(-L:h:L);
alphaOpt=(2/(1+(pi/N)));
V_old=zeros(N,N);


%{
Mathematicaly:
Coords(index) = (2*L*(index-1))/(N-1)-L
IndexN(coordinate) = (coordinate+L)*(N-1)/(2*L)+1

Como N = (2*L/h)+1
Coords(index) = h*(index-1) - L
IndexN(coordinate) = (coordinate+L)/h + 1
%}


for xIndex=1:N
    for yIndex = 1:N
        if xIndex == N || xIndex == 1
            V_old(xIndex,yIndex) = (h*(yIndex-1)-L)/L; %Or V_old(xIndex,yIndex) = h*(yIndex-((N+1)/2))/L;
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
for x = -L:h:L
    for y = -L:h:L
        if x == L || x == -L
            V_old((x+L)/h + 1,(y+L)/h + 1) = y/L;
        end
        
        if y == L
            V_old((x+L)/h + 1,(y+L)/h + 1) = 1;
        end
        if y == -L
            V_old((x+L)/h + 1,(y+L)/h + 1) = -1;
        end
    end
end
%}

%Method

V_new=V_old;

nIte = 0;
while true
    nIte = nIte + 1; %Incremento de n�mero de itera��es
    
    for xIndex = 2:N-1 %Fronteiras n�o alteradas
        for yIndex = 2:N-1
            %Condi��o extra: A fronteira do prisma tem valor inalter�vel (0)
            %Como resultado todo o campo dentro dessas fronteiras � 0
            if ~(((xIndex >= (N-1)/4 && xIndex <= 3*(N-1)/4) && (yIndex == (N-1)/4 || yIndex == 3*(N-1)/4)) || ((yIndex >= (N-1)/4 && yIndex <= 3*(N-1)/4) && (xIndex == (N-1)/4 || xIndex == 3*(N-1)/4)))
            
            %Por outro lado podemos assumir por m�todos te�ricos que a zona
            %dentro das fronteiras � intoc�vel e portanto temos a seguinte
            %condi��o alternativa
            %if ~((xIndex >= (N-1)/4 && xIndex <= 3*(N-1)/4) && (yIndex >= (N-1)/4 && yIndex <= 3*(N-1)/4))
            
                V_new(xIndex,yIndex)=(1-alphaOpt)*V_old(xIndex,yIndex)+alphaOpt*(V_new(xIndex,yIndex+1)+V_new(xIndex,yIndex-1)+V_new(xIndex+1,yIndex)+V_new(xIndex-1,yIndex))/4;
            end
        end
    end
    
    if (sqrt(sum(sum((V_new-V_old).^2)))/sqrt(sum(sum(V_new.^2)))) < tol %Condi��o de toler�ncia
        [X,Y] = meshgrid(-L:h:L);
        mesh(X,Y,V_new);
        break
    end
    
    V_old=V_new;
end

fprintf('N�mero de itera��es: %d.\n',nIte);