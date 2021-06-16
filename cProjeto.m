clear;
clc;

%An�lise do campo el�trico C)

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
            V_new(xIndex,yIndex)=(1-alphaOpt)*V_old(xIndex,yIndex)+alphaOpt*(V_new(xIndex,yIndex+1)+V_new(xIndex,yIndex-1)+V_new(xIndex+1,yIndex)+V_new(xIndex-1,yIndex))/4;
        end
    end
    
    if (sqrt(sum(sum((V_new-V_old).^2)))/sqrt(sum(sum(V_new.^2)))) < tol %Condi��o de toler�ncia
        [X,Y] = meshgrid(-L:h:L);
        mesh(X,Y,V_new);
        break
    end
    
    V_old=V_new;
end

[Ex,Ey]=gradient(V_new,h,h);
Ex=-Ex; Ey=-Ey; %Dado que os vectores t�m o sentido oposto

figure(2);
quiver(X,Y,Ex,Ey);
grid on;

E = norm([Ex(1,1),Ey(1,1)]); %O campo � uniforme podemos utilizar qualquer vetor;

fprintf('M�dulo do campo el�trico: %d.\n',E);
