clear;
clc;

%An?lise do campo el?trico C)

%Properties:
L = 1; %Width
h = 0.05; %Step
tol = 1e-7; %Tolerance

R = 1/2; %Cylinder's radius
P = [1/4 1/8]; %Cylinder's center [x y]
%........................................


N = (2*L/h)+1; % Or N = length(-L:h:L);
alphaOpt = (2/(1+(pi/N)));
V_old = zeros(N,N);


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

V_new = V_old;

nIte = 0;
while true
    nIte = nIte + 1; %Incremento de n?mero de itera??es
    
    for xIndex = 2:N-1 %Fronteiras n?o alteradas
        for yIndex = 2:N-1
            if (h*(xIndex-1) - L - P(1))^2+(h*(yIndex-1) - L - P(2))^2 <= R^2
                %Condi??o do cil?ndro (x^2+y^2<=R^2)
                V_new(xIndex,yIndex)=(1-alphaOpt)*V_old(xIndex,yIndex)+alphaOpt*(V_new(xIndex,yIndex+1)+V_new(xIndex,yIndex-1)+V_new(xIndex+1,yIndex)+V_new(xIndex-1,yIndex)+100*h^2)/4;
            else
                V_new(xIndex,yIndex)=(1-alphaOpt)*V_old(xIndex,yIndex)+alphaOpt*(V_new(xIndex,yIndex+1)+V_new(xIndex,yIndex-1)+V_new(xIndex+1,yIndex)+V_new(xIndex-1,yIndex))/4;
            end
        end
    end
    
    if (sqrt(sum(sum((V_new-V_old).^2)))/sqrt(sum(sum(V_new.^2)))) < tol %Condi??o de toler?ncia
        [X,Y] = meshgrid(-L:h:L);
        
        figure;
        mesh(X,Y,V_new);
        
        title('Potencial na superf?cie');
        xlabel('Eixo y');
        ylabel('Eixo x');
        zlabel('Diferen?a de potencial');
        break
    end
    
    V_old = V_new;
end

fprintf('N?mero de itera??es: %d.\n',nIte);

[Ex,Ey] = gradient(V_new,h,h);
Ex = -Ex; Ey = -Ey; %Dado que os vectores t?m o sentido oposto

figure;
quiver(X,Y,Ex,Ey,'Color','r');
grid on;
axis equal;

title('Campo el?trico');
xlabel('Eixo y');
ylabel('Eixo x');
