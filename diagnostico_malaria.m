clear;
clc;
clf;

%-----------------Longitudes y distancias------------------%
Lp=3.5;                      % Longitud de region de la carga positiva
Ln=2.5;                      % Longitud de region de la carga negativa
t=0.02;                      % t representa el grosor de tamaño que esta relacionado con la distancia de los electrodos
d=0.4;                       % distancia entre los electrodos
p=0.01;                      % ajusta posicion

% Definir las características de los electrodos.
ke=1/(4*pi*8.85*10^-12);     % K es la constante electrica de coulomb
Q=3e-5;                      % La carga de toda la barra
Nq=28;                       % Numero de cargas en la que se divide la carga en toda la barra

%----------Especificaciones de la gráfica------------%

xmin_original=-d/2-3*t;  xmax_original=-xmin_original;  % Puntos minimos y máximos en x originales
xmin=-d/2-3*t;  xmax=-xmin;                             % Vertices en x de la gráfica
ymin=2*(-Lp/2);   ymax=-ymin;                           % Vértices en y

%----Ajuste de la gráfica-------------%
if ymin <= -1  
    if xmin >= -0.5 && xmax <= 0.5
            xmin = -1.5;
            xmax = -xmin;
    end
end

%-----Crea vectores de la carga----------%
Ny=30;  Nx=Ny;
x=linspace(xmin, xmax, Nx); y=linspace(ymin, ymax, Ny);

%-----------Vertices de las placas-----------%

vertices2d=[[-d/2-t,Lp/2]    %1
    [-d/2,Lp/2]              %2
    [-d/2,-Lp/2]             %3  
    [-d/2-t,-Lp/2]           %4
    [d/2,Ln/2]               %5
    [d/2+t,Ln/2]             %6  
    [d/2+t,-Ln/2]            %7  
    [d/2,-Ln/2]];            %8  

%  Caras del rectangulo (placas)
facesP=[1 2 3 4 1];
facesN=[5 6 7 8 5];

% Color para las placas
colorP=[0.95,0,0];           
colorN=[0,0,0.7];            



%---------------------Comienza a posicionar las cargas------------------%
% Definir un diferencial de carga lineal
dq=Q/Nq;                     % Magnitud diferencial de carga

% Definir las posiciones de las cargas.
yp=linspace(-(1-p)*Lp/2,(1-p)*Lp/2,Nq); % Cargas positivas posiciones "Y"
xp(1:Nq)=-d/2-t/2;                      % Cargas positivas posiciones "X"
yn=linspace(-(1-p)*Ln/2,(1-p)*Ln/2,Nq);  % Cargas negativas posiciones "Y"
xn(1:Nq)=d/2+t/2;                       % Cargas negativas posiciones "X"


plot(xp, yp,'*') %Grafica para positivas
hold on
plot(xn,yn,'*') %Grafica cargas negativas.

%------Cálculo del campo eléctrico para cada punto XY (sin gradiente)------%
% % Inicializa los componentes del campo eléctrico.

 Ex = zeros(Nx, Ny);  % ¿Dónde deseas que se calcule este componente? Para cada punto de la carga
 Ey = zeros(Ny, Nx);
 
% % Calcula los componentes del campo eléctrico.
% Tres bucles for anidados comienzan aquí...

 for i= 1:Nx
     for j= 1:Ny
         for k= 1:Nq
             rp = sqrt((x(i)-xp(k)).^2 + (y(j)-yp(k)).^2);
             Ex(i,j) = Ex(i, j) + (ke .* dq ./rp.^2) .* ((x(i)-xp(k))./rp);
             Ey(i,j) = Ey(i, j) + (ke .* dq ./rp.^2) .* ((y(j)-yp(k))./rp);
             %negativo
             rn = sqrt((x(i)-xn(k)).^2 + (y(j)-yn(k)).^2);
             Ex(i,j) = Ex(i, j) - (ke .* dq ./rn.^2) .* ((x(i)-xn(k))./rn);
             Ey(i,j) = Ey(i, j) - (ke .* dq ./rn.^2) .* ((y(j)-yn(k))./rn);
         end
     end
 end 


% Datos del grafico de no gradiente
figure (1)
hold on
axis ([xmin xmax ymin ymax])
xlabel 'x position, mm'
ylabel 'y position, mm'
title 'Dielectrophoresis (No gradient)'
grid on

% Utiliza aquí el comando streamslice para trazar el campo E calculado...
% Líneas de campo eléctrico sin gradiente.
%Este es el streamslice del campo electrico sin gradiente
w = streamslice(x,y,Ex',Ey'); %Al usar la traspuesta (Ex' y Ey') convertimos los vectores de una dimensión en matrices de dos dimensiones
set (w, 'Color', 'b') %Le agrega el color a las flechas del campo
 

% Grafica las placas
patch('Faces',facesP,'Vertices',vertices2d,'FaceColor',colorP);
patch('Faces',facesN,'Vertices',vertices2d,'FaceColor',colorN);

%%%-- Grafica aquí ahora el campo E usando el gradiente. --%%%
% Usa líneas anteriores como plantilla, usa la misma estética. 

%Se repite el proceso de arriba pero ahora usando la formula E= -gradienteV

% Inicializa el potencial aquí
V(Nx,Ny)=0;

% Inicializa los componentes del campo eléctrico.
Ex = zeros(Nx, Ny);  
Ey = zeros(Nx, Ny);

% Calcula los componentes del campo eléctrico.

 for i= 1:Nx
     for j= 1:Ny
         for k= 1:Nq
             %positivo
             rp = sqrt((x(i)-xp(k)).^2 + (y(j)-yp(k)).^2);
             V(i,j)= V(i,j) +(ke.*dq./rp);
             %negativo
             rn = sqrt((x(i)-xn(k)).^2 + (y(j)-yn(k)).^2);
             V(i,j)= V(i,j) - (ke.*dq./rn);
         end
     end
 end 

% Calcula los componentes del campo eléctrico usando el gradiente de 
% potencial aquí.

[Ex, Ey] = gradient(-V'); % Se usa gradient para sacar el gradiente del potencial como lo explica la formula la formula E= -gradienteV y así poder graficar el campo

figure(2)
hold on
axis ([xmin xmax ymin ymax])
xlabel 'x position, mm'
ylabel 'y position, mm'
title 'Dielectrophoresis (gradient)'
grid on

% Establece la estética (para el gráfico de campo E) utilizando los 
% valores potenciales agregados recientemente (simplemente descomenta).
pcolor(x,y,V')                % Mapa de colores del voltaje.
colormap bone                 % Color

% Utiliza aquí el comando streamslice para trazar el campo E calculado...
streamslice(x,y,Ex,Ey)  

colorbar

% Dibuja las placas en el gráfico
patch('Faces',facesP,'Vertices',vertices2d,'FaceColor',colorP);
patch('Faces',facesN,'Vertices',vertices2d,'FaceColor',colorN);

contour(x,y,V',10,'-w','LineWidth',0.5)

%-----------------------Configuración de eritrocitos----------------------%
% Definir las características de los eritrocitos.

h=0.15;                                 % Parámetro ajustable (aquí entra biología) usando la letra "h". 
                                   % Considerar que h<0.03 se asociará con sangre sana y h>0.2 con sangre infectada 
                                   % (para un conjunto final de parámetros dados).

rad=(xmax_original-xmin_original)/40;  % Define el radio de los eritrocitos como una fracción de la escala del espacio x.

Ne=3;                       % Crea una variable para definir el No. de eritrocitos a modelar.
                            

dt=0.2;                             %   Variable de paso de tiempo "dt" en unidades arbitrarias (usa un valor de 0.2)

qe= 1e-6;                      % Crea una variable "qe" para definir la carga eléctrica de los eritrocitos. 
                               % Igualmente positivas y negativas - ¿por
                               % que?: No se modifica el signo de la carga ya que por su alto contenido en hierro tienden a moverse en dirección al campo electrico en general.
                               % utilizar un valor de 1 micro Coulombs.
                         
                            


%-----------------------Desplazamiento de eritrocitos---------------------%
% Inicializar el recuento de eritrocitos sanos e infectados.
healthy=0; infected=0;


% (3) Abre un bucle for que recorrerá el número total de eritrocitos (Ne) 
% definido por el usuario, y dentro de él inicializará la posición de los 
% eritrocitos (xe, ye) y las velocidades (Vx, Vy)

for ery= 1:Ne
    path=animatedline('LineWidth',2,'LineStyle',':','Color','#7E2F8E');
    dx=h*rand()*rad;       % Creemos que es la distancia que se ladea en x dependiendo de si es sano o no sano
    xe= 0;                 %xe es cero ya que los eritrocitos deben empezar en medio del plano para que pase por las placas
    ye= ymax;
    Vx= 0; 
    Vy= -1;

% (4) A continuación, dentro de este bucle y después de la inicialización 
% mencionada anteriormente, usa (define) un while que se ejecutará desde la 
% coordenada Y más alta del eritrocito, digamos "ye", hasta la coordenada Y
% más baja (más abajo) de el eritrocito (digamos "ye>ymin") 
% (5) Finalmente, dentro de este while abre otro bucle for que irá desde la
% carga 1 hasta la última carga en las placas (electrodos)..

    while ye>ymin
    % ¿Para qué es esto?
         addpoints(path,xe,ye); 
         head=scatter(xe,ye, 100,'filled','o','red');   %Scatter se utiliza para que pequeños circulos se grafiquen en la trayectoria del eritrocito
         drawnow  
          %--------------Cálculo de fuerza-------------------%
         % Inicializar fuerza eléctrica
         Fx=0;

        for i= 1:Nq
        
        %Basandonos en el plano cartersiano, las variables se le denominan
        %como p y n dependiendo de si la x y la y se encuentran en el lado
        %positivo o negativo del plano
             rnp=sqrt((xe-dx-xp(i))^2+(ye-yp(i))^2); %Placa roja x negativa y positiva
             rnn=sqrt((xe-dx-xn(i))^2+(ye-yn(i))^2); %Placa roja x negativa y ya negativa
             rpp=sqrt((xe+dx-xp(i))^2+(ye-yp(i))^2); %Placa azul, x positiva y positiva
             rpn=sqrt((xe+dx-xn(i))^2+(ye-yn(i))^2); %Placa azul, x positiva y negativa

             Fnn=ke*dq*qe/rnn^2;%Placa roja con x negativa y y negativas           
             Fnp=ke*dq*qe/rnp^2; %Placa roja con x negativas y y positiva
             Fpp=ke*dq*qe/rpp^2; %Placa azul con x podditiva y y positiva
             Fpn=ke*dq*qe/rpn^2; %Placa azul con x positiva y y negativa

             %Todo esto se hace en el scatter final, los valores de xe y de
             %ye se van actualizando con el for superior que se usa con los
             %valores de pp pn np nn, pero utilizando pitagoras sacando el
             %cateto
             Fx=Fx-Fpn*(xe+dx-xn(i))/rpn-Fpp*(xe+dx+xn(i))/rpp+Fnn*(xe-dx-xn(i))/rnn+Fnp*(xe-dx-xp(i))/rnp; 
    
        end
    
    m = 1;
    a = Fx/m;
    Vx =Vx + a*dt;
    xe = xe + Vx*dt + 0.5*a*dt.^2;
    ye = ye + dt*Vy;

    % ¿Para qué es esto?
        if ye>ymin + dx
        delete(head);  %Esto es para que vaya eliminando los puntos pasados y se observe el movimiento en tiempo real
        end
    end
    if xe<0.07
         healthy=healthy+1;
     else
         infected=infected+1;
     end
end

% Calcula la probabilidad de infección
infected_probability = 100*infected/(healthy+infected);

str = ['Infected probability of  ', num2str(infected_probability), ' %'];
annotation('textbox', [0.15, 0.8, 0.1, 0.1], 'String', str, 'FitBoxToText', 'on');

fprintf('Erythrocyte: %d has an associated position of : %f\n', ery, dx);

