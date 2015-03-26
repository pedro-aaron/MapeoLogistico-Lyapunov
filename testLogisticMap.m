%http://watermarkero.blogspot.mx/
%http://watermarkero.blogspot.com/2015/03/mapeo-logistico-y-calculo-del-exponente.html
%Mapeo Logístico y Cálculo del Exponente de Lyapunov
%Autor: Pedro Aaron Hernandez-Avalos
function testLogisticMap()
clc
close all

    %% graficar Mapeo logistico con diferentes valores de m
    figure;
    for m=0.4:0.4:4
        plotLogisticMap(m,'k',0);
        hold on;
    end
    title('Mapeo logistico, diferentes valores de \mu')
    ylabel('f(x)') % y-axis label
    xlabel('x') % x-axis label    
    hold off;

    
    %% graficar Mapeo Logistico con m fijo
    m = 4;
    figure;
    plotLogisticMap(m,'m',1);
    title('Diagrama de trayectoria')
    ylabel('f(x)') % y-axis label
    xlabel('x') % x-axis label    
    hold on;

    %% obtener trayectoria iniciando de punto inicial
    x0 = 0.1;
    nPoints = 200;
    trajectory = logisticMapTrajectory(m, x0, nPoints);
    
    %% graficar trayectoria
    plotTrajectory(trajectory)
    hold off;

    %% graficar la trayectoria como señal caótica
    x0 = 0.1;
    m = 4;
    nPoints = 1000;
    trajectory = logisticMapTrajectory(m, x0, nPoints);
    figure;
    plot(trajectory);

    %% Calcular el diagrama de bifurcación usando una trayectoria estable,
    % Se espera a que el sistema de estabilice, se eliminan los primeros 
    % k elementos    
    k = 300;
    step = 0.001;
    cont = 1;
    for m=0:step:4
        % el valor inicial es aleatorio
        x0 = rand(1);
        % se calcula la trayectoria
        trajectory = logisticMapTrajectory(m, x0, nPoints);
        %se eliminan los primeros k valores, para asegurar estabilización
        trayectoria_estable(cont,:) = trajectory((k+1):length(trajectory));
        mValues(cont,:) = m * ones(1,length(trajectory) - k);
        cont = cont + 1;
        disp(m);
    end

    %% graficar diagrama de bifurcación
    plotBifurcationDiagram(mValues,trayectoria_estable);
        
    %% calcular Lyapunov
    [n,~] = size(mValues);
    cont = 1;    
    for i=1:n
        if mValues(i,1) >3,
            lyap(cont) = lyapunov(trayectoria_estable(i,:),mValues(i,:));
            m(cont) = mValues(i,1);
            cont = cont + 1;   
        end
    end
    
    %% graficar Lyapunov
    figure;
    plot(m,lyap);        
    title('Exponente de Lyapunov')        
    ylabel('Exponente de Lyapunov') % y-axis label
    xlabel('\mu') % x-axis label    
end
%% Graficación del diagrama de bifurcación
function plotBifurcationDiagram(mValues,trayectoria)
    figure;
    title('Diagrama de bifurcación')    
    ylabel('f(x)') % y-axis label
    xlabel('\mu') % x-axis label    
    hold on;
    markerSize = 0.1;
    [n,~] = size(mValues);
    for i = 1:n
        scatter(mValues(i,:),trayectoria(i,:),markerSize,'Marker','.','MarkerEdgeColor','b','MarkerFaceColor','b')
    end
    hold off
end

%% cálculo del Exponente de Lyapunov
function lambda = lyapunov(trayectoria_estable,m)
    %calculo de la derivada absoluta de la trayectoria
    absDerivada = logisticMap_absDerivada(trayectoria_estable,m);
    % calculo del exponente
    lambda = sum(log(absDerivada))/length(trayectoria_estable);
    if lambda == -Inf, 
        lambda = -10; 
    end
end

%% Cálculo del absoluto de la derivada del Mapeo Logístico
function absDerivative = logisticMap_absDerivada(trayectoria_estable,m)
    % mapeo logistico: f(x) = mx(1 - x)
    % la derivada del mapeo logistico es: f'(x) = m - 2mx
    absDerivative = abs(m-2*m.*trayectoria_estable);
end

%% Cálculo de la función del Mapeo Logístico evaluada en los puntos x y con
% mu m
function fx = logisticMap(m, x)
    fx = m*x.*(1-x);
end

%% Graficación del Mapeo Logístico
function plotLogisticMap(m, color, showLeyend)
    X = 0:0.001:1;
    fx = logisticMap(m,X);        
    plot(X,fx,color);
    if showLeyend == 1,
        legend(['\mu = ' num2str(m)])
    end
end

%% Cálculo de las trayectorias del Mapeo Logístico, dado mu y x0
function trajectory = logisticMapTrajectory(m, Xo,nPoints)
    currentPoint = Xo;    
    for i = 1:nPoints
        trajectory(i) = currentPoint;        
        currentPoint = logisticMap(m,currentPoint);
    end
end

%% Graficación de los trayectorias
function plotTrajectory(trajectory)
    plot([0,1],[0,1],':r');
    axis([0,1,0,1]);
    x1 = trajectory(1);
    y1 = 0;
    x2 = trajectory(1);
    y2 = trajectory(2);
    plot([x1,x2],[y1,y2],'b','LineWidth',0.5);
    for i = 3:length(trajectory)
        x1 = trajectory(i-1);
        y1 = trajectory(i);
        x2 = trajectory(i-1);
        y2 = trajectory(i-1);
        x3 = trajectory(i-2);
        y3 = trajectory(i-1);
        
        plot([x1,x2,x3],[y1,y2,y3],'b','LineWidth',0.5);
    end
end