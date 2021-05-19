function [] = SensingModelTest()
    x = -2:0.005:2;
    y = -2:0.005:2;
    [X,Y] = meshgrid(x,y);
    Z = sqrt(X.^2 + Y.^2)  + sqrt(X.^2 + (100 - Y).^2) - 100;
    %surfc(X,Y,Z,'EdgeColor', 'none');
    contour(X,Y,Z,100);
    xlabel('x');
    ylabel('y');
    zlabel('delta d');
end