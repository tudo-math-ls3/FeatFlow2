function plotelem(celement,varargin)
%PLOTELEM finite element visualization.
% PLOT(celement) plots all basis functions on the reference element.
%
% PLOT(celement,n) plots all basis functions on the reference element and
% uses resolution n for visualization.
%
% PLOT(celement,n,u) plots the finite element solution u on the reference
% element and uses resolution n for visualization.

switch(nargin)
    case{2}
        n = varargin{1};
    case{3}
        n = varargin{1};
        u = varargin{2};
    otherwise
        n = 10;
end

switch(celement)
    case{'P1','E001'}
        % Linear finite element
        x=linspace(0,1,n);
        y=linspace(0,1,n);
        [X,Y]=meshgrid(x,y);
        tri=delaunay(X,Y);
        idx=find(X+Y>1);
        
        p1 = 1-X-Y;
        p2 = X;
        p3 = Y;      
                
        clf
        if exist('u'),
            p = u(1).*p1 + u(2).*p2 + u(3).*p3;
            p(idx)=nan;
            
            trisurf(tri,X,Y,p);
            title('(P1): FE-solution')
            xaxis('x-axis')
            yaxis('y-axis')
        else
            p1(idx)=nan;
            p2(idx)=nan;
            p3(idx)=nan;
            
            subplot(1,3,1)
            trisurf(tri,X,Y,p1);
            title('(P1): p_1(x,y)=1-x-y')
            xaxis('x-axis')
            yaxis('y-axis')
            
            subplot(1,3,2)
            trisurf(tri,X,Y,p2)
            title('(P1): p_2(x,y)=x')
            xaxis('x-axis')
            yaxis('y-axis')
            
            subplot(1,3,3)
            trisurf(tri,X,Y,p3)
            title('(P1): p_3(x,y)=y')
            xaxis('x-axis')
            yaxis('y-axis')
        end
        
    case{'P1T'}
        % Rotated linear finite element
        x=linspace(0,1,n);
        y=linspace(0,1,n);
        [X,Y]=meshgrid(x,y);
        tri=delaunay(X,Y);
        idx=find(X+Y>1);
        
        p1 =  1-2*Y;
        p2 = -1+2*X+2*Y;
        p3 =  1-2*X;      
                
        clf
        if exist('u'),
            p = u(1).*p1 + u(2).*p2 + u(3).*p3;
            p(idx)=nan;
            
            trisurf(tri,X,Y,p);
            title('(P1T): FE-solution')
            xaxis('x-axis')
            yaxis('y-axis')
        else
            p1(idx)=nan;
            p2(idx)=nan;
            p3(idx)=nan;
            
            subplot(1,3,1)
            trisurf(tri,X,Y,p1);
            title('(P1T): p_1(x,y)=1-2y')
            xaxis('x-axis')
            yaxis('y-axis')
            
            subplot(1,3,2)
            trisurf(tri,X,Y,p2)
            title('(P1T): p_2(x,y)=-1+2x+2y')
            xaxis('x-axis')
            yaxis('y-axis')
            
            subplot(1,3,3)
            trisurf(tri,X,Y,p3)
            title('(P1T): p_3(x,y)=1-2x')
            xaxis('x-axis')
            yaxis('y-axis')
        end
        
    case{'Q1','E011'}
        % Bilinear finite elemenet
        x=linspace(-1,1,n);
        y=linspace(-1,1,n);
        [X,Y]=meshgrid(x,y);
        
        p1 = 1/4*(1-X).*(1-Y);
        p2 = 1/4*(1+X).*(1-Y);
        p3 = 1/4*(1+X).*(1+Y);
        p4 = 1/4*(1-X).*(1+Y);
        
        clf
        if exist('u'),
            p = u(1).*p1 + u(2).*p2 +...
                u(3).*p3 + u(4).*p4;
            surf(X,Y,p);
            title('(Q1): FE-solution')
            xaxis('x-axis')
            yaxis('y-axis')
        else
            subplot(2,2,1)
            surf(X,Y,p1)
            title('(Q1): p_1(x,y)=1/4(1-x)(1-y)')
            xaxis('x-axis')
            yaxis('y-axis')
            
            subplot(2,2,2)
            surf(X,Y,p2)
            title('(Q1): p_2(x,y)=1/4(1+x)(1-y)')
            xaxis('x-axis')
            yaxis('y-axis')
            
            subplot(2,2,3)
            surf(X,Y,p3)
            title('(Q1): p_3(x,y)=1/4(1+x)(1+y)')
            xaxis('x-axis')
            yaxis('y-axis')
            
            subplot(2,2,4)
            surf(X,Y,p4)
            title('(Q1): p_4(x,y)=1/4(1-x)(1+y)')
            xaxis('x-axis')
            yaxis('y-axis')
        end
        
    case{'E030'}
        % Rotated bilinear finite elemenet (meanvalue-based variant)
        x=linspace(-1,1,n);
        y=linspace(-1,1,n);
        [X,Y]=meshgrid(x,y);
        
        p1 = -3/8*(X.^2-Y.^2) - Y/2 + 1/4;
        p2 =  3/8*(X.^2-Y.^2) + X/2 + 1/4;
        p3 = -3/8*(X.^2-Y.^2) + Y/2 + 1/4;
        p4 =  3/8*(X.^2-Y.^2) - X/2 + 1/4;
        
        clf
        if exist('u'),
            p = u(1).*p1 + u(2).*p2 +...
                u(3).*p3 + u(4).*p4;
            surf(X,Y,p);
            title('(E030): FE-solution')
            xaxis('x-axis')
            yaxis('y-axis')
        else
            subplot(2,2,1)
            surf(X,Y,p1)
            title('(E030): p_1(x,y)=-3/8(x^2-y^2)-1/2y+1/4')
            xaxis('x-axis')
            yaxis('y-axis')
            
            subplot(2,2,2)
            surf(X,Y,p2)
            title('(E030): p_2(x,y)=3/8(x^2-y^2)+1/2x+1/4')
            xaxis('x-axis')
            yaxis('y-axis')
            
            subplot(2,2,3)
            surf(X,Y,p3)
            title('(E030): p_3(x,y)=-3/8(x^2-y^2)+1/2y+1/4')
            xaxis('x-axis')
            yaxis('y-axis')
            
            subplot(2,2,4)
            surf(X,Y,p4)
            title('(E030): p_4(x,y)=3/8(x^2-y^2)-1/2x+1/4')
            xaxis('x-axis')
            yaxis('y-axis')
        end
        
    case{'E031'}
        % Rotated bilinear finite elemenet (midpoint-based variant)
        x=linspace(-1,1,n);
        y=linspace(-1,1,n);
        [X,Y]=meshgrid(x,y);
        
        p1 = -1/4*(X.^2-Y.^2) - Y/2 + 1/4;
        p2 =  1/4*(X.^2-Y.^2) + X/2 + 1/4;
        p3 = -1/4*(X.^2-Y.^2) + Y/2 + 1/4;
        p4 =  1/4*(X.^2-Y.^2) - X/2 + 1/4;
        
        clf
        if exist('u'),
            p = u(1).*p1 + u(2).*p2 +...
                u(3).*p3 + u(4).*p4;
            surf(X,Y,p);
            title('(E031): FE-solution')
            xaxis('x-axis')
            yaxis('y-axis')
        else
            subplot(2,2,1)
            surf(X,Y,p1)
            title('(E031): p_1(x,y)=-1/4(x^2-y^2)-1/2y+1/4')
            xaxis('x-axis')
            yaxis('y-axis')
            
            subplot(2,2,2)
            surf(X,Y,p2)
            title('(E031): p_2(x,y)=1/4(x^2-y^2)+1/2x+1/4')
            xaxis('x-axis')
            yaxis('y-axis')
            
            subplot(2,2,3)
            surf(X,Y,p3)
            title('(E031): p_3(x,y)=-1/4(x^2-y^2)+1/2y+1/4')
            xaxis('x-axis')
            yaxis('y-axis')
            
            subplot(2,2,4)
            surf(X,Y,p4)
            title('(E031): p_4(x,y)=1/4(x^2-y^2)-1/2x+1/4')
            xaxis('x-axis')
            yaxis('y-axis')
        end
        
    otherwise
        error('Unsupported element type!');
end
