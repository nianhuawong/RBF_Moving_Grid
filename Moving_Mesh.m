clc;clear;close all;
[WALL,Coord,Grid,wallNodes]  = read_grid('./naca0012-tri.cas', 0);
xCoord = Coord(:,1);
yCoord = Coord(:,2);
nNodes = size(xCoord,1);
nWallNodes = size(wallNodes,2);

fig = figure;
fig.Color = 'white'; hold on;
%% 参数
lamda = 1;      %波长
c = 0.1;        %波速
v = -0.2;       %游动速度
T = 2.0;        %周期
t = 0;          %起始时间
dt = 0.5;       %时间间隔
r0 = 10.0;      %紧支半径
basis = 11;      %基函数类型
while t < 10
    dy = zeros(nWallNodes,1);
    xCoord_new = xCoord;
    yCoord_new = yCoord;
    
    t = t + dt
%% 对物面点进行变形    
    for i = 1: nWallNodes
        wallIndex = wallNodes(i);
        xCoord_new(wallIndex) = xCoord(wallIndex) + v * t;
    end
    nose_x = min(xCoord_new(wallNodes));
    
    for i = 1: nWallNodes
        wallIndex = wallNodes(i);
        x = xCoord_new(wallIndex) - nose_x;
        y = yCoord(wallIndex);
        A = min( 1, t / T ) * ( 0.02 - 0.0825 * x + 0.1625 * x * x );
        dy(i) = A * sin( 2 * pi / lamda * ( x - c * t ));
        yCoord_new(wallIndex) = sign(yCoord(wallIndex)) * abs(yCoord(wallIndex)) + dy(i);
    end
     
%     PLOT(WALL, xCoord_new, yCoord_new)
%% 计算权重系数矩阵W   
    fai = zeros(nWallNodes,nWallNodes);
    for i = 1: nWallNodes
        wallIndex = wallNodes(i);
        x1 = xCoord(wallIndex);
        y1 = yCoord(wallIndex);
        for j = 1:nWallNodes
            wallIndex2 = wallNodes(j);
            x2 = xCoord(wallIndex2);
            y2 = yCoord(wallIndex2);
            dis = sqrt( ( x1 - x2 )^2 + ( y1 - y2 )^2 ) + 1e-40;
            fai(i,j) = RBF_func(dis, r0, basis);
        end
    end
    
    W = fai \ dy;
%% 利用W计算内场点的位移    
    fai = zeros(1,nWallNodes);
    for i = 1:nNodes
        xNode = xCoord(i);
        yNode = yCoord(i);
%         plot(xNode, yNode, 'b*');
        if( sum(i==wallNodes)~= 0 )
            continue;
        end
        for j = 1:nWallNodes
            wallIndex = wallNodes(j);
            xw = xCoord(wallIndex);
            yw = yCoord(wallIndex);
            
            dis = sqrt( ( xNode - xw )^2 + ( yNode - yw )^2 ) + 1e-40;
            
            fai(1,j) = RBF_func(dis, r0, basis);
        end
        
        dy = fai(1,:) * W;
        yCoord_new(i) = yCoord_new(i) + dy;
        
        dx = v * t;
        xCoord_new(i) = xCoord_new(i) + dx;
%         xNode = xCoord_new(i);
%         yNode = yCoord_new(i);
%         plot(xNode, yNode, 'ro');
%         kkk = 1;
    end
    
    PLOT(Grid, xCoord_new, yCoord_new, nose_x)
end

