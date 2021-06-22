function PLOT(AFT_stack, xCoord, yCoord, nose_x)
clf;
% fig = figure;
% fig.Color = 'white'; hold on;
for i = 1:size(AFT_stack,1)
    node1 = AFT_stack(i,1);
    node2 = AFT_stack(i,2);
   
    xx = [xCoord(node1),xCoord(node2)];
    yy = [yCoord(node1),yCoord(node2)];
    
    if AFT_stack(i,7)==3
        plot( xx, yy, '-k','LineWidth',1.5); hold on;
    else
        plot( xx, yy, '-r');hold on;
    end  
end

axis equal;
% axis([-1.2 1.2 -0.7 0.7])
axis([nose_x-0.5, nose_x+1.5, -0.7, 0.7])

pause(0.001)

