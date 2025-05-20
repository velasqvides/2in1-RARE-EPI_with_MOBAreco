function quiver2Dspokes(traj)
% based on: https://de.mathworks.com/matlabcentral/answers/360438-i-want-to-draw-an-arrow-in-vector
%plot2DSpokeS plots de spokes for a 2D radial k space trajectory
%   AT the moment it is assumed thta the trajectory is a 2D one,
%   dimensions: 3 x nSamples x nSpokes
% in bart the first file is  'y', then 'x' and then 'z'.


nSpokes = size(traj,3);
nKspaces = size(traj,6);
nSamples = size(traj,2);

colororder(summer(nSpokes)); hold on; axis equal; grid on;

for k = 1:nKspaces
    for i = 1:nSpokes
        QX=traj(2,[1,nSamples],i,1,1,k)';
        QY=traj(1,[1,nSamples],i,1,1,k)';
        Q=[QX QY];
        [~,UV] = gradient(Q);                                   
        UVX = [UV(1,1); 0];                                     
        UVY = [UV(1,2); 0];
        quiver(QX, QY, UVX, UVY, 0,'LineWidth',2.5,'MaxHeadSize',0.1)
        hold on;
         % Adding the number next to the tip of each vector
            tipX = QX(1) + UVX(1);
            tipY = QY(1) + UVY(1);
            offset = 0.1; % Adjust this offset for better visibility
            % text(tipX + offset, tipY + offset, num2str(i), 'FontSize', 8, 'Color', 'k');
    end
end
axis equal;
xlabel('K_{x}', 'FontSize', 14, 'FontWeight', 'bold');

% Set y-axis label with bold and larger font size
ylabel('K_{y}', 'FontSize', 14, 'FontWeight', 'bold');

end