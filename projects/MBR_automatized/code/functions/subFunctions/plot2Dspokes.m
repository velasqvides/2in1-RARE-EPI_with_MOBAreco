function plot2Dspokes(traj)

%plot2DSpokeS plots de spokes for a 2D radial k space trajectory
%   AT the moment it is assumed thta the trajectory is a 2D one,
%   dimensions: 3 x nSamples x nSpokes
% in bart the firt file is  'y', then 'x' and then 'z'.


nSpokes = size(traj,3);
nKspaces = size(traj,6);
colorI='r';
for k = 1:nKspaces
    for i = 1:nSpokes
%         if k==1 && ismember(i,firstSpokes)
%            colorI='k'; 
        
        plot(traj(2,:,i,1,1,k), traj(1,:,i,1,1,k),sprintf('.%s',colorI),'MarkerSize',10);
        hold on;
%         end
% colorI='g';
    end
end
axis equal;
xlabel('K_{x}'); ylabel('K_{y}');

end


    
