
%%%%%%%%%%%%% Checkpoint Setup %%%%%%%%%%%%%%%%%%%%%%%%%
% numIter = 10;
% startIter = 1;
% checkpointFilename = 'checkpoint.mat';
% 
%  if exist(checkpointFilename, 'file')
%     s = load(checkpointFilename);
%     startIter = s.i;
%     fprintf('Restarting from iteration %d\n', startIter);
%  end
% 
%   for i = startIter:numIter
%     fprintf('Starting iteration %d\n', i);
% %     expensiveComputation();
% %     save(checkpointFilename, 'i');
% %   end
%   % We succefully finished. Let's delete our checkpoint file
%   delete(checkpointFilename);
% 
%  % function expensiveComputation()
%     % Pretend to do lots of work!
%     pause(1);
%  % end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



s = load('./data/simres1215.mat')
s.simresults(1,end)