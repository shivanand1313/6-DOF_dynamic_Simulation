function ctrl = cont(t,~)

thrust = 756.4639; % N

%% trim condition
dele_trim = 0.1239; % rad
dela_trim = 0; % rad
delr_trim = 0; % rad

%% no control
% offset = 0; % rad

%% sinusoidal inputs
% offset = 0.02*sin(pi/2*t); % rad

%% doublet inputs
% if t < 1
%     offset = 0; % rad
% elseif t < 3
%     offset = 0.02; % rad    
% elseif t < 5
%     offset = -0.02; % rad
% else
%     offset = 0; % rad
% end

%% 3-2-1-1 inputs
% if t < 1
%     offset = 0; % rad
% elseif t < 4
%     offset = 0.02; % rad
% elseif t < 6
%     offset = -0.02; % rad
% elseif t < 7
%     offset = 0.02; % rad
% elseif t < 8
%     offset = -0.02; % rad
% else
%     offset = 0; % rad
% end

dele = dele_trim;
% dele = dele_trim + offset;
dela = dela_trim;
% dela = dela_trim + offset;
delr = delr_trim;
% delr = delr_trim + offset;

ctrl = [dele, dela, delr, thrust];
end