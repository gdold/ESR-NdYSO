function [init_rotm, mode] = setInitialAxes(magaxis,MWaxis)
% Sets initial magnetic field axis and MW axis 
% along given crystal axis denoted by x,y,z or D1,D2,b
% by outputting a rotation matrix and mode
% to transform the crystal frame to the
% lab x,y,z. 
%
% Using x,y,z to denote the crystal axes is to be avoided
% due to confusion with lab x,y,z but is done here for
% compatibility with other frame conventions, which can be added
% by editing the first sets of if/else statements
%
% Only works for MWaxis and MW axis that are along the
% crystal principal axes.

deg = pi/180;
rt = 1.0/sqrt(2);

% Sanitise input
if strcmp(magaxis,'x') || strcmp(magaxis,'D1')
    magaxis = 'x';
elseif strcmp(magaxis,'y') || strcmp(magaxis,'D2')
    magaxis = 'y';
elseif strcmp(magaxis,'z') || strcmp(magaxis,'b')
    magaxis = 'z';
else
    disp(['Magaxis "',magaxis,'" not understood, init_rotm is identity'])
    init_rotm = eye(3);
    mode = 'perpendicular';
    return
end
    
if strcmp(MWaxis,'x') || strcmp(MWaxis,'D1')
    MWaxis = 'x';
elseif strcmp(MWaxis,'y') || strcmp(MWaxis,'D2')
    MWaxis = 'y';
elseif strcmp(MWaxis,'z') || strcmp(MWaxis,'b')
    MWaxis = 'z';
else
    disp(['MWaxis "',MWaxis,'" not understood, init_rotm is identity'])
    init_rotm = eye(3);
    mode = 'perpendicular';
    return
end


switch magaxis
    case 'x'
        init_rotm = rotaxi2mat([rt,0,rt],180*deg);
        
        switch MWaxis
            case 'x'
                mode = 'parallel';
                return
            case 'y'
                mode = 'perpendicular';
                init_rotm = rotaxi2mat([0,0,1],90*deg)*init_rotm;
                return
            case 'z'
                mode = 'perpendicular';
                return
        end        

    case 'y'
        init_rotm = rotaxi2mat([0,rt,rt],180*deg);
        
        switch MWaxis
            case 'x'
                mode = 'perpendicular';
                init_rotm = rotaxi2mat([0,0,1],180*deg)*init_rotm;
                return
            case 'y'
                mode = 'parallel';
                return
            case 'z'
                mode = 'perpendicular';
                init_rotm = rotaxi2mat([0,0,1],-90*deg)*init_rotm;
                return
        end        

    
    case 'z'
        init_rotm = eye(3);
        
        switch MWaxis
            case 'x'
                mode = 'perpendicular';
                return
            case 'y'
                mode = 'perpendicular';
                init_rotm = rotaxi2mat([0,0,1],-90*deg)*init_rotm;
                return
            case 'z'
                mode = 'parallel';
                return
        end        
end