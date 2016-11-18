function [init_rotm, mode] = setInitialAxes(magaxis,MWaxis)
% Sets initial magnetic field axis and MW axis 
% along given crystal axis denoted by u,v,w or D1,D2,b
% by outputting a rotation matrix and mode
% to transform the crystal frame to the
% lab x,y,z. 
%
% Only works for MWaxis and MW axis that are along the
% crystal principal axes.

deg = pi/180;
rt = 1.0/sqrt(2);

% Sanitise input
if strcmp(magaxis,'u') || strcmp(magaxis,'D1')
    magaxis = 'u';
elseif strcmp(magaxis,'v') || strcmp(magaxis,'D2')
    magaxis = 'v';
elseif strcmp(magaxis,'w') || strcmp(magaxis,'b')
    magaxis = 'w';
else
    disp(['Magaxis "',magaxis,'" not understood, init_rotm is identity'])
    init_rotm = eye(3);
    mode = 'perpendicular';
    return
end
    
if strcmp(MWaxis,'u') || strcmp(MWaxis,'D1')
    MWaxis = 'u';
elseif strcmp(MWaxis,'v') || strcmp(MWaxis,'D2')
    MWaxis = 'v';
elseif strcmp(MWaxis,'w') || strcmp(MWaxis,'b')
    MWaxis = 'w';
else
    disp(['MWaxis "',MWaxis,'" not understood, init_rotm is identity'])
    init_rotm = eye(3);
    mode = 'perpendicular';
    return
end


switch magaxis
    case 'u'
        init_rotm = rotaxi2mat([rt,0,rt],180*deg);
        
        switch MWaxis
            case 'u'
                mode = 'parallel';
                return
            case 'v'
                mode = 'perpendicular';
                init_rotm = rotaxi2mat([0,0,1],90*deg)*init_rotm;
                return
            case 'w'
                mode = 'perpendicular';
                return
        end        

    case 'v'
        init_rotm = rotaxi2mat([0,rt,rt],180*deg);
        
        switch MWaxis
            case 'u'
                mode = 'perpendicular';
                init_rotm = rotaxi2mat([0,0,1],180*deg)*init_rotm;
                return
            case 'v'
                mode = 'parallel';
                return
            case 'w'
                mode = 'perpendicular';
                init_rotm = rotaxi2mat([0,0,1],-90*deg)*init_rotm;
                return
        end        

    
    case 'w'
        init_rotm = eye(3);
        
        switch MWaxis
            case 'u'
                mode = 'perpendicular';
                return
            case 'v'
                mode = 'perpendicular';
                init_rotm = rotaxi2mat([0,0,1],-90*deg)*init_rotm;
                return
            case 'w'
                mode = 'parallel';
                return
        end        
end