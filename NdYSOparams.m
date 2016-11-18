function Sys = NdYSOparams(Sys,param)
% Appends the chosen parameters to struct Sys
% 
% Allowed parameters:
% 'Maier-FlaigTensor'
% 'Maier-FlaigPrincipal'
% 'Wolfowicz'

Sys.Nucs = '145Nd';
Sys.I = 3.5;
Sys.S = 0.5;

switch param
    case 'Maier-FlaigTensor'
        % Site 1 is occupied to a much higher degree
        % than site 2 for Nd:YSO
        % Basis [D1 D2 b]
        % gyromagnetic ratio
        Sys.g = [1.30 0.62 0.22;
                 0.62 -2.07 1.62;
                 0.22 1.62 -2.86];
        Sys.A = [-37.1 -99.9 -83.4;
                 -99.9 -589.2 169.4;
                 -83.4 169.4 -678.4];
        if isfield(Sys,'gFrame')
            Sys = rmfield(Sys,'gFrame');
        end
        if isfield(Sys,'AFrame')
            Sys = rmfield(Sys,'AFrame');
        end
    
    case 'Maier-FlaigGavin'
        % Although Q*D*Q' gives the same tensor as above,
        % The results from these don't match M-Ftensor
        Sys.g = [-4.1407 1.4786 -0.9679];
        Sys.A = [-809.0 -0.183 -495.5];
        
        gFrame_zyz = [0.1708 0.9153 1.2319];
        AFrame_zyz = [-0.2127 0.9474 1.9135];
        
        Sys.gFrame = gFrame_zyz;
        Sys.AFrame = AFrame_zyz;      
        
        
    case 'Maier-FlaigPrincipal'
        Sys.g = [-0.96 1.48 -4.14];
        Sys.A = [-0.146 -495.5 -808.0]; %MHz
        
        gFrame_zyz = [12.4 38.6 -86.4]*(pi/180);
        AFrame_zyz = [106.4 37.5 -90.7]*(pi/180);
        
        % Sys.gFrame = eulang(inv(erot(gFrame_zyz))); % try inverting matrix
        % Sys.AFrame = eulang(inv(erot(AFrame_zyz)));
        
        Sys.gFrame = gFrame_zyz;
        Sys.AFrame = AFrame_zyz;
    
    case 'Wolfowicz'
        Sys.g = [-1.49 -0.98 -4.17];
        Sys.A = [398 0.1 827]; %MHz
        
        %gFrame_zxz = [192 39 183]*(pi/180); % Original euler convention
        %AFrame_zxz = [154 34 200]*(pi/180); % different to easyspin
        
        gFrame_zyz = [1.7802 -0.6807 -1.5184]; % Converted szxz->szyz
        AFrame_zyz = [1.1170 -0.5934 -1.2217]; % using euler.py
        
        %gFrame_zyz = [-1.5184 -0.6807 1.7802]; % Converted szxz->rzyz
        %AFrame_zyz = [-1.2217 -0.5934 1.1170]; % using euler.py
       
        %gFrame_zyz = [-1.3614 -0.6807 1.6232]; % Converted rzxz->rzyz
        %AFrame_zyz = [-2.0258 -0.5934 1.9199]; % using euler.py
        
        Sys.gFrame = gFrame_zyz;
        Sys.AFrame = AFrame_zyz;
end