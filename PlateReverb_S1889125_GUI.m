
classdef PlateReverb_S1889125_GUI < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                   matlab.ui.Figure
        PlaysoundButton            matlab.ui.control.Button
        ReverbtypeDropDownLabel    matlab.ui.control.Label
        ReverbtypeDropDown         matlab.ui.control.DropDown
        ReverbTimeKnobLabel        matlab.ui.control.Label
        ReverbTimeKnob             matlab.ui.control.Knob
        thicknessofplateKnobLabel  matlab.ui.control.Label
        thicknessofplateKnob       matlab.ui.control.Knob
        SwitchLabel                matlab.ui.control.Label
        Switch                     matlab.ui.control.Switch
        PLATEREVERBLabel           matlab.ui.control.Label
    end

    methods (Access = private)

        % Button pushed function: PlaysoundButton
        function PlaysoundButtonPushed(app, event)
           
            
            
            Reverbtype = app.ReverbtypeDropDown.Value;
            switch_auimp = app.Switch.Value;
            MinReverbTime = app.ReverbTimeKnob.Value;
            platethickness = app.thicknessofplateKnob.Value;
%             convolvecheck = app.convolveaudioSwitch.Value;
            
         
%-------------------------------------------------------------------------%
              % Selection for operating plate reverb algorithm
%-------------------------------------------------------------------------%

input_mode=switch_auimp;       %1:Audio input 
                               %2:Impulse input

input_material=Reverbtype;   %1:Steel (default). 
                             %2:Glass . (extra)
                             %3:Wood
                             %4:Ceramic
                    
operation_type=1;   %1:ndgrid mechanism use
                    %2:for-loop use
%    if input_mode=='2'                 
%  convolutionneeded=  convolvecheck;% 1:Off 2:On       
%    end         
%-------------------------------------------------------------------------%
              % physical  parameters 
%-------------------------------------------------------------------------%
% sample rate (Hz)
SR = 44100*1;              
% coordinate of excitation (normalised, 0-1)                     
xi = [0.2,0.6] ;
% coordinate of output (normalised, 0-1)
xo = [0.2,0.8;0.7,0.3]; 

% Steel Material
if Reverbtype=='1'
  % Length x
    Lx = 2; 
  % Lenght y
    Ly = 1; 
  % tension (N)
    T = 700; 
  % thickness (m)
    a = platethickness*10^-4; 
  % density (kg/m^3)
    rho = 8*10^3;
  % Minimum T60
    T60min = MinReverbTime; 
  % Maximum T60
    T60max = 10 ; 
  % Young's Modulus
    E = 2*10^11;
  %poisson ratio;
    pr= 0.3; 
 
   
 
    % Glass Plate
elseif Reverbtype=='2'
  % Length x
    Lx = .8;
  % Lenght y
    Ly = .7;
  % tension (N)
    T = 500;
  % thickness (m)
    a = platethickness*10^-4;
  % density (kg/m^3)
    rho = 2560;
  % Minimum T60
    T60min=MinReverbTime;
  % Maximum T60
    T60max=10;
  % Young's Modulus
    E=3*10^10;
  %poisson ratio;
    pr=0.25;
  
    

elseif Reverbtype=='3'  %wood
  % Length x
    Lx = .6;
  % Lenght y
    Ly = .3;
  % tension (N)
    T = 2000;
  % thickness (m)
    a = platethickness*10^-4;
  % density (kg/m^3)
    rho = 800;
  % Minimum T60
    T60min=MinReverbTime;
  % Maximum T60
    T60max=10;
  % Young's Modulus
    E=10*10^9;
  %poisson ratio;
    pr=0.45;

    
    
elseif Reverbtype=='4' %wood
  % Length x
    Lx = .7;
  % Lenght y
    Ly = .6;
  % tension (N)
    T = 900;
  % thickness (m)
    a = platethickness*10^-4;
  % density (kg/m^3)
    rho = 3800;
  % Minimum T60
    T60min=MinReverbTime;
  % Maximum T60
    T60max=10;
  % Young's Modulus
    E =4*10^11;
  %poisson ratio;
    pr=0.27;
  
end   
%-------------------------------------------------------------------------%
                        % Derived parameters 
%-------------------------------------------------------------------------%
  
 
  % Time Step
    k = 1/SR;  
  % wave speed
     c = sqrt(T/(rho*a)); 
  % Stiffness
    kappa = sqrt( (E*a^2)/(12*rho*(1-pr^2))); 
  % omega maximum valued squared
    wmax_sq = 4/k^2;   
  % Sigma(1)min
    sigma1 = 6*log(10)/T60max; 
  % Sigma(Q)max
    sigmaQ = 6*log(10)/T60min; 
  % Beta square at q=1
    beta1_sq = (1*pi/Lx)^2 + (1*pi/Ly)^2; 
  % Beta at max Q using quadratic solution
    betaQ_sq = ((-c^2)+sqrt(c^4+4*kappa^2*wmax_sq))/(2*kappa^2);
  % epsilon1
    ep1 = (sigma1-sigmaQ)/(beta1_sq-betaQ_sq); 
  % epsilon0
    ep0 = sigma1-(beta1_sq*((sigma1-sigmaQ)/(beta1_sq-betaQ_sq)));
    assert(ep1>=0)
    assert(ep0>=0)
  % Minimum grid spacing
    hmin = sqrt(0.5* (c.^2*k^2+sqrt(c.^4*k^4+16*kappa.^2.*k.^2)) ); 
  % Number of grid points at x  
    Nx = floor(Lx/hmin)-1 ; 
  % Number of grid points at y
    Ny=floor(Ly/hmin)-1;
  % Actual total geid spacing   
    h = sqrt(Lx*Ly)/sqrt(Nx*Ny); 
  % grid index of excitation at x and y coordinate
    li_x=floor(xi(1,1)*Nx);      
    li_y=floor(xi(1,2)*Ny);       
  % grid index of output at x and y coordinate at left channel
    lo_x1 = floor(xo(1,1)*Nx);    
    lo_y1=floor(xo(1,2)*Ny);
   % grid index of output at x and y coordinate at left channel  
    lo_x2=floor(xo(2,1)*Nx);
    lo_y2=floor(xo(2,2)*Ny);  
%-------------------------------------------------------------------------%
                        % Force
%-------------------------------------------------------------------------%

% force input is a audio signal
if input_mode=='1'
   %Input force
    [F,SR]=audioread('Cath_Very_Short2.wav' );
   % Checks if stereo and converts to mono
      if size(F,2) == 2   
    %convert to mono
        F = sum(F,2) / size(F,2);                
      end
      
% force input is an impulse        
elseif input_mode=='2' 
    F=zeros(10,1);
    F(1)=1;
   
end

%-------------------------------------------------------------------------%
                        %  Calculation 
%-------------------------------------------------------------------------%
% Waitbar initialising
f=waitbar(0,'Total time estimated is 54 sec..Plate reverb in progress....');
   pause(1)
 
% Maximum number of modes calculation
coef_BetaQ = sqrt((-c^2)+sqrt(c^4+4*kappa^2*wmax_sq)/(2*kappa^2)); 
 Qmax_x = floor((Lx/pi)*coef_BetaQ);
 Qmax_y = floor((Ly/pi)*coef_BetaQ);
 
% calculation of parameters using 'for loop'
if operation_type==1
for i=1:(Qmax_x)
    for j=1:(Qmax_y)
         % temporary beta square
            beta_sq_temp=((i*pi/Lx)^2)+((j*pi/Ly)^2);
         % temporary omega square   
            omega_sq_temp=c^2*beta_sq_temp + kappa^2*beta_sq_temp^2;
         % if omega square is less than equal to maximum omega square condition check   
        if omega_sq_temp<=wmax_sq
         % Beta square
            beta_sq(i,j)=beta_sq_temp;
         % Omega square   
            omega_sq(i,j)=omega_sq_temp;
         %Sigma q    
            sig_q(i,j)= ep0+ep1*beta_sq(i,j);
         % checks if sigma q is less than omega    
            assert((sig_q(i,j))<=sqrt(omega_sq(i,j)));
         % modal calculation at input    
            phi_in(i,j)=((2/sqrt(Lx*Ly)) .* sin((i)*i.*li_x.*pi/Lx).*sin(j*j.*li_y.*pi/Ly));
         % modal calculation at output left channel 
            phi_outl(i,j)=((2/sqrt(Lx*Ly)) .* sin((i).*lo_x1.*pi/Lx).*sin(j.*lo_y1.*pi/Ly));
         % modal calculation at output right channel   
            phi_outr(i,j)=((2/sqrt(Lx*Ly)) .* sin((i).*lo_x2.*pi/Lx).*sin(j.*lo_y2.*pi/Ly));
        end
    end
end

%calculation using 'NDGRID'
elseif operation_type==2
  % Initialise ndgrid for given number of modes 
   [qx,qy]=meshgrid(1:Qmax_x,1:Qmax_y);
  % Calculate beta square in a temporary variable  
    beta_sq_temp=((qx*pi/Lx).^2)+((qy*pi/Ly).^2);
  % Calculate omega square in a temporary variable  
    omega_sq_temp=c^2.*beta_sq_temp + kappa^2.*beta_sq_temp.^2;
  % Set the condition  
    M=omega_sq_temp<=wmax_sq;
  % Filter beta square with the given condition   
    beta_sq=(beta_sq_temp(M));
  % Filter omega square with the given condition  
    omega_sq=omega_sq_temp(M);
  % Calculate sigma q 
    sig_q= ep0+ep1*beta_sq;
  % Calculate  modal calculation at input and filter it with the condition 
    phi_in=((2/sqrt(Lx*Ly)) .* sin(qx.*li_x.*pi/Lx).*sin(qy.*li_y.*pi/Ly));
    phi_in=phi_in(M);
 % Calculate  modal calculation at output left channel  and filter it with the condition
    phi_outl=((2/sqrt(Lx*Ly)) .* sin(qx.*lo_x1.*pi/Lx).*sin(qy.*lo_y1.*pi/Ly));
    phi_outl=phi_outl(M);
 % Calculate  modal calculation at output right channel  and filter it with the condition  
   phi_outr=((2/sqrt(Lx*Ly)) .* sin(qx.*lo_x2.*pi/Lx).*sin(qy.*lo_y2.*pi/Ly));
   phi_outr=phi_outr(M);
end
%-------------------------------------------------------------------------%
                        %  Modal update/ Finite Difference
%-------------------------------------------------------------------------%
% Coefficient1
cf1 = (1+sig_q.*k);  
% Coefficient2
cf2 = (2-k^2.*omega_sq)./cf1; 
% Coefficient3
cf3 = (1-(sig_q.*k ))./cf1; 


%Total size of output if length of force + lenght of T60max
size_o = length(F)+T60max*SR;
% Initialise output vector for left  channel
out_l=zeros(size_o,1);
% Initialise output vector for left  channel
out_r=zeros(size_o,1);
%initial velocity
v0=0; 
%initial displacement
x0=0; 
%value at time step n=1
p2=x0;
%value at time step n=2
p1=x0+(k*v0); 
% Create force vector with addition of zeros at end to make the size of
% force equal to size_o
Fin=cat(1,F,zeros(T60max*SR,1));
%force coefficient
force_in=(Fin*k^2/(rho*a*h));
%Update waitbar
waitbar(.40,f,'Started calculating modal update');
pause(1)

waitbar(.55,f,'calculating Modal update/ FD..Please wait..')


 for i=1:size_o
    
     %Scheme
     p0 = p1.*cf2 - p2.*cf3;     
     %Input
     p0= p0+force_in(i).*phi_in;     
     %Output at left channel
     out_l(i)=sum(p0.*phi_outl,'all');
     %output at right channel
     out_r(i)=sum(p0.*phi_outr,'all');
     % update
     p2=p1;
     p1=p0;  

end
 
% Update wait bar
waitbar(.70,f,'calculation done');
pause(0.02)
waitbar(.75,f,'creating a stereo sound');
pause(0.02)

% Normalize output channel 
 out_l = -1 + (2)*(out_l-min(out_l)) / ( max(out_l)-min(out_l));
 out_r = -1 + (2)*(out_r-min(out_r)) / ( max(out_r)-min(out_r));
% create stereo  
 out=[out_l out_r];

% if force choosen is impulse then spectrogram initialisation

    
  
 %update wait bar
 waitbar(.90,f,'Sound ready to play'); 
 %update waitbar
 waitbar(.91,f,'playing reverb of choosen plate');
 % Play sound
 soundsc(out,SR)
 %Close update bar
 close(f)
        end
        % Value changed function: ReverbtypeDropDown
        function ReverbtypeDropDownValueChanged(app, event)
            value = app.ReverbtypeDropDown.Value;
            
        end

        % Value changed function: ReverbTimeKnob
        function ReverbTimeKnobValueChanged(app, event)
            value = app.ReverbTimeKnob.Value;
            
        end

        % Value changed function: thicknessofplateKnob
        function thicknessofplateKnobValueChanged(app, event)
            value = app.thicknessofplateKnob.Value;
            
        end

        % Value changed function: Switch
        function SwitchValueChanged(app, event)
            value = app.Switch.Value;
            
        end

%         % Value changed function: convolveaudioSwitch
%         function convolveaudioSwitchValueChanged(app, event)
%             value = app.convolveaudioSwitch.Value;
%             
%         end
    end

    % App initialization and construction
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure
            app.UIFigure = uifigure;
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'UI Figure';

            % Create PlaysoundButton
            app.PlaysoundButton = uibutton(app.UIFigure, 'push');
            app.PlaysoundButton.ButtonPushedFcn = createCallbackFcn(app, @PlaysoundButtonPushed, true);
            app.PlaysoundButton.Position = [277 42 100 22];
            app.PlaysoundButton.Text = 'Play sound';

            % Create ReverbtypeDropDownLabel
            app.ReverbtypeDropDownLabel = uilabel(app.UIFigure);
            app.ReverbtypeDropDownLabel.HorizontalAlignment = 'right';
            app.ReverbtypeDropDownLabel.Position = [47 353 70 22];
            app.ReverbtypeDropDownLabel.Text = 'Reverb type';

            % Create ReverbtypeDropDown
            app.ReverbtypeDropDown = uidropdown(app.UIFigure);
            app.ReverbtypeDropDown.Items = {'Steel plate', 'Glass plate', 'Wood', 'Ceramic'};
            app.ReverbtypeDropDown.ItemsData = {'1', '2', '3', '4'};
            app.ReverbtypeDropDown.ValueChangedFcn = createCallbackFcn(app, @ReverbtypeDropDownValueChanged, true);
            app.ReverbtypeDropDown.Position = [132 353 100 22];
            app.ReverbtypeDropDown.Value = '3';

            % Create ReverbTimeKnobLabel
            app.ReverbTimeKnobLabel = uilabel(app.UIFigure);
            app.ReverbTimeKnobLabel.HorizontalAlignment = 'center';
            app.ReverbTimeKnobLabel.Position = [153 124 73 22];
            app.ReverbTimeKnobLabel.Text = 'Reverb Time';

            % Create ReverbTimeKnob
            app.ReverbTimeKnob = uiknob(app.UIFigure, 'continuous');
            app.ReverbTimeKnob.Limits = [1 10];
            app.ReverbTimeKnob.ValueChangedFcn = createCallbackFcn(app, @ReverbTimeKnobValueChanged, true);
            app.ReverbTimeKnob.Position = [149 172 82 82];
            app.ReverbTimeKnob.Value = 5;

            % Create thicknessofplateKnobLabel
            app.thicknessofplateKnobLabel = uilabel(app.UIFigure);
            app.thicknessofplateKnobLabel.HorizontalAlignment = 'center';
            app.thicknessofplateKnobLabel.Position = [385 118 100 22];
            app.thicknessofplateKnobLabel.Text = 'thickness of plate';

            % Create thicknessofplateKnob
            app.thicknessofplateKnob = uiknob(app.UIFigure, 'continuous');
            app.thicknessofplateKnob.Limits = [1 10];
            app.thicknessofplateKnob.ValueChangedFcn = createCallbackFcn(app, @thicknessofplateKnobValueChanged, true);
            app.thicknessofplateKnob.Position = [395 174 80 80];
            app.thicknessofplateKnob.Value = 1;

            % Create SwitchLabel
            app.SwitchLabel = uilabel(app.UIFigure);
            app.SwitchLabel.HorizontalAlignment = 'center';
            app.SwitchLabel.Position = [399 374 42 22];
            app.SwitchLabel.Text = 'Switch';

            % Create Switch
            app.Switch = uiswitch(app.UIFigure, 'slider');
            app.Switch.Items = {'Audio output', 'Impulse response'};
            app.Switch.ItemsData = {'1', '2'};
            app.Switch.ValueChangedFcn = createCallbackFcn(app, @SwitchValueChanged, true);
            app.Switch.Position = [438 341 45 20];
            app.Switch.Value = '1';
                
            % Create PLATEREVERBLabel
            app.PLATEREVERBLabel = uilabel(app.UIFigure);
            app.PLATEREVERBLabel.FontName = 'Courier';
            app.PLATEREVERBLabel.FontSize = 30;
            app.PLATEREVERBLabel.FontWeight = 'bold';
            app.PLATEREVERBLabel.Position = [188 422 224 35];
            app.PLATEREVERBLabel.Text = 'PLATE REVERB';
        
        end
    end

    methods (Access = public)

        % Construct app
        function app = PlateReverb_S1889125_GUI

            % Create and configure components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end