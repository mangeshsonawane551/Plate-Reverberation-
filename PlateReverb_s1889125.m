%-------------------------------------------------------------------------%
 % Program Details: This program creates modal reverberator based on a 
 % physical model of a thin metal plate with simply supported boundary
 % conditions. The method used is the modal synthesis. Modes are calculated
 % and used in the modal update. Response of default plate, impulse
 % response and spectrogram can be extracted through this program.
 % Parameters for glass plate have also been used in this program.
 %
%-------------------------------------------------------------------------%

clc;
clear all;
tic
%-------------------------------------------------------------------------%
              % Selection for operating plate reverb algorithm
%-------------------------------------------------------------------------%
input_mode=1;       %1:Audio input 
                    %2:Impulse input

input_material=2;   %1:Steel (default). 
                    %2:Glass . (extra)
                    
operation_type=1;   %1:ndgrid mechanism use
                    %2:for-loop use
                    
 scheme=1;          %1: Accurate Scheme
                    %2: Exact Scheme
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
if input_material==1
  % Length x
    Lx = 2; 
  % Lenght y
    Ly = 1; 
  % tension (N)
    T = 700; 
  % thickness (m)
    a = 5e-4; 
  % density (kg/m^3)
    rho = 8e3;
  % Minimum T60
    T60min = 1; 
  % Maximum T60
    T60max =8 ; 
  % Young's Modulus
    E = 2e11;
  %poisson ratio;
    pr= 0.3; 
  % naming label input while saving the file
    material_type='default';
    
    % Glass Plate
elseif input_material==2
  % Length x
    Lx = .8;
  % Lenght y
    Ly = .7;
  % tension (N)
    T = 500;
  % thickness (m)
    a = 12e-3;
  % density (kg/m^3)
    rho = 2560;
  % Minimum T60
    T60min=4;
  % Maximum T60
    T60max=8;
  % Young's Modulus
    E=1e10;
  %poisson ratio;
    pr=0.25;
  % naming label input while saving the file
    material_type='extra';
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
if input_mode==1
   %Input force
    [F,SR]=audioread('Cath_Very_Short2.wav' );
   % Checks if stereo and converts to mono
      if size(F,2) == 2   
    %convert to mono
        F = sum(F,2) / size(F,2);                
      end
      % naming convention used for file
        audio_type='dryout';
% force input is an impulse        
elseif input_mode==2
    F=zeros(10,1);
    F(1)=1;
    audio_type= 'IR';
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
% Coefficient4 
cf4 = ((k^2)/(rho*a))./cf1;
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
force_in=Fin*k^2/(rho*a*h);
%Update waitbar
waitbar(.40,f,'Started calculating modal update');
pause(1)

waitbar(.55,f,'calculating modal update..Please wait..')

A=sqrt((sig_q.^2)-omega_sq).*k;
cf1_e = exp(-sig_q.*k).*(exp(A)+ exp(-A));
cf2_e = exp(-2.*sig_q*k);

 for i=1:size_o
     % Accurate scheme
    if scheme==1 
     %Scheme
     p0 = p1.*cf2 - p2.*cf3;  
     %Exact Scheme
    elseif scheme==2 
      p0=p1.*cf1_e - p2.*cf2_e;
    end
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
 if input_mode==2
   %update waitbar
    waitbar(.85,f,'converting to mono');
   % convert stereo to mono
    out = sum(out,2) / size(out,2);  
    NFFT=2*1024;
   %Window length
    window_length=floor(2*round(0.031*SR));
   %update waitbar 
   waitbar(.88,f,'spectrogram on its way!');
   % create a hann window 
       for n=1:(window_length)
            window(n)= 0.5*(1-cos(2*pi*n/window_length)); 
       end
   %number of windows samples without overlapping    
    overlap=floor(0.45*window_length); 
   %construct spectrogram
    [S,F,T] = spectrogram(out,window,window_length-overlap,NFFT,SR); 
   %Spectrogram properties 
    [Nf,Nw]=size(S);
    figtitle1 = 'Impulse response Spectrogram';
    figure('name',figtitle1)
    factor_colour=0.000001;
    S_one_sided=max(S(1:fix(length(F)/2),:),factor_colour); %keep only the positive frequency
    pcolor(T,F(1:fix(Nf/2)),10*log10(abs(S_one_sided)));
    shading interp;
    colormap('jet');
    colorbar;
    title('Impulse response');
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    image_filename=sprintf('PlateReverb_s1889125_%s_IRspec.png',material_type);
    saveas(gcf,image_filename);
 end
 
 %update wait bar
 waitbar(.90,f,'Sound ready to play');
 % File name for audio file to be saved
 filename=sprintf('PlateReverb_s1889125_%s_%s.wav',material_type,audio_type);
 %write audio 
 audiowrite(filename,out,SR);  
 %update waitbar
 waitbar(.91,f,'playing reverb of choosen plate');
 % Play sound
 soundsc(out,SR)
 %Close update bar
 close(f)
 %end Timer
 toc
 
 
 
 
 
 
 
 
 
            