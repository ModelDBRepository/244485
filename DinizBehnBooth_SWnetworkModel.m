function y = DinizBehnBooth_SWnetworkModel(t, x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Human Sleep-wake regulatory network model
% Based on: Gleit, Diniz Behn and Booth, J Biol Rhythms 28:339-355, 2013
% Slightly revised in: Booth, Xique and Diniz Behn, SIAM J Appl Dyn Sys
%                      16(2):1089-1112, 2017
% 
% Revisions include:
% - 3 sleep-wake state promoting neuronal populations
% - reduced equation form for neurotransmitter concentrations 
% - parameters for time constants of homeostatic sleep drive updated
%   based on data in Rusterholz etal, SLEEP 33:491-498, 2010.
% - parameters corrected in circadian oscillator model as in Serkh and
%   Forger, PLoS Comput Biol 10:e1003523, 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fLC = x(1);
  fVLPO = x(2);
  fR = x(3);
  fSCN = x(4);
  h = x(5);
  c  =  x(6);
  xc = x(7);
  n = x(8);
  
  
  % parameter to change number of REM bouts
  % -0.3 = no REM, -0.8 = 4 REM bouts
  betaR = -0.8;            

  %parmeters
  gACHLC=0.8; gGABALC=1.5;         
  gGABAVLPO=0;  gNEVLPO=1.5;       
  gACHR=2.2;  gNER=8;  gGABAR=1.0;  
  gGSCNLC = 0.5;  gGSCNVLPO = 0.1;  gGSCNR = 0.8;
  tauLC=23;  tauR=1;  tauVLPO=10; taufSCN=0.5;
  maxLC=6;  maxR=5;  maxVLPO=5;  maxSCN=7; 
  gammaAChR=2.5;  gammaNE=5;   gammaGABA=2.5;  gammaSCN=4;   
  betaLC=-0.4;  alphaLC=0.4;     
  alphaR=0.1;   
  alphaVLPO=0.2; 
  k1 = -0.1;        
  k2 = -0.0045;
  betaSCN=-0.1; alphaSCN=0.7; 
  thetaW=2;
  tauhw=946.8;  tauhs=202.2;  hmax=323.88;  hmin=0;
      
% feedback from sleep-wake to SCN
  gAChSCN=0;
  gNESCN=0;
 
  % Circadian oscillator parameters
  taux=24.2;  alph_0=0.05; beta=0.0075;
  G=33.75;  p=0.5;  k=0.55;  I_0=9500;  mu=0.23;
 
%Steady state neurotransmitter values
  cNE=tanh(fLC/gammaNE);
  cGABA=tanh(fVLPO/gammaGABA);
  cACh=tanh(fR/gammaAChR);
  cSCN=tanh(fSCN/gammaSCN);
  
% other functions
  sw=gAChSCN*cACh+gNESCN*cNE;
  betaVLPO=k2*h+k1; %homeostatically varying VLPO threshold
  
  %input to firing rate inf functions
  LC_INPUT = gACHLC*cACh+gGSCNLC*cSCN-gGABALC*cGABA;
  VLPO_INPUT = -gNEVLPO*cNE-gGSCNVLPO*cSCN-gGABAVLPO*cGABA;
  R_INPUT =gACHR*cACh-gNER*cNE+gGSCNR*cSCN-gGABAR*cGABA;
  SCN_INPUT = c + sw;

  %population steady state response functions
  LC_inf = maxLC*0.5*(1+tanh((LC_INPUT-betaLC)/alphaLC));
  R_inf = maxR*0.5*(1+tanh((R_INPUT-betaR)/alphaR));
  fSCN_inf = maxSCN*0.5*(1+tanh((SCN_INPUT-betaSCN)/alphaSCN));
  VLPO_inf=maxVLPO*0.5*(1+tanh((VLPO_INPUT-betaVLPO)/alphaVLPO));

% Light function, set for 10 h dark: 14 h light
    per = 1440;
    t1 = 600;
  if ( mod(t,per) >= t1 )    % 600 = 10 h dark
    I = 500;
  else
    I = 0;
  end

  %Circadian Functions

  alph = (alph_0)*((I^p)/((I_0)^p));
  B_hat = G*(1-n)*alph;
  B = B_hat*(1-0.4*c)*(1-0.4*xc);

  %ODES

  y = zeros(8 ,1);

  y(1) = (LC_inf-fLC)/tauLC;
  y(2) = (VLPO_inf-fVLPO)/tauVLPO;
  y(3) = (R_inf-fR)/tauR;
  y(4) = (fSCN_inf-fSCN)/taufSCN;

  y(5) = heaviside(fLC-thetaW)*(hmax-h)/tauhw + heaviside(thetaW-fLC)*(hmin-h)/tauhs;
  y(6) = (pi/720)*(xc+B);
  y(7) = (pi/720)*(mu*(xc-((4*xc^3)/3))-c*(((24/(0.99669*taux))^2)+k*B));
  y(8) = ((alph*(1-n))-(beta*n));

end