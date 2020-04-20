%% Cannot run this without LTE Toolbox
% https://www.mathworks.com/help/lte/examples/lte-parameterization-for-waveform-generation-and-simulation.html
%  The following example shows how to create a 20MHz, QPSK, 3/4 rate
%  waveform corresponding to transmission mode 8 ('Port7-8' transmission
%  scheme) with full allocation 

% 100 RBs.  20 MHz Signal has 100 RBs 
% 12 Subcarriers per RB
% So, 100 * 12 = 1200 Subcarriers

% 10 subframes
% 2 Slots per Subframe
% 7 Symbols per Slot
% So, 10 * 2 * 7 = 140 Symbols

% 1 ms per Subframe
% So we have 10 * 1 = 10 ms of time domain signal

% fs = 30,720,000
% number of samples = 307,200 in time domain signal
% So, yes we have 10 ms of time domain signal
% 307,200 / 140 ~= 2194.28 samples per symbol

dataStream = [];     % Define the input user data stream (was [1 0 0 1] )
params = struct();          % Initialize the parameter structure
params.NDLRB = 100;         % 20 MHz bandwidth
params.CellRefP = 1;        
params.PDSCH.PRBSet = (0:params.NDLRB-1)'; % Full allocation
params.PDSCH.TargetCodeRate = 3/4; % The target code rate
params.PDSCH.TxScheme = 'Port0'; % Transmission mode 0
params.PDSCH.NLayers = 1;          % 1 layer transmission
params.PDSCH.Modulation = 'QPSK';  % Modulation scheme
params.PDSCH.NSCID = 0;            % Scrambling identity
params.PDSCH.NTxAnts = 1;          % 1 transmit antennas
params.PDSCH.W = lteCSICodebook(params.PDSCH.NLayers,...
                        params.PDSCH.NTxAnts,0).'; % Precoding matrix

% Now use lteRMCDL to populate other parameter fields
fullParams = lteRMCDL(params);
% Generate the waveform using the full parameter set 'fullParams'
[rmcwaveform, rmcgrid, rmcconfig] = lteRMCDLTool(fullParams,dataStream);
% dlWaveform is the time domain waveform, dlGrid is the resource grid and
% dlParams is the full set of parameters used in the waveform generation.

% Populated resource grid, returned as a numeric 3-D array of resource elements
% for several subframes across all configured antenna ports
mesh(abs(rmcgrid))
view(2)
title('Resource Element Grid')
ylabel('Subcarriers')
xlabel('Symbols')

disp('LTE Downlink Configuration')
rmcconfig