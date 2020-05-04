clear all
close all
%clc

%warning ("off");

logging = true;     % log and print data
plotting = true;    % plot data, logging needs to be enabled
plotLength = 200;   % plot that contains this interval scans of data
plotFreq = 100;     % plot every this interval of scans


% system kinematics
sys = kinematics; % assign object

sys.joint1 = [-1,1,2]; 
sys.joint2 = [-1.2,-2,5];
sys.joint3 = [-2,-7,12];
sys.endEffectorRaw = [-0.5,-6,14];
sys.PivotAxis = [1,0,0]; 
sys.yawAxis = [0,0,1];
sys.naturalFreqX = 1.2;
sys.naturalFreqZ = 0.8;
sys.oscillationMagnitudeX = 0.001;
sys.oscillationMagnitudeZ = 0.001;
sys.loopFreq = 500;


sys = initKinematics(sys);

% tracking sensor reading and filtering
sensorUniform = 0.006; % in each read there could be such accuracy
sensorNormal = 0.01;

% imu reading and filtering % need bias
imuUniform = 0.001; % in each read there could be such uniform random in run bias
imuNormal = 50*10^-6*9.81*30; % Accelerometers noise density

%% Initialization

sensorFilter = filtering;
sensorFilter.loopFreq = sys.loopFreq;
sensorFilter.movingAvgNum = 25;
sensorFilter.lowpassCoeffX1 = [0.1 0];
sensorFilter.lowpassCoeffY1 = [0.1 0];
sensorFilter.lowpassCoeffZ1 = [0.1 0];
sensorFilter.sensorInterval = 1;
sensorFilter = initFilter(sensorFilter);


errRaw = [];
errMovingAvg = [];
errlowpass = [];
errKalman = [];

posTrue = [];
posReading = []; % sensor raw data
poslowpass = [];
posMovingAvg = [];
posKalman = [];
posKalmanPredict = [];

kalmanNorm = [];
accl = [];

lastendEffectorpos = zeros(3,3);
lastLastendEffectorpos = zeros(3,3);


% calculate next iteration
sys = calcEndEffector(sys);
    
%% Main Loop
while true

    tic; % start timer
    
    %% System Simulation
    sys.loopCount = sys.loopCount + 1;
    sys = calcEndEffector(sys);

    sensorFilter.loopCount = sensorFilter.loopCount + 1;

    % apply error on sensor readings, 3 endEffectors and their 3 axis exhibit 
    % same error characteristics
    % read value to multiple filters and filtering process
    
    sensorFilter.endEffectorReading = noiseApply(sys.endEffector1,sensorUniform,sensorNormal);
    
    sensorFilter.accl = ([sys.endEffector1]-lastendEffectorpos)-(lastendEffectorpos-lastLastendEffectorpos); % vel=pos*loopFreq;accl=vel*loopFreq;
    sensorFilter.accl = sensorFilter.accl*sys.loopFreq*sys.loopFreq; % two scans in between, the difference in velocity is one scan
    sensorFilter.accl = [noiseApply(sensorFilter.accl(1,:),imuUniform,imuNormal)];
                     
    lastLastendEffectorpos = lastendEffectorpos;
    lastendEffectorpos = [sys.endEffector1];
    
    
    %% Filter Comparison
    sensorFilter = movingAvg(sensorFilter);     % moving average filter
    sensorFilter = lowpassFilter(sensorFilter); % low pass filter
    
    
    %% Kalman Filter
    % Kalman filter prediction stage
    sensorFilter = kalmanPredict(sensorFilter); 
    
    % Kalman filter update stage
    sensorFilter = kalmanUpdateIMU(sensorFilter);
    % Kalman filter sensor update
    if mod(sys.loopCount,sensorFilter.sensorInterval)==0 % sensor update frequency
        sensorFilter = kalmanUpdatesensor(sensorFilter); 
    end

    %% Logging and Plotting
    if logging == true
        
        % endEffector 1 error of 3 axes combined log
        errRaw       = [errRaw       sum((sensorFilter.endEffectorReading-sys.endEffector1).^2)^0.5];       % raw error in this scan
        errMovingAvg = [errMovingAvg sum((sensorFilter.endEffector1AvgFiltered-sys.endEffector1).^2)^0.5];   % moving average error in this scan
        errlowpass   = [errlowpass sum((sensorFilter.endEffector1lowpassFiltered-sys.endEffector1).^2)^0.5]; 
        errKalman    = [errKalman sum((sensorFilter.endEffector1kalmanFiltered-sys.endEffector1).^2)^0.5]; 

        % position array log
        posTrue      = [posTrue sys.endEffector1(1)];
        posReading   = [posReading sensorFilter.endEffectorReading(1)];
        poslowpass   = [poslowpass sensorFilter.endEffector1lowpassFiltered(1)];
        posMovingAvg = [posMovingAvg sensorFilter.endEffector1AvgFiltered(1)];
        posKalman    = [posKalman sensorFilter.endEffector1kalmanFiltered(1)];
        posKalmanPredict = [posKalmanPredict sensorFilter.endEffector1kalmanFiltered(1)];
        
        kalmanNorm = [kalmanNorm norm(sensorFilter.P)];
        accl = [accl sensorFilter.accl(1)];
        
        fprintf("Actual: %.6f %.4f %.4f\n",sys.endEffector1(1),sys.endEffector1(2),sys.endEffector1(3))
        fprintf("Observ: %.6f %.4f %.4f\n",sensorFilter.endEffectorReading(1),sensorFilter.endEffectorReading(2),sensorFilter.endEffectorReading(3))
        fprintf("MovAvg: %.6f %.4f %.4f\n",sensorFilter.endEffector1AvgFiltered(1),sensorFilter.endEffector1AvgFiltered(2),sensorFilter.endEffector1AvgFiltered(3))
        fprintf("Lowpass: %.6f %.4f %.4f\n",sensorFilter.endEffector1lowpassFiltered(1),sensorFilter.endEffector1lowpassFiltered(2),sensorFilter.endEffector1lowpassFiltered(3))
        fprintf("Kalman: %.6f %.4f %.4f\n",sensorFilter.endEffector1kalmanFiltered(1),sensorFilter.endEffector1kalmanFiltered(2),sensorFilter.endEffector1kalmanFiltered(3))
        fprintf("Raw Err: %.5f  MovAvg Err: %.5f lowpass Err: %.5f Kalman Err: %.5f\n\n",errRaw(end),errMovingAvg(end),errlowpass(end),errKalman(end)) % error of each filter in current scan

    end

    % plot frequency controlled by plotFreq
    if plotting && length(errRaw)>plotLength && mod(sys.loopCount,plotFreq)==0
      
        % initialize figure
        if ~exist('figErr','var')
            figErr = figure(); % trend variable
            figPos = figure(); % true and reading position figure
            figKinematics = figure(); % spawn a new figure
            figNorm = figure(); % Kalman covariance            
            figAccl = figure();
        end
        
        % filter error
        figure(figErr);
        plot([1:plotLength+1]/sys.loopFreq,errRaw(end-plotLength:end));
        hold on;
        plot([1:plotLength+1]/sys.loopFreq,errMovingAvg(end-plotLength:end));
        plot([1:plotLength+1]/sys.loopFreq,errlowpass(end-plotLength:end));
        plot([1:plotLength+1]/sys.loopFreq,errKalman(end-plotLength:end));
		ylabel('Mean Squared Error');
        xlabel('Time (s)');
        title('Filter Error Comparison');
        legend('Raw Error', 'Moving Avg Error', 'Lowpass Error',['Kalman ' num2str(sys.loopFreq/sensorFilter.sensorInterval) 'Hz Error']);
        hold off;
        
        % position plot (single axis)
        figure(figPos);
        plot([1:plotLength+1]/sys.loopFreq,posReading(end-plotLength:end));
        hold on;
        plot([1:plotLength+1]/sys.loopFreq,posMovingAvg(end-plotLength:end));
        plot([1:plotLength+1]/sys.loopFreq,poslowpass(end-plotLength:end));
        plot([1:plotLength+1]/sys.loopFreq,posKalman(end-plotLength:end));
        plot([1:plotLength+1]/sys.loopFreq,posTrue(end-plotLength:end));
        ylabel('Position (m)');
        xlabel('Time (s)');
        title('Filter Estimated Position');
		legend('Sensor Reading','Moving Avg','Lowpass',['Kalman ' num2str(sys.loopFreq/sensorFilter.sensorInterval) 'Hz'],'True Pos');
        hold off;
        
        % Kalman covariance norm
        figure(figNorm);
        plot([1:plotLength+1]/sys.loopFreq,kalmanNorm(end-plotLength:end));
        hold on;
        ylabel('Norm');
        xlabel('Time (s)');
        title('Process Covariance Matrix Norm');
        legend('Kalman Covariance P');
        hold off;

        % Acceleration Reading
        figure(figAccl);
        plot([1:plotLength+1]/sys.loopFreq,accl(end-plotLength:end));
        hold on;
        ylabel('Accl (m/s^2)');
        xlabel('Time (s)');
        title('Acceleration Reading');
        legend('Acceleration Reading');
        hold off;      
        
        % draw kinematics      
        figure(figKinematics);
        vectarrow(sys.joing1Pos,sys.joing2Pos);
        vectarrow(sys.joing2Pos,sys.joing3Pos);
        vectarrow(sys.joing3Pos,sys.endEffector1);
        title('Kinematics');
        axis equal;
        drawnow;
        xlim([-16 16]);
        ylim([-16 16]);
        zlim([-5 16]);
        clf('reset');
        
    end
    
    if mod(sys.loopCount,1000)==0
        fprintf("\nTime: %.4f Time per scan: %.5f\n",sys.loopCount/sys.loopFreq,toc)
    end
end

    
