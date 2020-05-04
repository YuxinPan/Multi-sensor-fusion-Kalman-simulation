classdef filtering
    properties (Access = public)
        endEffectorReading;
        loopCount;
        loopFreq;
        
        % moving average filter
        movingAvgNum; % moving average array length
        movingAvgArrX1;
        movingAvgArrY1;
        movingAvgArrZ1;
        endEffector1AvgFiltered;

        % lowpass filter
        lowpassCoeffX1; % [b0 b1 b2 a1 a2]
        lowpassLastX1;  % [Xn-1 Xn-2 Yn-1 Yn-2]
        lowpassCoeffY1; % [b0 b1 b2 a1 a2]
        lowpassLastY1;  % [Xn-1 Xn-2 Yn-1 Yn-2]
        lowpassCoeffZ1; % [b0 b1 b2 a1 a2]
        lowpassLastZ1;  % [Xn-1 Xn-2 Yn-1 Yn-2]
        endEffector1lowpassFiltered;
        
        % Kalman filter
        H_IMU;          % IMU observation model
        H_sensor;       % tracking sensor observation model
        F;              % state transition matrix
        P;              % state covariance matrix
        Q;              % process motion noise
        R_IMU;          % IMU noise covariance
        R_sensor;       % tracking sensor noise covariance
        S;              
        W;              % Kalman gain
        inn;            % innovation matrix
        x_est;
        z_est;
        z;
        accl;
        sensorInterval; % between how many scans one sensor reading is available
        endEffector1kalmanFiltered;
        
    end
    methods
        function obj=initFilter(obj) % initialize filter
          
           obj.loopCount = 0;
          
           %%%% moving average init %%%%
           obj.movingAvgArrX1 = zeros(1,obj.movingAvgNum);
           obj.movingAvgArrY1 = zeros(1,obj.movingAvgNum);
           obj.movingAvgArrZ1 = zeros(1,obj.movingAvgNum);
           
           %%%% lowpass init %%%%
           obj.lowpassLastX1 = [0 0]; % last and last'last
           obj.lowpassLastY1 = [0 0];
           obj.lowpassLastZ1 = [0 0];
           
           %%%% Kalman init %%%%
           dt = 1/obj.loopFreq;
           
           % state transition matrix, kinematics model 
           obj.F = eye(9);
           
           transition = [1  dt 1/2*dt^2;
                         0  1  dt;
                         0  0  1];
           obj.F(1:3,1:3) = transition;
           obj.F(4:6,4:6) = transition;
           obj.F(7:9,7:9) = transition;
                
           % observation model 
           obj.H_sensor = zeros(3,9);
           obj.H_sensor(1,1) = 1;
           obj.H_sensor(2,4) = 1;
           obj.H_sensor(3,7) = 1;

           obj.H_IMU = zeros(3,9);
           obj.H_IMU(1,3) = 1;
           obj.H_IMU(2,6) = 1;
           obj.H_IMU(3,9) = 1;
           
           obj.P = zeros(9);
           
           % process motion noise
           obj.Q = 0.008*2 * eye(9);
           obj.Q(2,2) = 0.02*2;
           obj.Q(5,5) = 0.02*2;
           obj.Q(8,8) = 0.02*2;
           obj.Q(3,3) = 0.2*2;
           obj.Q(6,6) = 0.2*2;
           obj.Q(9,9) = 0.2*2;
           
           % sensor noise covariance
           obj.R_sensor = 8^2 * eye(3);
           obj.R_IMU = 0.001^2 * eye(3);
           
           % state initial
           obj.x_est = zeros(9,1);
           
        end
        function obj=kalmanPredict(obj) % Kalman predict

            obj.x_est = obj.F*obj.x_est; % F*x(k-1|k-1)

            obj.P = obj.F*obj.P*obj.F' + obj.Q; % F*P*F' + Q

            obj.endEffector1kalmanFiltered = [obj.x_est(1) obj.x_est(4) obj.x_est(7)]; % save positional data

        end
        function obj=kalmanUpdateIMU(obj) % Kalman update IMU
            
            obj.z = [obj.accl(1,1); % sensor reading % 3*1
                     obj.accl(1,2);
                     obj.accl(1,3)];
                 
            obj.z_est = obj.H_IMU*obj.x_est; % H*x(k|k-1)
            
            obj.inn = obj.z - obj.z_est; % innovation matrix
            
            obj.S = obj.H_IMU*obj.P*obj.H_IMU' + obj.R_IMU; % H*P*H' + R
            
            obj.W = obj.P*obj.H_IMU'*inv(obj.S); %P*H'*S^-1 % Kalman gain

            obj.x_est = obj.x_est + obj.W*obj.inn; % x(k|k-1) + W*inn
            
            obj.P = obj.P - obj.W*obj.S*obj.W'; % P-W*S*W'

            obj.endEffector1kalmanFiltered = [obj.x_est(1) obj.x_est(4) obj.x_est(7)]; % save positional data

        end  
        function obj=kalmanUpdatesensor(obj) % Kalman update tracking sensor
            
            obj.z = [obj.endEffectorReading(1); % sensor reading % 3*1
                     obj.endEffectorReading(2);
                     obj.endEffectorReading(3);];
                 
            obj.z_est = obj.H_sensor*obj.x_est; % H*x(k|k-1)
            
            obj.inn = obj.z - obj.z_est; % innovation matrix
            
            obj.S = obj.H_sensor*obj.P*obj.H_sensor' + obj.R_sensor; % H*P*H' + R
            
            obj.W = obj.P*obj.H_sensor'*inv(obj.S); %P*H'*S^-1 % Kalman gain

            obj.x_est = obj.x_est + obj.W*obj.inn; % x(k|k-1) + W*inn
            
            obj.P = obj.P - obj.W*obj.S*obj.W'; % P-W*S*W'

            obj.endEffector1kalmanFiltered = [obj.x_est(1) obj.x_est(4) obj.x_est(7)]; % save positional data

        end          
        function obj=lowpassFilter(obj) % lowpass filter
            obj.endEffector1lowpassFiltered(1) = obj.endEffectorReading(1)*obj.lowpassCoeffX1(1) + (1-obj.lowpassCoeffX1(1))*obj.lowpassLastX1(1) + obj.lowpassCoeffX1(2)*(obj.lowpassLastX1(1)-obj.lowpassLastX1(2));
            obj.lowpassLastX1(2) = obj.lowpassLastX1(1);
            obj.lowpassLastX1(1) = obj.endEffector1lowpassFiltered(1);

            obj.endEffector1lowpassFiltered(2) = obj.endEffectorReading(2)*obj.lowpassCoeffY1(1) + (1-obj.lowpassCoeffY1(1))*obj.lowpassLastY1(1) + obj.lowpassCoeffY1(2)*(obj.lowpassLastY1(1)-obj.lowpassLastY1(2));
            obj.lowpassLastY1(2) = obj.lowpassLastY1(1);
            obj.lowpassLastY1(1) = obj.endEffector1lowpassFiltered(2);

            obj.endEffector1lowpassFiltered(3) = obj.endEffectorReading(3)*obj.lowpassCoeffZ1(1) + (1-obj.lowpassCoeffZ1(1))*obj.lowpassLastZ1(1) + obj.lowpassCoeffZ1(2)*(obj.lowpassLastZ1(1)-obj.lowpassLastZ1(2));
            obj.lowpassLastZ1(2) = obj.lowpassLastZ1(1);
            obj.lowpassLastZ1(1) = obj.endEffector1lowpassFiltered(3);

        end
        function obj=movingAvg(obj) % moving average filter 
            obj.movingAvgArrX1(mod(obj.loopCount,obj.movingAvgNum)+1) = obj.endEffectorReading(1);
            obj.movingAvgArrY1(mod(obj.loopCount,obj.movingAvgNum)+1) = obj.endEffectorReading(2);
            obj.movingAvgArrZ1(mod(obj.loopCount,obj.movingAvgNum)+1) = obj.endEffectorReading(3);

            obj.endEffector1AvgFiltered = [mean(obj.movingAvgArrX1) mean(obj.movingAvgArrY1) mean(obj.movingAvgArrZ1)];

        end      
    end
end