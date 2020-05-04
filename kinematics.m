classdef kinematics
    properties (Access = public)
        joint1;
        joint2;
        joint3;
        endEffectorRaw;
        joing2PosV;
        joing3PosV;
        endEffector1V;
        PivotAxis;
        yawAxis;
        joint1Angle;    % angle of change in each iteration
        joint2Angle;    % angle of change in each iteration
        joint3Angle;    % angle of change in each iteration
        jointYawAngle;  % angle of slew, different axis applied here
        joing1Pos;      % pivot point coordinate
        joing2Pos;      % pivot point coordinate
        joing3Pos;      % pivot point coordinate
        endEffector1;
        naturalFreqX;
        naturalFreqZ;   % natural frequency on the slew axis
        oscillationMagnitudeX;
        oscillationMagnitudeZ;
        loopCount;
        loopFreq;
        figureKinematics;
        noiseCount;     % counter on the noise array
        impulse;
        impulseAxis;
    end
    methods
        function obj=initKinematics(obj)
            
            obj.loopCount = 1;
            obj.noiseCount = 1;
            
            obj.impulse = impulseGen(obj.loopFreq);
            
            % get raw vector from raw measurement
            obj.joing2PosV = obj.joint2-obj.joint1;
            obj.joing3PosV = obj.joint3-obj.joint2;
            obj.endEffector1V = obj.endEffectorRaw-obj.joint3;
                        
            % initial boom, stick, hear joint angle
            obj.joint1Angle = 30;
            obj.joint2Angle = 50;
            obj.joint3Angle = 12;
            obj.jointYawAngle = 0;
            
            % rotate vector to get rotated vector for each component
            obj.joing2PosV = rotVecAroundArbAxis(obj.joing2PosV,obj.PivotAxis,obj.joint1Angle);
            obj.joing3PosV = rotVecAroundArbAxis(obj.joing3PosV,obj.PivotAxis,obj.joint1Angle+obj.joint2Angle);
            obj.endEffector1V = rotVecAroundArbAxis(obj.endEffector1V,obj.PivotAxis,obj.joint1Angle+obj.joint2Angle+obj.joint3Angle);

            % assemble all vector back together
            obj.joing1Pos = obj.joint1;
            obj.joing2Pos = obj.joing1Pos+obj.joing2PosV;
            obj.joing3Pos = obj.joing2Pos+obj.joing3PosV;
            obj.endEffector1 = obj.joing3Pos+obj.endEffector1V;

        end
        function obj = calcEndEffector(obj)
                        
            % add natural frequency
            obj.joint1Angle = obj.oscillationMagnitudeX * sin(obj.loopCount/obj.loopFreq/obj.naturalFreqX);
            obj.joint2Angle = obj.oscillationMagnitudeX * sin(obj.loopCount/obj.loopFreq/obj.naturalFreqX);
            obj.joint3Angle = obj.oscillationMagnitudeX * sin(obj.loopCount/obj.loopFreq/obj.naturalFreqX);
            obj.jointYawAngle = obj.oscillationMagnitudeZ * sin(obj.loopCount/obj.loopFreq/obj.naturalFreqZ); 

            % rotate vector to get rotated vector for each component
            obj.joing2PosV = rotVecAroundArbAxis(obj.joing2PosV,obj.PivotAxis,obj.joint1Angle); % stick pose
            obj.joing3PosV  = rotVecAroundArbAxis(obj.joing3PosV, obj.PivotAxis,obj.joint1Angle+obj.joint2Angle); % head pose
            obj.endEffector1V       = rotVecAroundArbAxis(obj.endEffector1V,      obj.PivotAxis,obj.joint1Angle+obj.joint2Angle+obj.joint3Angle);

            % rotate boom slew for all vectors
            obj.joing2PosV = rotVecAroundArbAxis(obj.joing2PosV,obj.yawAxis,obj.jointYawAngle); % stick slew rotation
            obj.joing3PosV  = rotVecAroundArbAxis(obj.joing3PosV, obj.yawAxis,obj.jointYawAngle);
            obj.endEffector1V       = rotVecAroundArbAxis(obj.endEffector1V,      obj.yawAxis,obj.jointYawAngle);

            % assemble all vector back together
            obj.joing2Pos = obj.joing1Pos+obj.joing2PosV;
            obj.joing3Pos = obj.joing2Pos+obj.joing3PosV;
            obj.endEffector1 = obj.joing3Pos+obj.endEffector1V;

            % add random disturbance
            obj.noiseCount = obj.noiseCount + 1;
            sz = size(obj.impulse);
            if obj.noiseCount > sz(1)
                obj.impulse = impulseGen(obj.loopFreq);
                obj.impulseAxis = randi([1 3]);
                obj.noiseCount = 1;
            end
            
            obj.endEffector1(obj.impulseAxis) = obj.endEffector1(obj.impulseAxis) + obj.impulse(obj.noiseCount);
            
        end
    end
end