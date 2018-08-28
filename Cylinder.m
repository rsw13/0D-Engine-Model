classdef Cylinder < handle
    %CYLINDER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
                
        compressionRatio = 16.5;
        stroke = 135e-3;                             % (m)
        bore = 105e-3;                               % (m)
        clearance = stroke / (compressionRatio - 1); % (m)
        Vswept = stroke * pi * (bore/2)^2;           % (m^3)
        crankRadius = stroke/2;                      % (m)
        rodRatio = 1.6
        connectingRod = stroke * rodRatio;           % (m)
        
        theta = 0;        % crank angle (deg)
        stepSize          % crank angle step size (deg)
        volume            % Volume (m^3)
        ID                % Ignition delay
        N                 % Engine speed (rpm)
        
        Fsto = 14.01;     % Stoichiometric ratio
        hf = -212.0423    % Enthalpy of formation of n-dectane (KJ/Kg)
        FB = zeros(1,2);  % fraction of fuel burnt at (t-1) and (t)
        
        T = 200 + 273.15; % Intial temperature guess K
        F = 1;            % Intial equivalence ratio guess
        P = 1;            % Intial pressure guess bar
        m = 1;            % Intial mass guess kg
        
        dV                % Change in volume wrt crank angle
        dm                % Change in mass wrt crank angle
        dF                % Change in equivalence ratio wrt crank angle
        dT                % Change in temp. wrt crank angle
        
        massTrace = [];
        volumeTrace = [];
        temperatureTrace = [];
        pressureTrace = [];
        equivalenceTrace = [];
        thetaTrace = [];
        
        exhustValveOpens = 520;  % Exhaust
        exhustValveCloses = 560;
        inletValveOpens = 710;
        inletValveCloses = 180;
        combustionStarts = 360;
        
        % internal energy fitting constants
        k1 = [0.692, 39.17e-6, 52.9e-9, -228.62e-13, 227.58e-17];
        k2 = [3049.39, -5.7e-2, -9.5e-5, 21.53e-9, -200.26e-14];
        
        heatExchangeConstants = [0.13, 0, 0];
        
        Tref
        Pref
        Vref
        Tsf = 100 + 273.15     % Cylinder surface temp (K)
        
        % Connetions
        exhuast
        intake
        
    end
    
    methods
        function self = Cylinder(stepSize, dTheta)
            
           if nargin == 1
               
               self.stepSize = stepSize;
               
           elseif nargin == 2
               
               self.stepSize = stepSize;
               self.theta = dTheta;
               self.newVolume;
                              
           end
           
           
           
        end
        
        function newVolume(self)
            % to find the current volume of the cylinder 
            
            H = self.connectingRod + self.crankRadius;
            
            L = self.crankRadius * cos(self.theta) + ...
                sqrt(self.connectingRod^2 + (self.crankRadius * ...
                sin(self.theta)^2));
            
            self.volume = (pi() * (self.bore / 2) ^ 2) * ...
                (self. clearance + H - L);
           
        end
        
        function A = newSurfaceArea(self)
            % to find the current volume of the cylinder 
            
            H = self.connectingRod + self.crankRadius;
            
            L = self.crankRadius * cos(self.theta) + ...
                sqrt(self.connectingRod^2 + (self.crankRadius * ...
                sin(self.theta)^2));
            
            A = pi * self.bore * (self. clearance + H - L) + pi * ...
                (self.bore / 2)^2 * 2;
           
        end
        
        function changeInVol(self)
            % to find the rate of change of the volume of the cylinder wrt
            % crank angle
            dL = self.crankRadius * sin(self.theta) + ...
                (((self.connectingRod^2 - (self.crankRadius * ...
                sin(self.theta))^2)^(-1/2)) * self.crankRadius^2 * ...
                sin(self.theta) * cos(self.theta)); 
            
            self.dV = pi() * (self.bore / 2) ^ 2 * dL;
            
        end
        
        function changeInMass(self)
            if self.inletValveCloses <= self.theta < self.combustionStarts
                % valves are shut and no fuel is added, so no mass is added
                self.dm = 0;
                
                self.massTrace(end+1,:) = self.dm;
                
            elseif self.combustionStarts <= self.theta < ...
                    self.exhustValveOpens
                % valves are shut but fuel is added, so change in mass is 
                % equal to the rate of fuel burnt
                if self.theta == self.combustionStarts
                    
                    self.ID = 2.4 * self.F^(-0.2) * self.P^(-1.02) *...
                        exp(2100/self.T);
                    
                    self.FB(1) = 0;
                    
                end
                
                a1 = 2 + 1.25e-8 * (self.ID * self.N)^2.4;
                a2 = 5000;
                a3 = 14.2/(self.F ^ 0.644);
                a4 = 0.79 * a3 ^ 0.25;
                
                t = (self.theta - self.combustionStarts)/125;
                
                f1 = 1 - (1 - t^a1)^a2;
                f2 = 1 - exp(-a3 * t ^ a4);
                
                a = 0.85;
                b = 0.3;
                c = 0.4;
             
                beta = 1 - (a * self.F ^ b) / (self.ID ^ c);
                
                self.FB(2) = beta * f1 + (1 - beta) * f2;
                
                dmf = (self.FB(2) - self.FB(1))/self.stepSize;
                
                self.dm = dmf;
                
                self.massTrace(end+1,:) = self.dm;
                
            end
        end
        
        function changeInF(self)
           if self.inletValveCloses <= self.theta < self.combustionStarts
               % Compression
               
               self.dF = 0;
               self.equivalenceTrace(end+1,:) = self.dF;
               
               
           elseif self.combustionStarts <= self.theta < ...
                    self.exhustValveOpens
               % Combustion
               
               F1 = 1 + self.F * self.Fsto;
               self.dF = (F1 * self.dm) / (self.m * self.Fsto);
               
               self.equivalenceTrace = self.dF;
               
           elseif (self.exhustValveOpens <= self.theta < ...
                   self.exhustValveCloses) && (self.theta < ...
                   self.inletValveOpens)
                % Exhaust peroid 
                
                self.dF = 0;
                
           elseif self.inletValveOpens <= self.theta < ...
                    self.exhustValveCloses
                % Valve overlap
                
           elseif self.inletValveOpens <= self.theta < ...
                    (720 + self.inletValveCloses)
                % Inlet peroid    
           end
        end
           
        function changeInT(self)
            
           R = 0.287 + 0.02 * self.F;
               
           Tpowers = [1, self.T, self.T^2, self.T^3, self.T^4, self.T^5];
           
           u = self.k1 .* Tpowers(2:6) - self.k2 .* Tpowers(1:5) * self.F;
                      
           du_dT = self.k1 .* Tpowers(1,5) .* [1, 2, 3, 4, 5] - ...
                self.k2 .* Tpowers(1,4) .* [1, 2, 3, 4] * self.F;
            
           du_dF = - self.k2 .* Tpowers(1:5);
           
           Qloss = self.cylinderHeatLoss;
           
           if self.inletValveCloses <= self.theta < self.combustionStarts
                % Compression
                % NEED TO UPDATE DV AND VOL. AT START OF CYCLE
                
                self.dT = (Qloss * self.m - ((R * self.T * self.dV)/...
                    self.volume))/du_dT;
                
           elseif self.combustionStarts <= self.theta < ...
                    self.exhustValveOpens
                % Combustion
                
                self.dT = ((Qloss + self.dm * self.hf - u * self.dm) ...
                    * self.m - ((R * self.T * self.dV) / self.volume) ...
                    - du_dF * self.dF) / du_dT;
                
           elseif (self.exhustValveOpens <= self.theta < ...
                   self.exhustValveCloses) && (self.theta < ...
                   self.inletValveOpens)
                % Exhaust peroid 
                
           elseif self.inletValveOpens <= self.theta < ...
                    self.exhustValveCloses
                % Valve overlap
                
           elseif self.inletValveOpens <= self.theta < ...
                    (720 + self.inletValveCloses)
                % Inlet peroid
           
                
           end
        end
        
        function Qloss = cylinderHeatLoss(self)
            
           Cpis = 2 * self.N * self.stroke;
           cylinderArea = self.newSurfaceArea; 
           
           if self.inletValveCloses <= self.theta < self.combustionStarts
                % Compression
                
                self.heatExchangeConstants(2) = 2.28;
                self.heatExchangeConstants(2) = 0;
              
                              
           elseif self.combustionStarts <= self.theta < ...
                    self.exhustValveOpens
                % Combustion
                
                if self.theta == self.combustionStarts
                    
                    self.Tref = self.T;
                    self.Pref = self.P;
                    self.Vref = self.volume;
                end
               
                self.heatExchangeConstants(2) = 2.28;
                self.heatExchangeConstants(3) = 3.24e-3;
                
                Pmot = self.Pref * (self.Vref / self.volume)^(1.32);
              
           elseif (self.exhustValveOpens <= self.theta < ...
                   self.exhustValveCloses) && (self.theta < ...
                   self.inletValveOpens)
                % Exhaust peroid

                self.heatExchangeConstants(2) = 6.18;
                self.heatExchangeConstants(3) = 0;
               
           elseif self.inletValveOpens <= self.theta < ...
                    self.exhustValveCloses
                % Valve overlap

                self.heatExchangeConstants(2) = 6.18;
                self.heatExchangeConstants(3) = 0;
                
           elseif self.inletValveOpens <= self.theta < ...
                    (720 + self.inletValveCloses)
                % Inlet peroid

                self.heatExchangeConstants(2) = 6.18;
                self.heatExchangeConstants(3) = 0;
                 
               
           end
            
           htConstant = (self.heatExchangeConstants(1) * self.P ^ 0.8) /...
                   (self.bore^0.2 * self.T^0.53);
           
           
           if self.heatExchangeConstants(3) == 0
           
               ht = htConstant * (self.heatExchangeConstants(2) * Cpis)^0.8;
           else
           
               ht = htConstant * (self.heatExchangeConstants(2) * Cpis +...
                   (self.heatExchangeConstants(3) * self.Vswept * ...
                   self.Tref * (self.P - Pmot)) / (self.Pref * ...
                   self.Vref))^0.8;
           end
           
           Qloss = cylinderArea * ht * (self.T - self.Tsf);
           
        end
        
        function linkExhaust(self, exhaust)
            
            self.exhaust = exhaust;
            
        end
        
        function linkIntake(self, intake)
            
            self.intake = intake;
            
        end
        
        function linkCam(self, cam)
            
            self.cam = cam;
            
        end
        
        
    end
end

