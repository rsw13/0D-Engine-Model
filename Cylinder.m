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
        R 
        gamma             % specific heat ratio
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
        
        combustionStarts = 360; % Combustion start (deg) 
        
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
        cam
        
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
           
           self.updateConstants()
           
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
        
        function updateConstants(self)
           
            self.gamma = 1.4 - 0.16 * self.F;
            self.R = 0.287 + 0.02 * self.F;
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
            if self.cam.inlet_close <= self.theta < self.combustionStarts
                % valves are shut and no fuel is added, so no mass is added
                self.dm = 0;
                
                self.massTrace(end+1,:) = self.dm;
                
            elseif self.combustionStarts <= self.theta < ...
                    self.cam.exhaust_open
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
                
            elseif self.cam.exhaust_open <= self.theta < ...
                    self.cam.inlet_open
                % Exhaust peroid 
                
                valve_diameter = self.cam.exhaust_valve_diameter;
                lift, Cd = self.cam.exhaustLiftCd(self.theta);
               
                
                if self.P > self.exhuast.P
                                        
                    P1 = self.P;
                    P2 = self.exhuast.P;
                    T1 = self.T;
                                        
                    self.dm = - (self.valveFlow(valve_diameter, lift, ...
                        P1, P2, T1, Cd, self.gamma, self.dm));
                    
                elseif self.P < self.exhuast.P
                    
                    P1 = self.exhuast.P;
                    P2 = self.P;
                    T1 = self.exhuast.T;
                    
                    self.dm = self.valveFlow(valve_diameter, lift, P1, P2,...
                        T1, Cd, self.gamma, self.dm);
                    
                else
                    
                    self.dm = 0;
                    
                    
                end
                
            end
        end
        
        function changeInF(self)
           if self.cam.inlet_close <= self.theta < self.combustionStarts
               % Compression
               
               self.dF = 0;
               self.equivalenceTrace(end+1,:) = self.dF;
               
               
           elseif self.combustionStarts <= self.theta < ...
                    self.cam.exhaust_open
               % Combustion
               
               F1 = 1 + self.F * self.Fsto;
               self.dF = (F1 * self.dm) / (self.m * self.Fsto);
               
               self.equivalenceTrace = self.dF;
               
           elseif (self.cam.exhaust_open <= self.theta < ...
                   self.cam.exhaust_close) && (self.theta < ...
                   self.cam.inlet_open)
                % Exhaust peroid 
                
                self.dF = 0;
                
           elseif self.cam.inlet_open <= self.theta < ...
                    self.cam.exhaust_close
                % Valve overlap
                
           elseif self.cam.inlet_open <= self.theta < ...
                    (720 + self.cam.inlet_close)
                % Inlet peroid    
           end
        end
           
        function changeInT(self)
                           
           Tpowers = [1, self.T, self.T^2, self.T^3, self.T^4, self.T^5];
           
           u = self.k1 .* Tpowers(2:6) - self.k2 .* Tpowers(1:5) * self.F;
                      
           du_dT = self.k1 .* Tpowers(1,5) .* [1, 2, 3, 4, 5] - ...
                self.k2 .* Tpowers(1,4) .* [1, 2, 3, 4] * self.F;
            
           du_dF = - self.k2 .* Tpowers(1:5);
           
           Qloss = self.cylinderHeatLoss;
           
           if self.cam.inlet_close <= self.theta < self.combustionStarts
                % Compression
                % NEED TO UPDATE DV AND VOL. AT START OF CYCLE
                
                self.dT = (Qloss * self.m - ((self.R * self.T * self.dV)/...
                    self.volume))/du_dT;
                
           elseif self.combustionStarts <= self.theta < ...
                    self.cam.exhaust_open
                % Combustion
                
                self.dT = ((Qloss + self.dm * self.hf - u * self.dm) ...
                    * self.m - ((self.R * self.T * self.dV) / self.volume) ...
                    - du_dF * self.dF) / du_dT;
                
           elseif (self.cam.exhaust_open <= self.theta < ...
                   self.cam.exhaust_close) && (self.theta < ...
                   self.cam.inlet_open)
                % Exhaust peroid 
                
           elseif self.cam.inlet_open <= self.theta < ...
                    self.cam.exhaust_close
                % Valve overlap
                
           elseif self.cam.inlet_open <= self.theta < ...
                    (720 + self.cam.inlet_close)
                % Inlet peroid
           
                
           end
        end
        
        function Qloss = cylinderHeatLoss(self)
            
           Cpis = 2 * self.N * self.stroke;
           cylinderArea = self.newSurfaceArea; 
           
           if self.cam.inlet_close <= self.theta < self.combustionStarts
                % Compression
                
                self.heatExchangeConstants(2) = 2.28;
                self.heatExchangeConstants(2) = 0;
              
                              
           elseif self.combustionStarts <= self.theta < ...
                    self.cam.exhaust_open
                % Combustion
                
                if self.theta == self.combustionStarts
                    
                    self.Tref = self.T;
                    self.Pref = self.P;
                    self.Vref = self.volume;
                end
               
                self.heatExchangeConstants(2) = 2.28;
                self.heatExchangeConstants(3) = 3.24e-3;
                
                Pmot = self.Pref * (self.Vref / self.volume)^(1.32);
              
           elseif (self.cam.exhaust_open <= self.theta < ...
                   self.cam.exhaust_close) && (self.theta < ...
                   self.cam.inlet_open)
                % Exhaust peroid

                self.heatExchangeConstants(2) = 6.18;
                self.heatExchangeConstants(3) = 0;
               
           elseif self.cam.inlet_open <= self.theta < ...
                    self.cam.exhaust_close
                % Valve overlap

                self.heatExchangeConstants(2) = 6.18;
                self.heatExchangeConstants(3) = 0;
                
           elseif self.cam.inlet_open <= self.theta < ...
                    (720 + self.cam.inlet_close)
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
        
        function dm = valveFlow(valve_diameter, lift, P1, P2, T1, Cd, ...
                gamma, dm)
           
            
            A = lift * 2 * pi * valve_diameter;  % valve flow area

            rho = P1 / (self.R * T1);

            u = dm / (rho * A * Cd);

            P0 = P1 * (1 + ((gamma - 1) / 2) * (u / sqrt(gamma * self.R * T1)) ^ 2)^...
            (gamma / (gamma -1  ));   % stagnation pressure upstream of the valve

            critical_pressure_ratio = (2 / (gamma + 1) ) ^ (-gamma / (gamma - 1)); % pressure ratio at which the flow will choke

            pressure_ratio = P0 / P2;  % pressure ratio between upstream of the valve and downstream

            if pressure_ratio >= critical_pressure_ratio

                dm = Cd * A * P1 * sqrt((gamma / (self.R * T1)) * (2 / (gamma +1)) ^...
                ((gamma +1) / (gamma -1)));

            else

                dm = Cd * A * P1 * sqrt(((2 * gamma) / (gamma -1)) * (1 / (self.R * T1) ) * ...
                    ((P2 / P1)^(2/gamma) - (P2 / P1)^((gamma + 1)/gamma)));

            end
        end
    end
end

