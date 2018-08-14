classdef Cylinder < handle
    %CYLINDER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        stroke = 135;  % mm
        bore = 105;    % mm
        crankRadius = stroke/2; %mm
        rodRatio = 1.6
        connectingRod = stroke * rodRatio; % mm
        compressionRatio = 16.5;
        
        theta = 0;  % deg
        volume
        ID
        N % engine speed
        
        Fsto = 14.01; % Stoichiometric ratio
        
        T = 200 + 273.15; % Intial temperature guess K
        F = 1; % Intial equivalence ratio guess
        P = 1; % Intial pressure guess bar
        m = 1; %Intial mass guess kg
        
        dV
        
        volumeTrace = [];
        temperatureTrace = [];
        pressureTrace = [];
        equivalenceTrace = [];
        thetaTrace = [];
        
        exhustValveOpens = 520;
        exhustValveCloses = 560;
        inletValveOpens = 710;
        inletValveCloses = 180;
        combustionStarts = 360;
        
    end
    
    methods
        function self = Cylinder(dTheta)
            
           if nargin == 1
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
            self.volume = (pi() * (self.bore / 2) ^ 2) * (H - L);
           
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
                dm = 0;
                
            elseif self.combustionStarts <= self.theta < ...
                    self.exhustValveOpens
                if self.theta == self.combustionStarts
                    self.ID = 2.4 * self.F^(-0.2) * self.P^(-1.02) * exp(2100/self.T);
                    
                end
                
                k1 = 2 + 1.25e-8 * (self.ID * self.N)^2.4;
                k2 = 5000;
                k3 = 14.2/(self.F ^ 0.644);
                k4 = 0.79 * k3 ^ 0.25;
                
                t = (self.theta - self.combustionStarts)/125;
                
                f1 = 1 - (1 - t^k1)^k2;
                f2 = 1 - exp(-k3 * t ^ k4);
                
                a = 0.85;
                b = 0.3;
                c = 0.4;
             
                beta = 1 - (a * self.F ^ b) / (self.ID ^ c);
                
                FB = beta * f1 + (1 - beta) * f2;
                
                dm = dmf;
                
            end
    end
    
end

