classdef Cam < handle
    %CAM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        exhaust_open
        exhaust_close
        exhaust_lift_profile
        exhaust_max_lift
        
        inlet_open
        inlet_close
        inlet_lift_profile
        inlet_max_lift
        
    end
    
    methods
        function self = Cam(inlet_filename, exhaust_filename)
            
            exhaust_data = csvread(exhaust_filename, 1, 0);

            self.exhaust_open = exhaust_data(1, 3);
            self.exhaust_close = exhaust_data(1, 4);
            self.exhaust_max_lift = exhaust_data(1, 5);
            self.exhaust_lift_profile = exhaust_data(:, 1:2);


            inlet_data = csvread(inlet_filename, 1, 0);

            self.inlet_open = inlet_data(1, 3);
            self.inlet_close = inlet_data(1, 4);
            self.inlet_max_lift = inlet_data(1, 5);
            self.inlet_lift_profile = inlet_data(:, 1:2);
        end
        
        function lift = inletLift(self, theta)
            
            if mod(theta,2) == 0
            
                lift = self.inlet_lift_profile( (theta/2), 2);
                
            else
                                         
                l(1) = self.inlet_lift_profile( (theta/2 - 0.5), 2);
                l(2) = self.inlet_lift_profile( (theta/2 + 0.5), 2);
                
                lift = (l(1) + l(2)) / 2;
                
            end
            
            lift = lift * self.inlet_max_lift;
        end
        
        function lift = exhaustLift(self, theta)
            
            if mod(theta,2) == 0
            
                lift = self.exhaust_lift_profile( (theta/2), 2);
                
            else
                                         
                l(1) = self.exhaust_lift_profile( (theta/2 - 0.5), 2);
                l(2) = self.exhaust_lift_profile( (theta/2 + 0.5), 2);
                
                lift = (l(1) + l(2)) / 2;
                
            end
            
            lift = lift * self.exhaust_max_lift;
        end
        
    end
end
