classdef test1 < handle
    %TEST1 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        a
        linkObj
    end
    
    methods
        function self = test1(a)
            %TEST1 Construct an instance of this class
            %   Detailed explanation goes here
            self.a = a;
        end
        
        function link(self, obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            self.linkObj = obj;
        end
        
        function a = add(self)
            
            a = self.linkObj.a + self.a;
        end
    end
end

