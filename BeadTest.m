classdef BeadTest < matlab.unittest.TestCase
    

%     methods(TestClassSetup)
%         % Shared setup for the entire test class
%     end
%     
%     methods(TestMethodSetup)
%         % Setup for each test
% 
%     end


    
    methods(Test)
        % Test methods
        function testPositionforOriginalBead(testCase)
            bead = OriginalBead('position',[0 0],'speed',1);
            bead.velocity
            bead.velocity = [1 0];
            bead = bead.updatePosition(1);
            testCase.verifyEqual(bead.position, [1,0]);
        end


        function testPositionforStuckBead(testCase)
            bead = StuckBead();
            bead.updatePosition(0.01)
            testCase.verifyEqual(bead.position, [0,0]);
        end

    end
    
end