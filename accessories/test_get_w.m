classdef test_get_w < matlab.unittest.TestCase
    
    
    methods(Test)
        % Test methods
        
        function test_y(testCase)
            
            originalArray = ones(2,1,2);

            originalArray(1,:,:) = 6;
            originalArray(2,:,:) = 7;
            
            w = ones(4,1,4);
            w(1:2,:,:) = 6;
            w(3:4,:,:) = 7;


            expandedArray = get_w1(originalArray, 2);
            testCase.verifyEqual(w, expandedArray);
        end

        function test_one_w1(testCase)
            
            originalArray = ones(2,1,2);

            originalArray(1,:,:) = 6;
            originalArray(2,:,:) = 7;
            


            w = get_w1(originalArray, 1);
            testCase.verifyEqual(originalArray, w);
        end

        function test_one_w2(testCase)
            
            originalArray = ones(2,1,2);

            originalArray(1,:,:) = 6;
            originalArray(2,:,:) = 7;
            


            w = get_w2(originalArray, 1);
            testCase.verifyEqual(originalArray, w);
        end


    end
    
end