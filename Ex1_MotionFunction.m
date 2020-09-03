classdef Ex1_MotionFunction
    methods
        function x_dot = f_lossless(~, x)
            x_dot = [x(2); -x(1)];
        end
        function x_dot = f_damping(~, x, miu)
            x_dot = [x(2); -x(1)-x(2)*miu];
        end
    end
end

