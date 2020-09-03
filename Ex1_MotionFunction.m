classdef Ex1_MotionFunction
    properties
        miu_;
        m_;
        k_;
        u_;
    end
    methods
        function obj = Ex1_MotionFunction(m, k, miu, u)
            obj.m_ = m;
            obj.k_ = k;
            obj.miu_ = miu;
            obj.u_ = u;
        end
        function x_dot = f(obj, x)
           x_dot = [x(2); (obj.u_ - obj.k_*x(1) - obj.miu_*x(2))./obj.m_];
        end
        function x_dot = f_lossless(~, x)
            x_dot = [x(2); -x(1)];
        end
        function x_dot = f_damping(~, x, miu)
            x_dot = [x(2); -x(1)-x(2)*miu];
        end
    end
end

