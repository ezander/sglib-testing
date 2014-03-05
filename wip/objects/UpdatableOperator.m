classdef UpdatableOperator < Operator
    methods(Abstract)
        newop = update(op, varargin);
    end
end
