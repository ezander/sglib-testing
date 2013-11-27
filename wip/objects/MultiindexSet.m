classdef MultiindexSet
    properties
        m
        p
    end
    properties (SetAccess=protected)
        I
    end
    properties (Dependent)
        length
    end
    methods
        function obj=MultiindexSet(m,p)
            obj.m=m;
            obj.p=p;
            obj.I=multiindex(obj.m,obj.p);
        end
        function obj=set.m( obj, m )
            obj.m=m;
            obj.I=multiindex(obj.m,obj.p);
        end
        function obj=set.p( obj, p )
            obj.p=p;
            obj.I=multiindex(obj.m,obj.p);
        end
        function n=get.length( obj )
            n=size(obj.I,1);
        end
    end
end
