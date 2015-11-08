classdef FIRBandpassFilter < dsp.private.AbstractSampleRateEngine
    
    properties (Nontunable)
        CenterFrequency = 10e3;
        Bandwidth = 4e3;
        FilterOrder = 50;
        PassbandRipple = 0.1;
        StopbandAttenuation = 80;
    end
    properties (Logical,Nontunable)
        DesignForMinimumOrder = true;
    end
    properties (Access = protected)
        % pSampleRateDialog Default sample rate
        pSampleRateDialog = 44100;
    end
    properties (Access = protected)
       NumChannels = 1; 
       FilterObj
       IsFilterDesigned = false;
    end
    
    methods
        function obj = FIRBandpassFilter(varargin)
            setProperties(obj,nargin,varargin{:});
            obj.NumChannels = -1;
        end
        function set.CenterFrequency(obj,Fc)
            validateattributes(Fc,{'numeric'},{'finite','scalar','real',...
                'positive'},'');
            obj.CenterFrequency = Fc;
            needToDesignFilters(obj);
        end
        function set.Bandwidth(obj,bw)
            validateattributes(bw,{'numeric'},{'finite','scalar','real',...
                'positive'},'');
            obj.Bandwidth = bw;
            needToDesignFilters(obj);
        end
        function set.FilterOrder(obj,N)
            validateattributes(N,{'numeric'},{'finite','scalar','real',...
                'integer','positive'},'');
            obj.FilterOrder = N;
            needToDesignFilters(obj);
        end
        function set.PassbandRipple(obj,Apass)
            validateattributes(Apass,{'numeric'},{'finite','scalar','real',...
                'positive','<=',10},'');
            obj.PassbandRipple = Apass;
            needToDesignFilters(obj);
        end
        function set.StopbandAttenuation(obj,Astop)
            validateattributes(Astop,{'numeric'},{'finite','scalar','real',...
                'positive','<=',180},'');
            obj.StopbandAttenuation = Astop;
            needToDesignFilters(obj);
        end
        function F = getFilter(obj)
            if ~obj.IsFilterDesigned
                designFilter(obj);
            end
            F = clone(obj.FilterObj);
        end
        function fvtool(obj)
            FIR = getFilter(obj);
            fvtool(FIR,'Fs',obj.SampleRate);
        end
    end
    
    methods(Access = protected)
        function needToDesignFilters(obj)
            obj.IsFilterDesigned = false;
        end
        function designFilter(obj)
            Fs = obj.SampleRate;
            Ap = obj.PassbandRipple;
            Ast = obj.StopbandAttenuation;
            Fp = obj.Bandwidth/2;
            lpf = dsp.LowpassFilter('PassbandFrequency',Fp,...
                            'StopbandFrequency',Fp+1000,...
                            'PassbandRipple',Ap,...
                            'StopbandAttenuation',Ast,...
                            'SampleRate',Fs);
            FIR = getFilter(lpf);
            b = FIR.Numerator;
            wc = obj.CenterFrequency*2/Fs;
            L = size(b,2);
            b = 2*cos(pi*wc*(-(L-1)/2:(L-1)/2)).*b;
            FIR.Numerator = b;
            obj.FilterObj = FIR;
            
        end
        function setupImpl(obj,u)
            setupImpl@dsp.private.AbstractSampleRateEngine(obj,u);
            designFilter(obj);
            obj.NumChannels = getNumChannels(obj,u);
        end
        function y = stepImpl(obj,u)
            y = step(obj.FilterObj,u);
        end
        function resetImpl(obj)
            reset(obj.FilterObj);
        end
        function releaseImpl(obj)
            % Release the filter
            release(obj.FilterObj);
            obj.NumChannels = -1;
        end
        function validateInputsImpl(obj,u)
            % input cannot be unsigned complex fixed-point
            if isfi(u)
                coder.internal.errorIf(~issigned(u) && ~isreal(u), ...
                    'dsp:system:CICCompensationInterpolator:inputUnsignedComplexFixedPoint');
            end
            
            % Cache input data type for filter analysis
            if isempty(coder.target)
                if ~isempty(u)
                    cacheInputDataType(obj,u)
                end
            end
            
            if obj.NumChannels ~= -1
                coder.internal.errorIf(size(u,2) ~= obj.NumChannels, ['dsp:system',...
                    ':Shared:numChannels']);
            end
        end
        function s = saveObjectImpl(obj)
            % Default implementaion saves all public properties
            s = saveObjectImpl@dsp.private.AbstractSampleRateEngine(obj);
            % Private properties
            s.FilterObj = matlab.System.saveObject(obj.FilterObj);
            s.IsFilterDesigned = obj.IsFilterDesigned;
            s.pSR      = obj.pSR;
            if isLocked(obj)
                % All the following are set at setup
                s.NumChannels = obj.NumChannels;
            end
            s = saveFA(obj, s);   % Save filter analysis class properties
        end
        function s = loadObjectImpl(obj, s, wasLocked)
            if isfield(s, 'CoefficientsDataType')
                % backwards compatibility for objects saved before R2015b
                if strcmp(s.CoefficientsDataType.Signedness, 'Auto')
                    s.CoefficientsDataType.Signedness = 'Signed';
                end
            end
            if wasLocked
                % All the following were set at setup
                obj.NumChannels = s.NumChannels;
            end
            % Private properties
            obj.FilterObj = matlab.System.loadObject(s.FilterObj);
            
            obj.IsFilterDesigned = s.IsFilterDesigned;
            
            % Call base class method
            loadObjectImpl@dsp.private.AbstractSampleRateEngine(obj, s, wasLocked);
            
            if isfield(s, 'pSR')
                obj.pSR = s.pSR;
            else
                obj.pSR =  obj.SampleRate;
            end
            
            % Load filter analysis class properties
            loadFA(obj, s, wasLocked);
        end
        function flag = isInactivePropertyImpl(obj, prop)
            flag = isInactivePropertyImpl@dsp.private.AbstractSampleRateEngine(obj, prop);
            minord = obj.DesignForMinimumOrder;
            switch prop
                case {'FilterOrder'}
                    if minord
                        flag = true;
                    end
            end
        end
        
    end
    
end