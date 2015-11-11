classdef BandpassFilter < dsp.private.LPHPFilterBase
    %BandpassFilter FIR or IIR lowpass filter
    %   LPF = dsp.BandpassFilter returns a lowpass filter, LPF, which
    %   independently filters each channel of the input over time using the
    %   given design specifications.
    %
    %   LPF = dsp.BandpassFilter('PropertyName', PropertyValue, ...) returns
    %   a lowpass filter, LPF, with each specified property set to the
    %   specified value.
    %
    %   Step method syntax:
    %
    %   Y = step(LPF, X) filters the real or complex input signal X using
    %   the specified filter to produce the output Y. Each column of the
    %   input signal is filtered independently over time.
    %
    %   BandpassFilter methods:
    %
    %   step      - See above description for use of this method
    %   getFilter - Return underlying FIR or Biquad filter
    %   release   - Allow property value and input characteristics changes
    %   clone     - Create lowpass filter object with same property values
    %   isLocked  - Locked status (logical)
    %   reset     - Reset the internal states to initial conditions
    %
    %   BandpassFilter properties:
    %
    %   SampleRate             - Sample rate of input in Hz
    %   FilterType             - FIR or IIR (Biquad)
    %   DesignForMinimumOrder  - Option for minimum/specify order design
    %   FilterOrder            - Order of the filter
    %   CenterFrequency      - Passband frequency of lowpass filter
    %   Bandwidth      - Stopband frequency of lowpass filter
    %   PassbandRipple         - Allowed passband ripple of lowpass filter
    %   StopbandAttenuation    - Minimum stopband attenuation of the filter
    %
    %   This System object supports fixed-point operations. For more
    %   information, type dsp.BandpassFilter.helpFixedPoint.
    %
    %   % EXAMPLE: Design an equiripple lowpass FIR filter. Filter white
    %   % Gaussian noise with the resulting filter and estimate the
    %   % transfer function.
    %   % Initialize
    %   Fs = 96e3;
    %   LPF = dsp.BandpassFilter('SampleRate',Fs,'CenterFrequency',8e3);
    %   TFE = dsp.TransferFunctionEstimator('FrequencyRange','onesided',...
    %           'SpectralAverages',50);
    %   AP = dsp.ArrayPlot('PlotType','Line','YLimits',[-100 5],...
    %           'YLabel','Magnitude (dB)','XLabel','Frequency (Hz)',...
    %           'SampleIncrement',Fs/1024,'ShowLegend',true,...
    %           'Title',...
    %           'Magnitude Response; Exact (Ch 1), Estimate (Ch 2)');
    %   Htrue = freqz(LPF,513); % Compute exact transfer function
    %   % Stream processing loop
    %   Niter = 1000;
    %   for k = 1:Niter
    %       x = randn(1024,1);  % Input signal = white Gaussian noise
    %       y = step(LPF,x);   % Filter noise with FIR filter
    %       H = step(TFE,x,y); % Compute transfer function estimate
    %       step(AP,20*log10(abs([Htrue,H]))); % Plot estimate
    %   end
    %
    %   See also dsp.HighpassFilter, dsp.VariableBandwidthFIRFilter,
    %   dsp.VariableBandwidthIIRFilter, dsp.FIRFilter, dsp.BiquadFilter.
    
    % Copyright 2014-15 The MathWorks, Inc.
    
    %#codegen
    %#ok<*EMCA>
    
    properties (Nontunable)
        %CenterFrequency  Passband edge frequency (Hz)
        %   Specify the passband edge frequency of the filter in Hertz, as
        %   a positive real scalar smaller than half the sample rate. The
        %   default value of this property is 8 kHz.
        CenterFrequency = 10e3;
        %Bandwidth  Stopband edge frequency (Hz)
        %   Specify the stopband edge frequency of the filter in Hertz, as
        %   a positive real scalar smaller than half the sample rate. This
        %   property applies when you set the DesignForMinimumOrder
        %   property to true. The default value of this property is 12 kHz.
        Bandwidth = 4e3;
        
        StopbandFrequency1 = 3e3;
        PassbandFrequency1 = 5e3;
        PassbandFrequency2 = 8e3;
        StopbandFrequency2 = 12e3;
        StopbandAttenuation1 = 80;
        StopbandAttenuation2 = 60;
        
    end
    
    methods
        function obj = BandpassFilter(varargin)
            
            % Support name-value pair arguments
            setProperties(obj,nargin,varargin{:});
            obj.NumChannels = -1;
        end
        
        function set.CenterFrequency(obj,Fc)
            validateattributes(Fc,{'numeric'},...
            {'finite','scalar','real',...
                'positive'},'');
            obj.CenterFrequency = Fc;
            needToDesignFilters(obj);
        end
        
        function set.Bandwidth(obj,Fst)
            validateattributes(Fst,{'numeric'},{'finite','scalar','real',...
                'positive'},'');
            obj.Bandwidth = Fst;
            needToDesignFilters(obj);
        end
        
        function set.StopbandFrequency1(obj,Fst1)
            validateattributes(Fst1,{'numeric'},...
            {'finite','scalar','real',...
                'positive'},'');
            obj.StopbandFrequency1 = Fst1;
            needToDesignFilters(obj);
        end
        
        function set.PassbandFrequency1(obj,Fp1)
            validateattributes(Fp1,{'numeric'},...
            {'finite','scalar','real',...
                'positive'},'');
            obj.PassbandFrequency1 = Fp1;
            needToDesignFilters(obj);
        end
        
        function set.PassbandFrequency2(obj,Fp2)
            validateattributes(Fp2,{'numeric'},...
            {'finite','scalar','real',...
                'positive'},'');
            obj.PassbandFrequency2 = Fp2;
            needToDesignFilters(obj);
        end
        
        function set.StopbandFrequency2(obj,Fst2)
            validateattributes(Fst2,{'numeric'},...
            {'finite','scalar','real',...
                'positive'},'');
            obj.StopbandFrequency2 = Fst2;
            needToDesignFilters(obj);
        end
        
    end
    %----------------------------------------------------------------------
    methods (Access = protected)
        
        function b = computeFIRCoefficients(obj)
            Fs = getSampleRate(obj);
            
            Ap = obj.PassbandRipple;
            Rp = (10.^(Ap./20) - 1)./...
                (10.^(Ap./20) + 1);
            
            if obj.DesignForMinimumOrder
                Fst1 = obj.StopbandFrequency1./(Fs/2);
                Fp1 = obj.PassbandFrequency1./(Fs/2);
                Fp2 = obj.PassbandFrequency2./(Fs/2);
                Fst2 = obj.StopbandFrequency2./(Fs/2);
                % Adjust frequencies to maintain same transition width both
                % sides
                TW1 = Fp1 - Fst1;
                TW2 = Fst2 - Fp2;
                designTW = min(TW1,TW2);
                if TW1 ~= designTW
                    Fst1 = Fp1-designTW;
                elseif TW2 ~= designTW
                    Fst2 = designTW+Fp2;
                end
                Ast1 = obj.StopbandAttenuation1;
                Rs1 = 10.^(-Ast1./20);
                Ast2 = obj.StopbandAttenuation2;
                Rs2 = 10.^(-Ast2./20);
                b = firgr('minorder',[0 Fst1 Fp1 Fp2 Fst2 1],[0 0 1 1 0 0],...
                    [Rs1,Rp,Rs2]);
            else
                Ast = obj.StopbandAttenuation;
                Rs = 10.^(-Ast./20);
                Fp = obj.Bandwidth/2;
                Fpn = Fp./(Fs/2);
                N  = obj.FilterOrder;
                b = firceqrip(N,Fpn,[Rp,Rs],'passedge');
                wc = obj.CenterFrequency*2/Fs;
                L = size(b,2);
                b = 2*cos(pi*wc*(-(L-1)/2:(L-1)/2)).*b;
            end
        end
        
        function designFilter(obj,u)
            % Create filter System objects
            
            if isempty(coder.target)
                if obj.IsFilterDesigned
                    % This helps avoid duplicate construction of filters when
                    % setup() is called after getFilters().
                    % For codegen, this property is always false.
                    return
                end
            end
            
            coder.extrinsic('dsp.BandpassFilter.designBiquad');
            coder.extrinsic('dsp.private.LPHPFilterBase.firPrecision');
            coder.extrinsic('dsp.private.LPHPFilterBase.biquadPrecision');
            
            Fs = getSampleRate(obj);
            Fp = obj.CenterFrequency;
            Ap = obj.PassbandRipple;
            Ast = obj.StopbandAttenuation;

            switch obj.FilterType
                case 'FIR'
                    b = computeFIRCoefficients(obj);
                    if nargin  == 1
                        FIR = dsp.FIRFilter('Numerator', b);
                    elseif nargin > 1
                        % fixed-point input
                        w = fi(b, obj.CoefficientsDataType);  % vector of filter weights of type fi
                        
                        % compute output data type
                        if isinteger(u)
                            x = numerictype(class(u));
                        else
                            x = numerictype(u);                 % numerictype of input
                        end
                        [~,~,outputWL, outputFL] = coder.const(@dsp.private.LPHPFilterBase.firPrecision, ...
                            w,x);
                        outputDType = numerictype(true, outputWL, outputFL);
                        
                        FIR = dsp.FIRFilter('Numerator', b, ...
                            'CoefficientsDataType', 'Custom', ...
                            'CustomCoefficientsDataType', numerictype(true,w.WordLength,w.FractionLength), ...
                            'FullPrecisionOverride', false, ...
                            'OutputDataType', 'Custom', ...
                            'CustomOutputDataType', outputDType, ...
                            'RoundingMethod', obj.RoundingMethod);
                    end
                    
                    obj.FilterObj = FIR;
                    
%                     if isempty(coder.target)
%                         % Set metadata for FVTool
%                         if obj.DesignForMinimumOrder
%                             Fst = obj.Bandwidth;
%                             fdes = fdesign.lowpass(Fp,Fst,Ap,Ast,Fs);
%                             fmet = fdfmethod.eqriplp;
%                         else
%                             N  = obj.FilterOrder;
%                             fdes = fdesign.lowpass('N,Fp,Ap,Ast',N,Fp,Ap,Ast,Fs);
%                             fmet = fdfmethod.eqriplpfpass;
%                         end
%                         setMetaData(obj, ...
%                             fdes, fmet, ...
%                             [], 'equiripple');
%                         
%                         obj.IsFilterDesigned = true;
%                     end
                case 'IIR'
                    N  = obj.FilterOrder;
                    Fst = obj.Bandwidth;
                    
                    [SOS,SV] = coder.const(@obj.designBiquad,obj.DesignForMinimumOrder,N,Fp,Fst,Ap,Ast,Fs);
                    
                    if nargin == 1
                        % Floating point
                        IIR = dsp.BiquadFilter(...
                            'Structure', 'Direct form I',...
                            'SOSMatrix',SOS,...
                            'ScaleValues',SV);
                        
                    elseif nargin > 1
                        % fixed-point input
                        
                        % Determine DT
                        cTin = obj.CoefficientsDataType;
                        if isinteger(u)
                            x = numerictype(class(u));
                        else
                            x = numerictype(u);
                        end
                        [cT,secT,prodT,accT] = coder.const(@dsp.private.LPHPFilterBase.biquadPrecision,cTin,x);
                        
                        IIR = dsp.BiquadFilter(...
                            'Structure', 'Direct form I',...
                            'SOSMatrix',SOS,...
                            'ScaleValues',SV,...
                            'RoundingMethod', obj.RoundingMethod,...
                            'OverflowAction', 'Wrap',...
                            'SectionInputDataType', 'Custom',...
                            'CustomSectionInputDataType', secT,...
                            'SectionOutputDataType', 'Custom',...
                            'CustomSectionOutputDataType', secT,...
                            'NumeratorCoefficientsDataType','Custom',...
                            'CustomNumeratorCoefficientsDataType',cT,...
                            'DenominatorCoefficientsDataType', 'Custom',...
                            'CustomDenominatorCoefficientsDataType', cT,...
                            'ScaleValuesDataType', 'Custom',...
                            'CustomScaleValuesDataType',cT,...
                            'NumeratorProductDataType', 'Custom',...
                            'CustomNumeratorProductDataType', prodT, ...
                            'DenominatorProductDataType', 'Custom',...
                            'CustomDenominatorProductDataType', prodT,...
                            'NumeratorAccumulatorDataType','Custom',...
                            'CustomNumeratorAccumulatorDataType', accT, ...
                            'DenominatorAccumulatorDataType', 'Custom',...
                            'CustomDenominatorAccumulatorDataType', accT, ...
                            'OutputDataType', 'Custom',...
                            'CustomOutputDataType', secT ...
                            );
                    end
                    
                    obj.FilterObj = IIR;
                    
%                     if isempty(coder.target)
%                         % Set metadata for FVTool
%                         if obj.DesignForMinimumOrder
%                             fdes = fdesign.lowpass(Fp,Fst,Ap,Ast,Fs);
%                             fmet = fmethod.elliplpmin;
%                         else
%                             fdes = fdesign.lowpass('N,Fp,Ap,Ast',N,Fp,Ap,Ast,Fs);
%                             fmet = fmethod.elliplpastop;
%                         end
%                         setMetaData(obj, ...
%                             fdes, fmet, ...
%                             [], 'ellip');
%                         
%                         obj.IsFilterDesigned = true;
%                     end
            end
        end
        
        function validatePropertiesImpl(obj)
            
            Fs = obj.SampleRate;
            if (Fs == -1)
                return;
            end
            
            HalfFs = obj.SampleRate./2;
            Fc = obj.CenterFrequency;
            Ap = obj.PassbandRipple;
            Ast= obj.StopbandAttenuation;
            
            validateattributes(Ast,{'numeric'},{'finite','scalar','real',...
                'positive','>',Ap},'','StopbandAttenuation');
            
            bw = obj.Bandwidth;
            validateattributes(bw,{'numeric'},{'finite','scalar','real',...
                'positive','<',HalfFs},'','Bandwidth');
            
            validateattributes(Fc,{'numeric'},{'finite','scalar','real',...
                'positive','>',bw/2,'<',HalfFs-bw/2},'','CenterFrequency');
            
            
        end
        
        function flag = isInactivePropertyImpl(obj, prop)
            flag = isInactivePropertyImpl@dsp.private.LPHPFilterBase(obj, prop);
            minord = obj.DesignForMinimumOrder;
            switch prop
                case {'StopbandFrequency1','StopbandFrequency2',...
                      'PassbandFrequency1','PassbandFrequency2',...
                      'StopbandAttenuation1','StopbandAttenuation2'}
                    if ~minord
                        flag = true;
                    end
                case {'CenterFrequency','Bandwidth','StopbandAttenuation'}
                    if minord
                        flag = true;
                    end
            end
        end
         % Block Icon
        function icon = getIconImpl(~)
            iconStr = getString(message('dsp:system:HPLPFilter:LPBlockIcon'));
            icon = sprintf(iconStr);
        end
        
    end
    %---------------------------------------------------------------------
    methods (Static)
         function helpFixedPoint
            %helpFixedPoint Display BandpassFilter System
            %               object fixed-point information
            %   dsp.BandpassFilter.helpFixedPoint displays
            %   information about fixed-point properties and operations of
            %   the BandpassFilter System object.
            
            matlab.system.dispFixptHelp('dsp.BandpassFilter', ...
                dsp.BandpassFilter.getDisplayFixedPointPropertiesImpl);
         end
    end
    %---------------------------------------------------------------------
    methods (Static,Hidden)
        
        function props = getDisplayPropertiesImpl()
            props = {...
                
                'FilterType',...
                'DesignForMinimumOrder', ...
                'FilterOrder', ...
                'CenterFrequency', ...
                'Bandwidth', ...
                'StopbandFrequency1',...
                'PassbandFrequency1',...
                'PassbandFrequency2',...
                'StopbandFrequency2',...
                'StopbandAttenuation1',...
                'PassbandRipple', ...
                'StopbandAttenuation',...
                'StopbandAttenuation2',...
                'InheritSampleRate',...
                'SampleRate'};
        end
        
       
        
        function [SOS,SV] = designBiquad(DesignForMinimumOrder,N,Fp,Fst,Ap,Ast,Fs)
            
            if DesignForMinimumOrder
                
                fdes = fdesign.lowpass(Fp,Fst,Ap,Ast,Fs);
            else
                
                fdes = fdesign.lowpass('N,Fp,Ap,Ast',N,Fp,Ap,Ast,Fs);
            end
            fdo = fdopts.sosscaling;
            fdo.NumeratorConstraint='none';
            fdo.ScaleValueConstraint='unit';
            FilterObj = design(fdes,'ellip','SystemObject',true,...
                'FilterStructure','df2tsos','SOSScaleNorm','Linf',...
                'SOSScaleOpts',fdo);
            SOS = FilterObj.SOSMatrix;
            SV  = FilterObj.ScaleValues;
            SV(end) = 1; % Force an exact one
        end
    end
    %----------------------------------------------------------------------
    methods(Static,Access=protected)
        function groups = getPropertyGroupsImpl
            props = dsp.BandpassFilter.getDisplayPropertiesImpl;
            
            mainS = matlab.system.display.Section('Title', getString(message('dsp:system:Shared:Parameters')), ...
                'PropertyList', props,...
                'DependOnPrivatePropertyList',{'SampleRate'});
            filterVis = dsp.private.FilterVisualizationButton;
            mainS.Actions = filterVis.fvtoolAction;
            mainSG = matlab.system.display.SectionGroup('Title', getString(message('dsp:system:Shared:Main')), ...
                'Sections', mainS);

            dtSG = matlab.system.display.internal.DataTypesGroup(...
                'dsp.BandpassFilter');
            dtSG.PropertyList{2} = matlab.system.display.internal.DataTypeProperty('CoefficientsDataType', ...
                'Prefix', 'coeffs', ...
                'Description', getString(message('dsp:system:HPLPFilter:Coefficients')));
            groups = [mainSG, dtSG];
        end
        %------------------------------------------------------------------
        function header = getHeaderImpl
            % MATLAB System block header
            header = matlab.system.display.Header(...
                'dsp.BandpassFilter', ...
                'ShowSourceLink', true, ...
                'Title', getString(message('dsp:system:HPLPFilter:LPBlockName')),...
                'Text',  getString(message('dsp:system:HPLPFilter:LPBlockHeader')));
        end
    end

end
