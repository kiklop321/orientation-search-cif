classdef RAPIDutils
    %RAPIDUTILS Collection of utilities for manipulating RIGAKU R-AXIS data
    %   Collection of methods for manipulation with RIGAKU RAPID R-AXIS
    %   data.
    %
    % See: RAPIDread
    
properties
end
    
methods(Static,Access='public')
        
    function ret = getHeaderEntry(d,keyword)
    % r = getHeaderEntry(keyword) - return header information
    %     Extract and decode information from the headerInfo data field
    %     according to the specified keyword.
    %
    %     Examples:
    %         127.2600 = RAPIDutils.getHeaderEntry(d,'distance'), % (mm)
    %         [0.1,0.1] = RAPIDutils.getHeaderEntry(d,'pixSize'), % (mm)
    %         [,] = RAPIDutils.getHeaderEntry(d,'pixNumber'), % (pixels)
    %
    % See: RAPIDread,
    
    switch lower(keyword),
        case {'distance'}
            pattern = 'Distance: start = '; idx = strfind( d.headerInfo, pattern );
            ret = sscanf(d.headerInfo(idx+length(pattern):end),'%f',1); % (mm)
        case {'pixsize'}
            pattern = 'Pixel information'; idx = strfind( d.headerInfo, pattern );
            pattern = 'size: ('; idx1 = strfind( d.headerInfo(idx:end), pattern );
            ret = sscanf(d.headerInfo(idx+idx1+length(pattern)-1:end),'%f,%f',[1 2]); % (mm)
        case {'pixnumber'}
            pattern = 'Pixel information'; idx = strfind( d.headerInfo, pattern );
            pattern = 'number: ('; idx1 = strfind( d.headerInfo(idx:end), pattern );
            ret = sscanf(d.headerInfo(idx+idx1+length(pattern)-1:end),'%d,%d',[1 2]); % (pixels)
        otherwise
            error('RAPIDutils:Error','Unknown keyword.'),
    end
    end % getHeaderEntry
    
    function geom = getDiffractionGeometry(d)
    % g = getDifractionGeometry(d) - get diffraction experiment geometry
    %     Return a structure with geometrical parameters of the
    %     diffractometer, which can later be used to transform pixcel
    %     positions to angular and reciprocal space units and vice versa.
    %
    %     Returned structure:
    %                 distance: 127.2600 % (mm)
    %                  pixSize: [0.1000 0.1000] % (mm)
    %                pixNumber: [600 2560] % (pixels)
    %                     zero: [1285 3246] % (pixels)
    %                     tube: 'Cu'
    %                   Lambda: 1.5405 % (A)
    %                    muAir: 1.15e-3 % (1/mm)
    
    % See: RAPIDread, ...
        geom.distance = RAPIDutils.getHeaderEntry(d,'distance'); % (mm)
        geom.pixSize = RAPIDutils.getHeaderEntry(d,'pixSize'); % (mm)
        geom.pixNumber = RAPIDutils.getHeaderEntry(d,'pixNumber'); % (pixels)
        geom.zero = double( fliplr(d.zero) ); % (pixels)
        geom.tube = d.tube; 
        geom.Lambda = d.Lambda; % (A)
        geom.muAir = 1.15e-3; % (1/mm) % TODO:
        
    end % getDiffractionGeometry(d)
    
    function h = drawImgOrigin(geom,ha)
        % h = drawImgOrigin(geom,ha) - draw cross at origin of data image
        
        if ~exist('ha','var') || isempty(ha), ha = gca; end
        
        size = 40; % cross size in (pixels)
        
        axes(ha), %#ok<MAXES>
        h(1) = line(geom.zero(1)+[-size size]/2,geom.zero(2)+[0 0],'Color','r','LineWidth',2);
        h(2) = line(geom.zero(1)+[0 0],geom.zero(2)+[-size size]/2,'Color','r','LineWidth',2);
    end % drawImgOrigin
    
    function h = drawMarkAt(geom,tth,bt,ha)
        % h = drawMarkAt(geom,tth,bt,ha) - draw cross at (tth,bt) position
        
        if ~exist('ha','var') || isempty(ha), ha = gca; end
        
        size = 40; % cross size in (pixels)
        
        [x0,z0] = RAPIDutils.diff2ipPixels(geom,tth,bt);
        
        axes(ha), %#ok<MAXES>
        h(1) = line(x0+[-size size]/2,z0+[0 0],'Color','w','LineWidth',2);
        h(2) = line(x0+[0 0],z0+[-size size]/2,'Color','w','LineWidth',2);
    end % drawMarkAt
    
    function [i0,ei0,nn,xx,zz] = integIntOnDebyeRing(d,tth,dtth,varargin)
        % r = integIntOnDebyeRing(geom,tth,dtth) varargin       
        
        % check if any beta limits were set
        if ~any(cellfun( @(x)strcmpi('betaLim',x), varargin, 'UniformOutput', true)),
            % no beta limits set, use full available beta-range
            betaLim = [];
        else
            ind = find( cellfun( @(x)strcmpi('betaLim',x), varargin, 'UniformOutput', true) );
            betaLim = varargin{ind+1};
        end
    
        % get diffractometer geometry
        geom = RAPIDutils.getDiffractionGeometry(d);

        i0 = zeros(size(tth)); ei0 = i0; nn = i0;
        xx = cell(size(tth)); zz = cell(size(tth));
        for m=1:numel(tth),
            TTH = tth(m) + dtth(m,:);
            [minbeta,maxbeta] = RAPIDutils.betaLim(geom,TTH(:));
            X = []; Z = [];
            for k=1:numel(TTH),
                for l=find(~isnan(minbeta(k,:)) & ~isnan(maxbeta(k,:))),
                    % compare with user defined beta limits
                    if ~isempty(betaLim),
                        minbt = max(minbeta(k,l),betaLim(1));
                        maxbt = min(maxbeta(k,l),betaLim(2));
                    else
                        minbt = minbeta(k,l); maxbt = maxbeta(k,l);
                    end
                    if maxbt<=minbt, continue, end
                    % measure length of the Debye ring segment 
                    beta = linspace(minbt,maxbt,round((maxbt-minbt)/(2*pi)*50));
                    [x,z] = RAPIDutils.diff2ipPixels(geom,TTH(k),beta);
                    nbPoints = round( sum(sqrt(diff(x).^2+diff(z).^2)) );
                    % integrate along the ring segment
                    beta = linspace(minbt,maxbt,nbPoints);
                    [x,z] = RAPIDutils.diff2ipPixels(geom,TTH(k),beta);
                    x = round(x); z = round(z);
                    Lind = x>=1 & x<geom.pixNumber(1) & z>=1 & z<geom.pixNumber(2);
                    X = [X x(Lind)]; Z = [Z z(Lind)]; %#ok<AGROW>
                end
            end
            % convert sub-indexes to linear index
            IND = RAPIDutils.ipPixels2ImgInd(geom,X,Z);
            % find unique pixcels
            UIND = unique(IND);
            % get multiplicity
            MULT = histc(IND,UIND);
            % integrate over pixels
            i0(m) = sum(MULT.*double(d.image(UIND)));
            ei0(m) = sqrt( sum(MULT.^2.*double(d.image(UIND))) );
            nn(m) = numel(IND);
            xx{m} = X; zz{m} = Z;
        end
        
    end % integIntOnDebyeRing
    
    function [p,ep,beta,nn] = getProfileOnDebyeRing(d,tth,dtth,varargin)
        % 
        
        % check if any beta limits were set
        if ~any(cellfun( @(x)strcmpi('betaLim',x), varargin, 'UniformOutput', true)),
            % no beta limits set, use full available beta-range
            betaLim = [];
        else
            ind = find( cellfun( @(x)strcmpi('betaLim',x), varargin, 'UniformOutput', true) );
            betaLim = varargin{ind+1};
        end
        
        % check if also data-deviations should be calculated
        if ~any(cellfun( @(x)strcmpi('evalerr',x), varargin, 'UniformOutput', true)),
            evalerr = false;
        else
            ind = find( cellfun( @(x)strcmpi('evalerr',x), varargin, 'UniformOutput', true) );
            evalerr = varargin{ind+1};
        end
        
        % get diffractometer geometry
        geom = RAPIDutils.getDiffractionGeometry(d);
        % get beta limits for the given 2Theta
        [minbeta,maxbeta] = RAPIDutils.betaLim(geom,tth);
        
        beta = []; p = []; ep = []; nn = [];
        
        % loop over all segments in the following order
        order = [3 1 2 4];
        k = 1;
        for l=order,
            if isnan(minbeta(k,l)) || isnan(maxbeta(k,l)), continue, end
            % compare with user defined beta limits
            if ~isempty(betaLim),
                minbt = max(minbeta(k,l),betaLim(1));
                maxbt = min(maxbeta(k,l),betaLim(2));
            else
                minbt = minbeta(k,l); maxbt = maxbeta(k,l);
            end
            if maxbt<=minbt, continue, end
            % measure length of the Debye ring segment 
            bt = linspace(minbt,maxbt,round((maxbt-minbt)/(2*pi)*50));
            [x,z] = RAPIDutils.diff2ipPixels(geom,tth(k),bt);
            nbPoints = round( sum(sqrt(diff(x).^2+diff(z).^2)) );
            % calcualte profile
            beta0 = linspace(minbt,maxbt,nbPoints);
            [x,z] = RAPIDutils.diff2ipPixels(geom,repmat(tth(k)+dtth,numel(beta0),1),repmat(beta0.',1,numel(dtth)));
            Lind = x>=1 & x<geom.pixNumber(1) & z>=1 & z<geom.pixNumber(2);
            ind = nan(size(x)); ind(Lind) = RAPIDutils.ipPixels2ImgInd(geom,round(x(Lind)),round(z(Lind)));
            % pick-up data
            t = zeros(size(x)); w = t; t(Lind) = double(d.image(ind(Lind)));
            if numel(dtth)>1,
                w(:,1) = 0.5; w(:,2:end-1) = 1.0; w(:,end) = 0.5;
            else
                w(:,1) = 1.;
            end
            w(~Lind) = 0;
            % integrate 
            p0 = sum(t.*w,2).';
            nn0 = sum(w,2).'; % normalisation factor
            % calculate data error
            ep0 = nan(size(p0));
            if evalerr,
                for m=1:numel(beta0),
                    Lind = ~isnan(ind(m,:));
                    [uSind,uind] = unique(ind(m,Lind)); % uind is from ind(m,Lind)
                    tind = find(Lind); uind = tind(uind);
                    mult = histc(ind(m,tind),uSind);
                    ep0(m) = sqrt( sum(mult.^2.*w(uind).*t(uind)) );
                end
            end
            % add profile data
            if numel(beta)>1, % connecting data
                if abs(beta(end)-beta0(1))<1e-6, % old last, new first points identical
                    % get mean of two points
                    p(end) = (p(end)+p0(1))/2; ep(end) = (ep(end)+ep0(1))/2; nn(end) = (nn(end)+nn0(1))/2; %#ok<AGROW>
                else
                    % add first point normally
                    beta = [beta, beta0(1)]; p = [p p0(1)]; ep = [ep ep0(1)]; nn = [nn nn0(1)]; %#ok<AGROW>
                end
                % center is added normally
                beta = [beta, beta0(2:end-1)]; p = [p p0(2:end-1)]; ep = [ep ep0(2:end-1)]; nn = [nn nn0(2:end-1)]; %#ok<AGROW>
                if abs(mod(beta(1)-beta0(end),2*pi))<1e-6, % old first, new last point identical
                    % get mean of two points
                    p(1) = (p(1)+p0(end))/2; ep(1) = (ep(1)+ep0(end))/2; nn(1) = (nn(1)+nn0(end))/2; 
                else
                    % add last point normally
                    beta = [beta, beta0(end)]; p = [p p0(end)]; ep = [ep ep0(end)]; nn = [nn nn0(end)]; %#ok<AGROW>
                end    
            else
                beta = beta0; p = p0; ep = ep0; nn = nn0;
            end
        end
    end % getProfileOnDebyeRing
    
    function [p,tth,bt] = getTThBtIntensityMap(d,tth0,bt0,rtth0,rbt0)
        
        % constants
        deg2rad = pi/180.;
        
        % get diffractometer geometry
        geom = RAPIDutils.getDiffractionGeometry(d);
        
        % get local metrics
        [x,z] = RAPIDutils.diff2ipPixels(geom,tth0+[0.0 0.1]*deg2rad,bt0); dtth0 = 0.1*deg2rad/sqrt(diff(x)^2+diff(z)^2);
        [x,z] = RAPIDutils.diff2ipPixels(geom,tth0,bt0+[0.0 0.1]*deg2rad); dbt0  = 0.1*deg2rad/sqrt(diff(x)^2+diff(z)^2);

        % create (tth,bt) map grid
        tth = tth0 + (-rtth0/2:dtth0:rtth0/2); bt = bt0 + (-rbt0/2:dbt0:rbt0/2);
        [Xtth,Xbt] = meshgrid(tth,bt);

        % extract data
        [x,z] = RAPIDutils.diff2ipPixels(geom,Xtth,Xbt);
        x = round(x); z = round(z);
        Lind = x>=1 & x<geom.pixNumber(1) & z>=1 & z<geom.pixNumber(2);
        p = nan(size(x)); p(Lind) = d.image( RAPIDutils.ipPixels2ImgInd(geom,x(Lind),z(Lind)) );        
    end % getTThBtIntensityMap
    
    function h = drawImgDebyeRing(geom,tth,ha)
        % h = drawImgDebyeRing(tth,geom,ha) - draw Debye ring
        
        if ~exist('ha','var') || isempty(ha), ha = gca; end
        
        [minbeta,maxbeta] = RAPIDutils.betaLim(geom,tth);
        
        h = [];
        axes(ha), set(ha,'NextPlot','add'), %#ok<MAXES>
        %hold on,
        for k=1:numel(tth),
            for l=find(~isnan(minbeta(k,:)) & ~isnan(maxbeta(k,:))),
                beta = linspace(minbeta(k,l),maxbeta(k,l),round((maxbeta(k,l)-minbeta(k,l))/(2*pi)*500));
                [x,z] = RAPIDutils.diff2ipPixels(geom,tth(k),beta);
                h(end+1) = line(x,z); %#ok<AGROW>
                %plot([x(1) x(end)],[z(1) z(end)],'b*'),
            end
        end
%        h(1) = line(geom.zero(1)+[-size size]/2,geom.zero(2)+[0
%        0],'Color','b','LineWidth',2);
    end % drawImgDebyeRing
    
    function [eta,delta] = diff2ipAngles(tth,beta)
        % [eta,delta] = diff2ipAngles(tth,beta)
        
        % input unification
        if isscalar(beta) && ~isscalar(tth), beta = beta*ones(size(tth)); end
        if isscalar(tth) && ~isscalar(beta), tth = tth*ones(size(beta)); end

        % delta
        s = sin(tth);
        c = cos(tth);
        delta = asin( sin(beta).*s );
        % eta
        cb = cos(beta);
        Lindt = abs(cb)>1e-6;
        
        eta = nan(size(tth));
        eta(Lindt) = atan2( cb(Lindt).*s(Lindt), c(Lindt) );
        
        if ~all(Lindt), % beta = pi/2+n*pi
            ind = find(~Lindt);
            Lindt1 = tth(ind)<pi/2;
            eta(ind(Lindt1)) = 0.;
            Lindt1 = tth(ind)>pi/2;
            eta(ind(Lindt1)) = pi;
            Lindt1 = abs(tth(ind)-pi/2)<1e-6;
            eta(ind(Lindt1)) = pi/2;
        end
    end % diff2ipAngles
    
    function ind = ipPixels2ImgInd(geom,x,z)
        % ind = ipPixels2ImgInd(geom,x,z)
        
        ind = sub2ind(fliplr(geom.pixNumber),z,x);
    end % ipPixels2ImgInd
    
    function [x,z] = diff2ipPixels(geom,tth,beta)
        % [x,z] = diff2ipPixels(tth,beta)
        
        % transform to image-plate angles (eta,delta)
        [eta,delta] = RAPIDutils.diff2ipAngles(tth,beta);
        
        % calculate positions in pixels
        z = geom.zero(2) - geom.distance*tan(delta)/geom.pixSize(2);
        x = geom.zero(1) - geom.distance*eta/geom.pixSize(1);
    end % diff2ipPixels
    
    function jac = diff2ipAnglesJacobian(tth,beta)
        % jac = diff2ipAnglesJacobian(geom,tth,beta)
        
        jac = -2*sin(tth)./sqrt( 3+cos(2*tth)+2*cos(2*beta).*sin(tth).^2 );
    end % diff2ipAnglesJacobian
    
    function [tth,beta] = ipAngles2diff(eta,delta)
        % [tth,beta] = ipAngles2diff(eta,delta)
        % beta is returned in range (0, 2*pi), where beta=0 coincides with
        % -x axis
        
        beta = atan2( tan(delta), sin(eta) );
        tth = abs( acos( cos(delta) .* cos(eta) ) );
        Lind = beta<-pi/2; beta(Lind) = beta(Lind) + 2*pi;
    end % ipAngles2diff
    
    function [tth,beta] = ipPixels2diff(geom,x,z)
        % [tth,beta] = ipPixels2diff(geom,x,z)
        % beta is returned in range (-pi/2, 3/2*pi)
        
        % get image-plate angles (eta,delta)
        [eta,delta] = RAPIDutils.ipPixels2ipAngles(geom,x,z);
        % get diffraction angles
        [tth,beta] = RAPIDutils.ipAngles2diff(eta,delta);
    end % ipPixels2diff
    
    function [eta,delta] = ipPixels2ipAngles(geom,x,z)
        % [eta,delta] = ipPixels2ipAngles(geom,x,z)
        
        % calculate positions in ip-Angles
        delta = atan( (geom.zero(2)-z)*geom.pixSize(2)/geom.distance );
        eta = (geom.zero(1)-x)*geom.pixSize(1)/geom.distance;
    end % ipPixels2ipAngles
    
    function [minbeta,maxbeta] = betaLim(geom,tth)
        %
        % beta ... (-pi/2,3/2pi)
        
        tth = abs(tth);
        
        minbeta = nan(numel(tth),4);
        maxbeta = nan(numel(tth),4);
        
        maxdelta =  atan((geom.zero(2)-1)/geom.distance*geom.pixSize(2));
        mindelta = -atan((geom.pixNumber(2)-geom.zero(2))/geom.distance*geom.pixSize(2));
        
        maxeta =  (geom.zero(1)-1)/geom.distance*geom.pixSize(1);
        mineta = -(geom.pixNumber(1)-geom.zero(1))/geom.distance*geom.pixSize(1);
        
        s = sin(tth);
        smin = sin(mindelta);
        smax = sin(maxdelta);

        % left ip-side (single-part half-rings)
        Lind = tth<=maxeta & s<=abs(smin);
        minbeta(Lind,1) = -pi/2;
        Lind = tth<=maxeta & s>abs(smin);
        minbeta(Lind,1) = asin( smin./s(Lind) );
        Lind = tth<=maxeta & s<=smax;
        maxbeta(Lind,1) = pi/2;
        Lind = tth<=maxeta & s>smax;
        maxbeta(Lind,1) = asin( smax./s(Lind) );
        if maxeta<=pi/2,
            % (double-part half-rings - up)
            maxtth = acos( cos(maxdelta)*cos(maxeta) );
            Lind = tth>=maxeta & tth<=maxtth & s<=smax;
            maxbeta(Lind,1) = pi/2;
            Lind = tth>=maxeta & tth<=maxtth & s>smax;
            maxbeta(Lind,1) = asin( smax./s(Lind) );
            Lind = tth>=maxeta & tth<=maxtth;
            if any(Lind),
                [eta,delta] = RAPIDutils.diff2ipAngles(tth(Lind),maxbeta(Lind,1));
                minbeta(Lind,1) = acos( cos(eta).*cos(delta).*tan(maxeta)./s(Lind) );
            end
            % (double-part half-rings - bottom)
            maxtth = acos( cos(mindelta)*cos(maxeta) );
            Lind = tth>=maxeta & tth<=maxtth & s<=-smin;
            minbeta(Lind,3) = -pi/2;
            Lind = tth>=maxeta & tth<=maxtth & s>-smin;
            minbeta(Lind,3) = asin( smin./s(Lind) );
            Lind = tth>=maxeta & tth<=maxtth;
            if any(Lind),
                [eta,delta] = RAPIDutils.diff2ipAngles(tth(Lind),minbeta(Lind,3));
                maxbeta(Lind,3) = -acos( cos(eta).*cos(delta).*tan(maxeta)./s(Lind) );
            end
        else
            % correct also long-lines
            maxtth = acos( cos(maxdelta)*cos(maxeta) );
            Lind = tth<=maxeta & tth>=maxtth;
            maxbeta(Lind,1) = acos( tan(maxeta)./tan(tth(Lind)) );
            maxtth = acos( cos(mindelta)*cos(maxeta) );
            Lind = tth<=maxeta & tth>=maxtth;
            minbeta(Lind,1) = -acos( tan(maxeta)./tan(tth(Lind)) );
        end
        % right ip-side (single-part half-rings)
        Lind = tth<=abs(mineta) & s<=smax;
        minbeta(Lind,2) = pi/2;
        Lind = tth<=abs(mineta) & s>smax;
        minbeta(Lind,2) = pi - asin( smax./s(Lind) );
        Lind = tth<=abs(mineta) & s<=abs(smin);
        maxbeta(Lind,2) = 3/2*pi;
        Lind = tth<=abs(mineta) & s>abs(smin);
        maxbeta(Lind,2) = pi + asin( abs(smin)./s(Lind) );
        if abs(mineta)<=pi/2,
            % (double-part half-rings - up)
            maxtth = acos( cos(maxdelta)*cos(mineta) );
            Lind = tth>=abs(mineta) & tth<=maxtth & s<=smax;
            minbeta(Lind,2) = pi/2;
            Lind = tth>=abs(mineta) & tth<=maxtth & s>smax;
            minbeta(Lind,2) = pi - asin( smax./s(Lind) );
            Lind = tth>=abs(mineta) & tth<=maxtth;
            if any(Lind),
                [eta,delta] = RAPIDutils.diff2ipAngles(tth(Lind),maxbeta(Lind,1));
                maxbeta(Lind,2) = pi - acos( cos(eta).*cos(delta).*tan(abs(mineta))./s(Lind) );
            end
            % (double-part half-rings - bottom)
            maxtth = acos( cos(mindelta)*cos(mineta) );
            Lind = tth>=abs(mineta) & tth<=maxtth & s<=-smin;
            maxbeta(Lind,4) = 3/2*pi;
            Lind = tth>=abs(mineta) & tth<=maxtth & s>-smin;
            maxbeta(Lind,4) = pi - asin( smin./s(Lind) );
            Lind = tth>=abs(mineta) & tth<=maxtth;
            if any(Lind),
                [eta,delta] = RAPIDutils.diff2ipAngles(tth(Lind),maxbeta(Lind,4));
                minbeta(Lind,4) = pi + acos( cos(eta).*cos(delta).*tan(abs(mineta))./s(Lind) );
            end
        else
            % correct also long-lines
            maxtth = acos( cos(maxdelta)*cos(mineta) );
            Lind = tth<=abs(mineta) & tth>=maxtth;
            minbeta(Lind,2) = pi - acos( tan(abs(mineta))./tan(tth(Lind)) );
            maxtth = acos( cos(mindelta)*cos(mineta) );
            Lind = tth<=abs(mineta) & tth>=maxtth;
            maxbeta(Lind,2) = pi + acos( tan(abs(mineta))./tan(tth(Lind)) );
        end
    end % betaLim(geom,tth)
    
    function ha = image(d,ha)
    % h = image(d,hp) - iamge data
    %     Visualise data as an image. Uses HA axes handle if supplied.
    %     Returns iamge axes handle.
    %
    %
    % See: RAPIDread, ...
        % check input arguments
        if ~exist('ha','var'), ha = []; end
        
        % activate axes
        if isempty(ha), ha = findobj('Tag','RAPIDutils-image Axes'); end
        if isempty(ha),
            figure,
            ha = axes('Tag','RAPIDutils-image Axes'); % new axes
        else
            axes(ha), %#ok<MAXES>
            cla(ha),
        end
        hold on, axis ij, s = size(d.image); xlim([1 s(2)]), ylim([1,s(1)]),
        
        % get diffractometer geometry
        geom = RAPIDutils.getDiffractionGeometry(d);
        
        % plot data
        [x,y] = meshgrid(1:size(d.image,2),1:size(d.image,1)); Lind = (x-double(d.zero(2))).^2 + (y-double(d.zero(1))).^2 > 40^2; % exclude origin
        dd = sqrt(double(d.image)); n = ceil(maxmax(dd(Lind))); image(dd),
        colormap(hot(n)), hcb = colorbar;
        YTick = get(hcb,'YTick');
        YTickLabel = cell( numel(YTick) );
        for k=1:numel(YTick),
            YTickLabel{k} = sprintf('%d',YTick(k).^2);
        end
        set(hcb,'YTickLabel',YTickLabel),
        
        % draw a cross at origin
        RAPIDutils.drawImgOrigin(geom);
    end % ha = image(d,ha)
    
    function [dtth,dbeta] = getDiffPositionCorrZeroShift(geom,dx0,dz0,tth0,beta0)
        % [dtth,dbt] = getDiffPositionCorrZeroShift(geom,dx0,dz0,tth0,bt0)
        
        % get ip-pixcels positions
        [x,z] = RAPIDutils.diff2ipPixels(geom,tth0,beta0);
        % get corrected diffraction positions
        [tth,beta] = RAPIDutils.ipPixels2diff(geom,x-dx0,z-dz0);
        % calcualted correction
        dtth = cycleDistance(tth,tth0,2*pi); dbeta = cycleDistance(beta,beta0,2*pi); 
    end % getDiffPositionCorrZeroShift
    
    function [dtth,dbeta,displx] = getDiffPositionCorrSampleDispl(geom,displ,tth0,beta0,om,chi)
        % [dtth,dbeta,displx] = getDiffPositionCorrSampleDispl(geom,displ,om,chi,tth0,beta0)
        % displ ... sample displacement (in mm)
        
        if isscalar(om) && isscalar(chi),
            % calculate x-coordinate of sample displacement vector in the
            % diffractometer reference system
%             t = sample2difracMatrix(0.0,chi,om)*[0; 0; 1];
%             displx = t(1) * displ;
            displx = displ;

            % get delta ip-angle
            [eta0,delta0] = RAPIDutils.diff2ipAngles(tth0,beta0); %#ok<ASGLU>
                        
            % calculate correction
            dtth = asin( sin(tth0) .* displx./geom.distance .* cos(delta0) );
            dbeta = zeros(size(dtth));
        else
            dtth = nan(size(tth0)); dbeta = nan(size(tth0)); displx = nan(size(tth0)); 
            for k=1:numel(tth0),
                [dtth(k),dbeta(k),displx(k)] = RAPIDutils.getDiffPositionCorrSampleDispl(geom,displ,tth0(k),beta0(k),om(k),chi(k));
            end
        end
    end % getDiffPositionCorrSampleDispl
    
    function [dtth,dbeta] = getDiffPositionCorrIpTilt(tilt,tth0,beta0)
        % [dtth,dbeta] = getDiffPositionCorrIpTilt(geom,ddelta,tth0,beta0)
        
        [eta0,delta0] = RAPIDutils.diff2ipAngles(tth0,beta0);
        [tth,beta] = RAPIDutils.ipAngles2diff(eta0,delta0+tilt);
        
        dtth = tth - tth0;
        dbeta = beta - beta0;
    end % getDiffPositionCorrIpTilt
    
end % static, public methods

end % RAPIDutils

