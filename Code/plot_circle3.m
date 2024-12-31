function plot_circle3(thetaAngle, value, mark, varargin)
%% Circular diagram.
%% -------------------------------------------------------------------------------------------------------------------------------
%% Input
%  (1) thetaAngle: angles representing the filled value (0 - 360)  
%  (2) value: values corresponding to angles    
%  (3) mark: markers corresponding to angles (0 or 1)       
%  (4) optional inputs:
%                   >>> 'outerRadius', 'outR', outer circle radius (default, 10.5)
%                   >>> 'middleRadius', 'midR', middle circle radius (default, 10)
%                   >>> 'innerRadius', 'innR', inner circle radius (default, 5)
%                   >>> 'valueColor', color corresponding to 'value'
%                   >>> 'markColor', color corresponding to 'mark'  
%% Example
% thetaAngle = [0:90];
% value = [1:91];
% mark = zeros(1,91);
% mark(1:10) = 1;
% plot_circle3(thetaAngle, value, mark)
%---------------------------------------------------------------------------------------------------------------------------------------------------%
% - Z.K.X. 2024/04/02
%---------------------------------------------------------------------------------------------------------------------------------------------------%

%% default setting 
outerRadius = 10.5;
middleRadius = 10;  
innerRadius = 5;
centerX = 0;
centerY = 0;

valueColor = getColor;
markColor = [0.9 0.9 0.9];

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'outerRadius', 'outR'}
                outerRadius = varargin{i+1};
            case {'middleRadius', 'midR'}
                middleRadius = varargin{i+1};
            case {'innerRadius', 'innR'}
                innerRadius = varargin{i+1};
            case {'valueColor', 'vC'}
                valueColor = varargin{i+1};
            case {'markColor', 'mC'}
                markColor = varargin{i+1};                
        end
    end
end

if nargin < 3 || isempty(mark)
    mark = zeros(length(value),1);
end

value_color = color2value(valueColor,value,length(value));

%% plotting circles
rectangle('Position',[centerX-outerRadius, centerY-outerRadius, 2*outerRadius, 2*outerRadius],...
          'Curvature',[1,1], 'EdgeColor', 'k', 'LineWidth', 2);
hold on; axis equal;  

rectangle('Position',[centerX-middleRadius, centerY-middleRadius, 2*middleRadius, 2*middleRadius],...
          'Curvature',[1,1], 'EdgeColor', 'k', 'FaceColor', 'w', 'LineWidth', 2);

rectangle('Position',[centerX-innerRadius, centerY-innerRadius, 2*innerRadius, 2*innerRadius],...
          'Curvature',[1,1], 'EdgeColor', 'k', 'FaceColor', 'w', 'LineWidth', 2);

%% filling values
for i = 1:length(thetaAngle)
    startAngle = thetaAngle(i)-0.5;
    endAngle = thetaAngle(i)+0.5;
    
    startAngleRad = deg2rad(startAngle);
    endAngleRad = deg2rad(endAngle);
    
    theta = linspace(startAngleRad, endAngleRad);
    xOuterArc = centerX + middleRadius * cos(theta);
    yOuterArc = centerY + middleRadius * sin(theta);
    xInnerArc = centerX + innerRadius * cos(theta(end:-1:1));
    yInnerArc = centerY + innerRadius * sin(theta(end:-1:1));
    
    patch([xOuterArc xInnerArc], [yOuterArc yInnerArc], value_color(i,:),'EdgeColor', 'none');

    if mark(i) == 1
        xOuterArc2 = centerX + outerRadius * cos(theta);
        yOuterArc2 = centerY + outerRadius * sin(theta);
        xInnerArc2 = centerX + middleRadius * cos(theta(end:-1:1));
        yInnerArc2 = centerY + middleRadius * sin(theta(end:-1:1));

        patch([xOuterArc2 xInnerArc2], [yOuterArc2 yInnerArc2], markColor);
    end
end

axis off

end
%%
function corlor = getColor()

corlor = [0.368627450980392	0.309803921568628	0.635294117647059
0.332902962548228	0.357721391577151	0.662056358721288
0.297410156499750	0.403625043100944	0.685496577686323
0.262190498809226	0.448560808248045	0.706300992018127
0.227253292838168	0.493158051920390	0.725190871451502
0.192935898134102	0.537843600479998	0.742750968131035
0.194917463625426	0.586164422285955	0.742074164445842
0.238534540078604	0.637581524678404	0.719736058071074
0.294417778161078	0.687928657669767	0.687785439646511
0.352244600735673	0.732580099971520	0.658881765806444
0.412141038529684	0.766623612699109	0.646137917746237
0.471160671311987	0.792663020566297	0.643958573071747
0.527897955864135	0.814598814961382	0.644445362723125
0.582579112579655	0.834093012647500	0.645610790696595
0.635147044405028	0.852907942384945	0.645350435367962
0.685818965271568	0.872862404052358	0.641133513206310
0.739150480235605	0.894283510896355	0.628193274137407
0.792378256834256	0.915760968253365	0.610733399555541
0.841269745626912	0.935771299381019	0.595787791521905
0.882384514228659	0.952734869482390	0.590910411630564
0.913164829929038	0.965639048892307	0.604010373240172
0.938253226604679	0.978184038872055	0.637118155506859
0.959436484971286	0.989677093394554	0.679555277371867
0.977267622337847	0.998033622257167	0.719449292662263
0.992758944259555	1	0.745131436748527
1	0.996091982777108	0.744767771139698
1	0.979123685064294	0.714366179891001
1	0.953289252191960	0.665187646744921
1	0.922535647727882	0.609689085101958
0.997985926415667	0.890669761508102	0.560320574599493
0.994481038739963	0.859385127175737	0.524456694930908
0.994186645738301	0.823329455530606	0.488802813479127
0.995036348602311	0.783089803491765	0.453116664779688
0.995233286539529	0.740063303478644	0.419150215181151
0.993321362908456	0.695696526020084	0.388626398684196
0.989844055548408	0.649301715379478	0.361676483457733
0.987508070014149	0.596770827105191	0.335201127201337
0.983499347950103	0.540993028012321	0.309935250193174
0.975251549955637	0.485835612321005	0.286760910295060
0.960691141463968	0.436191051972805	0.266418034331519
0.941289731179181	0.394297645467696	0.256777872871590
0.921700868432271	0.354745948681751	0.266622680528681
0.899673120062352	0.316920673415453	0.285103654804381
0.872832685669552	0.280744705060415	0.302401525669884
0.839113445393712	0.246446678861223	0.309905351180865
0.799644534080146	0.212027150558578	0.306652708472718
0.757946714564209	0.173816603741874	0.299847504690544
0.714035973978013	0.130623234728162	0.289581176133295
0.667914767073927	0.0784367485504265	0.275887896758510
0.619607843137255	0.00392156862745120	0.258823529411765];

end

%%
function [value_color,n] = color2value(color,value,n)

% color: RGB N * 3 OR 1 * 3
% value: value to be mapped
% n: number of sections

if size(color,1) > 1
    n = size(color,1);
    a = [];
else
    if nargin < 3 | isempty(n)
        n = 10;
    end        
    a = 0:1/(n):1; a(1) = [];
end

if isempty(a)
    color = color;
else
    p = color;
    for i = 1:length(a)
        color(i,:) = p + ([1 1 1] - p)*( 1 - a(i)/1); 
    end
end

k = 100/(n);
k = [k:k:100];
for i = 1:numel(k)
    edge(i) = prctile(value,k(i));
end

edge = [0,edge];

for i = 1:length(edge) - 1
    f = find(value > (edge(i)) & value <= (edge(i+1)));        
    value_color(f,:) = repmat(color(i,:),numel(f),1);
end

f = find(mean(value,2)==0 | isnan(mean(value,2)));
value_color(f,:) = repmat([1 1 1],numel(f),1);
end 