%Sys = struct('g',2,'S',.5); % No hyperfine
%Sys = struct('g',2,'S',.5,'A',1,'Nucs','1H'); % hyperfine
Sys = struct('g',2,'S',.5,'A',1,'Nucs','1H'); % hyperfine
Exp = struct();
%Exp = struct('Range',[0,3]);

% Seems sum of the amplitudes for g=2 is always around 195.8949

field_steps = 50;
max_field = 50;

x = [];
y = [];
z = [];

for n = 0:field_steps
    mag_field = n*max_field/field_steps;
    fprintf('%3.f%% - %2.0f mT',100*n/field_steps,mag_field) % display percentage complete and current field
    %disp([int2str(mag_field),'mT'])
    Exp.Field = mag_field;
    [Freq,Amp,~,Trans] = resfreqs_matrix(Sys,Exp);
    sum(Amp)
    field = repelem(mag_field,length(Freq))';
    x = [x; field];
    y = [y; Freq];
    z = [z; Amp];
    
    %fprintf('%c%c%c%c%c%c%c%c%c%c%c%c',8,8,8,8,8,8,8,8,8,8,8,8)%delete last 12 characters
end
fprintf('\n')

%figure


scatter(x,y,[],z,'.')
colormap(flipud(hot))
cbar = colorbar();
cbar.Label.String = 'Amplitude';
%caxis([0.0,1.0])
%ylim(Exp.Range*1000)
xlabel('B (mT)')
ylabel('Transition frequency (MHz)')
title(['Spin 1/2'])
%title(['Mag axis: ',magaxis,'; MW axis: ',MWaxis])
%saveas(gcf,['figure',int2str(step),'.png'])
saveas(gcf,['spin-half.png']);
%text_label = {['Source: ',parameter_source],['Rotation axis: ',num2str(rot_axis')]};
%annotation('textbox',[.2 .5 .3 .3],'string',text_label,'FitBoxToText','on');
