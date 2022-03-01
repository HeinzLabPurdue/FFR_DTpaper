function add_stat_txt(mdl, params)
x_txt_val= params.x_txt_val;
y_txt_val= params.y_txt_val;
y_txt_gap= params.y_txt_gap;
fSize= params.fSize;
pValThresh= params.pValThresh;


if ~isempty(params.title)
    text(x_txt_val,y_txt_val,sprintf('%s', params.title), 'units', 'normalized', 'fontsize', fSize, 'FontName', 'Arial', 'FontAngle', 'italic');
    text(x_txt_val,y_txt_val-y_txt_gap,sprintf('%s', '\bf------------'), 'units', 'normalized', 'fontsize', fSize, 'FontName', 'Arial', 'Interpreter', 'latex', 'FontWeight', 'bold');
else 
    y_txt_val= y_txt_val+2*y_txt_gap;
end

% text(x_txt_val,y_txt_val-2*y_txt_gap,sprintf('$r=%.2f$', mdl.Coefficients.Estimate(2)), 'units', 'normalized', 'fontsize', fSize, 'FontName', 'Arial');
if mdl.Coefficients.pValue(2)>.05
    text(x_txt_val,y_txt_val-2*y_txt_gap,sprintf('p=%.2f', mdl.Coefficients.pValue(2)), 'units', 'normalized', 'fontsize', fSize, 'FontName', 'Arial', 'FontAngle', 'italic');
elseif mdl.Coefficients.pValue(2)>pValThresh
    text(x_txt_val,y_txt_val-2*y_txt_gap,sprintf('p=%.3f', mdl.Coefficients.pValue(2)), 'units', 'normalized', 'fontsize', fSize, 'FontName', 'Arial', 'FontAngle', 'italic');
else
    text(x_txt_val,y_txt_val-2*y_txt_gap,sprintf('p<%.3f', pValThresh), 'units', 'normalized', 'fontsize', fSize, 'FontName', 'Arial', 'FontAngle', 'italic');
end
text(x_txt_val,y_txt_val-3.5*y_txt_gap,sprintf('R^2=%.2f', mdl.Rsquared.Ordinary), 'units', 'normalized', 'fontsize', fSize, 'FontName', 'Arial', 'FontAngle', 'italic');