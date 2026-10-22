function [Weave, Tow] = WeaveAnalysisSilent(Wp)
s = evalc('[Weave, Tow] = RunWeave(Wp);');
end
