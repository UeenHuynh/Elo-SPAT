function meta = research_only_metadata()
% research_only_metadata  Shared disclaimer/citation metadata for exported artifacts.

meta = struct();
meta.disclaimer = ['Research-only output. This module is not clinically validated and must not be ' ...
    'used for diagnosis, treatment selection, or treatment delivery.'];
meta.machineModelNote = ['Machine timing and beam-delivery estimates are configurable proxy values, ' ...
    'not vendor-validated machine measurements.'];
meta.citationPlaceholder = 'TODO: add manuscript / software citation before external release.';
end
