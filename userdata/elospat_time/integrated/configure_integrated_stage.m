function cfg = configure_integrated_stage(cfg)
% configure_integrated_stage  Activate the requested incremental formulation.

stage = lower(char(cfg.integrated.componentStage));
switch stage
    case 'layer'
        cfg.lambdaB = 0;
        cfg.lambdaS = 0;
        cfg.lambdaES = 0;
        cfg.lambdaT = 0;
    case 'spot'
        cfg.lambdaB = 0;
        cfg.lambdaES = 0;
        cfg.lambdaT = 0;
    case 'full'
        % All configured terms are active.
    otherwise
        error('configure_integrated_stage:invalidStage', ...
            'componentStage must be ''layer'', ''spot'', or ''full''.');
end
end
