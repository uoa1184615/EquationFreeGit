function [env, versionString] = isMatlabEnv()
% Determine environment (MATLAB,Octave) and version. Based
% upon the matlab2tikz function getEnvironment().  
% But is it used easily in our folder structure??
% AJR, 31 Mar 2019
    persistent isMatlabEnvCache
    if isempty(isMatlabEnvCache)
        env = ~exist('OCTAVE_VERSION', 'builtin');
        if env,  vData = ver(env);
             versionString = vData.Version;
        else versionString = OCTAVE_VERSION;
        end
        % store in this functions's cache
        isMatlabEnvCache.env = env;
        isMatlabEnvCache.versionString = versionString;
    else % get from the cache
        env = isMatlabEnvCache.env;
        versionString = isMatlabEnvCache.versionString;
    end
end
