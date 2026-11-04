function [s, e] = alongTrackKm(lat, lon, latO, lonO, orient)
%ALONGTRACKKM  Signed distance (km) along a flight transect.
%
%   s = alongTrackKm(lat, lon, latO, lonO)
%   s = alongTrackKm(lat, lon, latO, lonO, orient)
%
%   lat, lon      transect points (degrees)
%   latO, lonO    origin: the point on the transect where the mean wind
%                 vector through the centre of Marion County crosses it
%                 (see makeOrigin below)
%   orient        'east' (default for E-W legs) or 'north' (for N-S legs):
%                 direction in which s increases.  If omitted, the dominant
%                 extent of the leg decides.
%
%   The transect is projected onto its own principal axis (PCA on the track
%   in km, as in Compare_WRF_Obs_v5), so the result is distance ALONG the
%   flight track rather than a longitude or latitude difference.  s = 0 at
%   the origin.  Returned e is the unit vector of the axis [east north].
%
%   Convention used in the revised manuscript (Figure 2 caption):
%   "distance along the flight track, with the origin at the intersection of
%    the transect with the mean wind vector through the center of Marion
%    County, and distance increasing toward the east".

    lat = lat(:); lon = lon(:);
    la0 = mean(lat, 'omitnan');
    kmE = (lon - lonO) * 111.32 * cosd(la0);     % km east of the origin
    kmN = (lat - latO) * 110.57;                 % km north of the origin
    ok  = ~isnan(kmE) & ~isnan(kmN);

    % principal axis of the track
    X = [kmE(ok) - mean(kmE(ok)), kmN(ok) - mean(kmN(ok))];
    [V, D] = eig(cov(X));
    [~, m] = max(diag(D));
    e = V(:, m)';                                 % [east north]

    if nargin < 5 || isempty(orient)
        if abs(e(1)) >= abs(e(2)), orient = 'east'; else, orient = 'north'; end
    end
    switch lower(orient)
        case 'east',  if e(1) < 0, e = -e; end
        case 'north', if e(2) < 0, e = -e; end
        otherwise, error('orient must be ''east'' or ''north''');
    end
    s = kmE * e(1) + kmN * e(2);
end
