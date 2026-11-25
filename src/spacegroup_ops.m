function [ops, info] = spacegroup_ops(varargin)
%SPACEGROUP_OPS  Symmetry operations of a space group from its symbol.
%
%   [ops, info] = spacegroup_ops('hall', '-I 4 2 3')
%   [ops, info] = spacegroup_ops('hm', 'I m -3 m', 'cell', cell)
%   [ops, info] = spacegroup_ops('number', 229, 'cell', cell)
%
%   Looks the space group up in the CCP4 symmetry library
%   (src/symmetry_data/syminfo.lib: all 230 groups in 541 settings) and
%   returns all its operations, primitive operations combined with the
%   centring translations, as text ('x,y,z', '-y+1/2,x,z+1/4', ...).
%   Several keys can be given; the first one that is found is used, in the
%   order hall, hm, number.
%
%   Settings that a plain H-M symbol or a number do not fix:
%     rhombohedral groups (R...)   hexagonal axes (:H) if gamma = 120 deg in
%                                  'cell', rhombohedral axes (:R) if
%                                  alpha = beta = gamma ~= 90 deg
%     two origin choices (:1, :2)  origin choice 2 (ITA standard) with a warning
%     number only                  the first (standard) setting of the library
%
%   info.number, info.hall, info.xhm, info.matched_by ('hall', 'hm', 'number')

ip = inputParser;
ip.addParameter('hall', '');
ip.addParameter('hm', '');
ip.addParameter('number', NaN);
ip.addParameter('cell', []);
ip.parse(varargin{:});
o = ip.Results;

lib = load_library();
k = [];
if ~isempty(o.hall)
    key = norm_symbol(o.hall);
    k = find(strcmp({lib.hall_n}, key), 1);
    by = 'hall';
end
if isempty(k) && ~isempty(o.hm)
    key = norm_symbol(regexprep(char(o.hm), '\s+[Ss]$', ''));   % 'P 32 2 1 S' -> 'p3221'
    k = match_hm(lib, key, o.cell);
    by = 'hm';
end
if isempty(k) && ~isnan(o.number)
    c = find([lib.number] == o.number);
    k = pick_setting(lib, c, o.cell);
    by = 'number';
    if ~isempty(k) && (~isempty(o.hall) || ~isempty(o.hm))
        warning('spacegroup_ops:setting', ['Symbol "%s" / "%s" is not in the library (non-standard ' ...
            'setting?); the standard setting %s of No. %d is used. Check that it matches the cell ' ...
            'and coordinates of the CIF.'], o.hall, o.hm, lib(k).xhm, o.number);
    end
end
if isempty(k)
    error('spacegroup_ops:notFound', 'Space group not found (hall "%s", hm "%s", number %g).', ...
        o.hall, o.hm, o.number);
end
e = lib(k);
ops = {};
for i = 1:numel(e.cenops)
    for j = 1:numel(e.symops)
        ops{end+1, 1} = combine(e.symops{j}, e.cenops{i}); %#ok<AGROW>
    end
end
info = struct('number', e.number, 'hall', e.hall, 'xhm', e.xhm, 'matched_by', by);
end

%% ------------------------------------------------------------------------
function lib = load_library()
persistent L
if ~isempty(L), lib = L; return; end
file = fullfile(fileparts(mfilename('fullpath')), 'symmetry_data', 'syminfo.lib');
lines = strtrim(regexp(fileread(file), '\r?\n', 'split'));
L = struct('number', {}, 'hall', {}, 'hall_n', {}, 'xhm', {}, 'xhm_n', {}, 'old_n', {}, ...
           'symops', {}, 'cenops', {});
e = [];
for i = 1:numel(lines)
    s = lines{i};
    if strcmp(s, 'begin_spacegroup')
        e = struct('number', NaN, 'hall', '', 'hall_n', '', 'xhm', '', 'xhm_n', '', 'old_n', {{}}, ...
                   'symops', {{}}, 'cenops', {{}});
    elseif strcmp(s, 'end_spacegroup')
        L(end+1) = e; %#ok<AGROW>
    elseif isempty(e)
        continue
    elseif startsWith(s, 'number')
        e.number = str2double(strtrim(s(7:end)));
    elseif startsWith(s, 'symbol Hall')
        e.hall = quoted(s, 1);  e.hall_n = norm_symbol(e.hall);
    elseif startsWith(s, 'symbol xHM')
        e.xhm = quoted(s, 1);   e.xhm_n = norm_symbol(e.xhm);
    elseif startsWith(s, 'symbol old')
        e.old_n = cellfun(@norm_symbol, {quoted(s, 1), quoted(s, 2)}, 'UniformOutput', false);
    elseif startsWith(s, 'symop')
        e.symops{end+1} = strtrim(s(6:end));
    elseif startsWith(s, 'cenop')
        e.cenops{end+1} = strtrim(s(6:end));
    end
end
lib = L;
end

function q = quoted(s, n)
t = regexp(s, '''([^'']*)''', 'tokens');
if numel(t) >= n, q = t{n}{1}; else, q = ''; end
end

function k = norm_symbol(s)
% 'P 1 21/c 1' -> 'p121/c1', 'F d -3 m :1' -> 'fd-3m:1'
k = lower(regexprep(char(s), '[\s_]', ''));
end

function k = match_hm(lib, key, cell)
xhm = {lib.xhm_n};
k = find(strcmp(xhm, key), 1);                          % exact, incl. :1 / :2 / :H / :R
if ~isempty(k), return; end
base = regexprep(xhm, ':.*$', '');
c = find(strcmp(base, key));                            % without setting suffix
if isempty(c)                                           % old CCP4 names / full symbols
    c = find(cellfun(@(o) any(strcmp(o, key)), {lib.old_n}));
end
k = pick_setting(lib, c, cell);
end

function k = pick_setting(lib, c, cell)
k = [];
if isempty(c), return; end
if numel(c) == 1, k = c; return; end
suffix = regexprep({lib(c).xhm_n}, '^[^:]*', '');
if any(strcmp(suffix, ':h')) && ~isempty(cell)
    if abs(cell(6) - 120) < 1e-3, k = c(strcmp(suffix, ':h')); else, k = c(strcmp(suffix, ':r')); end
    return
end
if any(strcmp(suffix, ':2'))
    k = c(strcmp(suffix, ':2'));
    warning('spacegroup_ops:origin', ['%s has two origin choices and the CIF does not say which; ' ...
        'origin choice 2 is used.'], lib(k).xhm);
    return
end
k = c(1);
end

function op = combine(symop, cenop)
% symop + centring translation, written as text
[W, w] = parse_op(symop);
[~, t] = parse_op(cenop);
w = w + t;
ax = 'xyz';
parts = cell(1, 3);
for r = 1:3
    s = '';
    for c = 1:3
        if W(r, c) == 1, s = [s '+' ax(c)]; elseif W(r, c) == -1, s = [s '-' ax(c)]; end %#ok<AGROW>
    end
    f = mod(w(r), 1);
    if f > 1e-9
        [nmr, dnm] = rat(f);
        s = sprintf('%s+%d/%d', s, nmr, dnm);
    end
    if startsWith(s, '+'), s = s(2:end); end
    parts{r} = s;
end
op = strjoin(parts, ',');
end

function [W, w] = parse_op(op)
parts = strsplit(lower(strrep(op, ' ', '')), ',');
W = zeros(3); w = zeros(3, 1);
for r = 1:3
    terms = strsplit(strrep(parts{r}, '-', '+-'), '+');
    for k = 1:numel(terms)
        tm = terms{k};
        if isempty(tm), continue; end
        ax = find(ismember('xyz', tm), 1);
        if isempty(ax)
            p = strsplit(tm, '/');
            if numel(p) == 2, w(r) = w(r) + str2double(p{1}) / str2double(p{2});
            else, w(r) = w(r) + str2double(tm); end
        else
            c = regexprep(tm, '[xyz]', '');
            if isempty(c), cv = 1; elseif strcmp(c, '-'), cv = -1; else, cv = str2double(c); end
            W(r, ax) = W(r, ax) + cv;
        end
    end
end
end
