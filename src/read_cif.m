function cif = read_cif(file)
%READ_CIF  Read the crystal structure from a CIF file.
%
%   cif = read_cif(file)
%
%   Reads the first data block of a Crystallographic Information File that
%   contains a cell (_cell_length_a; e.g. the structure block after a
%   'data_global' block with publication data) and returns
%     cif.name        _chemical_name_common / _chemical_formula_sum / file name
%     cif.cell        [a b c alpha beta gamma] (A, deg)
%     cif.space_group Hermann-Mauguin symbol (if given)
%     cif.sg_number   International Tables number (if given)
%     cif.ops         symmetry operations as text, e.g. {'x,y,z'; '-x+1/2,y,-z'}
%     cif.W, cif.w    the same as rotation matrices W (3x3xN) and
%                     translations w (3xN):  x' = W x + w  (fractional)
%     cif.atoms       table: label, element, x, y, z, occupancy, B (A^2)
%                     (B from _atom_site_B_iso_or_equiv, or 8 pi^2 U_iso;
%                      0 if not given); atoms without coordinates are skipped
%     cif.hall        Hall symbol (if given)
%     cif.symmetry_source  where the operations come from
%
%   Supported: _space_group_symop_operation_xyz or _symmetry_equiv_pos_as_xyz,
%   numbers with standard uncertainties (e.g. 3.1599(2)), quoted values,
%   ;-delimited text fields, element symbols with charges (O2-). If the file
%   lists no symmetry operations, they are taken from the Hall or
%   Hermann-Mauguin symbol or the space-group number (SPACEGROUP_OPS).

txt = fileread(file);
lines = regexp(txt, '\r?\n', 'split');

% --- choose the data block: the first one that contains a cell
starts = find(startsWith(lower(strtrim(lines)), 'data_'));
if isempty(starts), starts = 1; end
ends = [starts(2:end) - 1, numel(lines)];
block = 1;
for b = 1:numel(starts)
    if any(startsWith(lower(strtrim(lines(starts(b):ends(b)))), '_cell_length_a'))
        block = b; break
    end
end
lines = lines(starts(block):ends(block));

% --- tokenise: drop comments and text fields, keep track of loops
items = {};          % {tag, value} pairs outside loops
loops = {};          % each: struct('tags', {..}, 'values', {..})
in_text = false;
i = 0;
while i < numel(lines)
    i = i + 1;
    L = lines{i};
    if in_text                                   % inside ;...; text field
        if startsWith(L, ';'), in_text = false; end
        continue
    end
    if startsWith(L, ';'), in_text = true; continue; end
    L = strtrim(strip_comment(L));
    if isempty(L), continue; end
    if startsWith(lower(L), 'data_')
        continue                                 % (only the chosen block is left)
    end
    if strcmpi(L, 'loop_')
        tags = {};
        while i < numel(lines)                   % loop header
            Ln = strtrim(strip_comment(lines{i + 1}));
            if startsWith(Ln, '_'), tags{end+1} = lower(strtok(Ln)); i = i + 1; %#ok<AGROW>
            elseif isempty(Ln), i = i + 1;
            else, break; end
        end
        vals = {};
        while i < numel(lines)                   % loop body
            Ln = lines{i + 1};
            if startsWith(Ln, ';')               % text field inside a loop: skip
                i = i + 1;
                while i < numel(lines) && ~startsWith(lines{i + 1}, ';'), i = i + 1; end
                i = i + 1; vals{end+1} = ''; continue; %#ok<AGROW>
            end
            Ls = strtrim(strip_comment(Ln));
            if isempty(Ls), i = i + 1; continue; end
            if startsWith(Ls, '_') || strcmpi(Ls, 'loop_') || startsWith(lower(Ls), 'data_'), break; end
            vals = [vals, tokens(Ls)]; %#ok<AGROW>
            i = i + 1;
        end
        loops{end+1} = struct('tags', {tags}, 'values', {vals}); %#ok<AGROW>
        continue
    end
    if startsWith(L, '_')
        [tag, rest] = strtok(L);
        rest = strtrim(rest);
        if isempty(rest)                         % value on the next line
            while i < numel(lines)
                Ln = lines{i + 1};
                if startsWith(Ln, ';')           % text field
                    i = i + 1; buf = {strtrim(Ln(2:end))};
                    while i < numel(lines) && ~startsWith(lines{i + 1}, ';')
                        i = i + 1; buf{end+1} = strtrim(lines{i}); %#ok<AGROW>
                    end
                    i = i + 1; rest = strtrim(strjoin(buf, ' ')); break
                end
                Ls = strtrim(strip_comment(Ln));
                i = i + 1;
                if ~isempty(Ls), rest = Ls; break; end
            end
        end
        t = tokens(rest);
        if isempty(t), t = {''}; end
        items(end+1, :) = {lower(tag), strjoin(t, ' ')}; %#ok<AGROW>
    end
end

% --- cell
get = @(tag) item_value(items, tag);
cif.cell = [num(get('_cell_length_a')), num(get('_cell_length_b')), num(get('_cell_length_c')), ...
            num(get('_cell_angle_alpha')), num(get('_cell_angle_beta')), num(get('_cell_angle_gamma'))];
if any(isnan(cif.cell(1:3)))
    error('read_cif:cell', '%s: cell lengths not found.', file);
end
cif.cell(isnan(cif.cell)) = 90;
cif.space_group = first_nonempty({get('_space_group_name_h-m_alt'), get('_symmetry_space_group_name_h-m')});
cif.sg_number = num(first_nonempty({get('_space_group_it_number'), get('_symmetry_int_tables_number')}));
[~, fname] = fileparts(file);
cif.name = first_nonempty({get('_chemical_name_common'), get('_chemical_name_mineral'), ...
                           get('_chemical_formula_sum'), fname});

% --- symmetry operations
ops = loop_column(loops, {'_space_group_symop_operation_xyz', '_symmetry_equiv_pos_as_xyz'});
cif.hall = first_nonempty({get('_space_group_name_hall'), get('_symmetry_space_group_name_hall')});
if ~isempty(ops)
    cif.symmetry_source = 'operations listed in the CIF';
else
    % no operations in the file: take them from the space-group symbol
    if isempty(cif.hall) && isempty(cif.space_group) && isnan(cif.sg_number)
        error('read_cif:noSymmetry', '%s: neither symmetry operations nor a space-group symbol found.', file);
    end
    [ops, sg] = spacegroup_ops('hall', cif.hall, 'hm', cif.space_group, 'number', cif.sg_number, 'cell', cif.cell);
    cif.symmetry_source = sprintf('from the space-group %s: %s (No. %d)', sg.matched_by, sg.xhm, sg.number);
    if isempty(cif.space_group), cif.space_group = sg.xhm; end
    fprintf('  %s: no symmetry operations in the file, %s\n', file, cif.symmetry_source);
end
cif.ops = ops(:);
n = numel(ops);
cif.W = zeros(3, 3, n);
cif.w = zeros(3, n);
for k = 1:n
    [cif.W(:,:,k), cif.w(:,k)] = parse_symop(ops{k});
end

% --- atoms
lab = loop_column(loops, {'_atom_site_label'});
if isempty(lab), error('read_cif:atoms', '%s: no _atom_site_label loop.', file); end
x   = num(loop_column(loops, {'_atom_site_fract_x'}));
y   = num(loop_column(loops, {'_atom_site_fract_y'}));
z   = num(loop_column(loops, {'_atom_site_fract_z'}));
typ = loop_column(loops, {'_atom_site_type_symbol'});
occ = num(loop_column(loops, {'_atom_site_occupancy'}));
Bis = num(loop_column(loops, {'_atom_site_b_iso_or_equiv'}));
Uis = num(loop_column(loops, {'_atom_site_u_iso_or_equiv'}));
na = numel(lab);
if isempty(typ), typ = lab; end
if numel(occ) ~= na, occ = ones(na, 1); end
occ(isnan(occ)) = 1;
B = zeros(na, 1);
if numel(Bis) == na, B = Bis; elseif numel(Uis) == na, B = 8 * pi^2 * Uis; end
B(isnan(B)) = 0;
element = cellfun(@element_symbol, typ(:), 'UniformOutput', false);
cif.atoms = table(lab(:), element, x(:), y(:), z(:), occ(:), B(:), ...
    'VariableNames', {'label', 'element', 'x', 'y', 'z', 'occupancy', 'B'});
bad = isnan(cif.atoms.x) | isnan(cif.atoms.y) | isnan(cif.atoms.z);
if any(bad)
    warning('read_cif:noCoordinates', '%s: atoms without coordinates skipped: %s', ...
        file, strjoin(cif.atoms.label(bad)', ', '));
    cif.atoms(bad, :) = [];
end
end

%% ------------------------------------------------------------------------
function L = strip_comment(L)
% remove a # comment that is not inside quotes
q = false; qc = '';
for k = 1:numel(L)
    c = L(k);
    if ~q && (c == '''' || c == '"') && (k == 1 || isspace(L(k-1)))
        q = true; qc = c;
    elseif q && c == qc && (k == numel(L) || isspace(L(k+1)))
        q = false;
    elseif ~q && c == '#'
        L = L(1:k-1); return
    end
end
end

function t = tokens(L)
% split a line into values, respecting '...' and "..."
t = {};
k = 1; n = numel(L);
while k <= n
    while k <= n && isspace(L(k)), k = k + 1; end
    if k > n, break; end
    if L(k) == '''' || L(k) == '"'
        qc = L(k); e = k + 1;
        while e <= n && ~(L(e) == qc && (e == n || isspace(L(e+1)))), e = e + 1; end
        t{end+1} = L(k+1:min(e-1, n)); %#ok<AGROW>
        k = e + 1;
    else
        e = k;
        while e <= n && ~isspace(L(e)), e = e + 1; end
        t{end+1} = L(k:e-1); %#ok<AGROW>
        k = e;
    end
end
end

function v = item_value(items, tag)
v = '';
if isempty(items), return; end
k = find(strcmp(items(:,1), tag), 1);
if ~isempty(k), v = items{k, 2}; end
end

function col = loop_column(loops, tags)
col = {};
for k = 1:numel(loops)
    c = find(ismember(loops{k}.tags, tags), 1);
    if ~isempty(c)
        nt = numel(loops{k}.tags);
        vals = loops{k}.values;
        nrow = floor(numel(vals) / nt);
        col = vals(c:nt:nrow * nt)';
        return
    end
end
end

function x = num(s)
% numbers like 3.1599(2), '.', '?' -> NaN
if ischar(s) || isstring(s), s = {char(s)}; end
x = nan(numel(s), 1);
for k = 1:numel(s)
    v = regexprep(strtrim(s{k}), '\(.*\)$', '');
    if ~isempty(v) && ~any(strcmp(v, {'.', '?'}))
        x(k) = str2double(v);
    end
end
end

function v = first_nonempty(c)
v = '';
for k = 1:numel(c)
    if ~isempty(c{k}) && ~any(strcmp(c{k}, {'?', '.'})), v = c{k}; return; end
end
end

function el = element_symbol(s)
% 'Ti4+' -> 'Ti', 'TI1' -> 'Ti', 'Nb' -> 'Nb'
t = regexp(s, '^[A-Za-z]{1,2}', 'match', 'once');
if isempty(t), el = ''; return; end
el = [upper(t(1)), lower(t(2:end))];
if numel(el) == 2 && isempty(element_Z(el)), el = el(1); end   % e.g. 'O1' read as 'O'
end

function [W, w] = parse_symop(op)
% 'x-y+1/2, -z, 1/4+y' -> W (3x3), w (3x1)
parts = strsplit(lower(strrep(op, ' ', '')), ',');
if numel(parts) ~= 3, error('read_cif:symop', 'Cannot read symmetry operation "%s".', op); end
W = zeros(3); w = zeros(3, 1);
for r = 1:3
    e = strrep(parts{r}, '-', '+-');
    terms = strsplit(e, '+');
    for k = 1:numel(terms)
        tm = terms{k};
        if isempty(tm), continue; end
        ax = find(ismember('xyz', tm), 1);
        if isempty(ax)
            w(r) = w(r) + frac(tm);
        else
            c = regexprep(tm, '[xyz]', '');
            if isempty(c) || strcmp(c, '+'), cv = 1;
            elseif strcmp(c, '-'), cv = -1;
            else, cv = frac(strrep(c, '*', ''));
            end
            W(r, ax) = W(r, ax) + cv;
        end
    end
end
end

function v = frac(s)
if contains(s, '/')
    p = strsplit(s, '/');
    v = str2double(p{1}) / str2double(p{2});
else
    v = str2double(s);
end
end
