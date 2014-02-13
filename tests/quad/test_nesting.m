function test_nesting

show_nestings(true)

function show_nestings(nested)
rules = get_rules();
for i=1:size(rules,1)
    rule_func = rules{i,1};
    nmin = rules{i,2};
    if ~nested
        nmax = rules{i,3};
    else
        nmax = rules{i,4};
    end
    test_rule_nesting(rule_func, nmin, nmax, nested);
end

function rules = get_rules
rules = {...
    funcreate(@newton_cotes_rule, @funarg, 'open', true), 1, 20, 33; ...
    funcreate(@newton_cotes_rule, @funarg, 'open', false), 1, 20, 33; ...
    funcreate(@gauss_hermite_rule, @funarg), 1, 20, 50; ...
    funcreate(@gauss_legendre_rule, @funarg), 1, 9, 9; ...
    funcreate(@clenshaw_curtis_rule, @funarg, 'mode', 'cc'), 1, 20, 50; ...
    funcreate(@clenshaw_curtis_rule, @funarg, 'mode', 'cc'), 2, 20, 50; ...
    funcreate(@clenshaw_curtis_rule, @funarg, 'mode', 'fejer1'), 1, 27, 82; ...
    funcreate(@clenshaw_curtis_rule, @funarg, 'mode', 'fejer2'), 1, 20, 64;
    };

function test_rule_nesting(rule_func, nmin, nmax, nested)
underline(sprintf('Rule function: %s', disp_func(rule_func)), 'newlines', 1);
n1 = nmin;
while n1<nmax
    [x1, w1] = funcall(rule_func, n1);
    for n2=n1+1:nmax
        [x2, w2] = funcall(rule_func, n2);
        if all(ismember(x1, x2))
            fprintf('%d subsset of %d\n', n1, n2);
            if nested
                n1 = n2-1;
            end
            break;
        end
    end
    n1 = n1+1;
end

        
        
    