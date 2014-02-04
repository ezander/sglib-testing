function do_publish

clear
publishing_defaults

stylesheet = make_stylesheet;

publish_to_latex('demo_simplefem', true, 'imageFormat', 'png', ...
    'stylesheet', stylesheet, 'latex_filter', @tex_modify);



function stylesheet = make_stylesheet
xsl_dir = fileparts(which('publish'));
curr_dir = fileparts(mfilename);
ml_stylesheet = fullfile(xsl_dir,'private','mxdom2latex.xsl');
stylesheet = fullfile(curr_dir,'mxdom2latex.xsl');
modify_stylesheet(ml_stylesheet, stylesheet);



function tex_modify(filename)
endl = sprintf('\n');
str=readtextfile(filename);
str = regexprep(str, '.*\\begin{document}[\n ]*', '');
str = regexprep(str, '[\n ]*\\end{document}.*', '');


str = regexprep(str, '\\section\*', '\\section');
str = regexprep(str, '\\section\*', '\\section');
str = regexprep(str, '\\subsection\*', '\\subsection');

str = ['\begin{document}' endl str];
str = ['\input{preamble}' endl str];
str = ['\documentclass[12pt,a4paper]{scrartcl}' endl str];

str = [str endl '\end{document}' endl];

str = regexprep(str, '\\subsection.*{Contents}.*?\\end{itemize}', '\\tableofcontents');

% clc
% [str(1:1000) str(end-200:end)]


writetextfile(filename, str);


function modify_stylesheet(ml_stylesheet, stylesheet);
str=readtextfile(ml_stylesheet);
% <xsl:template match="mcode">\begin{verbatim}
% <xsl:value-of select="."/>
% \end{verbatim}
% </xsl:template>
str = regexprep(str, '(match="mcode".*?){verbatim}(.*?){verbatim}', '$1{lstlisting}$2{lstlisting}');

writetextfile(stylesheet, str);

function str=readtextfile(filename)
fid = fopen(filename, 'r');
if fid==-1
    error('sglib:publish_to_latex', 'Cannot open file (for reading): %s', filename);
end
[s, pos]=textscan(fid,  '%s', 'whitespace', '');
str=s{1}{1};
fclose(fid);

function writetextfile(filename, str)
fid = fopen(filename, 'w');
if fid==-1
    error('sglib:publish_to_latex', 'Cannot open file (for writing): %s', filename);
end
fprintf(fid, '%s', str);
fclose(fid);
