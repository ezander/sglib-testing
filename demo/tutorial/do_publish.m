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

% replace preamble until begin{document}
str = regexprep(str, '.*\\begin{document}[\n ]*', '');
str = ['\begin{document}' endl str];
str = ['\input{preamble}' endl str];
str = ['\documentclass[12pt,a4paper]{scrartcl}' endl str];

% replace contents section with latex table of contents
str = regexprep(str, '\\subsection.*{Contents}.*?\\end{itemize}', '\\tableofcontents');

% replace footer section
str = regexprep(str, '[\n ]*\\end{document}.*', '');
str = [str endl '\end{document}' endl];

% replace starred sectioning commands with unstarred ones
str = regexprep(str, '\\section\*', '\\section');
str = regexprep(str, '\\section\*', '\\section');
str = regexprep(str, '\\subsection\*', '\\subsection');

% make graphics centered
str = regexprep(str, '(\\includegraphics.*?})', '\\begin{center}$1\\end{center}');


writetextfile(filename, str);


function modify_stylesheet(ml_stylesheet, stylesheet);
str=readtextfile(ml_stylesheet);
str = regexprep(str, '(match="mcode".*?){verbatim}(.*?){verbatim}', '$1{lstlisting}[style=all,style=mcode]$2{lstlisting}');
str = regexprep(str, '(match="mcodeoutput".*?){verbatim}(.*?){verbatim}', '$1{lstlisting}[style=all,style=mcodeoutput]$2{lstlisting}');
str = regexprep(str, '({lstlisting}\[.*?\])(<)', sprintf('$1\n$2'));
str = regexprep(str, '{lightgray}', '{mcodeoutput}');
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
