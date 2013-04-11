#!/usr/bin/env python
#
# encoding:utf-8
#
# Fortrandoc
#
# Python-based parser for FEAST documentation.
#
# (c) 2006,2007 Thomas Rohkaemper <thomas.rohkaemper@uni-dortmund.de>
#
# Last change: November, 12th 2007
#
"""
Python-based parser for FEAST documentation.

usage: fortrandoc.py <output types> [options] <files>

options:
  -h, --help            show this help message and exit
  -d OUTPUTDIR, --output-dir=OUTPUTDIR
                        write all output files to the specified directory
  -s, --strip-extension
                        strip extension from the filename
  -e, --exit-on-error   exit on first error

  output types:
    -t, --tex           generate LaTeX output
    -H, --html          generate HTML output
    -x, --xml           generate XML output
    -g HEADERFILE, --generate-latex-definitions=HEADERFILE
                        generate LaTeX definitions

(c) 2006,2007 Thomas Rohkaemper
"""

# Workaround for different versions of Python
try:
    from re import Scanner
except ImportError:
    from sre import Scanner
    
import re

def error(error_code, msg):
    "Print error message and return if option --force was not set."

    print msg
    global options
    if not options.force:
        from sys import exit
        exit(error_code)

# =================================================================

class HTMLExporter:
    "Exports to a given HTML file."

    # Modify this to change the way the doc file is displayed as HTML.
    CSS = """<style type="text/css">
* {
    font-family: sans-serif;
    font-size: 10pt;
}

body {
    margin-left: 50px;
    margin-right: 50px;
}

.keyword {
    font-family: monospace;
}

.type {
    font-family: monospace;
    font-style: italic;
}

.variable, .error_id {
    font-family: monospace;
    font-weight: bold;
}

.description {
    margin-left: 25px;
    padding: 5px;
    border: 1px dotted lightgray;
    #border: 1px dotted black;
    background: #ffa;
    whitespace: pre;
    text-align: justify;
}

.subroutine {
    font-family: monospace;
    font-weight: bold;
}

.function {
    font-family: monospace;
    font-weight: bold;
}

.parameters {
    display: inline;
}

.subroutine_definition, .function_definition, .type_definition {
    border: 1px dotted black;
    margin: 12px;
    margin-bottom: 20px;
    padding: 5px;
}

.variable_definition {
    margin: 12px;
    margin-bottom: 18px;
    padding: 5px;
    border: 1px dotted lightgray;
}

.subroutinelist, .functionlist {
    list-style-type: square;
}

.go-to-top {
    display: block;
    margin: 0px;
    text-align: right;
}

.go-to-top a {
    font-size: xx-small;
}

h1 {
    font-size: 18pt;
}

h2, h3 {
    border-bottom: 1px dotted silver;
    #color: darkgray;
}
</style>
"""

    def export(self, module, f):
        "Write the given module as HTML to f."

        from os import path

        header =[
            '<html>\n<title>Documentation of Module %s (%s)</title>\n' % (module.name, path.basename(module.filename)),
            self.CSS,
            '<body>',
            '<h1>Module %s (%s)</h1>' % (module.name, path.basename(module.filename)),
            '<h2>Purpose</h2>\n<div class="description">\n%s\n</div>' % module.purpose
        ]
        f.writelines(header)

        f.writelines(['<a name="typedefs"><h2>Type definitions</h2></a>\n'])
        for typedef in module.types:
            self.write_type(typedef, f)

        f.writelines(['</div>\n<a name="constants"><h2>Constants</h2></a>\n'])
        for var in module.constants:
            self.write_variable(var, f)

        f.writelines(['</div>\n<a name="globalvars"><h2>Global variables</h2></a>\n'])
        for var in module.global_variables:
            self.write_variable(var, f)

        f.writelines(['</div>\n<a name="privatevars"><h2>Private variables</h2></a>\n'])
        for var in module.private_variables:
            self.write_variable(var, f)

        f.writelines(['<a name="subroutinelist"><h2>Subroutines</h2></a>\n<ul class="subroutinelist">\n'])
        for sub in module.subroutines:
            f.writelines(['<li><a href="#subroutine_%s">%s</a></li> ' % (sub.name, sub.name)])

        f.writelines(['</ul>\n<a name="functionlist"><h2>Functions</h2></a>\n<ul class="functionlist">\n'])
        for fn in module.functions:
            f.writelines(['<li><a href="#function_%s">%s</a></li> ' % (fn.name, fn.name)])
        f.writelines(['</ul>\n<h2>Descriptions</h2>\n'])
        for sub in module.subroutines:
            self.write_subroutine(sub, f)
            self.write_navigation(f)
        for fn in module.functions:
            self.write_function(fn, f)
            self.write_navigation(f)

        if module.remarks:
            f.writelines(['<h2>Remarks</h2>\n<div class="description">\n%s\n</div>\n' % module.remarks])
        footer = [
            '</body>\n</html>'
        ]
        f.writelines(footer)

    def write_navigation(self, f):
        f.writelines([
            '<div class="go-to-top">\n',
            '<a href="#subroutinelist">&rarr; List of subroutines</a> ',
            '<a href="#functionlist">&rarr; List of functions</a> ',
            '<a href="#">&uarr; Top</a>\n',
            '</div>\n'
        ])

    def write_type(self, typedef, f):
        lines = [
            '<div class="type_definition">',
            '<a name="type_%s">' % typedef.name,
            '<span class="keyword">type</span> <span class="type">%s</span>' % typedef.name,
            '</a>\n',
            '<div class="description">\n%s\n</div>\n' % typedef.description
        ]
        f.writelines(lines)
        for var in typedef.variables:
            self.write_variable(var, f, "type_%s_" % typedef.name)

        f.writelines(['</div>\n'])

    def write_subroutine(self, sub, f):
        lines = [
            '<div class="subroutine_definition">\n',
            '<a name="subroutine_%s">' % sub.name,
            '<h3><span class="keyword">subroutine</span> <span class="subroutine">%s</span>(<div class="parameters">' % sub.name
        ]
        f.writelines(lines)
        subs = []
        for var in sub.parameters:
            names = map(lambda var: var.name, sub.inputvars+sub.inoutputvars+sub.outputvars)
            if var in names:
                subs.append('<a href="#subroutine_%s_var_%s"><span class="variable">%s</span></a>' % (sub.name, var, var))
            else:
                subs.append('<span class="variable">%s</span>' % var)
        f.writelines([', '.join(subs)])
        f.writelines(['</div>)</h3></a>\n'])
        f.writelines(['<div class="description">\n%s\n</div>' % sub.description])

        if sub.references:
            f.writelines(['<h4>References</h4>\n'])
            for ref in sub.references:
                f.writelines(['<a href="#var_%s">%s</a>\n' % (ref, ref)])
            f.writelines(['\n\n'])

        if sub.inputvars:
            f.writelines(['<h4>Input variables</h4>'])
            for var in sub.inputvars:
                self.write_variable(var, f, 'subroutine_%s_' % sub.name)

        if sub.inoutputvars:
            f.writelines(['<h4>Input/Output variables</h4>'])
            for var in sub.inoutputvars:
                self.write_variable(var, f, 'subroutine_%s_' % sub.name)

        if sub.outputvars:
            f.writelines(['<h4>Output variables</h4>'])
            for var in sub.outputvars:
                self.write_variable(var, f, 'subroutine_%s_' % sub.name)

        f.writelines(['<h4>Error codes</h4>\n'])
        if sub.errors:
            f.writelines(['<table class="errortable">\n'])
            for (error, desc) in sub.errors:
                f.writelines(['<tr> <td class="error_id">%s</td> <td class="error_description">%s</td> </tr>\n' % (error, desc)])
            f.writelines(['</table>\n'])
        else:
            f.writelines(['<div class="error">No error codes defined.</div>\n'])

        if sub.remarks:
            f.writelines(['<h4>Remarks</h4>\n<div class="description">\n%s\n</div>\n' % sub.remarks])

        f.writelines(['</div>\n'])

    def write_function(self, fn, f):
        lines = [
            '<div class="function_definition">\n',
            '<a name="function_%s">' % fn.name,
            '<h3><span class="type">%s</span> <span class="keyword">function</span> <span class="subroutine">%s</span>(<div class="parameters">' % (fn.return_value, fn.name)
        ]
        f.writelines(lines)
        subs = []
        for var in fn.parameters:
            names = map(lambda var: var.name, fn.inputvars+fn.inoutputvars+fn.outputvars)
            if var in names:
                subs.append('<a href="#function_%s_var_%s"><span class="variable">%s</span></a>' % (fn.name, var, var))
            else:
                subs.append('<span class="variable">%s</span>' % var)
        f.writelines([', '.join(subs)])
        f.writelines(['</div>)</h3></a>\n'])
        f.writelines(['<div class="description">\n%s\n</div>' % fn.description])

        f.writelines([
            '<h4>Result</h4>\n',
            '<div class="result">\n%s\n</div>\n' % fn.result
        ])

        if fn.references:
            f.writelines(['<h4>References</h4>\n'])
            for ref in fn.references:
                f.writelines(['<a href="#var_%s">%s</a>\n' % (ref, ref)])
            f.writelines(['\n\n'])

        if fn.inputvars:
            f.writelines(['<h4>Input variables</h4>\n'])
            for var in fn.inputvars:
                self.write_variable(var, f, 'function_%s_' % fn.name)
                
        if fn.inoutputvars:
            f.writelines(['<h4>Input/Output variables</h4>\n'])
            for var in fn.inoutputvars:
                self.write_variable(var, f, 'function_%s_' % fn.name)

        if fn.outputvars:
            f.writelines(['<h4>Output variables</h4>\n'])
            for var in fn.outputvars:
                self.write_variable(var, f, 'function_%s_' % fn.name)

        f.writelines(['<h4>Error codes</h4>\n'])
        if fn.errors:
            f.writelines(['<table class="errortable">\n'])
            for (error, desc) in fn.errors:
                f.writelines(['<tr> <td class="error_id">%s</td> <td class="error_description">%s</td> </tr>\n' % (error, desc)])
            f.writelines(['</table>\n'])
        else:
            f.writelines(['<div class="error">No error codes defined.</div>\n'])

        if fn.remarks:
            f.writelines(['<h4>Remarks</h4>\n<div class="description">\n%s\n</div>\n' % fn.remarks])

        f.writelines(['</div>\n'])

    def write_variable(self, var, f, context=''):
        if not var:
            return

        if isinstance(var, Group): # it's a group!
            f.writelines([
                '<div class="group"><p>%s</p>\n' % group.description
            ])
            for v in var.variables:
                self.write_variable(v, f, context)
            f.writelines([
                '</div>\n' % group.description
            ])
        else: # it's a variable!
            lines = [
                '<div class="variable_definition">',
                '<div class="variable_name">\n',
                '<a name="%svar_%s">' % (context, var.name),
                '<span class="type">%s</span> ' % var.vartype,
                '<span class="variable">%s</span></a>' % var.name
            ]
            f.writelines(lines)
            if var.dimension:
                f.writelines([', %s' % var.dimension])
            if var.optional:
                f.writelines([' (optional)'])
            lines = [
                '\n</div>\n'
            ]
            f.writelines(lines)
            if var.description:
                lines = [
                    '<div class="description">',
                    '%s' % var.description,
                    '</div>\n'
                ]
                f.writelines(lines)
            f.writelines(['</div>\n'])

# =================================================================

class XMLExporter:
    "Export as XML to a given file."

    def export(self, module, f):
        from os import path

        f.writelines(['<module name="%s" filename="%s">\n    <purpose>%s</purpose>\n' % (module.name, path.basename(module.filename), module.purpose)])
        for typedef in module.types:
            self.write_type(typedef, f, 4)
        f.writelines(['    <constant>\n'])
        for var in module.constants:
            self.write_variable(var, f, 8)
        f.writelines(['    </constant>'])
        f.writelines(['    <global>\n'])
        for var in module.global_variables:
            self.write_variable(var, f, 8)
        f.writelines(['    </global>'])
        f.writelines(['    <private>\n'])
        for var in module.private_variables:
            self.write_variable(var, f, 8)
        f.writelines(['    </private>'])
        for sub in module.subroutines:
            self.write_subroutine(sub, f, 4)
        for fn in module.functions:
            self.write_function(fn, f, 4)
        if module.remarks:
            f.writelines(['    <remarks>%s</remarks>' % module.remarks])
        f.writelines(['</module>\n'])

    def write_variable(self, var, f, indent=0):
        prefix = ' ' * indent
        if isinstance(var, Group):
            f.writelines([
                prefix + '<group description="%s">\n' % var.description
            ])
            for v in var.variables:
                self.write_variable(v, f, indent+4)
            f.writelines([
                prefix + '</group>\n'
            ])
        else:
            f.writelines([prefix + '<variable name="%s" type="%s">%s</variable>\n' % (var.name, var.vartype, var.description or '')])

    def write_subroutine(self, sub, f, indent=0):
        prefix = ' ' * indent
        f.writelines([prefix + '<subroutine name="%s">\n' % sub.name])

        if sub.references:
            for ref in sub.references:
                f.writelines([prefix + '    <references>%s</references>\n' % ref])

        if sub.inputvars:
            f.writelines([prefix + '    <input>\n'])
            for inp in sub.inputvars:
                self.write_variable(inp, f, indent+8)
            f.writelines([prefix + '    </input>\n'])

        if sub.inoutputvars:
            f.writelines([prefix + '    <inoutput>\n'])
            for inout in sub.inoutputvars:
                self.write_variable(inout, f, indent+8)
            f.writelines([prefix + '    </inoutput>\n'])

        if sub.outputvars:
            f.writelines([prefix + '    <output>\n'])
            for out in sub.outputvars:
                self.write_variable(out, f, indent+8)
            f.writelines([prefix + '    </output>\n'])

        f.writelines([prefix + "    <errors>\n"])
        if sub.errors:
            for error, desc in sub.errors:
                f.writelines([prefix + '        <error code="%s">%s</error>\n' % (error, desc)])
        f.writelines([prefix + '    </errors>\n'])

        if sub.description:
            f.writelines([prefix + sub.description + "\n"])
        if sub.remarks:
            f.writelines(['    <remarks>%s</remarks>' % sub.remarks])
        f.writelines([prefix + '</subroutine>\n'])

    def write_function(self, fn, f, indent=0):
        prefix = ' ' * indent
        f.writelines([prefix + '<function name="%s" returns="%s">\n' % (fn.name, fn.return_value)])

        if fn.references:
            for ref in fn.references:
                f.writelines([prefix + '    <references>%s</references>\n' % ref])

        if fn.inputvars:
            f.writelines([prefix + '    <input>\n'])
            for inp in fn.inputvars:
                self.write_variable(inp, f, indent+8)
            f.writelines([prefix + '    </input>\n'])

        if fn.inoutputvars:
            f.writelines([prefix + '    <inoutput>\n'])
            for inout in fn.inoutputvars:
                self.write_variable(inout, f, indent+8)
            f.writelines([prefix + '    </inoutput>\n'])

        if fn.outputvars:
            f.writelines([prefix + '    <output>\n'])
            for out in fn.outputvars:
                self.write_variable(out, f, indent+8)
            f.writelines([prefix + '    </output>\n'])

        f.writelines([prefix + "    <result>%s</result>\n" % fn.result])
        f.writelines([prefix + "    <errors>\n"])
        if fn.errors:
            for error, desc in fn.errors:
                f.writelines([prefix + '        <error code="%s">%s</error>\n' % (error, desc)])
        f.writelines([prefix + '    </errors>\n'])
        f.writelines([prefix + fn.description + '\n'])
        if fn.remarks:
            f.writelines(['    <remarks>%s</remarks>' % fn.remarks])
        f.writelines([prefix + '</function>\n'])

    def write_type(self, typedef, f, indent=0):
        prefix = ' '*indent
        f.writelines([prefix + '<type name="%s">\n' % typedef.name])
        for var in typedef.variables:
            self.write_variable(var, f, indent+4)
        f.writelines([prefix + typedef.description])
        f.writelines([prefix + '</type>\n'])

# =================================================================

class LaTeXExporter:
    "Export as LaTeX to the given file."

    # Modify this to change the way the docs are displayed.
    HEADER = r"""
\documentclass[a4paper,draft]{article}
"""
    DEFINITIONS = r"""
\newcommand{\module}[2]{\chapter{Module #1 (#2)}} % name, filename
\newcommand{\modulepurpose}{\section{Purpose}}
\newcommand{\function}[1]{\textbf{function} #1} % name
\newcommand{\subroutine}[1]{\textbf{subroutine} #1} % name
\newcommand{\variable}[2]{\texttt{#1} \textbf{#2}} % type, name
\newcommand{\constant}[1]{\texttt{#1}} % name
\newcommand{\interface}[1]{\texttt{interface #1}}
\newcommand{\sectionvariables}[1]{\subsubsection{#1}} % name
\newcommand{\sectionerrors}{\subsubsection{Error Codes}}
\newcommand{\sectiontypes}{\section{Types}}
\newcommand{\sectionconstants}{\section{Constants}}
\newcommand{\sectionpublics}{\section{Public Module Interface}}
\newcommand{\sectionglobals}{\section{Public Variables}}
\newcommand{\sectionprivates}{\section{Private Variables}}
\newcommand{\sectionsubroutines}{\section{Subroutines}}
\newcommand{\sectionfunctions}{\section{Functions}}
\newcommand{\sectioninterfaces}{\section{Interfaces}}
\newcommand{\sectionremarks}{\section{Remarks}}
\newcommand{\result}[1]{\subsubsection*{Result}}
\newcommand{\returnvalue}[1]{\texttt{#1}}
\newenvironment{desc}{\list{}{\rightmargin0pt\leftmargin20pt}\item\relax}{\endlist}
\newenvironment{functiondef}[1]{\begin{description}\item[{}]\relax#1}{\end{description}}
\newenvironment{subroutinedef}[1]{\begin{description}\item[{}]\relax#1}{\end{description}}
\newenvironment{typedef}[1]{\begin{description}\item[{}]\relax Type #1}{\end{description}}
\newenvironment{variabledef}{}{}
\newenvironment{vargroup}[1]{\begin{description}\item[{}]\relax\textbf{#1}\\}{\end{description}}
\clubpenalty = 10000
\widowpenalty = 10000
\displaywidowpenalty = 10000
"""

    BEGIN_DOC = r"""
\begin{document}
% Uncomment the next line for Sans Serif output:
%\sf
"""

    FOOTER = """
\end{document}
"""
    def export(self, module, f):
        # f.writelines(self.HEADER+self.DEFINITIONS+self.BEGIN_DOC)
        from os import path

        f.writelines([
            self.cmd('module', (module.name, path.basename(module.filename))), '\n\n',
            self.cmd('modulepurpose'), '\n',
            self.begin('desc'),
            self.texify(module.purpose), '\n',
            self.end('desc'), '\n'
        ])

        if module.publics:
            f.writelines([self.cmd('sectionpublics'), '\n\n'])
            for var in module.publics:
                self.write_variable_or_group(var, f)

        if module.types:
            f.writelines(['\sectiontypes\n\n'])
            for typedef in module.types:
                self.write_type(typedef, f)

        if module.constants:
            f.writelines([self.cmd('sectionconstants'), '\n\n'])
            for var in module.constants:
                self.write_variable_or_group(var, f)

        if module.global_variables:
            f.writelines([self.cmd('sectionglobals'), '\n\n'])
            for var in module.global_variables:
                self.write_variable_or_group(var, f)

        if module.private_variables:
            f.writelines([self.cmd('sectionprivates'), '\n\n'])
            for var in module.private_variables:
                self.write_variable_or_group(var, f)

        if module.interfaces:
            f.writelines([self.cmd('sectioninterfaces'), '\n\n'])
            for interface, description in module.interfaces:
                self.write_interface(interface, description, f)

        if module.subroutines:
            f.writelines([self.cmd('sectionsubroutines'), '\n\n'])
            for sub in module.subroutines:
                self.write_subroutine(sub, f)

        if module.functions:
            f.writelines([self.cmd('sectionfunctions'), '\n\n'])
            for fn in module.functions:
                self.write_function(fn, f)

        if module.remarks:
            f.writelines([self.cmd('sectionremarks'), '\n\n', self.begin('desc'), self.texify(module.remarks), '\n', self.end('desc'), '\n'])

        # f.writelines(self.FOOTER)

    def write_definitions(self, f):
        f.writelines([self.DEFINITIONS])

    def write_type(self, typedef, f):
        f.writelines([
            self.begin('typedef', typedef.name)
        ])
        if typedef.description:
            f.writelines([
                self.begin('desc'),
                self.texify(typedef.description),
                self.end('desc')
            ])
        for var in typedef.variables:
            self.write_variable_or_group(var, f)
        f.writelines([ self.end('typedef'), '\n\n' ])

    def write_subroutine(self, sub, f):
        f.writelines([
            r'\begin{subroutinedef}{',
            self.cmd('subroutine', sub.name), '(', self.cmd('sloppy')
        ])
        if sub.parameters:
            params = ", ".join(map(lambda s: self.cmd('textbf', s), sub.parameters))
            f.writelines([params])
        f.writelines([')}\n\n'])

        if sub.description:
            f.writelines([self.env('desc', sub.description), '\n\n'])

        if sub.inputvars:
            f.writelines([self.cmd('sectionvariables', 'Input variables'), '\n\n'])
            for var in sub.inputvars:
                self.write_variable_or_group(var, f)

        if sub.inoutputvars:
            f.writelines([self.cmd('sectionvariables', 'Input/Output variables'), '\n\n'])
            for var in sub.inoutputvars:
                self.write_variable_or_group(var, f)

        if sub.outputvars:
            f.writelines([self.cmd('sectionvariables', 'Output variables'), '\n\n'])
            for var in sub.outputvars:
                self.write_variable_or_group(var, f)

        if sub.errors:
            f.writelines([self.cmd('sectionerrors'), '\n\n'])
            f.writelines([self.begin('description'), '\n'])
            for error, desc in sub.errors:
                f.writelines([r'\item[\constant{%s}] %s' % (self.texify(error), self.texify(desc)), '\n'])
            f.writelines([self.end('description'), '\n\n'])
        #else:
        #    f.writelines([self.begin('desc'), 'Error codes missing.\n', self.end('desc')])

        if sub.remarks:
            f.writelines([self.cmd('sectionvariables', 'Remarks'), '\n\n', self.begin('desc'), self.texify(sub.remarks), '\n', self.end('desc'), '\n'])
            
        f.writelines([
            self.end('subroutinedef'), '\n\n',
            self.cmd('vspace*', ['1em']),
            self.cmd('hrule'),
            self.cmd('bigskip'),
            r'\pagebreak[0]', '\n\n'
        ])

    def write_function(self, fn, f):
        f.writelines([
            r'\begin{functiondef}{',
            self.cmd('returnvalue', fn.return_value), ' ',
            self.cmd('function', fn.name), '(',
        ])
        if fn.parameters:
            params = ", ".join(map(lambda s: self.cmd('textbf', s), fn.parameters))
            f.writelines([params])
        f.writelines([')}\n\n'])

        if fn.description:
            f.writelines([self.env('desc', fn.description), '\n\n'])

        if fn.result:
            f.writelines([self.cmd('result'), '\n\n%s\n\n' % self.texify(fn.result)])

        if fn.inputvars:
            f.writelines([self.cmd('sectionvariables', 'Input variables'), '\n\n'])
            for var in fn.inputvars:
                self.write_variable_or_group(var, f)

        if fn.inoutputvars:
            f.writelines([self.cmd('sectionvariables', 'Input/Output variables'), '\n\n'])
            for var in fn.inoutputvars:
                self.write_variable_or_group(var, f)

        if fn.outputvars:
            f.writelines([self.cmd('sectionvariables', 'Output variables'), '\n\n'])
            for var in fn.outputvars:
                self.write_variable_or_group(var, f)

        if fn.errors:
            f.writelines([self.cmd('sectionerrors'), '\n\n'])
            f.writelines([self.begin('description'), '\n'])
            for error, desc in fn.errors:
                f.writelines([r'\item[\constant{%s}] %s' % (self.texify(error), self.texify(desc)), '\n'])
            f.writelines([self.end('description'), '\n\n'])
        #else:
        #    f.writelines([self.begin('desc'), 'Error codes missing.\n', self.end('desc')])

        if fn.remarks:
            f.writelines([self.cmd('sectionvariables', 'Remarks'), '\n\n', self.begin('desc'), self.texify(fn.remarks), '\n', self.end('desc'), '\n'])

        f.writelines([
            self.end('functiondef'), '\n\n',
            self.cmd('vspace*', ['1em']),
            self.cmd('hrule'),
            self.cmd('bigskip'),
            r'\pagebreak[0]', '\n\n'
        ])
                
    def write_variable_or_group(self, vg, f):
        if isinstance(vg, Group):
            f.writelines([self.begin('vargroup', vg.description), '\n'])
            for var in vg.variables:
                self.write_variable_in_list(var, f)
            f.writelines([self.end('vargroup'), '\n'])
        else:
            self.write_variable_in_list(vg, f)

    def write_variable_in_list(self, var, f):
        f.writelines(self.begin('variabledef'))
        if var.dimension:
            f.writelines([
                self.cmd('variable', [var.vartype, var.name]), ', ',
                self.cmd('texttt', var.dimension), '\n'
            ])
        else:
            f.writelines([
                self.cmd('variable', [var.vartype, var.name]), '\n'
            ])

        if var.optional:
            f.writelines([self.cmd('textit', ' optional'), '\n'])

        if var.description:
            f.writelines([ self.env('desc', var.description) ])

        f.writelines([self.end('variabledef'), '\n\n'])

    def write_interface(self, interface, description, f):
        f.writelines([self.cmd('interface', interface), '\n'])
        f.writelines([self.env('desc', description), '\n'])

    TEX = Scanner([
        (r'<tex>.+?</tex>', lambda s, t: t.replace('<tex>','').replace('</tex>', '')),
        (r'<verb>.+?</verb>', lambda s, t: t.replace('<verb>','\\begin{verbatim}').replace('</verb>', '\\end{verbatim}')),
        (r'<code>.+?</code>', lambda s, t: t.replace('<code>','\\begin{verbatim}').replace('</code>', '\\end{verbatim}')),
        (r'<bf>.+?</bf>', lambda s, t: t.replace('<bf>','\\textbf{').replace('</bf>', '}')),
        (r'<ul>.+?</ul>', lambda s, t: t.replace('<ul>','\\underline{').replace('</ul>', '}')),
        (r'<it>.+?</it>', lambda s, t: t.replace('<it>','\\textit{').replace('</it>', '}')),
        ('\\$', '\\$'),
        ('_', '\\_'),
        (r'\\', '\\\\'),
        (r'#', '\\#'),
        (r'%', '\\%'),
        (r'\s->\s', ' $\\to$ '),
        (r' -> ', ' $\\to$ '),
        (r' => ', ' $\\,\\Rightarrow\\,$ '),
        (re.escape(r'^'), '\\symbol{94}'),
        (re.escape(r'&'), '\\&'),
        ('.', lambda s, t: t)
    ], re.DOTALL)
        
    def texify(self, s):
        "Replace some characters with their tex equivalents"
        return ''.join(self.TEX.scan(s)[0])
        #return s.replace('\\', '\\\\').replace('$', '\\$').replace('_', '\\_').replace('#', '\\#').replace('<tex>', '').replace('<\tex>', '')

    def env(self, name, content):
        r"Create a string of the form \begin{name}content\end{name}"
        content = self.texify(content)
        return '\\begin{%s}\n%s\n\end{%s}\n' % (name, content, name)

    def cmd(self, name, content=None):
        r"Create a string of the form \name or \name{content}."
        if content:
            if hasattr(content, '__iter__'):
                mapped = map(self.texify, content)
                return r'\%s{%s}' % (name, "}{".join(mapped))
            content = self.texify(content)
            return r'\%s{%s}' % (name, content)
        return r'\%s' % name

    def begin(self, name, parameters=None):
        r"Create a string of the form \begin{name}."
        if parameters:
            if hasattr(parameters, '__iter__'):
                return "\\begin{%s}{%s}\n" % (name, "}{".join(map(self.texify, parameters)))
            return "\\begin{%s}{%s}\n" % (name, self.texify(parameters))
        return "\\begin{%s}\n" % name

    def end(self, name):
        r"Create a string of the formn \end{name}."
        return "\\end{%s}\n" % name

# =================================================================

class Function:
    "A function or subroutine. Functions without return value are subroutines."

    def __init__(self, name, parameters = [], return_value = None):
        self.name = name
        self.parameters = parameters
        self.return_value = return_value
        self.inputvars = []
        self.inoutputvars = []
        self.outputvars = []
        self.result = None
        self.description = None
        self.references = []
        self.errors = []
        self.remarks = "";

    def is_function(self):
        return self.return_value is not None

    def __str__(self):
        if self.return_value:
            return "%s function %s(%s) ==> %s ! %s" % (self.return_value, self.name, ", ".join(self.parameters), self.result, self.description)
        return "subroutine %s(%s) ==> %s ! %s" % (self.name, ", ".join(self.parameters), self.result, self.description)
    
# =================================================================

class Type:
    "A type."

    def __init__(self, name, description=None):
        self.name = name
        self.description = description
        self.variables = []

    def __str__(self):
        s = "type %s ! %s\n" % (self.name, self.description)
        if self.variables:
            for var in self.variables:
                s += "    " + str(var) + "\n"
        s += "end type %s" % self.name
        return s

# =================================================================

class Variable:
    "A variable."

    def __init__(self, name, vartype, description=None):
        self.name = name
        self.vartype = vartype
        self.description = description
        self.dimension = None
        self.optional = False
        
    def __str__(self):
        if self.description:
            return "%s %s ! %s" % (self.vartype, self.name, self.description)
        return "%s %s" % (self.vartype, self.name)
        
class Group:
    "A variable group."

    def __init__(self, description=None):
        self.description = description
        self.variables = []

    def __str__(self):
        s = "Group: %s\n" % self.description
        if self.variables:
            for var in self.variables:
                s += "    " + str(var) + "\n"
        return s
# =================================================================

class Node(list):
    "Represents an XML node."

    def __init__(self, name=None, text=None):
        self.name = name
        self.text = text
        self.attributes = {}
        self.body = None

    def set_text(self, text):
        self.text = text

    def get_text(self):
        return self.text

    def set_name(self, name):
        self.name = name

    def get_name(self):
        return self.name

    def __repr__(self):
        return "Node %s <%s>" % (self.name, list.__repr__(self))

# =================================================================

class FortranParser:
    "Find all XML tags in a fortran file and parse them."
        
    def __init__(self):
        self.indent = 0
        self.path = ""
        pass

    PATTERN = re.compile(r"\<\s*(?P<word>[^\>\s]+)([^\>]*?)\>(.*?)\</(?P=word)\s*\>", re.DOTALL)
    ATTRIBUTES = re.compile(r'\s*([^=]+?)\s*=\s*"([^"]*?)"')
    BODY = re.compile(r'(.*?)\<\s*(?P<word>[^\>]+)\>.*?\</(?P=word)\>\s', re.DOTALL)
        
    def tag(self, s, t):
        result =  self.PATTERN.search(t)
        if result:
            #print "Examining:", t
            #print "Found tag:", result.group(1)
            self.path = self.path + (' '*self.indent) + result.group(1) + '\n'
            #print "with body:", result.group(3)
            node = Node(result.group(1), result.group(3))
            #body = re.findall(self.BODY, result.group(3))
            #node.body = ''.join(map(lambda l: l[0], body))
            #print 'Body:', node.body.strip()
            attr = re.findall(self.ATTRIBUTES, result.group(2))
            for a in attr:
                node.attributes[a[0].strip()] = a[1]

            self.indent += 3
            subresult = self.parse(result.group(3))
            self.indent -= 3
            if subresult:
                for i in subresult[0]:
                    node.append(i)
            return node

    def unmatched_open(self, s, t):
        print self.path
        error(100, 'No closing tag found for %s.' % t)

    def unmatched_close(self, s, t):
        print self.path
        error(101, 'No opening tag found for %s.' % t)
            
    def parse(self, sourcecode):
        scanner = Scanner([
            (r'\<!--.*?-->', None), # ignore comments
            (r"\<\s*(?P<word>[^\>\s]+)[^\>]*\>.*?\</(?P=word)\s*\>", self.tag),
            #(r"\<(?P<word>[^\>]+)\>.*?\</(?P=word)\>", self.tag),
            (r'\</[a-zA-Z]+([ \t][^\>]*)\>', self.unmatched_close),
            (r'\<[a-zA-Z]+([ \t][^\>]*)\>', self.unmatched_open),
            (".", None)
        ], re.DOTALL)
        return scanner.scan(sourcecode)

# =================================================================

class Module:
    "Module contains all information about a fortran module as well as some functions to extract the information."

    def __init__(self):
        self.name = "Unknown"
        self.filename = None
        self.purpose = "Unknown"
        self.functions = []
        self.subroutines = []
        self.private_variables = []
        self.global_variables = []
        self.constants = []
        self.publics = []
        self.types = []
        self.interfaces = []
        self.remarks = "";

    def interpret(self, nodes):
        "Extract information from the nodes."

        for node in nodes:
            name = node.get_name()
            if name == 'name':
                self.name = node.get_text().strip()
            elif name == 'purpose':
                self.purpose = self.clean_description(node.get_text())
            elif name == 'remark':
                text = self.clean_description(node.get_text())
                self.remarks = text.strip()
            elif name == 'remarks':
                text = self.clean_description(node.get_text())
                self.remarks = text.strip()
            elif name == 'subroutine':
                sub = self.extract_subroutine(node, node.get_text())
            elif name == 'function':
                fn = self.extract_function(node, node.get_text())
            elif name == 'publicvars':
                self.extract_variables(node.get_text(), self.global_variables)
            elif name == 'privatevars':
                self.extract_variables(node.get_text(), self.private_variables)
            elif name == 'constants' or name == 'constantblock':
                self.extract_variables(node.get_text(), self.constants)
            elif name == 'publics':
                self.extract_variables(node.get_text(), self.publics)
            elif name == 'types':
                self.extract_types(node)
            elif name == 'interfaces':
                self.extract_interfaces(node.get_text())
            elif name == 'tex':
                print 'Warning: Skipping <tex> block in source code.'
            else:
                error(10, 'Unknown node type found: %s\n%s' % (name, node.get_text()))

    SUBROUTINE = re.compile("subroutine\s+(?P<name>.+?)\s*\((?P<params>.*?)\)\s*\n")

    ERR_DESCR = re.compile(r'(ERR_\w+)\s+(.*?)(?=ERR_|$)', re.DOTALL)

    def extract_subroutine(self, node, text):
        "Extract subroutine from source code."

        text = self.clean_code(text)

        m = re.search(self.SUBROUTINE, text)
        if m:
            params = map(lambda s: s.strip(), m.group("params").split(","))
            while '' in params:
                params.remove('')
            fn = Function(m.group("name"), params)
            self.subroutines.append(fn)
            
            for n in node:
                t = n.get_name()
                if t == 'input':
                    self.extract_variables(n.get_text(), fn.inputvars)
                elif t == 'inoutput' or t == 'inputoutput':
                    self.extract_variables(n.get_text(), fn.inoutputvars)
                elif t == 'output':
                    self.extract_variables(n.get_text(), fn.outputvars)
                elif t == 'result' or t == 'returns':
                    fn.result = self.clean_description(n.get_text())
                elif t == 'description':
                    text = self.clean_description(n.get_text())
                    fn.description = text.strip()
                elif t == 'remark':
                    text = self.clean_description(n.get_text())
                    fn.remarks = text.strip()
                elif t == 'remarks':
                    text = self.clean_description(n.get_text())
                    fn.remarks = text.strip()
                elif t == 'errors':
                    text = self.clean_description(n.get_text())
                    # filter out error tags with 'none'
                    if text.strip().find('none') >= 0:
                        fn.errors.append(('', 'none'))
                        continue
                    m = re.findall(self.ERR_DESCR, text)
                    if m:
                        for err_id, descr in m:
                            fn.errors.append((err_id, descr.strip()))
                elif t == 'globals':
                    m = re.findall(r'!\s*(\S+)', n.get_text())
                    if m:
                        fn.references.extend(m)
                else:
                    error(11, '<subroutine>: Unknown tag or tag in wrong environment: <%s>\n%s' % (n.get_name(), n.get_text()))
        else:
            error(12, 'Subroutine: Pattern did not match: %.120s...' % text)
                        
    FUNCTION = re.compile(r"\A\s*?(?P<retval>.+?)\s*?function\s+(?P<name>\w+)\s*\((?P<params>[^\)]*)\)\s*?")

    # regular expressions for code cleanup
    RE_XMLCOMMENT = re.compile(r'<!--.*?-->', re.DOTALL)
    RE_PREPROCESSOR = re.compile(r'^[ \t]*#\w+.*?$', re.MULTILINE)
    RE_MULTILINES = re.compile("\s*&\s*\n\s*")

    def clean_code(self, text):
        # remove XML comments
        text = re.sub(self.RE_XMLCOMMENT, '', text)
        # remove preprocessor directions
        text = re.sub(self.RE_PREPROCESSOR, '', text)
        # substitute multilines with one line
        text = re.sub(self.RE_MULTILINES, "", text)
        return text
        

    # regular expressions for description cleanup
    # RE_BEGIN_EXCLMARKS = re.compile(r'^[ \t]*![ \t]*', re.MULTILINE)
    RE_BEGIN_EXCLMARKS = re.compile(r'^[ \t]*![ \t]?', re.MULTILINE)
    RE_ONLY_HASH = re.compile(r'^[ \t]*#+[ \t]*$', re.MULTILINE)
    RE_BEGIN_HASH = re.compile(r'^[ \t]*#[ \t]?', re.MULTILINE)
    RE_END_HASH = re.compile(r'[ \t]*#[ \t]*$', re.MULTILINE)

    def clean_description(self, text):
        # remove XML comments
        text = re.sub(self.RE_XMLCOMMENT, '', text)
        # remove ! marks at the beginning of line
        text = re.sub(self.RE_BEGIN_EXCLMARKS, '', text)
        # remove # marks at the beginning and end of line
        text = re.sub(self.RE_ONLY_HASH, '', text)
        text = re.sub(self.RE_BEGIN_HASH, '', text)
        text = re.sub(self.RE_END_HASH, '', text)
        return text

    def extract_function(self, node, text):
        "Extract function from source code."

        text = self.clean_code(text)
        
        m = re.search(self.FUNCTION, text)
        if m:
            params = map(lambda s: s.strip(), m.group("params").split(","))
            while '' in params:
                params.remove('')
            fn = Function(m.group("name"), params, m.group("retval").strip())
            self.functions.append(fn)
            
            for n in node:
                t = n.get_name()
                if t == 'input':
                    self.extract_variables(n.get_text(), fn.inputvars)
                elif t == 'inoutput' or t == 'inputoutput':
                    self.extract_variables(n.get_text(), fn.inoutputvars)
                elif t == 'output':
                    self.extract_variables(n.get_text(), fn.outputvars)
                elif t == 'result' or t == 'returns':
                    fn.result = n.get_text().replace("!", "").strip()
                elif t == 'description':
                    text = self.clean_description(n.get_text())
                    fn.description = text.strip()
                elif t == 'remark':
                    text = self.clean_description(n.get_text())
                    fn.remarks = text.strip()
                elif t == 'remarks':
                    text = self.clean_description(n.get_text())
                    fn.remarks = text.strip()
                elif t == 'errors':
                    text = self.clean_description(n.get_text())
                    # filter out error tags with 'none'
                    if text.strip().find('none') >= 0:
                        fn.errors.append(('', 'none'))
                        continue
                    m = re.findall(self.ERR_DESCR, text)
                    if m:
                        for err_id, descr in m:
                            fn.errors.append((err_id, descr.strip()))
                elif t == 'globals':
                    m = re.findall(r'!\s*(\S+)', n.get_text())
                    if m:
                        fn.references.extend(m)
                else:
                    error(13, '<function>: Unknown tag or tag in wrong environment: <%s>\n%s' % (n.get_name(), n.get_text()))
        else:
            error(14, 'Function: Pattern did not match: %.120s...' % text)
      
    VARIABLE = re.compile(r"(?:^\s*!(.+?)$)|(?:^\s*(.+)\s*::\s*(.+?)\s*(?:=\s*(.+)\s*)?$)", re.MULTILINE)

    def extract_variables(self, text, l):
        "Extract variable and group information from source code."

        text = self.clean_code(text)

        DESCRIPTION = re.compile(r'description\s*=\s*"(.*?)"', re.DOTALL)
        self.active_group = None
        self.current_comment = ''

        # some sub functions for use with sre.scanner
        def begin_group(scanner, token):
            m = re.search(DESCRIPTION, token, re.DOTALL)
            if self.active_group:
                error(15, 'Found invalid <group>: %s' % token)
            if not m:
                error(16, 'Group does not have a valid description: %s' % token)
            self.active_group =  Group(m.group(1))
            l.append(self.active_group)

        def end_group(scanner, token):
            if not self.active_group.variables:
                error(15, 'No variable found in group: %s' % self.active_group.description)
            self.active_group = None
            self.current_comment = ''

        def commentblock(scanner, token):
            self.current_comment += token.strip()

        def blank_lines(scanner, token):
            self.current_comment = ''

        def variable_definition(scanner, token):
            match = re.search(self.VARIABLE, token)

            if match:
                m = match.groups()
                vartype = m[1].split(",")[0].strip()
                if self.current_comment:
                    self.current_comment = self.clean_description(self.current_comment)
                    v = Variable(m[2], vartype, self.current_comment)
                else:
                    v = Variable(m[2], vartype)
                dim = re.search(r"(dimension\(.*?\))", m[1])
                if dim:
                    v.dimension = dim.group(1)
                opt = re.search(r"optional", m[1])

                if opt:
                    v.optional = True
                if self.active_group:
                    self.active_group.variables.append(v)
                else:
                    l.append(v)
                self.current_comment = ''

        def dump(scanner, token):
            print "[%s]\n" % token,

        VARSCANNER = Scanner([
            (r'(?m)^[ \t]*?\n', blank_lines),
            (r'(?m)^[ \t]*![^\n]*?\<[ \t]*\/group[^\>]*\>[^\n$]*', end_group),
            (r'(?m)^[ \t]*![^\<]*\<[ \t]*group[ \t]+[^\>]*\>[^\n]*\n', begin_group),
            (r'(?m)([ \t]*![^\n]*\n)+', commentblock),
            (r'(?m)[^\n!]*::[^\n]*\n', variable_definition),
            #(r'(?m)^[^!\n]*?::.*?$', variable_definition),
            (r'(?m)\n', None),
            (r'.', None)
            #(r'(?m)\n', None),
        ])
        #], re.MULTILINE)
        VARSCANNER.scan(text)

    INTERFACE = re.compile(r'((?:^\s*![^\n]*\n)+)?(?:^\s*interface\s+(\w+)[ \t]*\n)', re.MULTILINE)

    def extract_interfaces(self, text):
        "Extract interfaces information"

        m = re.findall(self.INTERFACE, text)
        for comment, interface in m:
            comment = self.clean_description(comment)
            self.interfaces.append((interface, comment))

    def extract_types(self, node):
        "Extract information from <types>."
        for n in node:
            name = n.get_name()
            if name == 'type':
                self.extract_type(n.get_text())

    TYPE = re.compile(r'((?:\A\s*!.*\n)*)\s+type(?:,[ \t]*private)?(?:\s*\:\:)?\s+(?P<type>\S*?)\s*?\n\s+?(.*)\s+end\s+type\s+(?P=type)?', re.DOTALL)

    def extract_type(self, text):
        "Extract type information from source."

        text = self.clean_code(text)
        m = re.match(self.TYPE, text)
        if m:
            description = m.group(1)
            description = re.sub("^\s*!\s*", "", description)
            description = re.sub("\n\s*!\s*", "\n", description).strip()
            name = m.group("type")
            typedef = Type(name, description)
            self.extract_variables(m.group(3), typedef.variables)
            self.types.append(typedef)
        else:
            error(10, 'No type in <types>: %s' % text)

    def __str__(self):
        s = "Module %s\n" % self.name
        s += ("Purpose:\n%s\n" % self.purpose)
        for var in self.global_variables:
            s += "global " + str(var) + "\n"
        for var in self.private_variables:
            s += "private " + str(var) + "\n"
        for subroutine in self.subroutines:
            s += str(subroutine)+"\n"
            for var in subroutine.inputvars:
                s += "   IN: " + str(var) + "\n"
            for var in subroutine.inoutputvars:
                s += "INOUT: " + str(var) + "\n"
            for var in subroutine.outputvars:
                s += "  OUT: " + str(var) + "\n"
        for function in self.functions:
            s += str(function)+"\n"
            for var in subroutine.inputvars:
                s += "   IN: " + str(var) + "\n"
            for var in subroutine.inoutputvars:
                s += "INOUT: " + str(var) + "\n"
            for var in subroutine.outputvars:
                s += "  OUT: " + str(var) + "\n"
        return s

# =================================================================

def main():
    from optparse import OptionParser
    from sys import exit
    from sys import argv
    from os import path

    # define and parse command line options
    op = OptionParser()
    typeopts = op.add_option_group('output types')
    typeopts.add_option('-t', '--tex', action='store_true', dest='texoutput', default=False, help='generate LaTeX output')
    typeopts.add_option('-H', '--html', action='store_true', dest='htmloutput', default=False, help='generate HTML output')
    typeopts.add_option('-x', '--xml', action='store_true', dest='xmloutput', default=False, help='generate XML output')

    op.add_option('-g', '--generate-latex-definitions', dest='headerfile', default=None, help='generate LaTeX definitions')
    op.add_option('-d', '--output-dir', dest='outputdir', default=None, help='write all output files to the specified directory')
    op.add_option('-s', '--strip-extension', action='store_true', dest='stripextension', default=False, help='strip extension from the filename')
    op.add_option('-f', '--force', action='store_true', dest='force', default=False, help='continue after errors')
    op.set_usage('%prog <output types> [options] <files>')
    op.set_description('A parser for FEAST documentation.')

    # options are globally accessible
    global options
    (options, args) = op.parse_args()

    # any filenames given?
    if len(args) == 0 and not options.headerfile:
        print '%s: no input files given. Use -h for help.\n' % argv[0]
        exit(1)

    # does the output directory exist?
    if options.outputdir:
        if not path.exists(options.outputdir):
            error(2, 'Unable to access output directory: %s' % options.outputdir)

    # generate file with tex definitions
    if options.headerfile:
        print 'Writing tex definitions to %s.' % options.headerfile 
        try:
            try:
                f = open(options.headerfile, "w")
                latexexporter = LaTeXExporter()
                latexexporter.write_definitions(f)
            except e:
                error(3, 'Unable to write header file: %s' % e)
        finally:
            f.close()

    # for each filename: parse and output in the different formats
    for filename in args:
        # extract name and directory
        basedir, basefilename = path.split(filename)

        if options.outputdir:
            basedir = options.outputdir

        # extract information from file
        print 'Analysing %s...' % filename
        try:
            try:
                # read sourcecode
                f = open(filename)
                sourcecode = f.read()
            except IOError:
                error(4, 'Error while reading %s.' % filename)
        finally:
            if f: f.close()
    
        # invoke parser
        parser = FortranParser()
        nodes = parser.parse(sourcecode)[0]
        module = Module()
        module.filename = filename
        module.interpret(nodes)

        # compute path+basefilename
        if options.stripextension:
            basefilename = path.splitext(basefilename)[0] or basefilename
        fn = path.join(basedir, basefilename)

        # XML export
        if options.xmloutput:
            print "Generating %s.xml..." % fn
            try:
                try:
                    f = open(fn+".xml", "w")
                    xmlexporter = XMLExporter()
                    xmlexporter.export(module, f)
                except IOError, e:
                    error(5, 'Unable to write %s: %s' % (fn, e))
            finally:
                if f: f.close()

        # HTML export
        if options.htmloutput:
            print "Generating %s.html..." % fn
            try:
                try:
                    f = open(fn+".html", "w")
                    htmlexporter = HTMLExporter()
                    htmlexporter.export(module, f)
                except IOError, e:
                    error(6, 'Unable to write %s: %s' % (fn, e))
            finally:
                if f: f.close()
        
        # LaTeX export
        if options.texoutput:
            print "Generating %s.tex..." % fn
            try:
                try:
                    f = open(fn+".tex", "w")
                    latexexporter = LaTeXExporter()
                    latexexporter.export(module, f)
                except IOError, e:
                    error(7, 'Unable to write %s: %s' % (fn, e))
            finally:
                f.close()

        print

if __name__ == "__main__":
    main()

