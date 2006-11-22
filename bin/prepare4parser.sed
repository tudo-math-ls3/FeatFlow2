# Join continuation lines
# (i.e. lines that end with shell's continuation character '\')
# but ignore those lines that end with LaTeX force-linebreak: '\\'
#
# define label 'x' to jump back to it whenever a line has been joined in
# order to be able to catch also several subsequent continuation lines.
:x
# if line ends with continuation character...
/[^\][^\]\\$/ {
# ... then get next line from input and append it to pattern space, ...
N
# ... remove continuation character, Fortran 90 comment sign and following 
# white space.
s/[ 	]*\\\n[ 	]*![ 	]*/ /
# If a substitution has been made go back to label 'x'
tx
}

# Remove multi-line continuation character at the end of a line
# as they might be misinterpreted as column separator in a LaTeX 
# table environment.
s/&[ 	]*$//g;

# All remaining ampersands are converted into the according XML entity.
# (it would confuse the java XML parser as '&' starts a
#  entity reference in XML.)
s/\([ 	]*\)&\([ 	]*\)/\1\&amp;\2/g;

# A token looking like the start of a tag does actually not start one if it
# * does not start with "<" plus a letter
# * is not an XML comment: <!--
# * is not a closing tag: </...>
s/<\([^a-zA-Z?!\/]\)/\\ensuremath{\&lt;}\1/g;

# A token looking like the end of a tag does actually not end one if it
# * is preceded by white space, followed optionally by an equal sign and 
#   following by white space
s/\([ 	]\{1,\}\)>\([=]*\)\([ 	]\{1,\}\)/\1\\ensuremath{\&gt;\2}\3/;

# A token looking like the end of a tag does actually not end one if it
# * is not an XML end comment marker: -->
# * does not end with a letter plus ">"
s/\([^-][^"a-zA-Z?!\/]\)>/\1\&gt;/g