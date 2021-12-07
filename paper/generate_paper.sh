usage="usage: compiler"

compiler=${1:?"Error. You must specify pdflatex or latex."}
texfile=${2:?"Error. You must specify a TeX document."}

case "$compiler" in
  'latex')
    latex $texfile
    bibtex $texfile  
    latex $texfile
    latex $texfile
    ;;
  'pdflatex')
    pdflatex $texfile
    bibtex $texfile  
    pdflatex $texfile
    pdflatex $texfile
    ;;
esac
