% -*- LaTeX -*-

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{planckDoc}
\RequirePackage{ifthen,amssymb,graphicx,rotating,mathptm}
\newboolean{@color}
\setboolean{@color}{false}
\newboolean{@german}
\setboolean{@german}{false}
\newboolean{@status}
\setboolean{@status}{false}
\newboolean{@mpac}
\setboolean{@mpac}{false}
\newboolean{@pdf}
\setboolean{@pdf}{false}

\DeclareOption{color}{\setboolean{@color}{true}}
\DeclareOption{german}{\setboolean{@german}{true}}
\DeclareOption{status}{\setboolean{@status}{true}}
\DeclareOption{mpac}{\setboolean{@mpac}{true}}
\DeclareOption{pdf}{\setboolean{@pdf}{true}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions
\LoadClass[a4paper]{article}
\RequirePackage{a4wide}

\ifthenelse{\boolean{@pdf}}
  {\RequirePackage{color}
   \RequirePackage[pdftex]{hyperref}}
  {\RequirePackage[dvips]{color}}

\ifthenelse{\boolean{@german}}
  {\RequirePackage{german}
   \RequirePackage[iso]{umlaute}}{}

\renewcommand{\textfraction}{0}
\newcommand{\@iconfile}{}
\newcommand{\@iconpath}{.}
\renewcommand{\@title}{}
\newcommand{\@shorttitle}{}
\renewcommand{\@author}{}
\renewcommand{\@date}{}
\newcommand{\@authorised}{}
\newcommand{\@docid}{}
\newcommand{\@docrv}{}
\newcommand{\@abstract}{}

\renewcommand{\title}[1]{\renewcommand{\@title}{#1}}
\newcommand{\shorttitle}[1]{\renewcommand{\@shorttitle}{#1}}
\renewcommand{\author}[1]{\renewcommand{\@author}{#1}}
\renewcommand{\date}[1]{\renewcommand{\@date}{#1}}
\newcommand{\authorised}[1]{\renewcommand{\@authorised}{#1}}
\newcommand{\docid}[1]{\renewcommand{\@docid}{#1}}
\newcommand{\docrv}[1]{\renewcommand{\@docrv}{#1}}
\renewcommand{\abstract}[1]{\renewcommand{\@abstract}{#1}}

\newsavebox{\pbox}
\newsavebox{\appbox}
\newsavebox{\statbox}
\newsavebox{\changebox}

\newcommand{\approvalpage}[2]{%
  \sbox{\appbox}{\begin{minipage}{\textwidth}%
  \centerline{\textbf{\large Document Approval}}%
  \vspace{1in}\begin{center}%
  \begin{tabular}
    {|p{0.2\hsize}|p{0.2\hsize}|p{0.2\hsize}|p{0.2\hsize}|}%
  \hline%
  \textbf{Approved by} &%
  \textbf{Organisation} &%
  \textbf{Signature} &%
  \textbf{Date} \\%
  \hline#1\hline%
  \end{tabular}%
  \end{center}\vspace{1in}%
  \noindent\textbf{Distribution:}~#2%
  \end{minipage}}}

\newcommand{\statuspage}[1]{%
  \sbox{\statbox}{\begin{minipage}{\textwidth}%
  \centerline{\textbf{\large Document Status Sheet}}%
  \vspace{1in}%
  \begin{center}%
  \begin{tabular}
	  {|p{0.15\hsize}|p{0.15\hsize}|p{0.15\hsize}|p{0.4\hsize}|}%
  \hline%
  \multicolumn{4}{|l|}
    {\textbf{1.~Document Title:}~\@title}\\%
  \hline%
  \textbf{2.~Issue} &%
  \textbf{3.~Revision} &%
  \textbf{4.~Date} &%
  \textbf{5.~Reason for Change} \\%
  \hline#1\hline\end{tabular}%
  \end{center}\end{minipage}}}

\newcommand{\changerecord}[1]{%
  \sbox{\changebox}{%
  \begin{minipage}{\textwidth}%
  \centerline{\textbf{\large Document Change Record}}%
  \vspace{1in}%
  \begin{center}%
  \begin{tabular}{|p{0.2\hsize}|p{0.7\hsize}|}%
  \hline%
  \multicolumn{2}{|l|}
    {\textbf{1.~Document Title:}~\@title}\\%
  \multicolumn{2}{|l|}
    {\textbf{2.~Document Reference Number:}~\@docid}\\%
  \multicolumn{2}{|l|}
    {\textbf{3.~Document Issue/Revision Number:}~\@docrv}\\%
  \hline
  \textbf{4.~Section} & \textbf{5.~Reason for Change}\\%   
  \hline#1\hline\end{tabular}%
  \end{center}\end{minipage}}}

\newlength{\bigicon}
\setlength{\bigicon}{1in}
\newlength{\smallicon}
\ifthenelse{\boolean{@mpac}}
  {\setlength{\smallicon}{0.75\bigicon}
   \ifthenelse{\boolean{@pdf}}
     {\renewcommand{\@iconfile}{\@iconpath/mpacLogo.jpg}}
     {\renewcommand{\@iconfile}{\@iconpath/mpacLogo.eps}}}
  {\setlength{\smallicon}{0.5\bigicon}
   \ifthenelse{\boolean{@pdf}}
     {\renewcommand{\@iconfile}{\@iconpath/planckESA.jpg}}
     {\renewcommand{\@iconfile}{\@iconpath/planckESA.eps}}}
\setlength{\parindent}{0pt}
\setlength{\parskip}{0.75ex plus 0.2ex minus 0.1ex}
\addtolength{\headheight}{\smallicon}
\addtolength{\textheight}{-\smallicon}
\newlength{\iconsep}
\setlength{\iconsep}{6em}
\addtolength{\iconsep}{1ex}
\addtolength{\iconsep}{-\smallicon}

\newcounter{creq}
\newcounter{savecreq}
\newcommand{\creset}{\setcounter{creq}{\thesavecreq}}
\newenvironment{requirements}[1]{\begin{list}
  {\ifthenelse{\boolean{@color}}{\color{red}}{}
   \bfseries\upshape%
   #1-\thesubsection-\thecreq\hfill}
  {\usecounter{creq}%
   \setlength{\labelwidth}{6em}%
   \setlength{\labelsep}{1ex}%
   \setlength{\leftmargin}{\labelwidth}%
   \addtolength{\leftmargin}{\labelsep}
   \setlength{\parsep}{0.5ex plus 0.2ex minus 0.1ex}
   \setlength{\itemsep}{0ex plus 0.2ex}}}
  {\setcounter{savecreq}{\thecreq}\end{list}}

\pagestyle{myheadings}
\markboth{\usebox{\pbox}}{\usebox{\pbox}}

\newcommand{\frontpage}{%
  \sbox{\pbox}{\underline{%
    \includegraphics[width=\smallicon]{\@iconfile}%
    \hspace{\iconsep}%
    \ifthenelse{\equal{\@shorttitle}{}}
      {\large\rm\@title}
      {\large\rm\@shorttitle}}}
  \begin{titlepage}
  \begin{center}
    \Large{\textrm{\@title}}\\
    \vspace{3\baselineskip}
    \includegraphics[width=\bigicon]{\@iconfile}\\
    \large
    \vspace{3\baselineskip}
    \renewcommand{\arraystretch}{3}
    \begin{tabular}{@{}p{0.2\hsize}p{0.6\hsize}@{}}
    \ifthenelse{\boolean{@german}}
      {Dokument: & \@docid \\
       Version: & \@docrv; \@date \\
       Zusammen\-fassung: & \@abstract \\
       Autor: & \@author \\
       Autorisiert von: & \@authorised \\}
      {Document: & \@docid \\
       Revision: & \@docrv; \@date \\
       Abstract: & \@abstract \\
       Prepared by: & \@author \\
       Authorised by: & \@authorised \\}
    \end{tabular}
  \end{center}
  \end{titlepage}
  \clearpage
  \ifthenelse{\boolean{@status}}
    {\usebox{\appbox}\clearpage%
     \usebox{\statbox}\clearpage%
     \usebox{\changebox}\clearpage}{}
  \tableofcontents
  \clearpage}
