PKG:=docsrc

docsrc_idx: | $(DOCDIR)_mkdir
	cp $(SRCROOT)/docsrc/index_main.html $(DOCDIR)/index.html

docsrc_manual: docsrc_idx
	cd $(SRCROOT)/docsrc; \
	pdflatex simulation; \
	pdflatex simulation; \
	pdflatex simulation; \
	mv simulation.pdf $(DOCDIR)/; \
	pdflatex paper; \
	bibtex paper; \
	pdflatex paper; \
	pdflatex paper; \
	pdflatex paper; \
	mv paper.pdf $(DOCDIR)/; \
	rm -f *.aux *.out *.toc *.log *.tag *.bbl *.blg simulation.pdf paper.pdf

docsrc_code_doc: docsrc_idx
	rm -rf $(DOCDIR)/code_doc
	mkdir $(DOCDIR)/code_doc
	cd $(SRCROOT)/docsrc; \
	for i in c_utils libfftpack libsharp cxxsupport Healpix_cxx Modules; do \
	  doxygen $${i}.dox; \
	  cp planckESA.jpg htmldoc/; \
	  mv htmldoc $(DOCDIR)/code_doc/$${i}; \
	done; \
	rm *.tag; \
	cp index_code.html $(DOCDIR)/code_doc/index.html

docsrc_clean:
	cd $(SRCROOT)/docsrc; \
	rm -f *.aux *.out *.toc *.log *.tag *.bbl *.blg simulation.pdf paper.pdf
	cd $(SRCROOT)/docsrc; \
	rm -rf htmldoc

doc: docsrc_manual docsrc_code_doc
