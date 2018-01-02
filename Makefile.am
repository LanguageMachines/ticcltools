
ACLOCAL_AMFLAGS =-I m4 --install

SUBDIRS = include src m4 docs

EXTRA_DIST = bootstrap.sh AUTHORS TODO NEWS README.md

ChangeLog: NEWS
	git pull; git2cl > ChangeLog
