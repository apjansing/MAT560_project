#!/bin/bash
DESTDIR=`(python -m site --user-site)`
rm -r $DESTDIR/polypy
cp -r polypy $DESTDIR/
echo $DESTDIR/polypy
ls $DESTDIR/polypy
echo "polypy model installed."
