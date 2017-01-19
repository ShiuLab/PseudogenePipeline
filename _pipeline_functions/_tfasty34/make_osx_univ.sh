#!/bin/csh

make -f Makefile.os_x all
make -f Makefile.os_x install
make -f Makefile.os_x clean-up

make -f Makefile.os_x86 all
make -f Makefile.os_x86 install
make -f Makefile.os_x86 clean-up

foreach n ( ppc/* )
set f=$n:t
lipo -create ppc/$f i386/$f -output bin/$f
echo "Universal $f built"
end
echo "Done!"

