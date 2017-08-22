noinst_SCRIPTS = 1.cpu.lt.sh \
                 1.cpu.qt.sh \
                 2.cpu.lt.sh \
                 2.cpu.qt.sh \
                 4.cpu.lt.sh \
                 4.cpu.qt.sh \
                 8.cpu.lt.sh \
                 8.cpu.qt.sh \
                16.cpu.lt.sh \
                16.cpu.qt.sh 

TESTS = 1.cpu.lt.sh \
        1.cpu.qt.sh \
        2.cpu.lt.sh \
        2.cpu.qt.sh \
        4.cpu.lt.sh \
        4.cpu.qt.sh \
        8.cpu.lt.sh \
        8.cpu.qt.sh \
       16.cpu.lt.sh \
       16.cpu.qt.sh 

args=shear stress.out -SS -ms -disp -V

1.cpu.lt.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "../common.sh 1 LT ${args}" >> $@
	chmod +x $@

2.cpu.lt.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "../common.sh 2 LT ${args}" >> $@
	chmod +x $@

4.cpu.lt.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "../common.sh 4 LT ${args}" >> $@
	chmod +x $@

8.cpu.lt.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "../common.sh 8 LT ${args}" >> $@
	chmod +x $@

16.cpu.lt.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "../common.sh 16 LT ${args}" >> $@
	chmod +x $@

1.cpu.qt.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "../common.sh 1 QT ${args}" >> $@
	chmod +x $@

2.cpu.qt.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "../common.sh 2 QT ${args}" >> $@
	chmod +x $@

4.cpu.qt.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "../common.sh 4 QT ${args}" >> $@
	chmod +x $@

8.cpu.qt.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "../common.sh 8 QT ${args}" >> $@
	chmod +x $@

16.cpu.qt.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "../common.sh 16 QT ${args}" >> $@
	chmod +x $@
