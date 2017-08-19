noinst_SCRIPTS=1.cpu.lt.sh \
               2.cpu.lt.sh \
               4.cpu.lt.sh \
               8.cpu.lt.sh \
              16.cpu.lt.sh 

TESTS=1.cpu.lt.sh \
      2.cpu.lt.sh \
      4.cpu.lt.sh \
      8.cpu.lt.sh \
     16.cpu.lt.sh 

args=tension_pressure stress_tension_pressure.out -SS -disp -V

1.cpu.lt.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "../box.sh 1 LT ${args}" >> $@
	chmod +x $@

2.cpu.lt.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "../box.sh 2 LT ${args}" >> $@
	chmod +x $@

4.cpu.lt.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "../box.sh 4 LT ${args}" >> $@
	chmod +x $@

8.cpu.lt.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "../box.sh 8 LT ${args}" >> $@
	chmod +x $@

16.cpu.lt.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "../box.sh 16 LT ${args}" >> $@
	chmod +x $@
