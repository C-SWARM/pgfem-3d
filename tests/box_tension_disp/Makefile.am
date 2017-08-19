noinst_SCRIPTS=1.cpu.sh \
               2.cpu.sh \
               4.cpu.sh \
               8.cpu.sh \
              16.cpu.sh 

TESTS=1.cpu.sh \
      2.cpu.sh \
      4.cpu.sh \
      8.cpu.sh \
     16.cpu.sh 

args=tension_disp stress_tension_disp.out -SS -disp -V

1.cpu.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "../box.sh 1 ${args}" >> $@
	chmod +x $@

2.cpu.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "../box.sh 2 ${args}" >> $@
	chmod +x $@

4.cpu.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "../box.sh 4 ${args}" >> $@
	chmod +x $@

8.cpu.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "../box.sh 8 ${args}" >> $@
	chmod +x $@

16.cpu.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "../box.sh 16 ${args}" >> $@
	chmod +x $@
