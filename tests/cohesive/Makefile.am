noinst_SCRIPTS = 1_cpu.sh \
                 2_cpu.sh \
                 4_cpu.sh \
                 8_cpu.sh \
                16_cpu.sh

TESTS = 1_cpu.sh \
        2_cpu.sh \
        4_cpu.sh \
        8_cpu.sh \
       16_cpu.sh


# args=shear stress.out -SS -ms -disp -V

1_cpu.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "./common.sh 1" >> $@
	chmod +x $@

2_cpu.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "./common.sh 2" >> $@
	chmod +x $@

4_cpu.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "./common.sh 4" >> $@
	chmod +x $@

8_cpu.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "./common.sh 8" >> $@
	chmod +x $@

16_cpu.sh:
	echo "#!/usr/bin/env bash" > $@
	echo "./common.sh 16" >> $@
	chmod +x $@
